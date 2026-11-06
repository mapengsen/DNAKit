"""Key-name and pattern based pseudonymisation of metadata (shared implementation).

Used by ``dnakit.workflow_security.pseudonymize_metadata_fields`` and
``dnakit.knowledge_graph.pseudonymize_genomic_metadata_fields``.

This is NOT de-identification. It only (1) replaces/drops/hashes the values of keys whose normalised name is in a
documented list, at any nesting depth, and (2) optionally replaces e-mail addresses, phone numbers and SSN-like
patterns inside string values. Dates, ages, small geographic units, names inside free text and quasi-identifiers
are not handled. The input is never mutated, and neither the result's bookkeeping, the warnings, the error
messages nor the envelope summary contain any of the original values.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import hmac
import json
import re
from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError

REDACTED = "[REDACTED]"
MODES = ("redact", "drop", "hash")
MIN_SALT_CHARS = 8
MIN_DIGEST_CHARS = 8
MAX_DIGEST_CHARS = 64
DEFAULT_DIGEST_CHARS = 16
MAX_DEPTH = 64
MAX_REPORTED_PATHS = 10_000

# Normalised form (lower case, only [a-z0-9]) of the key names treated as sensitive by default.
DEFAULT_SENSITIVE_KEYS = frozenset(
    {
        "name",
        "firstname",
        "lastname",
        "fullname",
        "givenname",
        "familyname",
        "surname",
        "middlename",
        "patient",
        "patientid",
        "patientname",
        "mrn",
        "medicalrecordnumber",
        "ssn",
        "socialsecuritynumber",
        "email",
        "emailaddress",
        "phone",
        "phonenumber",
        "telephone",
        "mobile",
        "mobilephone",
        "fax",
        "address",
        "homeaddress",
        "streetaddress",
        "street",
        "city",
        "zip",
        "zipcode",
        "postalcode",
        "postcode",
        "dob",
        "dateofbirth",
        "birthdate",
        "birthday",
        "nationalid",
        "passport",
        "passportnumber",
        "insurance",
        "insuranceid",
        "insurancenumber",
        "accountnumber",
        "ipaddress",
        "deviceid",
        "photo",
        "photograph",
        "biometric",
        "biometrics",
        "licensenumber",
        "driverlicense",
    }
)

_PATTERNS = (
    ("email", re.compile(r"[A-Za-z0-9._%+\-]+@[A-Za-z0-9\-]+(?:\.[A-Za-z0-9\-]+)*\.[A-Za-z]{2,}")),
    ("ssn", re.compile(r"(?<![\w-])\d{3}[- ]\d{2}[- ]\d{4}(?![\w-])")),
    (
        "phone",
        re.compile(r"(?<![\w+])\+\d{1,3}[\s.\-]?\(?\d{1,4}\)?(?:[\s.\-]?\d{2,4}){2,4}(?!\w)"),
    ),
    (
        "phone",
        re.compile(
            r"(?<![\w+])(?:\+?1[\s.\-]?)?(?:\(\d{3}\)\s?|\d{3}[\s.\-]?)\d{3}[\s.\-]?\d{4}(?!\w)"
        ),
    ),
)
_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

WARNINGS = (
    "Pseudonymisation is key-name and pattern based only. It is NOT a guarantee of de-identification and does "
    "not by itself satisfy HIPAA or GDPR requirements.",
    "Dates are not shifted or generalised; ages over 89 and small geographic units are not generalised.",
    "Free text may still identify a person: only e-mail, phone and SSN-like patterns are scanned, names and other "
    "identifiers inside text are not detected; dictionary keys that themselves contain identifiers are not rewritten.",
    "Review the output manually before sharing it.",
)


def normalize_key(key: object) -> str:
    """Lower-case a key and strip everything except letters and digits (``Patient_ID`` -> ``patientid``)."""

    return re.sub(r"[^a-z0-9]", "", str(key).casefold())


def _opts(options: object) -> dict[str, object]:
    if options is None:
        return {}
    if not isinstance(options, Mapping):
        raise ConfigurationError(
            "options must be a mapping or None.",
            code="INVALID_PSEUDONYMIZE_OPTIONS",
            hint="Pass options={'mode': 'redact'} or leave options out.",
        )
    return dict(options)


def _key_set(options: Mapping[str, object], name: str) -> frozenset[str]:
    value = options.get(name)
    if value is None:
        return frozenset()
    items = [value] if isinstance(value, str) else value
    if not isinstance(items, (list, tuple, set, frozenset)) or not all(
        isinstance(item, str) for item in items
    ):
        raise ConfigurationError(
            f"options['{name}'] must be a string or a list of strings.",
            code="INVALID_PSEUDONYMIZE_KEYS",
            hint=f"Example: options={{'{name}': ['sample_owner', 'clinic']}}.",
        )
    normalized = {normalize_key(item) for item in items}
    if "" in normalized:
        raise ConfigurationError(
            f"options['{name}'] contains a key without letters or digits.",
            code="INVALID_PSEUDONYMIZE_KEYS",
            hint="Give key names such as 'clinic' or 'Sample_Owner'.",
        )
    return frozenset(normalized)


class _Settings:
    def __init__(self, options: Mapping[str, object]) -> None:
        mode = options.get("mode", "redact")
        if not isinstance(mode, str) or mode not in MODES:
            raise ConfigurationError(
                "options['mode'] must be one of 'redact', 'drop' or 'hash'.",
                code="INVALID_PSEUDONYMIZE_MODE",
                hint="'redact' replaces the value by \"[REDACTED]\"; 'hash' needs options['salt'].",
            )
        self.mode = mode
        self.salt: bytes | None = None
        self.digest_chars = DEFAULT_DIGEST_CHARS
        if mode == "hash":
            salt = options.get("salt")
            if not isinstance(salt, str) or len(salt) < MIN_SALT_CHARS or not salt.strip():
                raise ConfigurationError(
                    f"mode 'hash' requires options['salt'], a string of at least {MIN_SALT_CHARS} characters.",
                    code="PSEUDONYMIZE_SALT_REQUIRED",
                    hint="Use a secret random string (for example secrets.token_hex(16)) and keep it private.",
                )
            self.salt = salt.encode("utf-8")
            digest = options.get("digest_chars", DEFAULT_DIGEST_CHARS)
            if (
                isinstance(digest, bool)
                or not isinstance(digest, int)
                or not MIN_DIGEST_CHARS <= digest <= MAX_DIGEST_CHARS
            ):
                raise ConfigurationError(
                    f"options['digest_chars'] must be an integer from {MIN_DIGEST_CHARS} to {MAX_DIGEST_CHARS}.",
                    code="INVALID_DIGEST_CHARS",
                    hint=f"The default is {DEFAULT_DIGEST_CHARS}.",
                )
            self.digest_chars = digest
        scan = options.get("scan_text", True)
        if not isinstance(scan, bool):
            raise ConfigurationError(
                "options['scan_text'] must be true or false.",
                code="INVALID_PSEUDONYMIZE_OPTIONS",
                hint="Use scan_text=False only if string values are already known to be free of e-mail/phone/SSN.",
            )
        self.scan_text = scan
        self.keep = _key_set(options, "keep")
        self.sensitive = (DEFAULT_SENSITIVE_KEYS | _key_set(options, "keys")) - self.keep


def _load(value: object) -> tuple[object, int]:
    """Return (mapping-or-list, record_count) for a mapping, list of mappings or JSON text of those."""

    hint = "Pass a dict, a list of dicts, or JSON text of one of them."
    if isinstance(value, (str, bytes, bytearray)):
        try:
            text = value.decode("utf-8") if not isinstance(value, str) else value
            value = json.loads(text)
        except (ValueError, UnicodeDecodeError):
            raise InputFormatError(
                "The metadata text is not valid JSON.", code="INVALID_METADATA_JSON", hint=hint
            ) from None
    if isinstance(value, Mapping):
        return value, 1
    if isinstance(value, (list, tuple)) and all(isinstance(item, Mapping) for item in value):
        return value, len(value)
    raise InputFormatError(
        "Metadata must be a mapping, a list of mappings, or JSON text of those.",
        code="INVALID_METADATA_INPUT",
        hint=hint,
    )


class _Run:
    def __init__(self, settings: _Settings) -> None:
        self.s = settings
        self.n_keys = 0
        self.n_key_matches = 0
        self.pattern_counts = {"email": 0, "phone": 0, "ssn": 0}
        self.paths: list[str] = []
        self._seen: set[str] = set()
        self.paths_truncated = False

    def note(self, path: str) -> None:
        if path in self._seen:
            return
        if len(self.paths) >= MAX_REPORTED_PATHS:
            self.paths_truncated = True
            return
        self._seen.add(path)
        self.paths.append(path)

    def token(self, value: object) -> str:
        if isinstance(value, str):
            message = value
        else:
            message = json.dumps(
                value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), default=str
            )
        digest = hmac.new(self.s.salt or b"", message.encode("utf-8"), hashlib.sha256).hexdigest()
        return digest[: self.s.digest_chars]

    def scan(self, text: str, path: str) -> str:
        for label, pattern in _PATTERNS:
            text, count = pattern.subn(f"[REDACTED:{label}]", text)
            if count:
                self.pattern_counts[label] += count
                self.note(path)
        return text

    def walk(self, value: object, path: str, depth: int) -> object:
        if depth > MAX_DEPTH:
            raise InputFormatError(
                f"Metadata is nested deeper than {MAX_DEPTH} levels (or contains a cycle).",
                code="METADATA_TOO_DEEP",
                hint="Flatten the metadata or remove self-references.",
            )
        if isinstance(value, Mapping):
            out: dict[object, object] = {}
            for key, item in value.items():
                self.n_keys += 1
                child = path + (
                    f".{key}"
                    if _IDENTIFIER.fullmatch(str(key))
                    else "['" + str(key).replace("\\", "\\\\").replace("'", "\\'") + "']"
                )
                if normalize_key(key) in self.s.sensitive:
                    self.n_key_matches += 1
                    self.note(child)
                    if self.s.mode == "drop":
                        continue
                    out[key] = REDACTED if self.s.mode == "redact" else self.token(item)
                else:
                    out[key] = self.walk(item, child, depth + 1)
            return out
        if isinstance(value, (list, tuple)):
            return [
                self.walk(item, f"{path}[{index}]", depth + 1) for index, item in enumerate(value)
            ]
        if isinstance(value, str):
            return self.scan(value, path) if self.s.scan_text else value
        if isinstance(value, (_dt.datetime, _dt.date)):
            text = value.isoformat()
            return self.scan(text, path) if self.s.scan_text else text
        if value is None or isinstance(value, (bool, int, float)):
            return value
        raise InputFormatError(
            f"Unsupported metadata value of type {type(value).__name__} at {path}.",
            code="UNSUPPORTED_METADATA_VALUE",
            hint="Use JSON types only: dict, list, str, int, float, bool, None.",
        )


def pseudonymize_metadata(input: object, options: object = None) -> dict[str, object]:
    """Return the pseudonymised metadata and its bookkeeping (see ``pseudonymize_envelope``)."""

    settings = _Settings(_opts(options))
    data, record_count = _load(input)
    run = _Run(settings)
    metadata = run.walk(data, "$", 0)
    n_text = sum(run.pattern_counts.values())
    warnings = list(WARNINGS)
    if not settings.scan_text:
        warnings.append(
            "scan_text is disabled: string values were not scanned for e-mail, phone or SSN-like patterns."
        )
    if settings.mode == "hash":
        warnings.append(
            "mode 'hash' gives equal tokens for equal values, so records stay linkable and low-entropy values "
            "can be guessed by anyone who knows the salt; keep the salt secret and never publish it with the data."
        )
    return {
        "metadata": metadata,
        "mode": settings.mode,
        "n_redacted": run.n_key_matches + n_text,
        "n_key_matches": run.n_key_matches,
        "n_text_matches": n_text,
        "text_matches_by_pattern": dict(run.pattern_counts),
        "redacted_paths": run.paths,
        "redacted_paths_truncated": run.paths_truncated,
        "warnings": warnings,
        "_summary": {"record_count": record_count, "n_keys": run.n_keys},
    }


def pseudonymize_envelope(
    module: str, operation: str, input: object, options: object = None
) -> object:
    """Run the shared implementation and wrap it in the standard candidate envelope.

    The envelope's ``input_summary`` only has ``record_count`` and ``n_keys`` (no text length, no hash of the
    unprotected input).
    """

    result = pseudonymize_metadata(input, options)
    summary = result.pop("_summary")
    return candidate_ok(
        module,
        operation,
        input,
        result,
        options=_opts(options),
        warnings=result["warnings"],
        implementation="native",
        input_summary=summary,
    )
