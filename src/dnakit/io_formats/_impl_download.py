"""Resumable HTTP(S) download using Range requests (standard library only)."""

from __future__ import annotations

import hashlib
import http.client
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Mapping
from pathlib import Path

from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.io_formats._impl_common import Opts, bad_input, bad_option

OPTIONS = {
    "url": "http(s) URL (alternative to passing it as input)",
    "output": "destination file path (required unless input is a mapping with 'output')",
    "overwrite": "bool, default False: replace an existing destination file",
    "chunk_size": "int, default 65536: bytes read per step",
    "timeout": "number of seconds, default 30: socket timeout",
    "max_bytes": "int, default 1073741824: refuse to download more than this",
    "expected_sha256": "hex digest the finished file must have (a mismatch deletes the partial file)",
}
SEMANTICS = (
    "Download a file over HTTP(S) into 'output', resuming an interrupted transfer: bytes already stored in '<output>.part' are kept "
    "and requested again with 'Range: bytes=<n>-' (guarded by If-Range with the ETag/Last-Modified saved in '<output>.part.meta'). "
    "A server that ignores Range (HTTP 200) makes the download restart from byte 0; HTTP 416 with a matching size means the part "
    "file is already complete. The file is renamed to 'output' only after its size (and expected_sha256, if given) is verified; "
    "a transfer cut short keeps the '.part' file and raises DownloadError so that calling the function again resumes it."
)


def _build_opener() -> urllib.request.OpenerDirector:
    """An opener that speaks only http(s): no file:, ftp: or data: handlers, also for redirects."""

    opener = urllib.request.OpenerDirector()
    for handler in (
        urllib.request.HTTPHandler(),
        urllib.request.HTTPSHandler(),
        urllib.request.HTTPRedirectHandler(),
        urllib.request.HTTPDefaultErrorHandler(),
        urllib.request.HTTPErrorProcessor(),
        urllib.request.UnknownHandler(),
    ):
        opener.add_handler(handler)
    return opener


_OPENER = _build_opener()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def _total_from_range(header: str | None) -> int | None:
    match = re.fullmatch(r"\s*bytes\s+(?:(\d+)-(\d+)|\*)/(\d+|\*)\s*", header or "")
    return int(match.group(3)) if match and match.group(3).isdigit() else None


def resume_download(input: object, opts: Opts) -> dict[str, object]:
    raw: dict[str, object] = dict(input) if isinstance(input, Mapping) else {}
    url = raw.get("url") or opts.text("url") or (input if isinstance(input, str) else None)
    opts.text("url")
    out_value = raw.get("output") or opts.text("output")
    overwrite = opts.boolean("overwrite", False)
    chunk = opts.integer("chunk_size", 65536, 1024, 16 * 1024 * 1024) or 65536
    timeout = opts.number("timeout", 30.0, 0.1, 3600.0) or 30.0
    max_bytes = opts.integer("max_bytes", 1 << 30, 1) or (1 << 30)
    expected = opts.text("expected_sha256")
    if (
        not isinstance(url, str)
        or urllib.parse.urlsplit(url).scheme not in ("http", "https")
        or not urllib.parse.urlsplit(url).netloc
    ):
        raise bad_input(
            "resume_download needs an http:// or https:// URL (as input or in options['url']).",
            code="IO_FORMATS_BAD_URL",
        )
    if not out_value:
        raise bad_option("resume_download needs options['output']: the destination file path.")
    if expected is not None and not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
        raise bad_option("expected_sha256 must be a 64-digit hexadecimal string.")
    target = Path(str(out_value)).expanduser()
    if target.exists() and not overwrite:
        raise ConfigurationError(
            f"Refusing to overwrite existing file: {target}",
            code="IO_FORMATS_EXISTS",
            hint="Set overwrite=True or choose another output.",
        )
    if target.is_dir():
        raise ConfigurationError(
            f"Output path is a directory: {target}", code="IO_FORMATS_BAD_OUTPUT"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    meta_path = target.with_name(target.name + ".part.meta")

    def fail(message: str, code: str, hint: str | None = None) -> DownloadError:
        return DownloadError(message, code=code, context={"url": url}, hint=hint)

    attempts = 0
    restarted = False
    resumed_from = 0
    while True:
        attempts += 1
        offset = part.stat().st_size if part.exists() else 0
        meta: dict[str, str] = {}
        if offset and meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                meta = {}
        headers = {"Accept-Encoding": "identity", "User-Agent": "dnakit-io_formats"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
            validator = meta.get("etag") or meta.get("last_modified")
            if validator:
                headers["If-Range"] = validator
        request = urllib.request.Request(url, headers=headers)
        try:
            response = _OPENER.open(request, timeout=timeout)
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and offset and attempts == 1:
                total = _total_from_range(exc.headers.get("Content-Range"))
                exc.close()
                if total == offset:
                    resumed_from = offset
                    break  # the part file already holds the whole resource
                part.unlink(missing_ok=True)  # stale/oversized part file: start again
                restarted = True
                continue
            exc.close()
            raise fail(f"HTTP {exc.code} {exc.reason} for {url}", "DOWNLOAD_HTTP_ERROR") from exc
        except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
            raise fail(
                f"Cannot connect to {url}: {exc}",
                "DOWNLOAD_NETWORK_ERROR",
                "Check the URL and network; call again to resume.",
            ) from exc
        with response:
            status = response.status
            if status == 206:
                start_match = re.match(
                    r"\s*bytes\s+(\d+)-", response.headers.get("Content-Range", "")
                )
                if not start_match or int(start_match.group(1)) != offset:
                    raise fail(
                        "The server answered a Range request with an unexpected range.",
                        "DOWNLOAD_RANGE_MISMATCH",
                    )
                total = _total_from_range(response.headers.get("Content-Range"))
                resumed_from = offset
                mode = "ab"
            elif status == 200:
                length = response.headers.get("Content-Length")
                total = int(length) if length and length.isdigit() else None
                if offset:
                    restarted = True  # server ignored Range or the validator changed
                offset = 0
                resumed_from = 0
                mode = "wb"
            else:
                raise fail(f"Unexpected HTTP status {status} for {url}", "DOWNLOAD_HTTP_ERROR")
            if total is not None and total > max_bytes:
                raise fail(
                    f"The resource has {total} bytes, above max_bytes={max_bytes}.",
                    "DOWNLOAD_TOO_LARGE",
                    "Raise max_bytes if this is intended.",
                )
            validator = response.headers.get("ETag") or response.headers.get("Last-Modified")
            if validator:
                header_name = "etag" if response.headers.get("ETag") else "last_modified"
                meta_path.write_text(json.dumps({header_name: validator}), encoding="utf-8")
            written = offset
            try:
                with part.open(mode) as handle:
                    while True:
                        block = response.read(chunk)
                        if not block:
                            break
                        written += len(block)
                        if written > max_bytes:
                            raise fail(
                                f"The download exceeded max_bytes={max_bytes}.",
                                "DOWNLOAD_TOO_LARGE",
                            )
                        handle.write(block)
                        handle.flush()
            except (OSError, http.client.HTTPException) as exc:
                raise fail(
                    f"The transfer was interrupted after {written} bytes: {exc}",
                    "DOWNLOAD_NETWORK_ERROR",
                    "Call resume_download again with the same output to resume.",
                ) from exc
        if total is not None and written != total:
            raise fail(
                f"The transfer ended after {written} of {total} bytes.",
                "DOWNLOAD_LENGTH_MISMATCH",
                "Call resume_download again with the same output to resume.",
            )
        break
    size = part.stat().st_size
    digest = _sha256(part)
    if expected is not None and digest != expected.lower():
        part.unlink(missing_ok=True)
        meta_path.unlink(missing_ok=True)
        raise fail(
            f"SHA-256 mismatch: expected {expected.lower()}, got {digest}.",
            "DOWNLOAD_CHECKSUM_MISMATCH",
            "The partial file was deleted; download again.",
        )
    os.replace(part, target)
    meta_path.unlink(missing_ok=True)
    return {
        "output_path": str(target),
        "bytes": size,
        "sha256": digest,
        "restarted": restarted,
        "complete": True,
        "resumed_from": resumed_from,
    }
