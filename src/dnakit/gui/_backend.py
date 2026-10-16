"""Event handlers behind the DNAKit graphical interface.

Single-sequence and dataset tools run the ``dnakit`` command line in a subprocess, so the page
can show the equivalent command; the function explorer runs registry tools in ``ToolWorker``.
"""

from __future__ import annotations

import csv
import json
import os
import pprint
import re
import shutil
import subprocess
import sys
import time
import uuid
from importlib import resources
from pathlib import Path
from typing import Any

import gradio as gr

from dnakit.gui._settings import GuiSettings
from dnakit.gui._worker import ToolWorker
from dnakit.skills import read_skill_file, search_skills, skill_catalog_stats
from dnakit.tools.registry import default_tool_registry

PREVIEW_FIELDS = (
    "symbol_length",
    "strong_fraction",
    "cpg_observed_expected",
    "shannon_entropy_k2_bits",
    "linguistic_complexity_k3",
    "nn_tm_c",
)
BLOCKED_ARGUMENT = re.compile(
    r"(^|_)(url|uri|endpoint|host|token|api_key|password|secret)s?$", re.IGNORECASE
)
EXAMPLE_SEQUENCE = (
    "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    "CGCGTAACGCGCTTATCGGCGCGATCGCGCATTTCGCGCGGATCGTAGCTAGCATCGGGCGCGCTAA"
)
EFFECT_NOTE = {
    "local": "runs locally",
    "network": "queries a public database",
    "filesystem": "writes files",
    "model": "downloads or loads model weights",
    "external_process": "runs an external program",
}
SAFE_EFFECTS = frozenset({"local", "network"})
CLI_BOOTSTRAP = "import sys; from dnakit.cli.app import run; sys.argv[0] = 'dnakit'; run()"


def example_fasta() -> Path:
    """Return the bundled 26-record example FASTA file."""

    return Path(str(resources.files("dnakit.gui").joinpath("examples/example.fasta")))


def fmt(value: Any) -> str:
    """Format one value for a table cell."""

    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (dict, list)):
        return json.dumps(value)[:200]
    return str(value)


def shorten(value: Any) -> Any:
    """Abbreviate long strings and lists for the code snippets."""

    if isinstance(value, str) and len(value) > 70:
        return value[:32] + "..." + value[-12:]
    if isinstance(value, dict):
        return {k: shorten(v) for k, v in value.items()}
    if isinstance(value, list):
        return [shorten(v) for v in value[:20]] + (["..."] if len(value) > 20 else [])
    return value


def hide_svg(value: Any) -> Any:
    """Replace inline SVG text in the JSON view; the figure itself is drawn above it."""

    if isinstance(value, str) and "<svg" in value[:400]:
        return f"<svg figure, {len(value):,} characters, shown above>"
    if isinstance(value, dict):
        return {k: hide_svg(v) for k, v in value.items()}
    if isinstance(value, list):
        return [hide_svg(v) for v in value]
    return value


def find_svg(value: Any) -> str | None:
    """Return the first inline SVG document inside a result."""

    if isinstance(value, str) and "<svg" in value[:400]:
        return value
    items = list(value.values()) if isinstance(value, dict) else value
    if isinstance(items, list):
        for item in items:
            found = find_svg(item)
            if found:
                return found
    return None


def example_value(name: str, schema: dict[str, Any]) -> Any:
    """Return a starting value for one required argument of a tool schema."""

    if "default" in schema:
        return schema["default"]
    if schema.get("enum"):
        return schema["enum"][0]
    kinds = {schema.get("type")}
    for sub in schema.get("anyOf", []):
        if isinstance(sub, dict):
            kinds.add(sub.get("type"))
            kinds |= {s.get("type") for s in sub.get("anyOf", []) if isinstance(s, dict)}
    seq_like = re.search(
        r"(seq|value|target|left|right|query|dna|primer|template|read)", name, re.IGNORECASE
    )
    if "array" in kinds:
        return []
    if "object" in kinds and "string" not in kinds:
        return {}
    if "string" in kinds or seq_like:
        return EXAMPLE_SEQUENCE if seq_like else ""
    if "integer" in kinds:
        return 3 if name in {"k", "kmer", "kmer_size"} else 1
    if "number" in kinds:
        return 0.5
    if "boolean" in kinds:
        return False
    return None


def quote(value: str) -> str:
    """Quote one command-line argument for display, abbreviating long sequences."""

    if len(value) > 60:
        value = value[:28] + "..." + value[-28:]
    if re.fullmatch(r"[A-Za-z0-9_./=,:+-]+", value):
        return value
    return "'" + value.replace("'", "'\\''") + "'"


def code_block(cli: str, python: str) -> str:
    """Render the equivalent command line and Python code."""

    return f"**Command line**\n```bash\n{cli}\n```\n**Python**\n```python\n{python}\n```"


def fail(message: Any) -> str:
    """Render an error message."""

    return f"⚠️ {message}"


class GuiBackend:
    """State and event handlers shared by every visitor of one interface."""

    def __init__(self, settings: GuiSettings) -> None:
        self.settings = settings
        self.registry = default_tool_registry(verbose=False)
        self.stats = self.registry.stats()
        self.category_choices = [("All categories", "")] + [
            (f"{c['label_en']} ({c['agent_compatible']})", str(c["category"]))
            for c in self.registry.categories()
        ]
        self.skill_total = int(skill_catalog_stats()["count"])
        self.worker = ToolWorker()
        settings.workdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ folders and limits
    def session(self, state: str | None) -> tuple[str, Path]:
        """Return ``(state, folder)``: a private folder per visitor, or the local work folder."""

        if not self.settings.is_public:
            return state or "local", self.settings.workdir
        if not state:
            state = uuid.uuid4().hex
            self.cleanup_sessions()
        folder = self.settings.workdir / state
        folder.mkdir(parents=True, exist_ok=True)
        return state, folder

    def cleanup_sessions(self) -> None:
        """Delete visitor folders older than the configured age (public demo only)."""

        max_age = self.settings.session_max_age_s
        if max_age is None or not self.settings.workdir.exists():
            return
        now = time.time()
        for folder in self.settings.workdir.iterdir():
            try:
                if folder.is_dir() and now - folder.stat().st_mtime > max_age:
                    shutil.rmtree(folder, ignore_errors=True)
            except OSError:
                pass

    def stage_upload(self, upload: str | None, folder: Path) -> Path | None:
        """Copy an uploaded file into ``folder`` under a safe name."""

        if not upload:
            return None
        source = Path(upload)
        target = folder / re.sub(r"[^A-Za-z0-9._-]", "_", source.name)
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
        return target

    def clean_sequence(self, text: str) -> str:
        """Return the sequence from plain text or one FASTA record."""

        lines = [line for line in (text or "").splitlines() if not line.startswith(">")]
        sequence = re.sub(r"\s+", "", "".join(lines))
        if not sequence:
            raise gr.Error("Enter a DNA sequence (plain text or one FASTA record).")
        limit = self.settings.max_sequence_length
        if len(sequence) > limit:
            raise gr.Error(
                f"This tab accepts up to {limit:,} nt; use the Dataset tab for longer sequences."
            )
        return sequence

    def check_arguments(self, value: Any, folder: Path, key: str = "") -> None:
        """Reject path and credential arguments the public demo does not allow."""

        if not self.settings.confine_paths:
            return
        if key and BLOCKED_ARGUMENT.search(key):
            raise gr.Error(f"Argument '{key}' is disabled in the public demo.")
        if isinstance(value, dict):
            for k, v in value.items():
                self.check_arguments(v, folder, str(k))
        elif isinstance(value, list):
            for v in value:
                self.check_arguments(v, folder, key)
        elif isinstance(value, str) and ("/" in value or "\\" in value or value.startswith("~")):
            path = Path(value).expanduser()
            resolved = (path if path.is_absolute() else folder / path).resolve()
            if not resolved.is_relative_to(folder.resolve()):
                raise gr.Error(
                    "Paths must point to files in your session folder (upload a file first)."
                )

    def run_cli(self, args: list[str], cwd: Path) -> tuple[bool, Any, str]:
        """Run ``dnakit <args>`` in ``cwd``; return (ok, JSON or text, displayed command)."""

        shown = "dnakit " + " ".join(quote(a) for a in args)
        command = [sys.executable, "-c", CLI_BOOTSTRAP, *args]
        env = {**os.environ, "COLUMNS": "120", "NO_COLOR": "1", "PYTHONUTF8": "1"}
        timeout = self.settings.cli_timeout_s
        try:
            done = subprocess.run(
                command,
                cwd=cwd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout,
                env=env,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return False, f"Stopped after {timeout:.0f} s (time limit of this interface).", shown
        text = done.stdout.strip() or done.stderr.strip()
        if done.returncode != 0:
            return False, text, shown
        try:
            return True, json.loads(done.stdout), shown
        except json.JSONDecodeError:
            return True, text, shown

    # ------------------------------------------------------------------ Sequence tab
    def seq_descriptors(self, text: str, state: str | None) -> tuple[Any, ...]:
        state, folder = self.session(state)
        sequence = self.clean_sequence(text)
        ok, result, shown = self.run_cli(["describe", sequence], folder)
        code = code_block(
            shown,
            "from dnakit.descriptors import all_descriptors\n"
            "result = all_descriptors(sequence)  # 240 fixed fields",
        )
        if not ok:
            return state, fail(result), None, None, code
        values = result.get("values", {})
        missing = result.get("unavailable_reasons", {})
        rows = [[k, fmt(v)] for k, v in values.items()]
        rows += [[k, f"n/a ({r})"] for k, r in missing.items()]
        path = folder / "descriptors.csv"
        with open(path, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerows([["field", "value"], *rows])
        summary = f"**{len(values)} descriptor fields** for a {len(sequence):,} nt sequence."
        return state, summary, rows, str(path), code

    def seq_orfs(self, text: str, min_length: float, state: str | None) -> tuple[Any, ...]:
        state, folder = self.session(state)
        sequence = self.clean_sequence(text)
        length = int(min_length or 0)
        ok, result, shown = self.run_cli(["orfs", sequence, "--min-length", str(length)], folder)
        code = code_block(
            shown,
            "from dnakit.patterns import scan_orfs\n"
            f"result = scan_orfs(sequence, min_length={length})",
        )
        if not ok:
            return state, fail(result), None, code
        rows = [
            [
                h["strand"],
                h["frame"],
                h["coordinate_location"]["start"],
                h["coordinate_location"]["end"],
                h["nucleotide_length"],
                h.get("translation", ""),
            ]
            for h in result.get("hits", [])
        ]
        return (
            state,
            f"**{len(rows)} ORFs** in six reading frames (0-based, half-open).",
            rows,
            code,
        )

    def seq_motif(self, text: str, pattern: str, mode: str, state: str | None) -> tuple[Any, ...]:
        state, folder = self.session(state)
        sequence = self.clean_sequence(text)
        motif = (pattern or "").strip()
        if not motif:
            raise gr.Error("Enter a motif, e.g. TATAWAW (IUPAC) or GAATTC (exact).")
        ok, result, shown = self.run_cli(["search", sequence, motif, "--mode", mode], folder)
        code = code_block(
            shown,
            "from dnakit.patterns import scan_motif\n"
            f"result = scan_motif(sequence, {motif!r}, mode={mode!r}, strand='both')",
        )
        if not ok:
            return state, fail(result), None, code
        rows = [
            [
                h["strand"],
                h["coordinate_location"]["start"],
                h["coordinate_location"]["end"],
                h["matched_sequence"],
            ]
            for h in result.get("hits", [])
        ]
        return state, f"**{len(rows)} matches** on both strands.", rows, code

    def seq_compare(
        self, text: str, other: str, method: str, k: float, state: str | None
    ) -> tuple[Any, ...]:
        state, folder = self.session(state)
        left, right = self.clean_sequence(text), self.clean_sequence(other)
        args = ["compare", left, right, "--method", method]
        k_arg = ""
        if method == "kmer_jaccard":
            args += ["--k", str(int(k or 3))]
            k_arg = f", k={int(k or 3)}"
        ok, result, shown = self.run_cli(args, folder)
        code = code_block(
            shown,
            "from dnakit.similarity import compare\n"
            f"result = compare(left, right, method={method!r}{k_arg})",
        )
        if not ok:
            return state, fail(result), None, code
        value = result.get("distance", result.get("similarity", result.get("value")))
        head = f"**{method}**: {fmt(value)}" if value is not None else f"**{method}** finished."
        return state, head, result, code

    def seq_validate(self, text: str, state: str | None) -> tuple[Any, ...]:
        state, folder = self.session(state)
        sequence = self.clean_sequence(text)
        ok, result, shown = self.run_cli(["validate", sequence], folder)
        code = code_block(
            shown, "from dnakit.standardize import validate\n\nresult = validate(sequence)"
        )
        if not ok:
            return state, fail(result), None, code
        verdict = "valid" if result.get("valid", result.get("ok", True)) else "has issues"
        return state, f"Sequence **{verdict}** ({len(sequence):,} symbols).", result, code

    # ------------------------------------------------------------------ Dataset tab
    def dataset_input(self, upload: str | None, use_example: bool, folder: Path) -> Path:
        """Return the uploaded file, or a copy of the example, inside ``folder``."""

        staged = self.stage_upload(upload, folder)
        if staged is not None:
            return staged
        if use_example:
            target = folder / "example.fasta"
            shutil.copyfile(example_fasta(), target)
            return target
        raise gr.Error("Upload a FASTA/FASTQ/CSV file or tick 'Use the example file'.")

    def ds_deduplicate(
        self, upload: str | None, use_example: bool, equivalence: str, state: str | None
    ) -> tuple[Any, ...]:
        state, folder = self.session(state)
        source = self.dataset_input(upload, use_example, folder)
        out = folder / f"{source.stem}.dedup.fasta"
        ok, result, shown = self.run_cli(
            [
                "deduplicate",
                source.name,
                out.name,
                "--equivalence",
                equivalence,
                "--output-format",
                "fasta",
                "--overwrite",
            ],
            folder,
        )
        code = code_block(
            shown,
            "from dnakit.datasets import deduplicate\nfrom dnakit.io import read_set\n\n"
            f"result = deduplicate(read_set({source.name!r}), equivalence={equivalence!r})",
        )
        if not ok:
            return state, fail(result), None, None, code
        groups = [
            [
                g["representative_id"],
                len(g["member_ids"]),
                ", ".join(g["member_ids"])[:300],
                ", ".join(g.get("orientations", []))[:120],
            ]
            for g in result.get("groups", [])
            if len(g["member_ids"]) > 1
        ]
        head = (
            f"Removed **{result.get('duplicate_count', 0)}** duplicates; "
            f"**{len(result.get('groups', []))}** unique sequences kept."
        )
        return state, head, groups, [str(out)], code

    def ds_split(
        self,
        upload: str | None,
        use_example: bool,
        method: str,
        ratios: str,
        seed: float,
        state: str | None,
    ) -> tuple[Any, ...]:
        state, folder = self.session(state)
        source = self.dataset_input(upload, use_example, folder)
        out = folder / f"{source.stem}_split"
        shutil.rmtree(out, ignore_errors=True)
        seed_value = int(seed or 0)
        ok, result, shown = self.run_cli(
            [
                "split",
                source.name,
                out.name,
                "--method",
                method,
                "--ratios",
                ratios,
                "--seed",
                str(seed_value),
                "--output-format",
                "fasta",
                "--overwrite",
            ],
            folder,
        )
        code = code_block(
            shown,
            "from dnakit.datasets import SplitConfig, split\nfrom dnakit.io import read_set\n\n"
            f"records = read_set({source.name!r})\n"
            f"result = split(records, config=SplitConfig(method={method!r}, ratios=..., "
            f"seed={seed_value}))",
        )
        if not ok:
            return state, fail(result), None, None, code
        counts = result.get("counts", {})
        files = sorted(str(p) for p in out.iterdir()) if out.exists() else []
        head = "Subset sizes: " + ", ".join(f"**{k}** {v}" for k, v in counts.items())
        rows = [[a["record_id"], a["split"]] for a in result.get("assignments", [])]
        return state, head, rows, files, code

    def ds_convert(
        self, upload: str | None, use_example: bool, target_format: str, state: str | None
    ) -> tuple[Any, ...]:
        state, folder = self.session(state)
        source = self.dataset_input(upload, use_example, folder)
        suffix = {"genbank": "gb"}.get(target_format, target_format)
        out = folder / f"{source.stem}.converted.{suffix}"
        ok, result, shown = self.run_cli(
            [
                "convert",
                source.name,
                out.name,
                "--output-format",
                target_format,
                "--overwrite",
                "--no-progress",
            ],
            folder,
        )
        code = code_block(
            shown,
            "from dnakit.io import read_set, write\n\n"
            f"write(read_set({source.name!r}), {out.name!r}, format={target_format!r})",
        )
        if not ok:
            return state, fail(result), None, code
        count = result.get("record_count", "?")
        return state, f"Wrote **{count}** records as {target_format}.", [str(out)], code

    def ds_descriptor_table(
        self,
        upload: str | None,
        use_example: bool,
        state: str | None,
        progress: gr.Progress = gr.Progress(),  # noqa: B008 - Gradio injects the tracker
    ) -> tuple[Any, ...]:
        state, folder = self.session(state)
        source = self.dataset_input(upload, use_example, folder)
        table = folder / f"{source.stem}.records.csv"
        ok, result, _ = self.run_cli(
            [
                "convert",
                source.name,
                table.name,
                "--output-format",
                "csv",
                "--overwrite",
                "--no-progress",
            ],
            folder,
        )
        if not ok:
            return state, fail(result), None, None, ""
        with open(table, newline="", encoding="utf-8") as fh:
            records = [(r["id"], r["sequence"]) for r in csv.DictReader(fh)]
        max_records = self.settings.max_table_records
        if max_records is not None and len(records) > max_records:
            raise gr.Error(f"This interface computes descriptors for up to {max_records} records.")
        max_length = self.settings.max_record_length
        too_long = [rid for rid, seq in records if max_length is not None and len(seq) > max_length]
        if too_long:
            raise gr.Error(
                f"{len(too_long)} records are longer than {max_length:,} nt "
                f"(first: {too_long[0]}); install DNAKit locally for long sequences."
            )
        rows: list[list[Any]] = []
        fields: list[str] = []
        for rid, seq in progress.tqdm(records, desc="descriptors"):
            ok, res = self.worker.call(
                "dnakit.descriptors.all_descriptors",
                {"value": seq},
                folder,
                timeout=self.settings.call_timeout_s,
            )
            if not ok:
                return state, fail(f"{rid}: {res}"), None, None, ""
            values = res.get("values", {})
            fields = fields or list(values)
            rows.append([rid, *[values.get(f) for f in fields]])
        out = folder / f"{source.stem}.descriptors.csv"
        with open(out, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerows([["id", *fields], *rows])
        shown = [f for f in PREVIEW_FIELDS if f in fields] or fields[:6]
        preview = [[r[0], *[fmt(r[1 + fields.index(f)]) for f in shown]] for r in rows[:50]]
        head = (
            f"Feature matrix: **{len(rows)} records x {len(fields)} descriptors** "
            f"(preview shows {len(shown)} columns; download the CSV for all of them)."
        )
        code = code_block(
            "dnakit describe <sequence>   # one record",
            "from dnakit.descriptors import all_descriptors\n"
            "table = [all_descriptors(r.sequence).values for r in records]",
        )
        return state, head, gr.update(value=preview, headers=["id", *shown]), [str(out)], code

    # ------------------------------------------------------------------ Function explorer
    def runnable(self, effect: str) -> bool:
        """Whether a tool with this effect may run in this interface."""

        return effect in SAFE_EFFECTS or self.settings.allow_side_effects

    def fx_search(self, query: str, category: str) -> tuple[Any, ...]:
        if not query and not category:
            matches: tuple[Any, ...] = ()
        else:
            matches = self.registry.search(
                query or "", category=category or None, limit=40, mode="bm25"
            )
        rows = [
            [t.name, t.category, EFFECT_NOTE.get(t.effect, t.effect), t.description[:160]]
            for t in matches
        ]
        names = [t.name for t in matches]
        hint = "" if rows else " (try e.g. 'cpg island', 'primer', 'kmer')"
        return (
            f"{len(rows)} functions{hint}",
            rows,
            gr.update(choices=names, value=names[0] if names else None),
        )

    def fx_describe(self, name: str) -> tuple[Any, ...]:
        if not name:
            return "", "{}", {}, gr.update(interactive=False), gr.update(visible=False, value=False)
        tool: dict[str, Any] = self.registry.describe(name)
        schema: dict[str, Any] = tool.get("input_schema", {})
        required = schema.get("required", [])
        properties = schema.get("properties", {})
        example = {p: example_value(p, properties.get(p, {})) for p in required}
        effect = str(tool["effect"])
        needs_consent = effect not in SAFE_EFFECTS
        if not self.runnable(effect):
            badge = f"⛔ {EFFECT_NOTE.get(effect)}: install DNAKit locally to run it"
        elif needs_consent:
            badge = f"⚠️ {EFFECT_NOTE.get(effect)}: tick the confirmation box below to run it"
        else:
            badge = "✅ runs here"
        doc = tool.get("full_description") or tool["description"]
        md = (
            f"### `{tool['name']}`\n{badge} · category **{tool['category']}**\n\n"
            f"```python\n{tool['function']}{str(tool['signature'])[:600]}\n```\n\n{doc[:2500]}"
        )
        consent = gr.update(visible=needs_consent and self.runnable(effect), value=False)
        return (
            md,
            json.dumps(example, indent=2),
            schema,
            gr.update(interactive=self.runnable(effect)),
            consent,
        )

    def fx_upload(self, upload: str | None, state: str | None) -> tuple[Any, ...]:
        state, folder = self.session(state)
        staged = self.stage_upload(upload, folder)
        if staged is None:
            return state, ""
        return state, f"Saved as `{staged.name}`; use this file name in path arguments."

    def fx_run(
        self, name: str, arguments: str, consent: bool, state: str | None
    ) -> tuple[Any, ...]:
        state, folder = self.session(state)
        if not name:
            raise gr.Error("Pick a function first.")
        tool = self.registry.resolve(name)
        if not self.runnable(tool.effect):
            raise gr.Error("This function is not executed in the public demo.")
        side_effects = tool.effect not in SAFE_EFFECTS
        if side_effects and not consent:
            raise gr.Error(f"This function {EFFECT_NOTE.get(tool.effect)}; tick the box to allow.")
        try:
            args = json.loads(arguments or "{}")
        except json.JSONDecodeError as exc:
            raise gr.Error(f"Arguments are not valid JSON: {exc}") from exc
        if not isinstance(args, dict):
            raise gr.Error("Arguments must be a JSON object.")
        self.check_arguments(args, folder)
        before = {p: p.stat().st_mtime for p in folder.rglob("*") if p.is_file()}
        started = time.time()
        ok, result = self.worker.call(
            tool.name,
            args,
            folder,
            timeout=self.settings.call_timeout_s,
            allow_side_effects=side_effects,
        )
        elapsed = time.time() - started
        module, function = tool.name.rsplit(".", 1)
        call = pprint.pformat(shorten(args), width=80, sort_dicts=False)
        py = f"from {module} import {function}\n\nresult = {function}(**{call})"
        mcp_arguments: dict[str, Any] = {"tool_name": tool.name, "arguments": shorten(args)}
        if side_effects:
            mcp_arguments["allow_side_effects"] = True
        mcp = json.dumps({"tool": "call_dnakit_tool", "arguments": mcp_arguments}, indent=2)
        code = (
            f"**Python**\n```python\n{py}\n```\n"
            f"**MCP (what an agent sends)**\n```json\n{mcp[:3000]}\n```"
        )
        if not ok:
            return state, fail(result), None, "", gr.update(value=None, visible=False), code
        new_files = [
            str(p) for p in folder.rglob("*") if p.is_file() and before.get(p) != p.stat().st_mtime
        ]
        svg = find_svg(result)
        html = f'<div class="dk-plot">{svg}</div>' if svg else ""
        return (
            state,
            f"✅ Finished in {elapsed:.1f} s.",
            hide_svg(result),
            html,
            gr.update(value=new_files or None, visible=bool(new_files)),
            code,
        )

    def fx_stop(self) -> str:
        self.worker.stop()
        return fail("Stopped by the user.")

    # ------------------------------------------------------------------ Skills tab
    def sk_search(self, query: str) -> tuple[Any, ...]:
        if not (query or "").strip():
            return "Enter a search term.", [], gr.update(choices=[], value=None)
        found = search_skills(query, limit=40)
        items = found.get("skills", [])
        rows = [
            [s["id"], s.get("category", ""), s.get("source_id", ""), s.get("description", "")[:200]]
            for s in items
        ]
        names = [s["id"] for s in items]
        return (
            f"{found.get('total', len(rows))} matching Skills (showing {len(rows)})",
            rows,
            gr.update(choices=names, value=names[0] if names else None),
        )

    def sk_read(self, skill_id: str) -> str:
        if not skill_id:
            return ""
        try:
            doc = read_skill_file(skill_id, limit=24000)
        except Exception as exc:  # unknown or unreadable Skill
            return fail(exc)
        text = re.sub(r"\A---\n.*?\n---\n", "", str(doc.get("content", "")), flags=re.S)
        meta = doc.get("skill", {})
        more = (
            "\n\n*(truncated; install the Skill to read the rest)*"
            if doc.get("next_offset")
            else ""
        )
        source = doc.get("attribution", meta.get("source_id", ""))
        head = (
            f"**{skill_id}** · {meta.get('category', '')} · license {meta.get('license', '?')} · "
            f"source [{source}]({doc.get('upstream_url', '')})\n\n"
            f"> Install for an agent: `dnakit skills install {skill_id} --agent claude-code` "
            f"(or `--agent codex`, `--agent cursor`)\n\n---\n"
        )
        return head + text + more


def pick_first_column(evt: gr.SelectData, table: Any) -> Any:
    """Return a dropdown update selecting the first cell of the clicked table row."""

    try:
        row = evt.index[0]
        value = table.iloc[row, 0] if hasattr(table, "iloc") else table[row][0]
        return gr.update(value=value)
    except (IndexError, KeyError, TypeError):
        return gr.update()
