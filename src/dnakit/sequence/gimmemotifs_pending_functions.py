"""GimmeMotifs wrappers for gimme motifs, gimme prediction and the maelstrom report.

Every public function takes explicit local files and plain arguments. Nothing is looked up in a
genomepy installation, downloaded or sent to a service. The gimme motifs wrapper scans with known
motifs (a PFM file) and, when de novo tools are requested, runs those tools; each tool must be
configured for gimmemotifs and installed, otherwise BackendUnavailableError is raised with an
install hint. Upstream gimmemotifs functions are called with their own parameters except where a
docstring says otherwise.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require, workspace
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_EXTRA = "gimmemotifs"
_LOG_TAIL = 2000
_ANALYSES = ("small", "medium", "large", "xl")
_TOOL_HINT = (
    "Install the {name} program, then add it to the gimmemotifs configuration file "
    "(~/.config/gimmemotifs/gimmemotifs.cfg) with its executable path; the gimmemotifs "
    "documentation lists the section for each tool."
)


def _require(module: str) -> Any:
    return require(module, _EXTRA)


def _version() -> str | None:
    try:
        return importlib.metadata.version("gimmemotifs")
    except importlib.metadata.PackageNotFoundError:
        return None


def _file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} does not exist: {path}")
    return path


def _folder(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a folder path.")
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} does not exist: {path}")
    return path


def _fresh(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new path: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _number(value: object, name: str, *, low: float, high: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number.")
    if value < low or (high is not None and value > high):
        bound = f"at least {low}" if high is None else f"between {low} and {high}"
        raise ConfigurationError(f"{name} must be {bound}.")
    return float(value)


def _names(value: object, name: str) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of names.")
    items: list[str] = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise ConfigurationError(f"{name} must hold non-empty strings.")
        items.append(item.strip())
    return items


def _tool(name: object) -> Any:
    """Return an installed, configured gimmemotifs tool program or raise BackendUnavailableError."""
    tools = _require("gimmemotifs.tools")
    if not isinstance(name, str) or name.lower() not in tools.__tools__:
        known = ", ".join(sorted(tools.__tools__))
        raise ConfigurationError(f"Unknown gimmemotifs tool {name!r}; choose from: {known}.")
    hint = _TOOL_HINT.format(name=name)
    try:
        program = tools.get_tool(name)
    except ValueError as exc:
        raise BackendUnavailableError(
            f"gimmemotifs has no configuration for the {name} program: {exc}", hint=hint
        ) from exc
    if not program.is_installed():
        raise BackendUnavailableError(
            f"The {name} program is not installed or not executable.", hint=hint
        )
    return program


def _tail(text: object) -> str:
    return str(text or "")[-_LOG_TAIL:]


def _result(operation: str, outputs: Mapping[str, object], **extra: Any) -> dict[str, Any]:
    return {
        "status": "complete",
        "operation": operation,
        "backend": "gimmemotifs",
        "version": _version(),
        "outputs": {key: str(value) for key, value in outputs.items()},
        **extra,
    }


def _written(path: Path, what: str) -> Path:
    if not path.exists():
        raise BackendExecutionError(f"gimmemotifs did not write {what}: {path}")
    return path


def gimmemotifs_gimme_motifs(
    sample_fasta: str,
    pfm_file: str,
    output_dir: str,
    *,
    background_fasta: str | None = None,
    known: bool = True,
    denovo_tools: Sequence[str] | None = None,
    analysis: str = "xl",
    fraction: float = 0.2,
    single: bool = False,
    report: bool = True,
    ncpus: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run ``gimme motifs`` on an explicit sample FASTA file.

    Known motifs from ``pfm_file`` are scanned on the sample and on a background, and the ROC,
    precision-recall and enrichment statistics are written to ``gimme.roc.report.txt``. A
    non-redundant subset is selected, its scan results are written as BED files, and the HTML
    reports ``gimme.motifs.html`` and ``gimme.motifs.redundant.html`` are written when ``report`` is
    true. ``background_fasta`` is a FASTA file of background sequences; when it is None, gimmemotifs
    generates a random background from the sample, as it does without a genome.

    ``denovo_tools`` names de novo tools (for example ``["MEME", "Homer", "BioProspector"]``, the
    upstream default set); they run only when each is configured and installed, otherwise
    BackendUnavailableError is raised before any work starts. The default is known motifs only.

    Differences from the upstream command line: the genome is not an input, so z-score and GC
    normalisation of the scores are off, as they are upstream when no genome is given. Upstream
    normalisation needs a genomepy-installed genome (a FASTA path gives an empty background), and
    its motif cache fails under xxhash 4.x (``Motif.hash`` raises TypeError), so it is not offered.
    ``output_dir`` must not exist yet.
    """
    if analysis not in _ANALYSES:
        raise ConfigurationError(f"analysis must be one of: {', '.join(_ANALYSES)}.")
    if not known and not denovo_tools:
        raise ConfigurationError("Choose known=True, give denovo_tools, or both.")
    sample = _file(sample_fasta, "sample_fasta")
    pfm = _file(pfm_file, "pfm_file")
    background = None if background_fasta is None else _file(background_fasta, "background_fasta")
    frac = _number(fraction, "fraction", low=0.0, high=1.0)
    if frac == 0.0:
        raise ConfigurationError("fraction must be above 0.")
    denovo = None
    if denovo_tools:
        denovo = _names(denovo_tools, "denovo_tools")
        for name in denovo:
            _tool(name)
    folder = workspace(output_dir)

    motif_command = _require("gimmemotifs.commands.motifs")
    args = SimpleNamespace(
        outdir=str(folder),
        sample=str(sample),
        size=0,
        genome=None,
        background=None if background is None else str(background),
        pfmfile=str(pfm),
        known=bool(known),
        denovo=denovo is not None,
        tools=None if denovo is None else ",".join(denovo),
        analysis=analysis,
        fraction=frac,
        single=bool(single),
        zscore=False,
        gc=False,
        ncpus=int(_number(ncpus, "ncpus", low=1)),
        noprogress=True,
        report=bool(report),
    )
    try:
        motif_command.motifs(args)
    except SystemExit as exc:
        raise BackendExecutionError(
            "gimme motifs found no motifs to report; add known motifs or de novo tools."
        ) from exc

    roc_report = _written(folder / "gimme.roc.report.txt", "the ROC statistics table")
    outputs: dict[str, object] = {
        "roc_report": roc_report,
        "motif_scan_results": folder / "motif_scan_results",
    }
    for name in ("gimme.motifs.html", "gimme.motifs.redundant.html", "combined.pfm"):
        if (folder / name).exists():
            outputs[name] = folder / name
    with roc_report.open(encoding="utf-8") as handle:
        scored = max(sum(1 for _ in handle) - 1, 0)
    return _result(
        "gimme_motifs",
        outputs,
        known=bool(known),
        denovo_tools=denovo or [],
        motifs_scored=scored,
    )


def gimmemotifs_predict_motifs(
    tool: str,
    fasta: str,
    output_pfm: str,
    *,
    params: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run one de novo motif tool on a FASTA file (``gimme prediction``).

    ``tool`` is a gimmemotifs tool name such as ``MEME`` or ``Homer``; it must be configured and
    installed, otherwise BackendUnavailableError is raised. ``params`` are the tool parameters
    (the mapping that ``gimme prediction`` reads from its ``--paramfile``). The predicted motifs
    are written to ``output_pfm`` in PFM format, one motif after another. The tool's standard
    output and error are returned as bounded tails.
    """
    program = _tool(tool)
    infile = _file(fasta, "fasta")
    target = _fresh(output_pfm, "output_pfm")
    if params is None:
        settings: dict[str, Any] = {}
    elif isinstance(params, Mapping):
        settings = dict(params)
    else:
        raise ConfigurationError("params must be a mapping of tool parameters.")
    motifs, stdout, stderr = program.run(str(infile), settings)
    with target.open("w", encoding="utf-8") as handle:
        for motif in motifs:
            handle.write(f"{motif.to_pfm()}\n")
    return _result(
        "predict_motifs",
        {"pfm": target},
        tool=str(tool),
        motifs=len(motifs),
        stdout_tail=_tail(stdout),
        stderr_tail=_tail(stderr),
    )


def gimmemotifs_maelstrom_html_report(
    run_folder: str,
    pfm_file: str,
    *,
    threshold: float = 3,
    plot_all_motifs: bool = False,
    plot_no_motifs: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the HTML report of a finished maelstrom run folder.

    ``run_folder`` is a maelstrom output folder that holds ``final.out.txt`` (the motif-by-method
    table). ``pfm_file`` is the motif database used by the run, for the logos and factor names.
    Motifs whose absolute value is below ``threshold`` in every method column are left out of the
    report. The report is written to ``run_folder/gimme.maelstrom.report.html``, and an existing
    report in that folder is replaced, as upstream does. ``plot_all_motifs`` draws a logo for every
    motif; ``plot_no_motifs`` leaves the logo column out.
    """
    folder = _folder(run_folder, "run_folder")
    table = _file(folder / "final.out.txt", "final.out.txt in run_folder")
    pfm = _file(pfm_file, "pfm_file")
    cutoff = _number(threshold, "threshold", low=0.0)
    pd = _require("pandas")
    df = pd.read_table(table, index_col=0)
    value_cols = [
        str(col) for col in df.columns if "corr" not in str(col) and "% with motif" not in str(col)
    ]
    shown = int((df[value_cols].abs() >= cutoff).any(axis=1).sum())
    report = _require("gimmemotifs.report")
    report.maelstrom_html_report(
        str(folder),
        str(table),
        str(pfm),
        threshold=cutoff,
        plot_all_motifs=bool(plot_all_motifs),
        plot_no_motifs=bool(plot_no_motifs),
    )
    html = _written(folder / "gimme.maelstrom.report.html", "the maelstrom HTML report")
    outputs: dict[str, object] = {"report": html}
    if (folder / "logos").is_dir():
        outputs["logos"] = folder / "logos"
    return _result(
        "maelstrom_html_report",
        outputs,
        threshold=cutoff,
        motifs_in_report=shown,
        methods=len(value_cols),
    )


__all__ = [
    "gimmemotifs_gimme_motifs",
    "gimmemotifs_maelstrom_html_report",
    "gimmemotifs_predict_motifs",
]

_register_public_verbose(globals())
