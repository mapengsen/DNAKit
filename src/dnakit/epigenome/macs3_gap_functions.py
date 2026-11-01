"""MACS3 peak-calling and signal-track operations as typed DNAKit facades.

The user installs MACS3 (``python -m pip install macs3``, or ``conda install -c bioconda macs3``). Each function runs one
``macs3`` subcommand with an explicit typed signature mapped to fixed flags; no shell is used and the program is resolved
from an explicit path, then PATH. Outputs are new explicit paths or fresh directories, so an existing result is never
overwritten. Results are JSON-serialisable dictionaries with status, backend, version, command, return code, existing
outputs and a bounded log tail. Defaults are the defaults of the MACS3 command-line options they map to.

Facades do not call the ``MACS3.Commands.*_cmd.run`` functions in-process; those take optparse option objects. Options that
MACS3 marks obsolete (``--to-large``, ``--ratio``, ``--shiftsize``) are not exposed. Not every option of a subcommand is
exposed; the ones that change the algorithm's inputs and outputs are.
"""

from __future__ import annotations

import re
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import ConfigurationError, DNAKitError

_HINT = "Install MACS3 3.0.5 or later: python -m pip install macs3 (or conda install -c bioconda macs3)."
_ALIGNMENT_FORMATS = ("AUTO", "BAM", "SAM", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "BOWTIE")
_PAIRED_FORMATS = ("BAMPE", "BEDPE")
_TAG_FORMATS = (*_ALIGNMENT_FORMATS, *_PAIRED_FORMATS)
_CALLPEAK_FORMATS = (*_TAG_FORMATS, "FRAG")
_GENOME_SHORTCUTS = ("hs", "mm", "ce", "dm")
_NAME = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]*")
_SINGLE_END_D = re.compile(r"predicted fragment length is (\d+) bps")
_PAIRED_END_D = re.compile(r"Average insertion length of all pairs is (\d+) bps")
_R_VECTOR = re.compile(r"^(p|m|ycorr|xcorr|altd)\s*<-\s*c\((.*)\)\s*$", re.MULTILINE)


# ---------------------------------------------------------------------------
# Shared option mapping
# ---------------------------------------------------------------------------


def _name(value: str, name: str) -> str:
    """A file-name prefix or file name: letters, digits, underscore, dot and dash; no path separators."""
    text = _cli.text(value, name)
    if _NAME.fullmatch(text) is None:
        raise ConfigurationError(
            f"{name} may contain only letters, digits, '_', '.' and '-': {text!r}"
        )
    return text


def _inputs(paths: Sequence[str], name: str) -> list[str]:
    """One or more existing files, as absolute paths (so they can never be read as flags)."""
    if isinstance(paths, str) or not paths:
        raise ConfigurationError(f"{name} must be a non-empty list of existing files.")
    return [str(_cli.existing_file(path, name)) for path in paths]


def _optional_file(path: str | None, name: str) -> str | None:
    return None if path is None else str(_cli.existing_file(path, name))


def _genome_size(value: str | float) -> str:
    """A MACS3 shortcut (hs, mm, ce, dm) or an effective genome size in bp."""
    if isinstance(value, str) and value in _GENOME_SHORTCUTS:
        return value
    try:
        size = 0.0 if isinstance(value, bool) else float(value)
    except (TypeError, ValueError):
        size = 0.0
    if not size > 0:
        raise ConfigurationError(
            f"genome_size must be one of {', '.join(_GENOME_SHORTCUTS)} or a positive number."
        )
    return format(size, ".12g")


def _optional_int(value: int | None, name: str, *, minimum: int = 1) -> str | None:
    return None if value is None else str(_cli.integer(value, name, minimum=minimum))


def _keep_duplicates(value: str | int, name: str = "keep_duplicates") -> str:
    """auto, all or a positive integer: how many tags at one position and strand are kept."""
    if isinstance(value, int) and not isinstance(value, bool):
        return str(_cli.integer(value, name, minimum=1))
    if value in ("auto", "all"):
        return str(value)
    if isinstance(value, str) and value.isdigit() and int(value) >= 1:
        return value
    raise ConfigurationError(f"{name} must be 'auto', 'all' or a positive integer.")


def _mfold(value: Sequence[int]) -> tuple[str, str]:
    """The lower and upper fold-enrichment bounds for building the shifting model (lower < upper)."""
    if len(value) != 2:
        raise ConfigurationError("mfold must be two integers: (lower, upper).")
    low = _cli.integer(value[0], "mfold lower", minimum=1)
    high = _cli.integer(value[1], "mfold upper", minimum=1)
    if low >= high:
        raise ConfigurationError("mfold lower bound must be smaller than the upper bound.")
    return str(low), str(high)


def _switch(flag: str, on: bool) -> list[str]:
    return [flag] if on else []


def _valued(flag: str, value: object | None) -> list[str]:
    return [] if value is None else [flag, str(value)]


def _new_directory_path(path: str | Path, name: str) -> Path:
    """Where a fresh output directory will be created. Refuses an existing path; creates nothing (see _make_directory)."""
    out = Path(path).expanduser().resolve()
    if out.exists():
        raise ConfigurationError(f"{name} must not exist yet: {out}")
    return out


def _make_directory(out: Path) -> Path:
    out.mkdir(parents=True)
    return out


def _output_file(path: str | Path, name: str) -> Path:
    return _cli.new_file(path, name)


def _executable(macs3_path: str | None) -> str:
    """The macs3 program (explicit path, then PATH); raised before any output directory is created."""
    return _cli.program("macs3", explicit=macs3_path, hint=_HINT)


def _macs3(
    subcommand: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    executable: str,
) -> dict[str, Any]:
    """Run ``macs3 <subcommand> <arguments>``; completed only when the exit code is 0 and every expected file exists."""
    return _cli.run(
        "macs3",
        [subcommand, *arguments],
        cwd=cwd,
        expected=list(expected),
        explicit=executable,
        hint=_HINT,
    )


def _existing(path: Path) -> str | None:
    return str(path) if path.is_file() else None


# ---------------------------------------------------------------------------
# callpeak: call peaks from alignment files
# ---------------------------------------------------------------------------


def call_chip_peaks(
    treatment: Sequence[str],
    output_dir: str,
    *,
    control: Sequence[str] | None = None,
    name: str = "NA",
    file_format: str = "AUTO",
    genome_size: str | float = "hs",
    tag_size: int | None = None,
    keep_duplicates: str | int = "1",
    max_count: int | None = None,
    barcodes: str | None = None,
    save_bedgraph: bool = False,
    spmr: bool = False,
    trackline: bool = False,
    nomodel: bool = False,
    shift: int = 0,
    extsize: int = 200,
    bandwidth: int = 300,
    min_fragment: int = 20,
    mfold: Sequence[int] = (5, 50),
    fix_bimodal: bool = False,
    qvalue: float | None = 0.05,
    pvalue: float | None = None,
    scale_to: str | None = None,
    down_sample: bool = False,
    seed: int | None = None,
    nolambda: bool = False,
    small_local: int = 1000,
    large_local: int = 10000,
    max_gap: int | None = None,
    min_length: int | None = None,
    broad: bool = False,
    broad_cutoff: float = 0.1,
    cutoff_analysis: bool = False,
    call_summits: bool = False,
    fe_cutoff: float = 1.0,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call peaks from ChIP-seq or ATAC-seq alignments (macs3 callpeak).

    ``treatment`` (and optional ``control``) are alignment files in ``file_format``. Outputs go to a new ``output_dir``
    as ``<name>_peaks.narrowPeak`` (or ``.broadPeak`` and ``.gappedPeak`` with ``broad``) and ``<name>_peaks.xls``; more
    files follow the options (summits, bedGraph tracks, cutoff analysis). ``qvalue`` and ``pvalue`` are mutually
    exclusive: pass ``qvalue=None`` to use ``pvalue``. Unless ``nomodel`` is set, MACS3 builds a fragment-size model
    and writes ``<name>_model.r``; if the model cannot be built, the subcommand fails and the status is failed.
    Result keys: ``peak_files``, ``xls_file``, ``model_script`` (path or None).
    """
    treat = _inputs(treatment, "treatment")
    ctrl = None if not control else _inputs(control, "control")
    out = _new_directory_path(output_dir, "output_dir")
    prefix = _name(name, "name")
    fmt = _cli.choice(file_format, _CALLPEAK_FORMATS, "file_format")
    if qvalue is not None and pvalue is not None:
        raise ConfigurationError(
            "qvalue and pvalue are mutually exclusive: pass qvalue=None to use pvalue."
        )
    if pvalue is None and qvalue is None:
        raise ConfigurationError("Either qvalue or pvalue must be given.")
    low, high = _mfold(mfold)
    keep = _keep_duplicates(keep_duplicates)
    max_gap_value = _optional_int(max_gap, "max_gap")
    min_length_value = _optional_int(min_length, "min_length")
    if scale_to is not None:
        _cli.choice(scale_to, ("large", "small"), "scale_to")
    seed_value = None if seed is None else str(_cli.integer(seed, "seed", minimum=0))
    args: list[str] = ["-t", *treat]
    if ctrl:
        args += ["-c", *ctrl]
    args += [
        "-f",
        fmt,
        "-g",
        _genome_size(genome_size),
        "--keep-dup",
        keep,
        "--outdir",
        str(out),
        "-n",
        prefix,
        "--shift",
        str(_cli.integer(shift, "shift")),
        "--extsize",
        str(_cli.integer(extsize, "extsize", minimum=1)),
        "--bw",
        str(_cli.integer(bandwidth, "bandwidth", minimum=1)),
        "--d-min",
        str(_cli.integer(min_fragment, "min_fragment", minimum=1)),
        "-m",
        low,
        high,
        "--slocal",
        str(_cli.integer(small_local, "small_local", minimum=0)),
        "--llocal",
        str(_cli.integer(large_local, "large_local", minimum=0)),
        "--fe-cutoff",
        str(_cli.number(fe_cutoff, "fe_cutoff", minimum=0.0)),
    ]
    args += _valued("-s", _optional_int(tag_size, "tag_size"))
    args += _valued(
        "--max-count",
        None if max_count is None else _cli.integer(max_count, "max_count", minimum=0),
    )
    args += _valued("--barcodes", _optional_file(barcodes, "barcodes"))
    args += _valued("--scale-to", scale_to)
    args += _valued("--seed", seed_value)
    args += _valued("--max-gap", max_gap_value)
    args += _valued("--min-length", min_length_value)
    if pvalue is not None:
        args += ["-p", str(_cli.number(pvalue, "pvalue", minimum=0.0, maximum=1.0))]
    else:
        args += ["-q", str(_cli.number(qvalue, "qvalue", minimum=0.0, maximum=1.0))]
    if broad:
        args += [
            "--broad",
            "--broad-cutoff",
            str(_cli.number(broad_cutoff, "broad_cutoff", minimum=0.0, maximum=1.0)),
        ]
    args += (
        _switch("-B", save_bedgraph) + _switch("--SPMR", spmr) + _switch("--trackline", trackline)
    )
    args += _switch("--nomodel", nomodel) + _switch("--fix-bimodal", fix_bimodal)
    args += _switch("--down-sample", down_sample) + _switch("--nolambda", nolambda)
    args += _switch("--cutoff-analysis", cutoff_analysis) + _switch("--call-summits", call_summits)

    stem = out / f"{prefix}"
    peak_files = (
        [Path(f"{stem}_peaks.broadPeak"), Path(f"{stem}_peaks.gappedPeak")]
        if broad
        else [Path(f"{stem}_peaks.narrowPeak")]
    )
    expected = [*peak_files, Path(f"{stem}_peaks.xls")]
    if call_summits and not broad:
        expected.append(Path(f"{stem}_summits.bed"))
    if save_bedgraph:
        expected += [Path(f"{stem}_treat_pileup.bdg"), Path(f"{stem}_control_lambda.bdg")]
    if cutoff_analysis:
        expected.append(Path(f"{stem}_cutoff_analysis.txt"))
    executable = _executable(macs3_path)
    _make_directory(out)
    result = _macs3("callpeak", args, cwd=out, expected=expected, executable=executable)
    result["output_dir"] = str(out)
    result["peak_files"] = [str(p) for p in peak_files if p.is_file()]
    result["xls_file"] = _existing(Path(f"{stem}_peaks.xls"))
    result["model_script"] = None if nomodel else _existing(Path(f"{stem}_model.r"))
    return result


# ---------------------------------------------------------------------------
# predictd and the fragment-model plot
# ---------------------------------------------------------------------------


def _predict(
    inputs: list[str],
    output_dir: str,
    *,
    file_format: str,
    genome_size: str | float,
    tag_size: int | None,
    bandwidth: int,
    min_fragment: int,
    mfold: Sequence[int],
    script_name: str,
    macs3_path: str | None,
) -> dict[str, Any]:
    fmt = _cli.choice(file_format, _TAG_FORMATS, "file_format")
    low, high = _mfold(mfold)
    gsize = _genome_size(genome_size)
    bw = _cli.integer(bandwidth, "bandwidth", minimum=1)
    dmin = _cli.integer(min_fragment, "min_fragment", minimum=1)
    tsize = _optional_int(tag_size, "tag_size")
    out = _new_directory_path(output_dir, "output_dir")
    script = out / script_name
    args = [
        "-i",
        *inputs,
        "-f",
        fmt,
        "-g",
        gsize,
        "--bw",
        str(bw),
        "--d-min",
        str(dmin),
        "-m",
        low,
        high,
        "--outdir",
        str(out),
        "--rfile",
        script_name,
    ]
    args += _valued("-s", tsize)
    # Paired-end input only reports the average insertion length; no R script is written in that mode.
    expected = [] if fmt in _PAIRED_FORMATS else [script]
    executable = _executable(macs3_path)
    _make_directory(out)
    result = _macs3("predictd", args, cwd=out, expected=expected, executable=executable)
    text = result.get("log_tail", "")
    single = _SINGLE_END_D.search(text)
    paired = _PAIRED_END_D.search(text)
    if single is not None:
        result["fragment_length"], result["fragment_length_source"] = (
            int(single.group(1)),
            "single_end_model",
        )
    elif paired is not None:
        result["fragment_length"], result["fragment_length_source"] = (
            int(paired.group(1)),
            "paired_end_average",
        )
    else:
        result["fragment_length"], result["fragment_length_source"] = None, None
    result["output_dir"] = str(out)
    result["model_script"] = _existing(script)
    return result


def predict_fragment_size(
    ifile: Sequence[str],
    output_dir: str,
    *,
    file_format: str = "AUTO",
    genome_size: str | float = "hs",
    tag_size: int | None = None,
    bandwidth: int = 300,
    min_fragment: int = 20,
    mfold: Sequence[int] = (5, 50),
    script_name: str = "predictd_model.R",
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Predict the fragment (d) size from alignments (macs3 predictd) and write the model R script.

    Single-end input builds the peak model: the predicted size is reported as ``fragment_length`` and the R script
    ``script_name`` (containing the model's forward and reverse tag profiles and cross-correlation) is written to the new
    ``output_dir``. Paired-end input (BAMPE, BEDPE) reports the average insertion length instead and writes no script.
    A single-end file without enough strand-symmetric peaks fails with the upstream warning in ``log_tail``.
    """
    return _predict(
        _inputs(ifile, "ifile"),
        output_dir,
        file_format=file_format,
        genome_size=genome_size,
        tag_size=tag_size,
        bandwidth=bandwidth,
        min_fragment=min_fragment,
        mfold=mfold,
        script_name=_name(script_name, "script_name"),
        macs3_path=macs3_path,
    )


def _parse_model_script(path: Path) -> dict[str, Any]:
    """The vectors of a model R script written by macs3: percentages per distance, cross-correlation and alternatives."""
    vectors = {key: body for key, body in _R_VECTOR.findall(path.read_text(encoding="utf-8"))}
    if not {"p", "m", "xcorr", "ycorr", "altd"} <= vectors.keys():
        raise DNAKitError(f"No model vectors found in the MACS3 R script: {path}")

    def numbers(key: str) -> list[float]:
        return [float(item) for item in vectors[key].split(",") if item.strip()]

    forward, reverse = numbers("p"), numbers("m")
    width = len(forward)
    return {
        "distance_from_middle": [i - (width - 1) / 2 for i in range(width)],
        "forward_percent": forward,
        "reverse_percent": reverse,
        "xcorr_lag": numbers("xcorr"),
        "xcorr_correlation": numbers("ycorr"),
        "alternative_fragment_lengths": [round(value) for value in numbers("altd")],
    }


def plot_fragment_model(
    ifile: Sequence[str],
    output_dir: str,
    *,
    file_format: str = "AUTO",
    genome_size: str | float = "hs",
    tag_size: int | None = None,
    bandwidth: int = 300,
    min_fragment: int = 20,
    mfold: Sequence[int] = (5, 50),
    script_name: str = "predictd_model.R",
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the fragment-length model as an R plotting script and return its data and, when possible, the image.

    macs3 predictd writes the R script (the function behind the MACS3 ``model2r_script`` writer); the model values are
    parsed from it: ``forward_percent``/``reverse_percent`` (tag profiles, in percent per distance from the middle),
    ``xcorr_lag``/``xcorr_correlation`` and ``alternative_fragment_lengths``. The image is rendered only when ``Rscript``
    is on PATH; otherwise ``image`` is None and ``image_note`` says so. The image is ``<script_name>_model.pdf`` in
    ``output_dir`` (the name the MACS3 script gives its PDF). The default script name is the predictd default.
    """
    script_file = _name(script_name, "script_name")
    result = _predict(
        _inputs(ifile, "ifile"),
        output_dir,
        file_format=file_format,
        genome_size=genome_size,
        tag_size=tag_size,
        bandwidth=bandwidth,
        min_fragment=min_fragment,
        mfold=mfold,
        script_name=script_file,
        macs3_path=macs3_path,
    )
    out = Path(result["output_dir"])
    script = out / script_file
    result["script"] = str(script) if script.is_file() else None
    result["data"] = None
    result["image"] = None
    result["image_note"] = None
    if result["status"] != "completed" or not script.is_file():
        return result
    result["data"] = _parse_model_script(script)
    pdf = out / f"{script_file}_model.pdf"
    if shutil.which("Rscript") is None:
        result["image_note"] = (
            "Rscript is not on PATH: the R script was written but no image was rendered. "
            f"Run `Rscript {script}` in {out} to draw {pdf.name}."
        )
        return result
    rendered = _cli.run(
        "Rscript",
        [str(script)],
        cwd=out,
        expected=[pdf],
        explicit=None,
        hint="Install R (for example conda install -c conda-forge r-base) to render the fragment model figure.",
    )
    result["image"] = str(pdf) if rendered["status"] == "completed" else None
    result["image_note"] = None if rendered["status"] == "completed" else rendered["log_tail"]
    return result


# ---------------------------------------------------------------------------
# bedGraph and track operations: bdgcmp, bdgopt, bdgdiff, cmbreps, bdgpeakcall, bdgbroadcall
# ---------------------------------------------------------------------------

_BDGCMP_METHODS = ("ppois", "qpois", "subtract", "logFE", "FE", "logLR", "slogLR", "max")
_BDGOPT_METHODS = ("multiply", "add", "p2q", "max", "min")
_CMBREPS_METHODS = ("fisher", "max", "mean")
_CALLVAR_FERMI = ("auto", "on", "off")


def _one_file(path: str, name: str) -> str:
    return str(_cli.existing_file(path, name))


def _choices(values: Sequence[str], allowed: Sequence[str], name: str) -> list[str]:
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list.")
    return [_cli.choice(value, allowed, name) for value in values]


def _bdgcmp_scores(
    treatment: str,
    control: str,
    output_dir: str,
    *,
    name: str,
    methods: Sequence[str],
    scaling_factor: float,
    pseudocount: float,
    macs3_path: str | None,
) -> dict[str, Any]:
    tfile = _one_file(treatment, "treatment")
    cfile = _one_file(control, "control")
    chosen = _choices(methods, _BDGCMP_METHODS, "methods")
    scale = _cli.number(scaling_factor, "scaling_factor", minimum=0.0)
    pseudo = _cli.number(pseudocount, "pseudocount", minimum=0.0)
    prefix = _name(name, "name")
    out = _new_directory_path(output_dir, "output_dir")
    expected = [out / f"{prefix}_{method}.bdg" for method in chosen]
    executable = _executable(macs3_path)
    _make_directory(out)
    args = [
        "-t",
        tfile,
        "-c",
        cfile,
        "-S",
        str(scale),
        "-p",
        str(pseudo),
        "-m",
        *chosen,
        "--outdir",
        str(out),
        "--o-prefix",
        prefix,
    ]
    result = _macs3("bdgcmp", args, cwd=out, expected=expected, executable=executable)
    result["output_dir"] = str(out)
    result["score_files"] = [str(path) for path in expected if path.is_file()]
    return result


def compare_bedgraph_tracks(
    treatment: str,
    control: str,
    output_dir: str,
    *,
    name: str = "NA",
    methods: Sequence[str] = ("ppois",),
    scaling_factor: float = 1.0,
    pseudocount: float = 0.0,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score a treatment bedGraph against a control bedGraph, one score track per method (macs3 bdgcmp).

    ``treatment`` is a pileup track and ``control`` a control lambda track (both written by ``callpeak`` with
    ``save_bedgraph``). Each method writes ``<output_dir>/<name>_<method>.bdg``; ``ppois`` is the Poisson p-value as
    -log10, ``qpois`` its q-value, ``FE`` the fold enrichment (pseudocount applied after depth scaling).
    """
    return _bdgcmp_scores(
        treatment,
        control,
        output_dir,
        name=name,
        methods=methods,
        scaling_factor=scaling_factor,
        pseudocount=pseudocount,
        macs3_path=macs3_path,
    )


def operate_bedgraph_scores(
    bedgraph: str,
    output_file: str,
    *,
    method: str = "p2q",
    extra_params: Sequence[float] | None = None,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transform the score column of a MACS3 bedGraph (macs3 bdgopt).

    ``multiply``, ``add``, ``max`` and ``min`` take exactly one extra parameter (``extra_params=[x]``; use ``add`` with
    ``-x`` to subtract). ``p2q`` converts -log10 p-value scores to q-values by Benjamini-Hochberg and takes none; it
    needs a bedGraph covering the whole genome, as MACS3 documents.
    """
    source = _one_file(bedgraph, "bedgraph")
    chosen = _cli.choice(method, _BDGOPT_METHODS, "method")
    params: list[str] = []
    if chosen == "p2q":
        if extra_params:
            raise ConfigurationError("method 'p2q' takes no extra parameter.")
    else:
        if not extra_params or len(extra_params) != 1:
            raise ConfigurationError(f"method {chosen!r} needs exactly one extra parameter.")
        params = [str(_cli.number(extra_params[0], "extra_params", minimum=None))]
    out = _cli.new_file(output_file, "output_file")
    executable = _executable(macs3_path)
    args = ["-i", source, "-m", chosen, *(["-p", *params] if params else []), "-o", str(out)]
    return _macs3("bdgopt", args, cwd=out.parent, expected=[out], executable=executable)


def differential_bedgraph_peaks(
    treatment1: str,
    treatment2: str,
    control1: str,
    control2: str,
    output_condition1: str,
    output_condition2: str,
    output_common: str,
    *,
    cutoff: float = 3.0,
    min_length: int = 200,
    max_gap: int = 100,
    depth1: float = 1.0,
    depth2: float = 1.0,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call regions that differ between two conditions from their pileup and control tracks (macs3 bdgdiff).

    The four inputs are the pileup (``treatment``) and control lambda (``control``) bedGraphs of each condition. The three
    new output files receive the regions specific to condition 1, specific to condition 2 and common to both. ``cutoff``
    is the log10 likelihood-ratio cutoff; ``depth1``/``depth2`` are the sequencing depths in millions of non-redundant reads.
    """
    files = [
        _one_file(path, name)
        for path, name in (
            (treatment1, "treatment1"),
            (treatment2, "treatment2"),
            (control1, "control1"),
            (control2, "control2"),
        )
    ]
    cut = _cli.number(cutoff, "cutoff", minimum=0.0)
    gap = _cli.integer(max_gap, "max_gap", minimum=1)
    length = _cli.integer(min_length, "min_length", minimum=1)
    depth_a = _cli.number(depth1, "depth1", minimum=0.0)
    depth_b = _cli.number(depth2, "depth2", minimum=0.0)
    outputs = [
        _cli.new_file(output_condition1, "output_condition1"),
        _cli.new_file(output_condition2, "output_condition2"),
        _cli.new_file(output_common, "output_common"),
    ]
    executable = _executable(macs3_path)
    args = [
        "--t1",
        files[0],
        "--t2",
        files[1],
        "--c1",
        files[2],
        "--c2",
        files[3],
        "-C",
        str(cut),
        "-l",
        str(length),
        "-g",
        str(gap),
        "--d1",
        str(depth_a),
        "--d2",
        str(depth_b),
        "-o",
        *(str(path) for path in outputs),
    ]
    return _macs3("bdgdiff", args, cwd=outputs[0].parent, expected=outputs, executable=executable)


def combine_replicate_bedgraphs(
    inputs: Sequence[str],
    output_file: str,
    *,
    method: str = "fisher",
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine the score tracks of two or more replicates into one bedGraph (macs3 cmbreps).

    ``fisher`` (default) combines -log10 p-values from bdgcmp ``ppois`` tracks with Fisher's method; ``max`` and ``mean``
    take the maximum or the average of the scores as they are.
    """
    files = [_one_file(path, "inputs") for path in inputs] if not isinstance(inputs, str) else []
    if len(files) < 2:
        raise ConfigurationError(
            "combine_replicate_bedgraphs needs at least two replicate bedGraph files."
        )
    chosen = _cli.choice(method, _CMBREPS_METHODS, "method")
    out = _cli.new_file(output_file, "output_file")
    executable = _executable(macs3_path)
    args = ["-i", *files, "-m", chosen, "-o", str(out)]
    return _macs3("cmbreps", args, cwd=out.parent, expected=[out], executable=executable)


def _bdg_peaks(
    subcommand: str,
    arguments: list[str],
    *,
    bedgraph: str,
    output_file: str,
    macs3_path: str | None,
    trackline: bool,
    call_summits: bool,
) -> dict[str, Any]:
    source = _one_file(bedgraph, "bedgraph")
    out = _cli.new_file(output_file, "output_file")
    executable = _executable(macs3_path)
    args = ["-i", source, *arguments, "-o", str(out)]
    args += _switch("--call-summits", call_summits) + _switch("--no-trackline", not trackline)
    return _macs3(subcommand, args, cwd=out.parent, expected=[out], executable=executable)


def call_peaks_from_bedgraph(
    bedgraph: str,
    output_file: str,
    *,
    cutoff: float = 5.0,
    min_length: int = 200,
    max_gap: int = 30,
    call_summits: bool = False,
    trackline: bool = True,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call narrow peaks from a score bedGraph such as a bdgcmp track (macs3 bdgpeakcall).

    Regions whose score is at least ``cutoff`` are merged across gaps up to ``max_gap`` and kept when at least
    ``min_length`` bp long. ``trackline`` (default on, as in MACS3) writes the UCSC track line.
    """
    numbers = [
        "-c",
        str(_cli.number(cutoff, "cutoff")),
        "-l",
        str(_cli.integer(min_length, "min_length", minimum=1)),
        "-g",
        str(_cli.integer(max_gap, "max_gap", minimum=0)),
    ]
    return _bdg_peaks(
        "bdgpeakcall",
        numbers,
        bedgraph=bedgraph,
        output_file=output_file,
        macs3_path=macs3_path,
        trackline=trackline,
        call_summits=call_summits,
    )


def call_broad_peaks_from_bedgraph(
    bedgraph: str,
    output_file: str,
    *,
    cutoff_peak: float = 2.0,
    cutoff_link: float = 1.0,
    min_length: int = 200,
    lvl1_max_gap: int = 30,
    lvl2_max_gap: int = 800,
    trackline: bool = True,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call nested broad peaks from a score bedGraph (macs3 bdgbroadcall), written as a gappedPeak table (15 columns).

    ``cutoff_peak`` selects the broad regions, ``cutoff_link`` (lower) links the nearby cores inside them; gaps up to
    ``lvl1_max_gap`` merge cores and up to ``lvl2_max_gap`` merge broad regions.
    """
    numbers = [
        "-c",
        str(_cli.number(cutoff_peak, "cutoff_peak")),
        "-C",
        str(_cli.number(cutoff_link, "cutoff_link")),
        "-l",
        str(_cli.integer(min_length, "min_length", minimum=1)),
        "-g",
        str(_cli.integer(lvl1_max_gap, "lvl1_max_gap", minimum=0)),
        "-G",
        str(_cli.integer(lvl2_max_gap, "lvl2_max_gap", minimum=0)),
    ]
    return _bdg_peaks(
        "bdgbroadcall",
        numbers,
        bedgraph=bedgraph,
        output_file=output_file,
        macs3_path=macs3_path,
        trackline=trackline,
        call_summits=False,
    )


# ---------------------------------------------------------------------------
# alignment operations: pileup, filterdup, randsample, refinepeak, callvar, hmmratac
# ---------------------------------------------------------------------------


def pileup_reads_to_bedgraph(
    ifile: Sequence[str],
    output_file: str,
    *,
    file_format: str = "AUTO",
    barcodes: str | None = None,
    max_count: int | None = None,
    both_direction: bool = False,
    extsize: int = 200,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pile up alignments (single-end reads extended to ``extsize`` bp, or fragments) into a bedGraph (macs3 pileup).

    Reads are extended downstream by ``extsize`` bp, or in both directions with ``both_direction``. For paired-end
    formats (BAMPE, BEDPE) and FRAG the fragments are piled up as they are; ``barcodes`` and ``max_count`` apply to FRAG.
    This is the implementation behind ``macs3 pileup``, which the console routes to ``pileup_v2_cmd.run``.
    """
    inputs = _inputs(ifile, "ifile")
    fmt = _cli.choice(file_format, _CALLPEAK_FORMATS, "file_format")
    extend = _cli.integer(extsize, "extsize", minimum=1)
    count = None if max_count is None else _cli.integer(max_count, "max_count", minimum=0)
    barcode_file = _optional_file(barcodes, "barcodes")
    out = _cli.new_file(output_file, "output_file")
    executable = _executable(macs3_path)
    args = ["-i", *inputs, "-o", str(out), "-f", fmt, "--extsize", str(extend)]
    args += (
        _valued("--barcodes", barcode_file)
        + _valued("--max-count", count)
        + _switch("-B", both_direction)
    )
    return _macs3("pileup", args, cwd=out.parent, expected=[out], executable=executable)


def filter_duplicate_reads(
    ifile: Sequence[str],
    output_file: str,
    *,
    file_format: str = "AUTO",
    genome_size: str | float = "hs",
    tag_size: int | None = None,
    keep_duplicates: str | int = "auto",
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove duplicate reads and write the rest as BED (BEDPE for paired-end input) (macs3 filterdup).

    ``keep_duplicates``: ``auto`` (default) keeps the number of tags at one position and strand that a binomial test
    allows for the genome size; ``all`` keeps every tag; an integer keeps at most that many.
    """
    inputs = _inputs(ifile, "ifile")
    fmt = _cli.choice(file_format, _TAG_FORMATS, "file_format")
    gsize = _genome_size(genome_size)
    tsize = _optional_int(tag_size, "tag_size")
    keep = _keep_duplicates(keep_duplicates)
    out = _cli.new_file(output_file, "output_file")
    executable = _executable(macs3_path)
    args = ["-i", *inputs, "-f", fmt, "-g", gsize, "--keep-dup", keep, "-o", str(out)]
    args += _valued("-s", tsize)
    return _macs3("filterdup", args, cwd=out.parent, expected=[out], executable=executable)


def random_sample_reads(
    ifile: Sequence[str],
    output_file: str,
    *,
    number: int | None = None,
    percentage: float | None = None,
    seed: int | None = None,
    file_format: str = "AUTO",
    tag_size: int | None = None,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly keep a number or a percentage of the alignments and write them as BED (macs3 randsample).

    Give exactly one of ``number`` (approximately that many tags) or ``percentage`` (80.0 keeps 80%). ``seed`` is passed
    to MACS3, but repeated runs with one seed did not always give the same sample when tested, so the sample is not
    guaranteed to be identical between runs.
    """
    inputs = _inputs(ifile, "ifile")
    if (number is None) == (percentage is None):
        raise ConfigurationError("Give exactly one of number or percentage.")
    fmt = _cli.choice(file_format, _TAG_FORMATS, "file_format")
    if number is not None:
        amount = ["-n", str(_cli.integer(number, "number", minimum=1))]
    else:
        amount = ["-p", str(_cli.number(percentage, "percentage", minimum=0.0, maximum=100.0))]
    seed_value = None if seed is None else str(_cli.integer(seed, "seed", minimum=0))
    tsize = _optional_int(tag_size, "tag_size")
    out = _cli.new_file(output_file, "output_file")
    executable = _executable(macs3_path)
    args = ["-i", *inputs, *amount, "-o", str(out), "-f", fmt]
    args += _valued("--seed", seed_value) + _valued("-s", tsize)
    return _macs3("randsample", args, cwd=out.parent, expected=[out], executable=executable)


def refine_peak_summits(
    peaks_bed: str,
    ifile: Sequence[str],
    output_file: str,
    *,
    file_format: str = "AUTO",
    cutoff: float = 5.0,
    window_size: int = 200,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Refine the summits of candidate peaks from the alignments (macs3 refinepeak, SPP-style).

    ``peaks_bed`` lists the candidate regions; each summit is re-placed within ``window_size`` bp on either side using the
    tag positions of ``ifile``, and regions below ``cutoff`` (the weighted score) are dropped. The help text of MACS3 says
    100 bp for the window; the parser default that is applied is 200, which this facade uses.
    """
    bed = _one_file(peaks_bed, "peaks_bed")
    inputs = _inputs(ifile, "ifile")
    fmt = _cli.choice(file_format, _ALIGNMENT_FORMATS, "file_format")
    out = _cli.new_file(output_file, "output_file")
    executable = _executable(macs3_path)
    args = [
        "-b",
        bed,
        "-i",
        *inputs,
        "-f",
        fmt,
        "-c",
        str(_cli.number(cutoff, "cutoff")),
        "-w",
        str(_cli.integer(window_size, "window_size", minimum=1)),
        "-o",
        str(out),
    ]
    return _macs3("refinepeak", args, cwd=out.parent, expected=[out], executable=executable)


def call_variants_in_peaks(
    peaks_bed: str,
    treatment_bam: str,
    output_vcf: str,
    *,
    control_bam: str | None = None,
    min_base_quality: int = 20,
    max_duplicates: int = 1,
    gq_hetero: int = 0,
    gq_homo: int = 0,
    fermi: str = "auto",
    fermi_min_overlap: int = 30,
    top2alleles_min_ratio: float = 0.8,
    altallele_min_count: int = 2,
    max_allowed_ar: float = 0.95,
    processes: int = 1,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call SNVs and small indels inside peak regions from a BAM with MD tags (macs3 callvar) and write a VCF.

    The BAM must be sorted and indexed and carry MD tags (``samtools calmd``). MACS3 reports a region with no readable reads
    as skipped, with exit code 0; the result's ``warnings`` lists this when it appears in the log tail. ``fermi`` is
    ``auto`` (assemble only where needed), ``on`` or ``off``.
    """
    bed = _one_file(peaks_bed, "peaks_bed")
    treat = _one_file(treatment_bam, "treatment_bam")
    control = _optional_file(control_bam, "control_bam")
    base_quality = _cli.integer(min_base_quality, "min_base_quality", minimum=0)
    duplicates = _cli.integer(max_duplicates, "max_duplicates", minimum=1)
    hetero = _cli.integer(gq_hetero, "gq_hetero", minimum=0)
    homo = _cli.integer(gq_homo, "gq_homo", minimum=0)
    assembly = _cli.choice(fermi, _CALLVAR_FERMI, "fermi")
    overlap = _cli.integer(fermi_min_overlap, "fermi_min_overlap", minimum=1)
    ratio = _cli.number(top2alleles_min_ratio, "top2alleles_min_ratio", minimum=0.0, maximum=1.0)
    alt_count = _cli.integer(altallele_min_count, "altallele_min_count", minimum=1)
    max_ar = _cli.number(max_allowed_ar, "max_allowed_ar", minimum=0.0, maximum=1.0)
    workers = _cli.integer(processes, "processes", minimum=1)
    out = _cli.new_file(output_vcf, "output_vcf")
    executable = _executable(macs3_path)
    args = [
        "-b",
        bed,
        "-t",
        treat,
        "-o",
        str(out),
        "-Q",
        str(base_quality),
        "-D",
        str(duplicates),
        "-g",
        str(hetero),
        "-G",
        str(homo),
        "-F",
        assembly,
        "--fermi-overlap",
        str(overlap),
        "--top2alleles-mratio",
        str(ratio),
        "--altallele-count",
        str(alt_count),
        "--max-ar",
        str(max_ar),
        "-m",
        str(workers),
    ]
    args += _valued("-c", control)
    result = _macs3("callvar", args, cwd=out.parent, expected=[out], executable=executable)
    skipped = result.get("log_tail", "").count("No reads found in this peak")
    result["warnings"] = (
        [
            f"{skipped} peak(s) in the log tail had no reads MACS3 could read; check the BAM has MD tags and an index."
        ]
        if skipped
        else []
    )
    return result


def call_atac_peaks_hmmratac(
    ifile: Sequence[str],
    output_dir: str,
    *,
    name: str = "NA",
    file_format: str = "BAMPE",
    lower_cutoff: int = 10,
    upper_cutoff: int = 20,
    prescan_cutoff: float = 1.2,
    min_region: int = 100,
    max_training: int = 1000,
    random_seed: int = 10151,
    hmm_type: str = "gaussian",
    training_regions: str | None = None,
    save_states: bool = False,
    macs3_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call ATAC-seq accessible regions from paired-end fragments with the HMMRATAC hidden Markov model (macs3 hmmratac).

    ``ifile`` holds fragments in BAMPE, BEDPE or FRAG format. Outputs go to a new ``output_dir``:
    ``<name>_accessible_regions.narrowPeak``, ``<name>_cutoff_analysis.tsv`` and ``<name>_model.json``, plus
    ``<name>_states.bed`` with ``save_states``. Training regions are the regions whose fold change over the average lies
    between ``lower_cutoff`` and ``upper_cutoff``; when none qualify MACS3 stops with "Not enough training regions".
    ``training_regions`` (a BED) replaces the automatic choice.
    """
    inputs = _inputs(ifile, "ifile")
    prefix = _name(name, "name")
    fmt = _cli.choice(file_format, ("BAMPE", "BEDPE", "FRAG"), "file_format")
    low = _cli.integer(lower_cutoff, "lower_cutoff", minimum=0)
    high = _cli.integer(upper_cutoff, "upper_cutoff", minimum=1)
    if low >= high:
        raise ConfigurationError("lower_cutoff must be smaller than upper_cutoff.")
    prescan = _cli.number(prescan_cutoff, "prescan_cutoff", minimum=0.0)
    min_len = _cli.integer(min_region, "min_region", minimum=1)
    max_train = _cli.integer(max_training, "max_training", minimum=1)
    seed = _cli.integer(random_seed, "random_seed", minimum=0)
    model_type = _cli.choice(hmm_type, ("gaussian", "poisson"), "hmm_type")
    regions = _optional_file(training_regions, "training_regions")
    out = _new_directory_path(output_dir, "output_dir")
    stem = out / prefix
    expected = [
        Path(f"{stem}_accessible_regions.narrowPeak"),
        Path(f"{stem}_cutoff_analysis.tsv"),
        Path(f"{stem}_model.json"),
    ]
    if save_states:
        expected.append(Path(f"{stem}_states.bed"))
    executable = _executable(macs3_path)
    _make_directory(out)
    args = [
        "-i",
        *inputs,
        "-f",
        fmt,
        "-n",
        prefix,
        "--outdir",
        str(out),
        "-l",
        str(low),
        "-u",
        str(high),
        "-c",
        str(prescan),
        "--minlen",
        str(min_len),
        "--maxTrain",
        str(max_train),
        "--randomSeed",
        str(seed),
        "--hmm-type",
        model_type,
    ]
    args += _valued("-t", regions) + _switch("--save-states", save_states)
    result = _macs3("hmmratac", args, cwd=out, expected=expected, executable=executable)
    result["output_dir"] = str(out)
    result["accessible_regions"] = _existing(Path(f"{stem}_accessible_regions.narrowPeak"))
    return result


__all__ = [
    "call_atac_peaks_hmmratac",
    "call_broad_peaks_from_bedgraph",
    "call_chip_peaks",
    "call_peaks_from_bedgraph",
    "call_variants_in_peaks",
    "combine_replicate_bedgraphs",
    "compare_bedgraph_tracks",
    "differential_bedgraph_peaks",
    "filter_duplicate_reads",
    "operate_bedgraph_scores",
    "pileup_reads_to_bedgraph",
    "plot_fragment_model",
    "predict_fragment_size",
    "random_sample_reads",
    "refine_peak_summits",
]
_register_public_verbose(globals())
