"""pysam wrapper for the bcftools genotype-export formats that the existing wrappers do not reach.

``export_vcf_to_phasing_format`` runs the bcftools build bundled with pysam. The package is
imported only when the function runs. A missing package raises ``BackendUnavailableError`` with the
``bio-reads`` install hint.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError
from dnakit.variants.pysam_functions import (
    _call,
    _choice,
    _input,
    _integer,
    _provenance,
    _pysam,
    _submodule,
)

_FORMATS = {
    "gensample": ("--gensample", (".gen.gz", ".samples"), ".gen.gz"),
    "hapsample": ("--hapsample", (".hap.gz", ".samples"), ".hap.gz"),
    "haplegendsample": (
        "--haplegendsample",
        (".hap.gz", ".legend.gz", ".samples"),
        ".hap.gz",
    ),
}
_MAX_ROWS = 100_000


def _prefix(value: object) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError("output_prefix must be a path prefix.")
    text = str(value)
    if text.endswith(("/", "\\")) or any(char in text for char in "\x00\r\n,"):
        raise ConfigurationError("output_prefix must end with a file-name prefix.")
    path = Path(text).expanduser().resolve()
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for output_prefix does not exist: {path.parent}")
    return path


def export_vcf_to_phasing_format(
    source: str | Path,
    output_prefix: str | Path,
    *,
    fmt: str = "gensample",
    max_rows: int = _MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Export VCF genotypes to the IMPUTE2 or SHAPEIT text formats with ``bcftools convert``.

    ``fmt="gensample"`` writes ``<prefix>.gen.gz`` (three genotype probabilities per sample and
    site) and ``<prefix>.samples``; ``"hapsample"`` writes ``<prefix>.hap.gz`` (two haplotype
    alleles per sample, so phased genotypes are expected) and ``<prefix>.samples``;
    ``"haplegendsample"`` writes ``<prefix>.hap.gz``, ``<prefix>.legend.gz`` and
    ``<prefix>.samples``. Sites are named ``chrom:pos_ref_alt``. The output files must not exist
    yet. The result lists the files, the sample IDs and the first ``max_rows`` lines of the main
    data file as text.
    """

    path = _input(source, "source")
    selected = _choice(fmt, "fmt", tuple(_FORMATS))
    limit = _integer(max_rows, "max_rows", minimum=1, maximum=_MAX_ROWS)
    prefix = _prefix(output_prefix)
    flag, suffixes, main_suffix = _FORMATS[selected]
    targets = [Path(str(prefix) + suffix) for suffix in suffixes]
    for target in targets:
        if target.exists():
            raise ConfigurationError(f"{target.name} already exists and will not be overwritten.")
    pysam = _pysam()
    bcftools = _submodule("pysam.bcftools")
    try:
        # catch_stdout=False lets bcftools write the files; the default capture writes none.
        _call(
            pysam,
            "bcftools convert",
            bcftools.convert,
            flag,
            str(prefix),
            str(path),
            catch_stdout=False,
        )
    except BaseException:
        for target in targets:
            target.unlink(missing_ok=True)
        raise
    missing = [target.name for target in targets if not target.is_file()]
    if missing:
        raise ConfigurationError(f"bcftools convert did not write {', '.join(missing)}.")
    main = Path(str(prefix) + main_suffix)
    with gzip.open(main, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    samples_file = Path(str(prefix) + ".samples")
    sample_lines = samples_file.read_text(encoding="utf-8").splitlines()
    # gensample and hapsample files have a header and a "0 0 0" line; haplegendsample has a header
    skip = 1 if selected == "haplegendsample" else 2
    sample_ids = [line.split()[0] for line in sample_lines[skip:] if line.strip()]
    return {
        "status": "completed",
        "operation": "export_vcf_to_phasing_format",
        "source": str(path),
        "format": selected,
        "files": [str(target) for target in targets],
        "site_count": len(lines),
        "sample_ids": sample_ids,
        "rows": lines[:limit],
        "rows_truncated": len(lines) > limit,
        "provenance": _provenance(
            "bcftools convert; sites are named chrom:pos_ref_alt; the main file is "
            f"{main.name} and rows are its lines"
        ),
    }


__all__ = ["export_vcf_to_phasing_format"]

_register_public_verbose(globals(), tuple(__all__))
