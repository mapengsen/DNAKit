"""cooler operations that DNAKit did not wrap (Fig. 2 coverage gap, batch 0 verify).

``cooler_resolve_bins`` implements the BINS argument of ``cooler load`` and ``cooler cload``. The argument is
either ``<chromsizes file>:<binsize in bp>``, which gives fixed-size genomic bins, or the path of a
three-column BED file that defines the bins. The function writes the resolved bin table to a new BED file.
It follows ``cooler.cli._util.parse_bins`` (cooler 0.10.4). One difference: a BED file that fails to parse
(for example a file whose only row has two fields) raises ``ConfigurationError``. Upstream raises
``AttributeError`` there, because its error handler refers to ``pd.parser.CParserError`` and pandas 2.x has no
``pd.parser``.

The function needs cooler, installed with the ``dnakit[bio-cooler]`` extra.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_BED_COLUMNS = ["chrom", "start", "end"]


def _new_output(path: object, name: str) -> Path:
    if not isinstance(path, (str, os.PathLike)) or not str(path):
        raise ConfigurationError(f"{name} must be a file path.")
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _bins_from_chromsizes(text: str, util: Any) -> tuple[Any, int]:
    chromsizes_text, binsize_text = text.rsplit(":", 1)
    chromsizes_path = Path(chromsizes_text).expanduser()
    if not chromsizes_path.is_file():
        raise ConfigurationError(f'Chromsizes file "{chromsizes_text}" not found.')
    try:
        binsize = int(binsize_text)
    except ValueError as exc:
        raise ConfigurationError(
            f'Expected integer binsize argument (bp), got "{binsize_text}".'
        ) from exc
    if binsize <= 0:
        raise ConfigurationError("The binsize must be a positive number of bases.")
    try:
        chromsizes = util.read_chromsizes(str(chromsizes_path), all_names=True)
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f'Could not read chromsizes file "{chromsizes_text}": {exc}') from exc
    return util.binnify(chromsizes, binsize), binsize


def _bins_from_bed(path: Path, pandas: Any) -> Any:
    try:
        table = pandas.read_csv(
            path,
            sep="\t",
            names=_BED_COLUMNS,
            usecols=[0, 1, 2],
            dtype={"chrom": str},
        )
    except (pandas.errors.ParserError, pandas.errors.EmptyDataError, UnicodeDecodeError, ValueError) as exc:
        raise ConfigurationError(f'Failed to parse bins file "{path}": {exc}') from exc
    if table.isna().any().any():
        raise ConfigurationError(f'Bins file "{path}" has rows with missing chrom, start or end fields.')
    try:
        table["start"] = pandas.to_numeric(table["start"], errors="raise").astype("int64")
        table["end"] = pandas.to_numeric(table["end"], errors="raise").astype("int64")
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f'Bins file "{path}" has non-integer start or end fields: {exc}') from exc
    return table[_BED_COLUMNS]


def cooler_resolve_bins(bins: str, output_bed: str, *, verbose: bool = True) -> dict[str, Any]:
    """Resolve a cooler BINS argument and write the bin table to a new BED file.

    ``bins`` is ``<chromsizes file>:<binsize in bp>`` (fixed-size bins from ``cooler.util.binnify``) or the
    path of a three-column BED file (chrom, start, end) that defines the bins. ``output_bed`` must not exist
    yet. Returns the source kind, the number of bins, the chromosomes in order of first appearance, the
    binsize (``None`` for a BED source) and the output path.
    """

    if not isinstance(bins, str) or not bins.strip():
        raise ConfigurationError("bins must be a non-empty string.")
    target = _new_output(output_bed, "output_bed")
    cooler = require("cooler", _COOLER_EXTRA)
    util = require("cooler.util", _COOLER_EXTRA)
    pandas = require("pandas", _COOLER_EXTRA)
    text = bins.strip()
    _, without_drive = os.path.splitdrive(text)
    if ":" in without_drive:
        table, binsize = _bins_from_chromsizes(text, util)
        source = "chromsizes"
    else:
        bed_path = Path(text).expanduser()
        if not bed_path.is_file():
            raise ConfigurationError(
                "Expected bins to be either <path to bins BED file> or "
                "<path to chromsizes file>:<binsize in bp>; got an unknown path."
            )
        table = _bins_from_bed(bed_path, pandas)
        binsize = None
        source = "bed"
    table.to_csv(target, sep="\t", header=False, index=False)
    chromosomes = list(dict.fromkeys(str(name) for name in table["chrom"]))
    return {
        "operation": "resolve_bins",
        "status": "completed",
        "backend": "cooler",
        "package": f"cooler {cooler.__version__}",
        "source": source,
        "binsize": binsize,
        "bin_count": int(len(table)),
        "chromosomes": chromosomes,
        "output_bed": str(target),
        "coordinates": "0-based half-open BED",
    }


__all__ = ["cooler_resolve_bins"]

_register_public_verbose(globals())
