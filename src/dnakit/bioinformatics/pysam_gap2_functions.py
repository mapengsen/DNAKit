"""pysam operations not covered by pysam_gap_functions: quality-string conversion, tabix iteration and
samtools-style utilities (bam2fq, bamshuf, cram_size, depad, pad2unpad, targetcut).

Each function is a thin, typed wrapper over pysam and returns plain Python values or the text that
pysam returns. Generic command dispatchers and global configuration helpers are not exposed here.
"""

from __future__ import annotations

from array import array
from collections.abc import Iterable, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require as _require
from dnakit.exceptions import ConfigurationError


class _LazyPysam:
    """Import pysam on first use: it is an optional dependency (extra ``bio-reads``)."""

    def __getattr__(self, name: str) -> Any:
        return getattr(_require("pysam", "bio-reads"), name)


pysam = _LazyPysam()


def _args(arguments: Sequence[str], name: str) -> list[str]:
    if isinstance(arguments, str) or not all(isinstance(item, str) for item in arguments):
        raise ConfigurationError(f"{name} must be a sequence of strings.")
    return list(arguments)


def qualitystring_to_array(quality: str, *, offset: int = 33, verbose: bool = True) -> list[int]:
    """Convert a Phred quality string to integer quality values (pysam.qualitystring_to_array)."""
    if not isinstance(quality, str):
        raise ConfigurationError("quality must be a string.")
    return [int(value) for value in pysam.qualitystring_to_array(quality, offset=offset)]


def array_to_qualitystring(
    qualities: Iterable[int], *, offset: int = 33, verbose: bool = True
) -> str:
    """Convert integer quality values to a Phred quality string (pysam.array_to_qualitystring)."""
    return str(pysam.array_to_qualitystring(array("B", [int(q) for q in qualities]), offset=offset))


def qualities_to_qualitystring(
    qualities: Iterable[int], *, offset: int = 33, verbose: bool = True
) -> str:
    """Convert a list of quality scores to a quality string (pysam.qualities_to_qualitystring)."""
    return str(
        pysam.qualities_to_qualitystring(array("B", [int(q) for q in qualities]), offset=offset)
    )


def tabix_records(path: str, parser: Any, *, verbose: bool = True) -> list[Any]:
    """Read all records from a tabix-indexed file with a pysam parser (pysam.tabix_iterator)."""
    with open(path, "rb") as handle:
        return list(pysam.tabix_iterator(handle, parser))


def bam_to_fastq(arguments: Sequence[str], output_path: str, *, verbose: bool = True) -> str:
    """Convert alignments to FASTQ with samtools fastq (pysam.bam2fq) and return the FASTQ text.

    pysam writes samtools' output to the process stdout, so the result is written to ``output_path``
    with ``-0`` and read back. The output path must not exist yet.
    """
    from pathlib import Path

    target = Path(output_path)
    if target.exists():
        raise ConfigurationError("output_path must not exist yet.")
    pysam.bam2fq(*_args(arguments, "arguments"), "-0", str(target))
    return target.read_text(encoding="utf-8")


def shuffle_bam(arguments: Sequence[str], *, verbose: bool = True) -> str | list[str]:
    """Run bamshuf with explicit arguments through pysam and return its output."""
    return pysam.bamshuf(*_args(arguments, "arguments"))


def cram_container_sizes(arguments: Sequence[str], *, verbose: bool = True) -> str | list[str]:
    """Report CRAM container sizes with explicit arguments through pysam (cram_size)."""
    return pysam.cram_size(*_args(arguments, "arguments"))


def remove_padding(arguments: Sequence[str], *, verbose: bool = True) -> str | list[str]:
    """Remove alignment padding from a padded alignment with explicit arguments (depad)."""
    return pysam.depad(*_args(arguments, "arguments"))


def pad_to_unpadded(arguments: Sequence[str], *, verbose: bool = True) -> str | list[str]:
    """Convert a padded alignment to unpadded coordinates with explicit arguments (pad2unpad)."""
    return pysam.pad2unpad(*_args(arguments, "arguments"))


def target_cut(arguments: Sequence[str], *, verbose: bool = True) -> str | list[str]:
    """Cut reads around target regions with explicit arguments through pysam (targetcut)."""
    return pysam.targetcut(*_args(arguments, "arguments"))


__all__ = [
    "array_to_qualitystring",
    "bam_to_fastq",
    "cram_container_sizes",
    "pad_to_unpadded",
    "qualities_to_qualitystring",
    "qualitystring_to_array",
    "remove_padding",
    "shuffle_bam",
    "tabix_records",
    "target_cut",
]

_register_public_verbose(globals())
