"""deepTools computeGCBias and correctGCBias as typed DNAKit facades.

Both need a 2bit genome (``genome``) and an effective genome size. computeGCBias writes the GC-bias
frequency table; correctGCBias uses that table to write a GC-corrected BAM (or bedGraph-style output as
given by deepTools). Same contract as the other deepTools facades.
"""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.epigenome.deeptools_cli_facades import _existing, _new, _run
from dnakit.exceptions import ConfigurationError


def compute_gc_bias(
    bam: str,
    genome_2bit: str,
    effective_genome_size: int,
    frequencies_output: str,
    *,
    fragment_length: int | None = None,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the GC-bias frequency table of a BAM file against a 2bit genome (deepTools computeGCBias)."""
    if effective_genome_size < 1:
        raise ConfigurationError("effective_genome_size must be a positive integer.")
    bam_path = _existing(bam, "bam")
    genome = _existing(genome_2bit, "genome_2bit")
    out = _new(frequencies_output, "frequencies_output")
    args = ["-b", str(bam_path), "-g", str(genome), "--effectiveGenomeSize", str(effective_genome_size),
            "-freq", str(out), "-p", str(number_of_processors)]
    if fragment_length is not None:
        if fragment_length < 1:
            raise ConfigurationError("fragment_length must be a positive integer.")
        args += ["-l", str(fragment_length)]
    return _run("computeGCBias", args, [out], deeptools_path)


def correct_gc_bias(
    bam: str,
    genome_2bit: str,
    effective_genome_size: int,
    frequencies: str,
    output_bam: str,
    *,
    number_of_processors: int = 1,
    deeptools_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a GC-corrected BAM using a frequency table from compute_gc_bias (deepTools correctGCBias)."""
    if effective_genome_size < 1:
        raise ConfigurationError("effective_genome_size must be a positive integer.")
    bam_path = _existing(bam, "bam")
    genome = _existing(genome_2bit, "genome_2bit")
    freq = _existing(frequencies, "frequencies")
    out = _new(output_bam, "output_bam")
    args = ["-b", str(bam_path), "-g", str(genome), "--effectiveGenomeSize", str(effective_genome_size),
            "-freq", str(freq), "-o", str(out), "-p", str(number_of_processors)]
    return _run("correctGCBias", args, [out], deeptools_path)


__all__ = ["compute_gc_bias", "correct_gc_bias"]
_register_public_verbose(globals())
