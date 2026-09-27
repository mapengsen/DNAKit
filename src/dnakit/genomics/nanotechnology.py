"""DNA nanostructure construction and oxDNA export via optional scadnano.

Cite Doty, Lee and Sterin (2020), https://doi.org/10.4230/LIPIcs.DNA.2020.9.
Creating a geometry is not molecular simulation or experimental validation.
"""

from __future__ import annotations

import copy
import json
import math
import re
from importlib.metadata import version
from typing import Any, cast

from dnakit._progress import Progress
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def design_dna_duplex(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Create a two-strand scadnano design with explicitly assigned complementary DNA.

    Requires ``dnakit[dna-nanotechnology]``. Return a scadnano JSON object that
    can be saved as a .sc file or passed to ``export_oxdna_design``. This simple
    duplex constructor does not optimize origami folding or mechanical stability.
    """
    if (
        not isinstance(sequence, str)
        or not 1 <= len(sequence) <= 100_000
        or re.fullmatch("[ACGTacgt]+", sequence) is None
    ):
        raise ConfigurationError("sequence must contain 1..100000 canonical DNA bases.")
    sc = require("scadnano", "dna-nanotechnology")
    length = len(sequence)
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Construct DNA geometry", total=2)
        forward = sc.Strand([sc.Domain(helix=0, forward=True, start=0, end=length)])
        reverse = sc.Strand([sc.Domain(helix=0, forward=False, start=0, end=length)])
        design = sc.Design(
            helices=[sc.Helix(max_offset=length)], strands=[forward, reverse], grid=sc.Grid.square
        )
        progress.advance(task)
        design.assign_dna(forward, sequence.upper())
        result = json.loads(design.to_json())
        progress.advance(task)
    return cast(dict[str, Any], result)


def export_oxdna_design(design: dict[str, Any], *, verbose: bool = True) -> dict[str, Any]:
    """Export caller-supplied scadnano geometry to oxDNA configuration and topology text.

    No files are written and no oxDNA process is run. All strands must have
    complete canonical sequences; wildcard bases cannot become arbitrary bases.
    Only linear strands are supported: scadnano's exporter does not preserve
    circular bonds. Non-finite geometry and invalid backend output are rejected.
    Relaxation, simulation settings and physical interpretation are separate tasks.
    """
    if not isinstance(design, dict):
        raise ConfigurationError("design must be a scadnano JSON mapping.")
    try:
        json.dumps(design, allow_nan=False)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(
            "design must contain JSON values and finite geometry numbers."
        ) from exc
    sc = require("scadnano", "dna-nanotechnology")
    try:
        model = sc.Design.from_scadnano_json_map(copy.deepcopy(design))
    except (ValueError, TypeError, KeyError, IndexError, ArithmeticError, AssertionError) as exc:
        raise ConfigurationError(f"Invalid scadnano design: {exc}") from exc
    if not model.strands:
        raise ConfigurationError("The design must contain at least one strand.")
    nucleotides = 0
    for strand in model.strands:
        if strand.circular:
            raise ConfigurationError(
                "Circular strands are unsupported by the scadnano oxDNA exporter; "
                "export a deliberately linearized design."
            )
        dna = strand.dna_sequence
        if (
            not isinstance(dna, str)
            or re.fullmatch("[ACGT]+", dna) is None
            or len(dna) != strand.dna_length()
        ):
            raise ConfigurationError(
                "Assign complete canonical DNA to every strand before exporting."
            )
        nucleotides += len(dna)
    if nucleotides > 200_000:
        raise ConfigurationError("Export is limited to 200000 nucleotides.")
    with Progress(disable=not verbose) as progress:
        task = progress.add_task("Export oxDNA geometry", total=1)
        try:
            configuration, topology = model.to_oxdna_format(warn_duplicate_strand_names=False)
        except (
            ValueError,
            TypeError,
            KeyError,
            IndexError,
            ArithmeticError,
            AssertionError,
        ) as exc:
            raise ConfigurationError(
                f"scadnano cannot export this geometry to oxDNA: {exc}"
            ) from exc
        configuration_lines = configuration.splitlines()
        topology_lines = topology.splitlines()
        if (
            len(configuration_lines) != nucleotides + 3
            or len(topology_lines) != nucleotides + 1
            or topology_lines[0].split() != [str(nucleotides), str(len(model.strands))]
        ):
            raise ConfigurationError(
                "scadnano exported inconsistent oxDNA nucleotide or strand counts."
            )
        try:
            for line in configuration_lines:
                fields = line.split("=", 1)[-1].split()
                if not fields or not all(math.isfinite(float(value)) for value in fields):
                    raise ValueError("non-finite geometry")
        except (ValueError, OverflowError) as exc:
            raise ConfigurationError(
                "scadnano exported invalid or non-finite oxDNA geometry."
            ) from exc
        progress.advance(task)
    return {
        "configuration": configuration,
        "topology": topology,
        "strand_count": len(model.strands),
        "nucleotide_count": nucleotides,
        "backend_version": version("scadnano"),
        "simulation_executed": False,
        "citations": ["https://doi.org/10.4230/LIPIcs.DNA.2020.9"],
    }
