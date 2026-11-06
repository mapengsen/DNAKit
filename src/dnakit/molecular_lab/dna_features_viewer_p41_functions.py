"""File-based feature-map operations using dna_features_viewer (P41)."""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError


def read_gff_feature_record(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Read the first GFF record using the upstream loader, retaining feature coordinates."""
    require("BCBio.GFF", "bio-dna-features-viewer")
    record = require("dna_features_viewer.biotools", "bio-dna-features-viewer").load_record(
        str(h.existing(path, "path")), filetype="gff"
    )
    return {
        "status": "completed",
        "backend": "dna_features_viewer",
        "id": record.id,
        "length": len(record),
        "features": [
            {
                "type": f.type,
                "start": int(f.location.start),
                "end": int(f.location.end),
                "strand": f.location.strand,
                "qualifiers": dict(f.qualifiers),
            }
            for f in record.features
        ],
    }


def _record(sequence: str, features: list[dict[str, Any]]) -> Any:
    if not isinstance(sequence, str) or not sequence or len(sequence) > 1000000:
        raise ConfigurationError("sequence must contain 1 to 1,000,000 characters.")
    module = require("dna_features_viewer", "bio-dna-features-viewer")
    rows = []
    for row in features:
        start = h.integer(row["start"], "start", minimum=0, maximum=len(sequence))
        end = h.integer(row["end"], "end", minimum=start + 1, maximum=len(sequence))
        strand = row.get("strand", 1)
        if strand not in (-1, 0, 1):
            raise ConfigurationError("strand must be -1, 0 or 1.")
        rows.append(
            module.GraphicFeature(
                start=start,
                end=end,
                strand=strand,
                label=row.get("label"),
                color=row.get("color", "#7245dc"),
            )
        )
    return module.GraphicRecord(sequence=sequence, features=rows)


def export_feature_map_record(
    sequence: str,
    features: list[dict[str, Any]],
    output: str,
    *,
    record_id: str = "sequence",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a GraphicRecord into an annotated DNA GenBank record and write it."""
    target = h.new_file(output, "output")
    record = _record(sequence, features).to_biopython_record(sequence)
    record.id = h.text(record_id, "record_id")
    record.annotations["molecule_type"] = "DNA"
    require("Bio.SeqIO", "bio-dna-features-viewer").write(record, str(target), "genbank")
    return {
        "status": "completed",
        "backend": "dna_features_viewer",
        "feature_count": len(record.features),
        "outputs": [str(target)],
    }


def plot_feature_map(
    sequence: str,
    features: list[dict[str, Any]],
    output: str,
    *,
    layout: str = "nucleotide",
    nucleotides_per_line: int = 60,
    lines_per_page: int = 5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw nucleotide/translation tracks, multiple lines or PDF pages.

    layout=nucleotide/translation/multiline/multipage selects the native method.
    Returns the real image/PDF path and the supplied underlying feature data.
    """
    h.choice(layout, ("nucleotide", "translation", "multiline", "multipage"), "layout")
    width = h.integer(nucleotides_per_line, "nucleotides_per_line", minimum=1, maximum=10000)
    pages = h.integer(lines_per_page, "lines_per_page", minimum=1, maximum=100)
    target = h.new_file(output, "output")
    if target.suffix.lower() not in (
        (".pdf",) if layout == "multipage" else (".png", ".svg", ".pdf")
    ):
        raise ConfigurationError("Use PDF for multipage, or PNG/SVG/PDF for other layouts.")
    record = _record(sequence, features)
    plt = require("matplotlib.pyplot", "bio-dna-features-viewer")
    before = set(plt.get_fignums())
    try:
        if layout == "multipage":
            record.plot_on_multiple_pages(str(target), nucl_per_line=width, lines_per_page=pages)
        elif layout == "multiline":
            figure, _ = record.plot_on_multiple_lines(nucl_per_line=width, plot_sequence=True)
            figure.savefig(target)
        else:
            ax, _ = record.plot()
            if layout == "nucleotide":
                record.plot_sequence(ax)
            else:
                record.plot_translation(ax, location=(0, len(sequence) - len(sequence) % 3))
            ax.figure.savefig(target, bbox_inches="tight")
    finally:
        for number in set(plt.get_fignums()) - before:
            plt.close(number)
    return {
        "status": "completed",
        "backend": "dna_features_viewer",
        "layout": layout,
        "sequence": sequence,
        "features": features,
        "outputs": [str(target)],
    }


__all__ = ["export_feature_map_record", "plot_feature_map", "read_gff_feature_record"]
_register_public_verbose(globals(), tuple(__all__))
