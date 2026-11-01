"""Dependency-free, deterministic SVG visualizations for DNAKit results."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.visualization._candidates import *  # noqa: F403
from dnakit.visualization._candidates import __all__ as _candidate_all
from dnakit.visualization.advanced import (
    plot_alignment,
    plot_circular_map,
    plot_linear_map,
)
from dnakit.visualization.config import (
    DistributionPlotConfig,
    HeatmapConfig,
    ImageExportConfig,
    LimitPolicy,
    SaveConfig,
    SequencePlotConfig,
    SVGTheme,
)
from dnakit.visualization.confusion import plot_confusion_matrix_csv
from dnakit.visualization.distribution import (
    DistributionTask,
    plot_csv_distributions,
    plot_distribution,
    plot_label_comparison,
)
from dnakit.visualization.embedding_errors import plot_embedding_errors
from dnakit.visualization.export import (
    ImageExportFormat,
    PathLike,
    save_html_report,
    save_image,
    save_svg,
)
from dnakit.visualization.figure import save_figure
from dnakit.visualization.learning import plot_learning_curves_csv
from dnakit.visualization.matrix import plot_similarity_matrix
from dnakit.visualization.metric import plot_metric_comparison_csv
from dnakit.visualization.report import build_html_report
from dnakit.visualization.results import HTMLReportArtifact, SVGArtifact, VisualizationSaveResult
from dnakit.visualization.sequence import Highlight, plot_sequence

__all__ = [
    "DistributionPlotConfig",
    "DistributionTask",
    "HTMLReportArtifact",
    "HeatmapConfig",
    "Highlight",
    "ImageExportConfig",
    "ImageExportFormat",
    "LimitPolicy",
    "PathLike",
    "SVGArtifact",
    "SVGTheme",
    "SaveConfig",
    "SequencePlotConfig",
    "VisualizationSaveResult",
    "build_html_report",
    "plot_alignment",
    "plot_circular_map",
    "plot_confusion_matrix_csv",
    "plot_csv_distributions",
    "plot_distribution",
    "plot_embedding_errors",
    "plot_label_comparison",
    "plot_learning_curves_csv",
    "plot_linear_map",
    "plot_metric_comparison_csv",
    "plot_sequence",
    "plot_similarity_matrix",
    "save_figure",
    "save_html_report",
    "save_image",
    "save_svg",
    *_candidate_all,
]


# Coverage-gap exports (see paper/write_paper/figure/fig2/task_coverage/coverage_gap/waves/CANONICAL_MODULES.md)
from .pygenomeviz_gap_functions import (
    compare_genomes_mmseqs,
    compare_genomes_mummer,
    compare_genomes_progressive_mauve,
    load_example_dataset,
    load_example_gff,
    plot_genome_features,
)

__all__ = [*__all__, "compare_genomes_mmseqs", "compare_genomes_mummer", "compare_genomes_progressive_mauve", "load_example_dataset", "load_example_gff", "plot_genome_features"]

# Public integration of the reviewed P00–P45 operation packs.
from .biopython_p19_graphics_functions import (
    chromosome_ideogram,
    genome_diagram_crosslinks,
)
from .pygenomeviz_p37_functions import (
    align_genomes_blast,
    align_genomes_last,
    calc_genbank_gc_content,
    draw_genome_segment_features,
    export_genome_html_viewer,
    extract_genome_features,
    extract_gff_exon_features,
    filter_alignment_coordinates,
    list_contig_features,
    parse_alignment_coordinates,
    read_alignment_coordinates_tsv,
    read_contig_sequences,
    read_gff_contig_sizes,
    select_features_within_range,
    write_alignment_coordinates_tsv,
    write_genbank_cds_fasta,
    write_genome_fasta,
)
from .pygenomeviz_p37_promoter_functions import (
    draw_genome_promoters,
)

__all__ = [
    *__all__,
    "align_genomes_blast",
    "align_genomes_last",
    "calc_genbank_gc_content",
    "chromosome_ideogram",
    "draw_genome_promoters",
    "draw_genome_segment_features",
    "export_genome_html_viewer",
    "extract_genome_features",
    "extract_gff_exon_features",
    "filter_alignment_coordinates",
    "genome_diagram_crosslinks",
    "list_contig_features",
    "parse_alignment_coordinates",
    "read_alignment_coordinates_tsv",
    "read_contig_sequences",
    "read_gff_contig_sizes",
    "select_features_within_range",
    "write_alignment_coordinates_tsv",
    "write_genbank_cds_fasta",
    "write_genome_fasta",
]

_register_public_verbose(globals())
