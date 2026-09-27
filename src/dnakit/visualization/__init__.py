"""Dependency-free, deterministic SVG visualizations for DNAKit results."""

from dnakit._verbosity import register_public as _register_public_verbose
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
]

_register_public_verbose(globals())
