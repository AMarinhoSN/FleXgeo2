from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class DistanceResult:
    long_df: object
    summary_df: object
    reference_label: str


@dataclass(slots=True)
class ResidueClusteringResult:
    assignments_df: object
    summary_df: object


@dataclass(slots=True)
class ResidueRangeClusteringResult:
    assignments_df: object
    summary_df: object


@dataclass(slots=True)
class OutputArtifacts:
    raw_csv: Path | None = None
    residue_summary_csv: Path | None = None
    model_summary_csv: Path | None = None
    overall_model_summary_csv: Path | None = None
    overview_plot: Path | None = None
    chains_dir: Path | None = None
    distance_long_csv: Path | None = None
    distance_summary_csv: Path | None = None
    distance_heatmap: Path | None = None
    distance_matrix_dir: Path | None = None
    cluster_assignments_csv: Path | None = None
    cluster_summary_csv: Path | None = None
    cluster_plots_dir: Path | None = None
    range_cluster_assignments_csv: Path | None = None
    range_cluster_summary_csv: Path | None = None
    range_cluster_plots_dir: Path | None = None


@dataclass(slots=True)
class AnalysisResult:
    pdb_file: Path
    raw_df: object
    residue_summary_df: object
    model_summary_df: object
    overall_model_summary_df: object
    distance_result: DistanceResult | None = None
    residue_clustering: ResidueClusteringResult | None = None
    residue_range_clustering: ResidueRangeClusteringResult | None = None
    outputs: OutputArtifacts | None = None
