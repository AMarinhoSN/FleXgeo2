from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

from flexgeo2.config import AnalysisConfig, PathLike


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
    distance_long_csv: Path | None = None
    distance_summary_csv: Path | None = None
    distance_heatmap: Path | None = None
    distance_matrix_dir: Path | None = None
    cluster_assignments_csv: Path | None = None
    cluster_summary_csv: Path | None = None
    cluster_map_plot: Path | None = None
    residue_plots_dir: Path | None = None
    range_cluster_assignments_csv: Path | None = None
    range_cluster_summary_csv: Path | None = None
    range_cluster_plots_dir: Path | None = None
    readme: Path | None = None
    run_manifest: Path | None = None
    # Above MAX_CHAINS_PER_FIGURE chains, the overview, heatmap and cluster map are
    # written per chain (overview_A.png, ...): listed here, with the fields above None.
    per_chain_plots: list[Path] = field(default_factory=list)


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
    config: AnalysisConfig | None = None
    outputs: OutputArtifacts | None = None

    def save(self, output_dir: PathLike, overwrite: bool = False) -> OutputArtifacts:
        """Write the output folder for this result: tables, figures, README.md, run.json.

        Uses the run's output settings (plot format, residue plots, distance matrices)
        with ``output_dir``. Like a run, it refuses a folder that already has files unless
        ``overwrite`` is true. Returns the paths written, also kept in ``outputs``.
        """
        from flexgeo2.outputs import OutputWriter, check_output_dir

        config = self.config or AnalysisConfig(pdb_file=self.pdb_file)
        output = replace(config.output, output_dir=output_dir, overwrite=overwrite)
        check_output_dir(output)
        # run.json and README.md describe the configuration the files were written with.
        self.config = replace(config, output=output)
        self.outputs = OutputWriter(output).write(
            self,
            max_models_in_plot=config.max_models_in_plot,
            hide_model_traces=config.hide_model_traces,
        )
        return self.outputs
