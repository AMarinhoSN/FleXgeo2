from __future__ import annotations

from pathlib import Path

from flexgeo2.config import OutputConfig
from flexgeo2.models import AnalysisResult, OutputArtifacts
from flexgeo2.plotting import (
    ChainGeometryPlotter,
    DistanceHeatmapPlotter,
    OverviewPlotter,
    ResidueClusterPlotter,
    ResidueRangeClusterPlotter,
    sanitize_chain_id,
)


class OutputWriter:
    """Write CSVs and plots for an analysis result."""

    def __init__(
        self,
        config: OutputConfig,
        overview_plotter: OverviewPlotter | None = None,
        chain_plotter: ChainGeometryPlotter | None = None,
        distance_plotter: DistanceHeatmapPlotter | None = None,
        residue_cluster_plotter: ResidueClusterPlotter | None = None,
        residue_range_cluster_plotter: ResidueRangeClusterPlotter | None = None,
    ) -> None:
        self.config = config
        self.overview_plotter = overview_plotter or OverviewPlotter()
        self.chain_plotter = chain_plotter or ChainGeometryPlotter()
        self.distance_plotter = distance_plotter or DistanceHeatmapPlotter()
        self.residue_cluster_plotter = residue_cluster_plotter or ResidueClusterPlotter()
        self.residue_range_cluster_plotter = (
            residue_range_cluster_plotter or ResidueRangeClusterPlotter()
        )

    @staticmethod
    def write_distance_matrix_csv(distance_long_df, output_path: str | Path) -> None:
        matrix_df = (
            distance_long_df.pivot(
                index="model", columns="residue_label", values="distance_to_reference"
            ).sort_index()
        )
        matrix_df.to_csv(output_path)

    def write(self, result: AnalysisResult, max_models_in_plot: int, hide_model_traces: bool):
        if not self.config.write_files:
            return OutputArtifacts()

        if self.config.output_dir is None:
            raise ValueError("OutputConfig.output_dir must be set when write_files=True.")

        output_dir = Path(self.config.output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        chains_dir = output_dir / "chains" if self.config.verbose else None
        if chains_dir is not None:
            chains_dir.mkdir(parents=True, exist_ok=True)

        artifacts = OutputArtifacts(
            raw_csv=output_dir / "geometry_descriptors.csv",
            residue_summary_csv=output_dir / "residue_summary.csv",
            model_summary_csv=(
                output_dir / "model_summary_by_chain.csv" if self.config.verbose else None
            ),
            overall_model_summary_csv=output_dir / "model_summary_overall.csv",
            overview_plot=plots_dir / "ensemble_overview.png",
            chains_dir=chains_dir,
            distance_long_csv=(
                output_dir / "distance_to_reference_long.csv"
                if result.distance_result is not None and self.config.verbose
                else None
            ),
            distance_summary_csv=(
                output_dir / "distance_to_reference_summary.csv"
                if result.distance_result is not None
                else None
            ),
            distance_heatmap=(
                plots_dir / "distance_to_reference_heatmap.png"
                if result.distance_result is not None
                else None
            ),
            distance_matrix_dir=(
                output_dir / "distance_matrices"
                if result.distance_result is not None and self.config.verbose
                else None
            ),
            cluster_assignments_csv=(
                output_dir / "residue_cluster_assignments.csv"
                if result.residue_clustering is not None and self.config.verbose
                else None
            ),
            cluster_summary_csv=(
                output_dir / "residue_cluster_summary.csv"
                if result.residue_clustering is not None
                else None
            ),
            cluster_plots_dir=(
                output_dir / "cluster_plots"
                if result.residue_clustering is not None
                else None
            ),
            range_cluster_assignments_csv=(
                output_dir / "residue_range_cluster_assignments.csv"
                if result.residue_range_clustering is not None and self.config.verbose
                else None
            ),
            range_cluster_summary_csv=(
                output_dir / "residue_range_cluster_summary.csv"
                if result.residue_range_clustering is not None
                else None
            ),
            range_cluster_plots_dir=(
                output_dir / "range_cluster_plots"
                if result.residue_range_clustering is not None
                else None
            ),
        )

        result.raw_df.to_csv(artifacts.raw_csv, index=False)
        result.residue_summary_df.to_csv(artifacts.residue_summary_csv, index=False)
        result.overall_model_summary_df.to_csv(artifacts.overall_model_summary_csv, index=False)
        self.overview_plotter.plot(result.residue_summary_df, artifacts.overview_plot)

        if artifacts.model_summary_csv is not None:
            result.model_summary_df.to_csv(artifacts.model_summary_csv, index=False)

        if result.distance_result is not None:
            result.distance_result.summary_df.to_csv(artifacts.distance_summary_csv, index=False)
            self.distance_plotter.plot(
                result.distance_result.long_df,
                artifacts.distance_heatmap,
                f"Distance to reference: {result.distance_result.reference_label}",
            )
            if artifacts.distance_matrix_dir is not None:
                artifacts.distance_matrix_dir.mkdir(parents=True, exist_ok=True)
                result.distance_result.long_df.to_csv(artifacts.distance_long_csv, index=False)
                for chain in result.distance_result.long_df["chain"].drop_duplicates():
                    chain_key = sanitize_chain_id(chain)
                    chain_distance_long_df = result.distance_result.long_df[
                        result.distance_result.long_df["chain"] == chain
                    ].copy()
                    self.write_distance_matrix_csv(
                        chain_distance_long_df,
                        artifacts.distance_matrix_dir / f"{chain_key}_distance_matrix.csv",
                    )

        if result.residue_clustering is not None:
            artifacts.cluster_plots_dir.mkdir(parents=True, exist_ok=True)
            result.residue_clustering.summary_df.to_csv(artifacts.cluster_summary_csv, index=False)
            if artifacts.cluster_assignments_csv is not None:
                result.residue_clustering.assignments_df.to_csv(
                    artifacts.cluster_assignments_csv, index=False
                )
            for (chain, residue_label), residue_cluster_df in result.residue_clustering.assignments_df.groupby(
                ["chain", "residue_label"], dropna=False
            ):
                chain_key = sanitize_chain_id(chain)
                self.residue_cluster_plotter.plot(
                    residue_cluster_df,
                    artifacts.cluster_plots_dir / f"{chain_key}_{residue_label}_clusters.png",
                )

        if result.residue_range_clustering is not None:
            artifacts.range_cluster_plots_dir.mkdir(parents=True, exist_ok=True)
            result.residue_range_clustering.summary_df.to_csv(
                artifacts.range_cluster_summary_csv, index=False
            )
            if artifacts.range_cluster_assignments_csv is not None:
                result.residue_range_clustering.assignments_df.to_csv(
                    artifacts.range_cluster_assignments_csv, index=False
                )
            for (chain, range_label), range_cluster_df in result.residue_range_clustering.assignments_df.groupby(
                ["chain", "range_label"], dropna=False
            ):
                chain_key = sanitize_chain_id(chain)
                safe_range_label = str(range_label).replace("/", "_")
                self.residue_range_cluster_plotter.plot(
                    range_cluster_df,
                    artifacts.range_cluster_plots_dir
                    / f"{chain_key}_{safe_range_label}_clusters.png",
                )

        if chains_dir is not None:
            self._write_verbose_chain_outputs(
                result=result,
                chains_dir=chains_dir,
                max_models_in_plot=max_models_in_plot,
                hide_model_traces=hide_model_traces,
            )

        return artifacts

    def _write_verbose_chain_outputs(
        self,
        result: AnalysisResult,
        chains_dir: Path,
        max_models_in_plot: int,
        hide_model_traces: bool,
    ) -> None:
        for chain in result.residue_summary_df["chain"].drop_duplicates():
            chain_key = sanitize_chain_id(chain)
            chain_output_dir = chains_dir / chain_key
            chain_output_dir.mkdir(parents=True, exist_ok=True)

            chain_raw_df = result.raw_df[result.raw_df["chain"] == chain].copy()
            chain_summary_df = result.residue_summary_df[
                result.residue_summary_df["chain"] == chain
            ].copy()
            chain_model_summary_df = result.model_summary_df[
                result.model_summary_df["chain"] == chain
            ].copy()

            chain_raw_df.to_csv(chain_output_dir / "geometry_descriptors.csv", index=False)
            chain_summary_df.to_csv(chain_output_dir / "residue_summary.csv", index=False)
            chain_model_summary_df.to_csv(chain_output_dir / "model_summary.csv", index=False)
            self.chain_plotter.plot(
                chain_raw_df=chain_raw_df,
                chain_summary_df=chain_summary_df,
                output_path=chain_output_dir / "curvature_torsion.png",
                show_model_traces=not hide_model_traces,
                max_models_in_plot=max_models_in_plot,
            )

            if result.distance_result is not None:
                chain_distance_long_df = result.distance_result.long_df[
                    result.distance_result.long_df["chain"] == chain
                ].copy()
                chain_distance_summary_df = result.distance_result.summary_df[
                    result.distance_result.summary_df["chain"] == chain
                ].copy()
                if not chain_distance_long_df.empty:
                    chain_distance_long_df.to_csv(
                        chain_output_dir / "distance_to_reference_long.csv", index=False
                    )
                    chain_distance_summary_df.to_csv(
                        chain_output_dir / "distance_to_reference_summary.csv", index=False
                    )
                    self.write_distance_matrix_csv(
                        chain_distance_long_df,
                        chain_output_dir / "distance_to_reference_matrix.csv",
                    )
                    self.distance_plotter.plot(
                        chain_distance_long_df,
                        chain_output_dir / "distance_to_reference_heatmap.png",
                        f"Chain {chain}: distance to {result.distance_result.reference_label}",
                    )

            if result.residue_clustering is not None:
                chain_cluster_df = result.residue_clustering.assignments_df[
                    result.residue_clustering.assignments_df["chain"] == chain
                ].copy()
                chain_cluster_summary_df = result.residue_clustering.summary_df[
                    result.residue_clustering.summary_df["chain"] == chain
                ].copy()
                if not chain_cluster_df.empty:
                    chain_cluster_plot_dir = chain_output_dir / "cluster_plots"
                    chain_cluster_plot_dir.mkdir(parents=True, exist_ok=True)
                    chain_cluster_df.to_csv(
                        chain_output_dir / "residue_cluster_assignments.csv", index=False
                    )
                    chain_cluster_summary_df.to_csv(
                        chain_output_dir / "residue_cluster_summary.csv", index=False
                    )
                    for residue_label, residue_cluster_df in chain_cluster_df.groupby(
                        "residue_label", dropna=False
                    ):
                        self.residue_cluster_plotter.plot(
                            residue_cluster_df,
                            chain_cluster_plot_dir / f"{residue_label}_clusters.png",
                        )

            if result.residue_range_clustering is not None:
                chain_range_cluster_df = result.residue_range_clustering.assignments_df[
                    result.residue_range_clustering.assignments_df["chain"] == chain
                ].copy()
                chain_range_cluster_summary_df = result.residue_range_clustering.summary_df[
                    result.residue_range_clustering.summary_df["chain"] == chain
                ].copy()
                if not chain_range_cluster_df.empty:
                    chain_range_cluster_plot_dir = chain_output_dir / "range_cluster_plots"
                    chain_range_cluster_plot_dir.mkdir(parents=True, exist_ok=True)
                    chain_range_cluster_df.to_csv(
                        chain_output_dir / "residue_range_cluster_assignments.csv",
                        index=False,
                    )
                    chain_range_cluster_summary_df.to_csv(
                        chain_output_dir / "residue_range_cluster_summary.csv", index=False
                    )
                    for range_label, range_cluster_df in chain_range_cluster_df.groupby(
                        "range_label", dropna=False
                    ):
                        safe_range_label = str(range_label).replace("/", "_")
                        self.residue_range_cluster_plotter.plot(
                            range_cluster_df,
                            chain_range_cluster_plot_dir / f"{safe_range_label}_clusters.png",
                        )
