from __future__ import annotations

from pathlib import Path

from flexgeo2.config import OutputConfig
from flexgeo2.distances import DistanceService
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
        DistanceService.to_matrix(distance_long_df).to_csv(output_path)

    def write(self, result: AnalysisResult, max_models_in_plot: int, hide_model_traces: bool):
        if not self.config.write_files:
            return OutputArtifacts()

        if self.config.output_dir is None:
            raise ValueError("OutputConfig.output_dir must be set when write_files=True.")

        verbose = self.config.verbose
        output_dir = Path(self.config.output_dir).resolve()
        geometry_dir = output_dir / "geometry"
        reference_dir = output_dir / "reference" if result.distance_result is not None else None
        clusters_dir = output_dir / "clusters" if result.residue_clustering is not None else None
        range_clusters_dir = (
            output_dir / "range_clusters" if result.residue_range_clustering is not None else None
        )
        chains_dir = output_dir / "chains" if verbose else None
        for directory in (geometry_dir, reference_dir, clusters_dir, range_clusters_dir):
            if directory is not None:
                directory.mkdir(parents=True, exist_ok=True)

        artifacts = OutputArtifacts(
            raw_csv=geometry_dir / "descriptors.csv",
            residue_summary_csv=geometry_dir / "residues.csv",
            model_summary_csv=geometry_dir / "models_by_chain.csv" if verbose else None,
            overall_model_summary_csv=geometry_dir / "models.csv",
            overview_plot=output_dir / "overview.png",
            chains_dir=chains_dir,
            distance_long_csv=(
                reference_dir / "distances.csv" if reference_dir is not None and verbose else None
            ),
            distance_summary_csv=reference_dir / "residues.csv" if reference_dir else None,
            distance_heatmap=reference_dir / "heatmap.png" if reference_dir else None,
            distance_matrix_dir=(
                reference_dir / "matrices" if reference_dir is not None and verbose else None
            ),
            cluster_assignments_csv=(
                clusters_dir / "assignments.csv" if clusters_dir is not None and verbose else None
            ),
            cluster_summary_csv=clusters_dir / "residues.csv" if clusters_dir else None,
            cluster_plots_dir=clusters_dir / "residue_plots" if clusters_dir else None,
            range_cluster_assignments_csv=(
                range_clusters_dir / "assignments.csv"
                if range_clusters_dir is not None and verbose
                else None
            ),
            range_cluster_summary_csv=(
                range_clusters_dir / "ranges.csv" if range_clusters_dir else None
            ),
            range_cluster_plots_dir=range_clusters_dir,
        )

        result.raw_df.to_csv(artifacts.raw_csv, index=False)
        result.residue_summary_df.to_csv(artifacts.residue_summary_csv, index=False)
        result.overall_model_summary_df.to_csv(artifacts.overall_model_summary_csv, index=False)
        self.overview_plotter.plot(
            result.residue_summary_df,
            artifacts.overview_plot,
            raw_df=result.raw_df,
            show_model_traces=not hide_model_traces,
            max_models_in_plot=max_models_in_plot,
        )

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
                    chain_distance_long_df = result.distance_result.long_df[
                        result.distance_result.long_df["chain"] == chain
                    ].copy()
                    self.write_distance_matrix_csv(
                        chain_distance_long_df,
                        artifacts.distance_matrix_dir / f"{sanitize_chain_id(chain)}.csv",
                    )

        if result.residue_clustering is not None:
            artifacts.cluster_plots_dir.mkdir(parents=True, exist_ok=True)
            result.residue_clustering.summary_df.to_csv(artifacts.cluster_summary_csv, index=False)
            if artifacts.cluster_assignments_csv is not None:
                result.residue_clustering.assignments_df.to_csv(
                    artifacts.cluster_assignments_csv, index=False
                )
            self._plot_residue_clusters(
                result.residue_clustering.assignments_df,
                artifacts.cluster_plots_dir,
                with_chain_prefix=True,
            )

        if result.residue_range_clustering is not None:
            result.residue_range_clustering.summary_df.to_csv(
                artifacts.range_cluster_summary_csv, index=False
            )
            if artifacts.range_cluster_assignments_csv is not None:
                result.residue_range_clustering.assignments_df.to_csv(
                    artifacts.range_cluster_assignments_csv, index=False
                )
            self._plot_range_clusters(
                result.residue_range_clustering.assignments_df,
                artifacts.range_cluster_plots_dir,
                with_chain_prefix=True,
            )

        if chains_dir is not None:
            self._write_verbose_chain_outputs(
                result=result,
                chains_dir=chains_dir,
                max_models_in_plot=max_models_in_plot,
                hide_model_traces=hide_model_traces,
            )

        return artifacts

    @staticmethod
    def residue_plot_name(chain, order: int, name: str, with_chain_prefix: bool) -> str:
        """Zero-padded residue number first, so files sort in sequence order."""
        stem = f"{int(order):04d}_{name}"
        return f"{sanitize_chain_id(chain)}_{stem}.png" if with_chain_prefix else f"{stem}.png"

    def _plot_residue_clusters(
        self, assignments_df, plots_dir: Path, with_chain_prefix: bool
    ) -> None:
        for (chain, order, name), residue_cluster_df in assignments_df.groupby(
            ["chain", "order", "name"], dropna=False
        ):
            self.residue_cluster_plotter.plot(
                residue_cluster_df,
                plots_dir / self.residue_plot_name(chain, order, name, with_chain_prefix),
            )

    def _plot_range_clusters(
        self, assignments_df, plots_dir: Path, with_chain_prefix: bool
    ) -> None:
        for (chain, range_label), range_cluster_df in assignments_df.groupby(
            ["chain", "range_label"], dropna=False
        ):
            stem = str(range_label).replace("/", "_")
            if with_chain_prefix:
                stem = f"{sanitize_chain_id(chain)}_{stem}"
            self.residue_range_cluster_plotter.plot(range_cluster_df, plots_dir / f"{stem}.png")

    def _write_verbose_chain_outputs(
        self,
        result: AnalysisResult,
        chains_dir: Path,
        max_models_in_plot: int,
        hide_model_traces: bool,
    ) -> None:
        """Mirror the top-level layout under chains/<chain>/ for each chain."""
        for chain in result.residue_summary_df["chain"].drop_duplicates():
            chain_output_dir = chains_dir / sanitize_chain_id(chain)
            geometry_dir = chain_output_dir / "geometry"
            geometry_dir.mkdir(parents=True, exist_ok=True)

            chain_raw_df = result.raw_df[result.raw_df["chain"] == chain].copy()
            chain_summary_df = result.residue_summary_df[
                result.residue_summary_df["chain"] == chain
            ].copy()
            chain_model_summary_df = result.model_summary_df[
                result.model_summary_df["chain"] == chain
            ].copy()

            chain_raw_df.to_csv(geometry_dir / "descriptors.csv", index=False)
            chain_summary_df.to_csv(geometry_dir / "residues.csv", index=False)
            chain_model_summary_df.to_csv(geometry_dir / "models.csv", index=False)
            self.chain_plotter.plot(
                chain_raw_df=chain_raw_df,
                chain_summary_df=chain_summary_df,
                output_path=chain_output_dir / "overview.png",
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
                    reference_dir = chain_output_dir / "reference"
                    reference_dir.mkdir(parents=True, exist_ok=True)
                    chain_distance_long_df.to_csv(reference_dir / "distances.csv", index=False)
                    chain_distance_summary_df.to_csv(reference_dir / "residues.csv", index=False)
                    self.write_distance_matrix_csv(
                        chain_distance_long_df, reference_dir / "matrix.csv"
                    )
                    self.distance_plotter.plot(
                        chain_distance_long_df,
                        reference_dir / "heatmap.png",
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
                    clusters_dir = chain_output_dir / "clusters"
                    (clusters_dir / "residue_plots").mkdir(parents=True, exist_ok=True)
                    chain_cluster_df.to_csv(clusters_dir / "assignments.csv", index=False)
                    chain_cluster_summary_df.to_csv(clusters_dir / "residues.csv", index=False)
                    self._plot_residue_clusters(
                        chain_cluster_df, clusters_dir / "residue_plots", with_chain_prefix=False
                    )

            if result.residue_range_clustering is not None:
                chain_range_cluster_df = result.residue_range_clustering.assignments_df[
                    result.residue_range_clustering.assignments_df["chain"] == chain
                ].copy()
                chain_range_cluster_summary_df = result.residue_range_clustering.summary_df[
                    result.residue_range_clustering.summary_df["chain"] == chain
                ].copy()
                if not chain_range_cluster_df.empty:
                    range_clusters_dir = chain_output_dir / "range_clusters"
                    range_clusters_dir.mkdir(parents=True, exist_ok=True)
                    chain_range_cluster_df.to_csv(
                        range_clusters_dir / "assignments.csv", index=False
                    )
                    chain_range_cluster_summary_df.to_csv(
                        range_clusters_dir / "ranges.csv", index=False
                    )
                    self._plot_range_clusters(
                        chain_range_cluster_df, range_clusters_dir, with_chain_prefix=False
                    )
