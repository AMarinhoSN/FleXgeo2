from __future__ import annotations

from pathlib import Path

from flexgeo2.config import OutputConfig
from flexgeo2.distances import DistanceService
from flexgeo2.models import AnalysisResult, OutputArtifacts
from flexgeo2.plotting import (
    ClusterMapPlotter,
    DistanceHeatmapPlotter,
    OverviewPlotter,
    ResidueClusterPlotter,
    ResidueRangeClusterPlotter,
    sanitize_chain_id,
)
from flexgeo2.report import write_report, written_files

# Folders FleXgeo2 creates inside output_dir, deepest first so they can be removed in order.
OUTPUT_SUBDIRS = (
    "reference/matrices",
    "clusters/residue_plots",
    "geometry",
    "reference",
    "clusters",
    "range_clusters",
)


class OutputDirectoryNotEmptyError(FileExistsError):
    """The output folder already has files and overwriting was not requested."""

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir
        super().__init__(
            f"Output folder {output_dir} is not empty. Set OutputConfig(overwrite=True) to "
            "replace the outputs of an earlier run, or choose another output_dir."
        )


def check_output_dir(config: OutputConfig) -> None:
    """Fail if writing would mix new outputs with files already in the output folder.

    Hidden files (e.g. ``.DS_Store``) do not count. With ``overwrite=True`` a folder
    holding an earlier run is accepted; ``remove_previous_outputs`` clears it later.
    """
    if not config.write_files:
        return
    if config.output_dir is None:
        raise ValueError("OutputConfig.output_dir must be set when write_files=True.")

    output_dir = Path(config.output_dir).resolve()
    if not output_dir.exists():
        return
    if not output_dir.is_dir():
        raise ValueError(f"Output path {output_dir} exists and is not a folder.")
    if not any(not entry.name.startswith(".") for entry in output_dir.iterdir()):
        return
    if not config.overwrite:
        raise OutputDirectoryNotEmptyError(output_dir)
    # Every run writes README.md and run.json together, so a README.md without run.json
    # belongs to someone else (e.g. --output-dir pointing at a project folder).
    if (output_dir / "README.md").exists() and not (output_dir / "run.json").exists():
        raise ValueError(
            f"{output_dir} has a README.md that FleXgeo2 did not write; overwriting would "
            "replace it. Choose another output folder."
        )


def remove_previous_outputs(output_dir: Path) -> None:
    """Delete the files an earlier run wrote, so none are left stale; keep all others."""
    for name in written_files(output_dir):
        (output_dir / name).unlink(missing_ok=True)
    for name in OUTPUT_SUBDIRS:
        directory = output_dir / name
        if directory.is_dir() and not any(directory.iterdir()):
            directory.rmdir()


class OutputWriter:
    """Write CSVs and plots for an analysis result."""

    def __init__(
        self,
        config: OutputConfig,
        overview_plotter: OverviewPlotter | None = None,
        distance_plotter: DistanceHeatmapPlotter | None = None,
        residue_cluster_plotter: ResidueClusterPlotter | None = None,
        residue_range_cluster_plotter: ResidueRangeClusterPlotter | None = None,
        cluster_map_plotter: ClusterMapPlotter | None = None,
    ) -> None:
        self.config = config
        self.overview_plotter = overview_plotter or OverviewPlotter()
        self.distance_plotter = distance_plotter or DistanceHeatmapPlotter()
        self.residue_cluster_plotter = residue_cluster_plotter or ResidueClusterPlotter()
        self.residue_range_cluster_plotter = (
            residue_range_cluster_plotter or ResidueRangeClusterPlotter()
        )
        self.cluster_map_plotter = cluster_map_plotter or ClusterMapPlotter()

    @staticmethod
    def write_distance_matrix_csv(distance_long_df, output_path: str | Path) -> None:
        DistanceService.to_matrix(distance_long_df).to_csv(output_path)

    def write(self, result: AnalysisResult, max_models_in_plot: int, hide_model_traces: bool):
        if not self.config.write_files:
            return OutputArtifacts()

        check_output_dir(self.config)
        verbose = self.config.verbose
        output_dir = Path(self.config.output_dir).resolve()
        if self.config.overwrite:
            remove_previous_outputs(output_dir)
        geometry_dir = output_dir / "geometry"
        reference_dir = output_dir / "reference" if result.distance_result is not None else None
        clusters_dir = output_dir / "clusters" if result.residue_clustering is not None else None
        range_clusters_dir = (
            output_dir / "range_clusters" if result.residue_range_clustering is not None else None
        )
        for directory in (geometry_dir, reference_dir, clusters_dir, range_clusters_dir):
            if directory is not None:
                directory.mkdir(parents=True, exist_ok=True)

        artifacts = OutputArtifacts(
            raw_csv=geometry_dir / "descriptors.csv",
            residue_summary_csv=geometry_dir / "residues.csv",
            model_summary_csv=geometry_dir / "models_by_chain.csv" if verbose else None,
            overall_model_summary_csv=geometry_dir / "models.csv",
            overview_plot=output_dir / "overview.png",
            distance_long_csv=reference_dir / "distances.csv" if reference_dir else None,
            distance_summary_csv=reference_dir / "residues.csv" if reference_dir else None,
            distance_heatmap=reference_dir / "heatmap.png" if reference_dir else None,
            distance_matrix_dir=(
                reference_dir / "matrices" if reference_dir is not None and verbose else None
            ),
            cluster_assignments_csv=clusters_dir / "assignments.csv" if clusters_dir else None,
            cluster_summary_csv=clusters_dir / "residues.csv" if clusters_dir else None,
            cluster_map_plot=clusters_dir / "clusters.png" if clusters_dir else None,
            cluster_plots_dir=(
                clusters_dir / "residue_plots" if clusters_dir is not None and verbose else None
            ),
            range_cluster_assignments_csv=(
                range_clusters_dir / "assignments.csv" if range_clusters_dir else None
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
            cluster_summary_df=(
                result.residue_clustering.summary_df
                if result.residue_clustering is not None
                else None
            ),
            distance_summary_df=(
                result.distance_result.summary_df if result.distance_result is not None else None
            ),
        )

        if artifacts.model_summary_csv is not None:
            result.model_summary_df.to_csv(artifacts.model_summary_csv, index=False)

        if result.distance_result is not None:
            result.distance_result.long_df.to_csv(artifacts.distance_long_csv, index=False)
            result.distance_result.summary_df.to_csv(artifacts.distance_summary_csv, index=False)
            self.distance_plotter.plot(
                result.distance_result.long_df,
                artifacts.distance_heatmap,
                f"Distance to reference: {result.distance_result.reference_label}",
            )
            if artifacts.distance_matrix_dir is not None:
                artifacts.distance_matrix_dir.mkdir(parents=True, exist_ok=True)
                for chain in result.distance_result.long_df["chain"].drop_duplicates():
                    chain_distance_long_df = result.distance_result.long_df[
                        result.distance_result.long_df["chain"] == chain
                    ].copy()
                    self.write_distance_matrix_csv(
                        chain_distance_long_df,
                        artifacts.distance_matrix_dir / f"{sanitize_chain_id(chain)}.csv",
                    )

        if result.residue_clustering is not None:
            result.residue_clustering.assignments_df.to_csv(
                artifacts.cluster_assignments_csv, index=False
            )
            result.residue_clustering.summary_df.to_csv(artifacts.cluster_summary_csv, index=False)
            self.cluster_map_plotter.plot(
                result.residue_clustering.summary_df,
                artifacts.cluster_map_plot,
                result.residue_clustering.assignments_df,
            )
            if artifacts.cluster_plots_dir is not None:
                artifacts.cluster_plots_dir.mkdir(parents=True, exist_ok=True)
                self._plot_residue_clusters(
                    result.residue_clustering.assignments_df, artifacts.cluster_plots_dir
                )

        if result.residue_range_clustering is not None:
            result.residue_range_clustering.assignments_df.to_csv(
                artifacts.range_cluster_assignments_csv, index=False
            )
            result.residue_range_clustering.summary_df.to_csv(
                artifacts.range_cluster_summary_csv, index=False
            )
            self._plot_range_clusters(
                result.residue_range_clustering.assignments_df, artifacts.range_cluster_plots_dir
            )

        # Last, so the guide and manifest describe the files that now exist.
        artifacts.readme, artifacts.run_manifest = write_report(result, output_dir)
        return artifacts

    @staticmethod
    def residue_plot_name(chain, order: int, name: str) -> str:
        """Zero-padded residue number first, so files sort in sequence order."""
        return f"{sanitize_chain_id(chain)}_{int(order):04d}_{name}.png"

    def _plot_residue_clusters(self, assignments_df, plots_dir: Path) -> None:
        for (chain, order, name), residue_cluster_df in assignments_df.groupby(
            ["chain", "order", "name"], dropna=False
        ):
            self.residue_cluster_plotter.plot(
                residue_cluster_df, plots_dir / self.residue_plot_name(chain, order, name)
            )

    def _plot_range_clusters(self, assignments_df, plots_dir: Path) -> None:
        for (chain, range_label), range_cluster_df in assignments_df.groupby(
            ["chain", "range_label"], dropna=False
        ):
            stem = f"{sanitize_chain_id(chain)}_{str(range_label).replace('/', '_')}"
            self.residue_range_cluster_plotter.plot(range_cluster_df, plots_dir / f"{stem}.png")
