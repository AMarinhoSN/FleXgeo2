from __future__ import annotations

from collections import Counter
from pathlib import Path

from flexgeo2.config import MAX_CHAINS_PER_FIGURE, OutputConfig
from flexgeo2.distances import DistanceService
from flexgeo2.geometry import DMAX_DETAIL_COLUMNS
from flexgeo2.models import AnalysisResult, OutputArtifacts
from flexgeo2.plotting import (
    ClusterMapPlotter,
    DistanceHeatmapPlotter,
    OverviewPlotter,
    ResiduePlotter,
    ResidueRangeClusterPlotter,
    sanitize_chain_id,
)
from flexgeo2.report import write_report, written_files
from flexgeo2.selection import parse_residue_selections, select_residues

# Folders FleXgeo2 creates inside output_dir, deepest first so they can be removed in order.
OUTPUT_SUBDIRS = (
    "reference/matrices",
    "geometry",
    "reference",
    "clusters",
    "range_clusters",
    "residue_plots",
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
        residue_range_cluster_plotter: ResidueRangeClusterPlotter | None = None,
        cluster_map_plotter: ClusterMapPlotter | None = None,
        residue_plotter: ResiduePlotter | None = None,
    ) -> None:
        self.config = config
        self.overview_plotter = overview_plotter or OverviewPlotter()
        self.distance_plotter = distance_plotter or DistanceHeatmapPlotter()
        self.residue_range_cluster_plotter = (
            residue_range_cluster_plotter or ResidueRangeClusterPlotter()
        )
        self.cluster_map_plotter = cluster_map_plotter or ClusterMapPlotter()
        self.residue_plotter = residue_plotter or ResiduePlotter()

    @staticmethod
    def write_distance_matrix_csv(distance_long_df, output_path: str | Path) -> None:
        DistanceService.to_matrix(distance_long_df).to_csv(output_path)

    def write(self, result: AnalysisResult, max_models_in_plot: int, hide_model_traces: bool):
        if not self.config.write_files:
            return OutputArtifacts()

        check_output_dir(self.config)
        output_dir = Path(self.config.output_dir).resolve()
        ext = self.config.plot_format
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
        chains = list(result.residue_summary_df["chain"].drop_duplicates())
        per_chain = len(chains) > MAX_CHAINS_PER_FIGURE
        names = chain_file_names(chains)

        artifacts = OutputArtifacts(
            raw_csv=geometry_dir / "descriptors.csv",
            residue_summary_csv=geometry_dir / "residues.csv",
            # Only a separate answer when there is more than one chain.
            model_summary_csv=(
                geometry_dir / "models_by_chain.csv"
                if result.raw_df["chain"].nunique(dropna=False) > 1
                else None
            ),
            overall_model_summary_csv=geometry_dir / "models.csv",
            overview_plot=None if per_chain else output_dir / f"overview.{ext}",
            distance_long_csv=reference_dir / "distances.csv" if reference_dir else None,
            distance_summary_csv=reference_dir / "residues.csv" if reference_dir else None,
            distance_heatmap=(
                reference_dir / f"heatmap.{ext}" if reference_dir and not per_chain else None
            ),
            distance_matrix_dir=(
                reference_dir / "matrices"
                if reference_dir is not None and self.config.distance_matrices
                else None
            ),
            cluster_assignments_csv=clusters_dir / "assignments.csv" if clusters_dir else None,
            cluster_summary_csv=clusters_dir / "residues.csv" if clusters_dir else None,
            cluster_map_plot=(
                clusters_dir / f"clusters.{ext}" if clusters_dir and not per_chain else None
            ),
            range_cluster_assignments_csv=(
                range_clusters_dir / "assignments.csv" if range_clusters_dir else None
            ),
            range_cluster_summary_csv=(
                range_clusters_dir / "ranges.csv" if range_clusters_dir else None
            ),
            range_cluster_plots_dir=range_clusters_dir,
            residue_plots_dir=output_dir / "residue_plots" if self.config.plot_residues else None,
        )

        result.raw_df.to_csv(artifacts.raw_csv, index=False)
        residue_columns = [
            column
            for column in result.residue_summary_df.columns
            if column not in DMAX_DETAIL_COLUMNS
        ]
        result.residue_summary_df[residue_columns].to_csv(
            artifacts.residue_summary_csv, index=False
        )
        result.overall_model_summary_df.to_csv(artifacts.overall_model_summary_csv, index=False)
        for figure_chains, path in self._figures(output_dir, "overview", names, per_chain):
            self.overview_plotter.plot(
                _in_chains(result.residue_summary_df, figure_chains),
                path,
                raw_df=result.raw_df,
                show_model_traces=not hide_model_traces,
                max_models_in_plot=max_models_in_plot,
                cluster_summary_df=(
                    result.residue_clustering.summary_df
                    if result.residue_clustering is not None
                    else None
                ),
                distance_summary_df=(
                    result.distance_result.summary_df
                    if result.distance_result is not None
                    else None
                ),
            )
            if per_chain:
                artifacts.per_chain_plots.append(path)

        if artifacts.model_summary_csv is not None:
            result.model_summary_df.to_csv(artifacts.model_summary_csv, index=False)

        if result.distance_result is not None:
            result.distance_result.long_df.to_csv(artifacts.distance_long_csv, index=False)
            result.distance_result.summary_df.to_csv(artifacts.distance_summary_csv, index=False)
            long_df = result.distance_result.long_df
            compared = {chain: names[chain] for chain in long_df["chain"].drop_duplicates()}
            for figure_chains, path in self._figures(reference_dir, "heatmap", compared, per_chain):
                self.distance_plotter.plot(
                    _in_chains(long_df, figure_chains),
                    path,
                    f"Distance to reference: {result.distance_result.reference_label}",
                )
                if per_chain:
                    artifacts.per_chain_plots.append(path)
            if artifacts.distance_matrix_dir is not None:
                artifacts.distance_matrix_dir.mkdir(parents=True, exist_ok=True)
                for chain in result.distance_result.long_df["chain"].drop_duplicates():
                    chain_distance_long_df = result.distance_result.long_df[
                        result.distance_result.long_df["chain"] == chain
                    ].copy()
                    self.write_distance_matrix_csv(
                        chain_distance_long_df,
                        artifacts.distance_matrix_dir / f"{names[chain]}.csv",
                    )

        if result.residue_clustering is not None:
            result.residue_clustering.assignments_df.to_csv(
                artifacts.cluster_assignments_csv, index=False
            )
            result.residue_clustering.summary_df.to_csv(artifacts.cluster_summary_csv, index=False)
            clustering = result.residue_clustering
            clustered = {chain: names[chain] for chain in clustering.summary_df["chain"].unique()}
            for figure_chains, path in self._figures(
                clusters_dir, "clusters", clustered, per_chain
            ):
                self.cluster_map_plotter.plot(
                    _in_chains(clustering.summary_df, figure_chains),
                    path,
                    _in_chains(clustering.assignments_df, figure_chains),
                )
                if per_chain:
                    artifacts.per_chain_plots.append(path)

        if result.residue_range_clustering is not None:
            result.residue_range_clustering.assignments_df.to_csv(
                artifacts.range_cluster_assignments_csv, index=False
            )
            result.residue_range_clustering.summary_df.to_csv(
                artifacts.range_cluster_summary_csv, index=False
            )
            self._plot_range_clusters(
                result.residue_range_clustering.assignments_df,
                artifacts.range_cluster_plots_dir,
                names,
            )

        if artifacts.residue_plots_dir is not None:
            artifacts.residue_plots_dir.mkdir(parents=True, exist_ok=True)
            self._plot_chosen_residues(result, artifacts.residue_plots_dir, names)

        # Last, so the guide and manifest describe the files that now exist.
        artifacts.readme, artifacts.run_manifest = write_report(result, output_dir)
        return artifacts

    def _figures(self, directory: Path, stem: str, names: dict, per_chain: bool):
        """``(chains, path)`` of each figure: all chains in one, or one file per chain.

        ``names`` maps each chain in the figure to its file-name part.
        """
        ext = self.config.plot_format
        if not per_chain:
            return [(list(names), directory / f"{stem}.{ext}")]
        return [([chain], directory / f"{stem}_{name}.{ext}") for chain, name in names.items()]

    @staticmethod
    def residue_plot_name(chain_name: str, order: int, name: str, ext: str = "png") -> str:
        """Zero-padded residue number first, so files sort in sequence order.

        ``chain_name`` is the chain's file-name part (see ``chain_file_names``).
        """
        return f"{chain_name}_{int(order):04d}_{name}.{ext}"

    def _plot_chosen_residues(self, result: AnalysisResult, plots_dir: Path, names: dict) -> None:
        summary = result.residue_summary_df
        residues_by_chain = {
            chain: set(group["order"]) for chain, group in summary.groupby("chain", dropna=False)
        }
        chosen = select_residues(
            parse_residue_selections(self.config.plot_residues), residues_by_chain
        )
        # Cluster assignments carry the same curvature and torsion plus the cluster label.
        points = (
            result.residue_clustering.assignments_df
            if result.residue_clustering is not None
            else result.raw_df
        )
        reference = result.distance_result.long_df if result.distance_result else None

        for chain, order in chosen:
            residue = summary[(summary["chain"] == chain) & (summary["order"] == order)].iloc[0]
            reference_point = None
            if reference is not None:
                rows = reference[(reference["chain"] == chain) & (reference["order"] == order)]
                if not rows.empty:
                    reference_point = (
                        rows["reference_curvature"].iloc[0],
                        rows["reference_torsion"].iloc[0],
                    )
            self.residue_plotter.plot(
                points[(points["chain"] == chain) & (points["order"] == order)],
                plots_dir
                / self.residue_plot_name(
                    names[chain], order, residue["name"], self.config.plot_format
                ),
                dmax=residue["dmax"],
                reference=reference_point,
            )

    def _plot_range_clusters(self, assignments_df, plots_dir: Path, names: dict) -> None:
        for (chain, range_label), range_cluster_df in assignments_df.groupby(
            ["chain", "range_label"], dropna=False
        ):
            stem = f"{names[chain]}_{str(range_label).replace('/', '_')}"
            self.residue_range_cluster_plotter.plot(
                range_cluster_df, plots_dir / f"{stem}.{self.config.plot_format}"
            )


def chain_file_names(chains) -> dict:
    """File-name part for each chain ID, unique even on case-insensitive file systems.

    IDs are used as they are ("/" replaced, blank as "unassigned") unless two differ only
    by case, as in large assemblies with chains A and a (macOS and Windows would treat
    overview_A.png and overview_a.png as the same file). Then the IDs that are not all
    upper case get a "_lower" suffix: overview_A.png and overview_a_lower.png. Names that
    still clash (identical after "/" is replaced, or multi-letter IDs) get a number.
    """
    names = {chain: sanitize_chain_id(chain) for chain in chains}
    # Distinct spellings per case-folded name; more than one means a case clash.
    spellings = Counter(name.casefold() for name in set(names.values()))
    for chain, name in names.items():
        if spellings[name.casefold()] > 1 and name != name.upper():
            names[chain] = f"{name}_lower"
    taken: set[str] = set()
    for chain, name in names.items():
        unique, number = name, 1
        while unique.casefold() in taken:
            number += 1
            unique = f"{name}_{number}"
        names[chain] = unique
        taken.add(unique.casefold())
    return names


def _in_chains(df, chains: list):
    """Rows of ``df`` whose chain is one of ``chains``."""
    return df[df["chain"].isin(chains)]
