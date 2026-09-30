from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING

from flexgeo2.config import AnalysisConfig, PathLike

if TYPE_CHECKING:
    # Imported for annotations only: pandas takes most of a second to import, and
    # `import flexgeo2` (and `flexgeo2 --help`) should stay fast.
    import pandas as pd
    from matplotlib.figure import Figure


@dataclass(slots=True)
class DistanceResult:
    """Distance of every model to the reference in (curvature, torsion) space.

    - ``long_df``: one row per model and residue, with the model's and the reference's
      curvature and torsion and ``distance_to_reference`` (``reference/distances.csv``).
    - ``summary_df``: one row per residue: mean, SD, min and max distance
      (``reference/residues.csv``).
    - ``reference_label``: the reference, e.g. ``"input model 1"``.
    """

    long_df: pd.DataFrame
    summary_df: pd.DataFrame
    reference_label: str

    @property
    def heatmap_title(self) -> str:
        return f"Distance to reference: {self.reference_label}"


@dataclass(slots=True)
class ResidueClusteringResult:
    """HDBSCAN clusters of each residue's (curvature, torsion) values across models.

    - ``assignments_df``: one row per model and residue, with ``cluster`` (-1 is noise)
      and ``cluster_probability`` (``clusters/assignments.csv``).
    - ``summary_df``: one row per residue: ``n_clusters`` and ``noise_fraction``
      (``clusters/residues.csv``).
    """

    assignments_df: pd.DataFrame
    summary_df: pd.DataFrame


@dataclass(slots=True)
class ResidueRangeClusteringResult:
    """HDBSCAN clusters of models by the geometry of whole residue ranges.

    - ``assignments_df``: one row per model and range, with ``cluster`` (-1 is noise),
      ``cluster_probability`` and the PCA coordinates ``pc1`` and ``pc2``
      (``range_clusters/assignments.csv``).
    - ``summary_df``: one row per range: ``n_clusters`` and ``noise_fraction``
      (``range_clusters/ranges.csv``).
    """

    assignments_df: pd.DataFrame
    summary_df: pd.DataFrame


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
    """Everything a FleXgeo2 run computed; ``save()`` writes it to an output folder.

    - ``raw_df``: Melodia descriptors, one row per model and residue
      (``geometry/descriptors.csv``).
    - ``residue_summary_df``: one row per residue: curvature and torsion statistics
      across models and ``dmax`` (``geometry/residues.csv``, which leaves out the
      ``*_dmax_*`` columns that show how ``dmax`` was computed).
    - ``model_summary_df``: one row per model and chain (``geometry/models_by_chain.csv``,
      written when the input has more than one chain).
    - ``overall_model_summary_df``: one row per model, all chains together
      (``geometry/models.csv``).
    - ``distance_result``, ``residue_clustering``, ``residue_range_clustering``: results
      of the optional analyses, ``None`` when they did not run.
    - ``config``: the configuration of the run; ``outputs``: the files written, if any.
    """

    pdb_file: Path
    raw_df: pd.DataFrame
    residue_summary_df: pd.DataFrame
    model_summary_df: pd.DataFrame
    overall_model_summary_df: pd.DataFrame
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

    # Figures. Each returns a matplotlib Figure drawn in the FleXgeo2 style; save it with
    # flexgeo2.save_figure() to get the same editable fonts as the files FleXgeo2 writes.

    def plot_overview(self, chains: str | list[str] | None = None) -> Figure:
        """Curvature, torsion, dmax and, when they ran, clusters per residue and distance to
        the reference along the sequence (``overview.png``), for all or the given chains."""
        from flexgeo2.plotting import OverviewPlotter

        config = self._run_config()
        return OverviewPlotter().render(
            self._in_chains(self.residue_summary_df, chains),
            raw_df=self.raw_df,
            show_model_traces=not config.hide_model_traces,
            max_models_in_plot=config.max_models_in_plot,
            **self._overview_extras(),
        )

    def plot_distance_heatmap(self, chains: str | list[str] | None = None) -> Figure:
        """Distance of every model to the reference at every residue
        (``reference/heatmap.png``). Needs a run with a reference."""
        from flexgeo2.plotting import DistanceHeatmapPlotter

        distance = self._needs(self.distance_result, "a reference (reference=ReferenceConfig)")
        return DistanceHeatmapPlotter().render(
            self._in_chains(distance.long_df, chains), distance.heatmap_title
        )

    def plot_cluster_map(self, chains: str | list[str] | None = None) -> Figure:
        """Cluster of every model at every residue (``clusters/clusters.png``). Needs a run
        with per-residue clustering."""
        from flexgeo2.plotting import ClusterMapPlotter

        clustering = self._needs(
            self.residue_clustering, "per-residue clustering (cluster_residues=True)"
        )
        return ClusterMapPlotter().render(
            self._in_chains(clustering.summary_df, chains),
            self._in_chains(clustering.assignments_df, chains),
        )

    def plot_residue(self, residue: str | int) -> Figure:
        """Curvature vs torsion of one residue, one point per model (as ``residue_plots/``).

        ``residue`` is a residue number, or ``"A:45"`` to name the chain; points are
        coloured by cluster when per-residue clustering ran.
        """
        from flexgeo2.plotting import ResiduePlotter
        from flexgeo2.selection import parse_residue_selections, select_residues

        residues_by_chain = {
            chain: set(group["order"])
            for chain, group in self.residue_summary_df.groupby("chain", dropna=False)
        }
        chosen = select_residues(parse_residue_selections([str(residue)]), residues_by_chain)
        if len(chosen) != 1:
            found = ", ".join(f"{chain}:{order}" for chain, order in chosen)
            raise ValueError(
                f"Residue '{residue}' matches {len(chosen)} residues ({found}); give one "
                "residue, with its chain if needed (e.g. 'A:45')."
            )
        points, dmax, reference = self._residue_plot_inputs(*chosen[0])
        return ResiduePlotter().render(points, dmax=dmax, reference=reference)

    def plot_residue_range(self, residues: str) -> Figure:
        """Models of one clustered residue range on its first two principal components
        (``range_clusters/<chain>_<start-end>.png``), e.g. ``"10-20"`` or ``"A:10-20"``.
        Needs a run that clustered that range."""
        from flexgeo2.plotting import ResidueRangeClusterPlotter
        from flexgeo2.selection import parse_residue_selections

        clustering = self._needs(
            self.residue_range_clustering,
            "range clustering (cluster_residue_ranges=[...])",
        )
        assignments = clustering.assignments_df
        selections = parse_residue_selections([residues])
        if len(selections) != 1:
            raise ValueError(
                f"Give one residue range, e.g. '10-20' or 'A:10-20', not '{residues}'."
            )
        [selection] = selections
        matches = assignments[
            (assignments["range_start"] == selection.start)
            & (assignments["range_end"] == selection.end)
            & (selection.chain is None or assignments["chain"] == selection.chain)
        ]
        ranges = matches[["chain", "range_label"]].drop_duplicates()
        if len(ranges) != 1:
            clustered = assignments[["chain", "range_label"]].drop_duplicates()
            listed = ", ".join(f"{row.chain}:{row.range_label}" for row in clustered.itertuples())
            problem = "is not a clustered range" if ranges.empty else "is in several chains"
            raise ValueError(f"Residue range '{residues}' {problem}. Clustered: {listed}.")
        return ResidueRangeClusterPlotter().render(matches)

    def _run_config(self) -> AnalysisConfig:
        return self.config or AnalysisConfig(pdb_file=self.pdb_file)

    @staticmethod
    def _needs(analysis, what: str):
        if analysis is None:
            raise ValueError(f"This figure needs a run with {what}.")
        return analysis

    @staticmethod
    def _in_chains(df: pd.DataFrame, chains: str | list[str] | None) -> pd.DataFrame:
        if chains is None:
            return df
        chains = [chains] if isinstance(chains, str) else list(chains)
        available = list(df["chain"].drop_duplicates())
        missing = [chain for chain in chains if chain not in available]
        if missing:
            raise ValueError(
                f"Chain(s) {', '.join(map(str, missing))} not in this result "
                f"({', '.join(map(str, available))})."
            )
        return df[df["chain"].isin(chains)]

    def _overview_extras(self) -> dict:
        """The optional overview panels: clusters per residue, distance to the reference."""
        return {
            "cluster_summary_df": (
                self.residue_clustering.summary_df if self.residue_clustering else None
            ),
            "distance_summary_df": (
                self.distance_result.summary_df if self.distance_result else None
            ),
        }

    def _residue_plot_inputs(self, chain, order: int):
        """``(points, dmax, reference)`` for the curvature vs torsion plot of a residue."""
        # Cluster assignments carry the same curvature and torsion plus the cluster label.
        points = self.residue_clustering.assignments_df if self.residue_clustering else self.raw_df
        summary = self.residue_summary_df
        dmax = summary[(summary["chain"] == chain) & (summary["order"] == order)]["dmax"].iloc[0]
        reference = None
        if self.distance_result is not None:
            rows = self.distance_result.long_df
            rows = rows[(rows["chain"] == chain) & (rows["order"] == order)]
            if not rows.empty:
                reference = (
                    rows["reference_curvature"].iloc[0],
                    rows["reference_torsion"].iloc[0],
                )
        return points[(points["chain"] == chain) & (points["order"] == order)], dmax, reference
