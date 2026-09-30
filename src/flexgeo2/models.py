from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING

from flexgeo2.config import AnalysisConfig, PathLike

if TYPE_CHECKING:
    # Imported for annotations only: pandas takes most of a second to import, and
    # `import flexgeo2` (and `flexgeo2 --help`) should stay fast.
    import pandas as pd


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
