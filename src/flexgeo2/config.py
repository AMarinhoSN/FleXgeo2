from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PathLike = str | Path

# File formats for figures; matplotlib picks the writer from the extension.
PLOT_FORMATS = ("png", "pdf", "svg")

# The overview, distance heatmap and cluster map stack one panel per chain. Above this
# many chains they are written as one file per chain instead (e.g. overview_A.png).
MAX_CHAINS_PER_FIGURE = 4


@dataclass(slots=True)
class ReferenceConfig:
    """Reference-state comparison options."""

    model_id: str | None = None
    pdb_file: PathLike | None = None
    pdb_model_id: str | None = None


@dataclass(slots=True)
class ClusteringConfig:
    """Clustering options for per-residue and residue-range analysis."""

    cluster_residues: bool = False
    cluster_residue_ranges: list[str] = field(default_factory=list)
    min_cluster_size: int = 5
    min_samples: int | None = None


@dataclass(slots=True)
class OutputConfig:
    """Output policy for file writing."""

    output_dir: PathLike | None = Path("results")
    write_files: bool = True
    # Also write the distances to the reference as one models x residues CSV per chain.
    distance_matrices: bool = False
    # Residues to plot (curvature vs torsion), e.g. ["45", "A:50-52"]; see selection.py.
    plot_residues: list[str] = field(default_factory=list)
    # File format of every figure: one of PLOT_FORMATS.
    plot_format: str = "png"
    # Allow writing into a non-empty output_dir, replacing an earlier run's outputs.
    overwrite: bool = False


@dataclass(slots=True)
class AnalysisConfig:
    """Configuration for a full FleXgeo2 analysis run."""

    pdb_file: PathLike
    chains: list[str] | None = None
    n_jobs: int = 1
    max_models_in_plot: int = 12
    hide_model_traces: bool = False
    dmax_outlier_fraction: float = 0.01
    reference: ReferenceConfig | None = None
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
