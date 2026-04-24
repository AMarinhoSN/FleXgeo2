from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PathLike = str | Path


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
    verbose: bool = False
    write_files: bool = True


@dataclass(slots=True)
class AnalysisConfig:
    """Configuration for a full FleXgeo2 analysis run."""

    pdb_file: PathLike
    chains: list[str] | None = None
    n_jobs: int = 1
    max_models_in_plot: int = 12
    hide_model_traces: bool = False
    reference: ReferenceConfig | None = None
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
