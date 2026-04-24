# Python Library Guide

FleXgeo2 can be used as a Python library when you want to run analyses inside notebooks, inspect data frames in memory, or integrate geometry descriptors into another workflow.

## High-level API

The main entry point is `FlexGeo2App`.

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App

result = FlexGeo2App().run(
    AnalysisConfig(pdb_file="path/to/ensemble.pdb")
)
```

The top-level configuration objects are:

- `AnalysisConfig`: full analysis settings
- `OutputConfig`: output directory and file-writing policy
- `ReferenceConfig`: reference-state comparison settings
- `ClusteringConfig`: per-residue and residue-range clustering settings

## AnalysisConfig

```python
from flexgeo2 import AnalysisConfig, ClusteringConfig, OutputConfig, ReferenceConfig

config = AnalysisConfig(
    pdb_file="path/to/ensemble.pdb",
    chains=["A"],
    n_jobs=4,
    max_models_in_plot=12,
    hide_model_traces=False,
    dmax_outlier_fraction=0.01,
    reference=ReferenceConfig(model_id="1"),
    clustering=ClusteringConfig(
        cluster_residues=True,
        cluster_residue_ranges=["45-54"],
        min_cluster_size=8,
        min_samples=4,
    ),
    output=OutputConfig(output_dir="path/to/results", verbose=True),
)
```

## Running without writing files

Use `OutputConfig(write_files=False)` when you only want the in-memory result data frames.

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App, OutputConfig

result = FlexGeo2App().run(
    AnalysisConfig(
        pdb_file="path/to/ensemble.pdb",
        output=OutputConfig(write_files=False),
    )
)
```

## Result objects

`FlexGeo2App().run(config)` returns an `AnalysisResult` with:

- `pdb_file`: resolved input path
- `raw_df`: normalized Melodia descriptor table
- `residue_summary_df`: per-residue curvature, torsion, and `dmax` summary
- `model_summary_df`: per-chain model summary
- `overall_model_summary_df`: model summary across selected chains
- `distance_result`: optional `DistanceResult`
- `residue_clustering`: optional `ResidueClusteringResult`
- `residue_range_clustering`: optional `ResidueRangeClusteringResult`
- `outputs`: `OutputArtifacts` with paths to written files

`DistanceResult` contains `long_df`, `summary_df`, and `reference_label`.

`ResidueClusteringResult` and `ResidueRangeClusteringResult` each contain `assignments_df` and `summary_df`.

## CLI to Python examples

Basic run:

```bash
flexgeo2 path/to/ensemble.pdb
```

```python
AnalysisConfig(pdb_file="path/to/ensemble.pdb")
```

Chain filtering:

```bash
flexgeo2 path/to/ensemble.pdb --chain A --n-jobs 4
```

```python
AnalysisConfig(
    pdb_file="path/to/ensemble.pdb",
    chains=["A"],
    n_jobs=4,
)
```

Reference model comparison:

```bash
flexgeo2 path/to/ensemble.pdb --reference-model 1
```

```python
AnalysisConfig(
    pdb_file="path/to/ensemble.pdb",
    reference=ReferenceConfig(model_id="1"),
)
```

External reference PDB:

```bash
flexgeo2 path/to/ensemble.pdb \
  --reference-pdb path/to/reference.pdb \
  --reference-pdb-model 1
```

```python
AnalysisConfig(
    pdb_file="path/to/ensemble.pdb",
    reference=ReferenceConfig(
        pdb_file="path/to/reference.pdb",
        pdb_model_id="1",
    ),
)
```

Residue clustering:

```bash
flexgeo2 path/to/ensemble.pdb \
  --cluster-residues \
  --cluster-min-size 8 \
  --cluster-min-samples 4
```

```python
AnalysisConfig(
    pdb_file="path/to/ensemble.pdb",
    clustering=ClusteringConfig(
        cluster_residues=True,
        min_cluster_size=8,
        min_samples=4,
    ),
)
```

Residue-range clustering:

```bash
flexgeo2 path/to/ensemble.pdb \
  --chain A \
  --cluster-residue-range 45-54 \
  --cluster-residue-range 90-99
```

```python
AnalysisConfig(
    pdb_file="path/to/ensemble.pdb",
    chains=["A"],
    clustering=ClusteringConfig(
        cluster_residue_ranges=["45-54", "90-99"],
    ),
)
```

For a notebook dedicated to these mappings, see `notebooks/05_cli_to_python_mapping.ipynb`.

## Service-level API

Advanced users can call services directly.

```python
from flexgeo2.geometry import GeometryService

geometry = GeometryService()
raw_df = geometry.load_structure("path/to/ensemble.pdb", n_jobs=4)
raw_df = geometry.filter_chains(raw_df, ["A"])
raw_df = geometry.normalize(raw_df)
summary_df = geometry.summarize(raw_df, dmax_outlier_fraction=0.01)
```

```python
from flexgeo2.clustering import ClusteringService

clusters = ClusteringService()
assignments_df, cluster_summary_df = clusters.cluster_residues(
    raw_df,
    min_cluster_size=5,
    min_samples=None,
)
```

```python
from flexgeo2.config import OutputConfig
from flexgeo2.outputs import OutputWriter

writer = OutputWriter(OutputConfig(output_dir="path/to/results", verbose=False))
artifacts = writer.write(
    result,
    max_models_in_plot=12,
    hide_model_traces=False,
)
```
