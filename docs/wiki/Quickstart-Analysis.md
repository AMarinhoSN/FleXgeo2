# Quickstart Analysis

This page shows the smallest useful FleXgeo2 workflows: one command-line run and one Python run. Both compute Melodia descriptors, summarize residue-level geometry, and write the default output files.

## Command-line run

```bash
flexgeo2 path/to/structure.pdb
```

By default, outputs are written to `results/`.

To choose a different output directory:

```bash
flexgeo2 path/to/structure.pdb --output-dir path/to/results
```

To keep one chain:

```bash
flexgeo2 path/to/structure.pdb --chain A
```

## Python run

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App

config = AnalysisConfig(pdb_file="path/to/structure.pdb")
result = FlexGeo2App().run(config)
```

To filter a chain and write to a specific directory:

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App, OutputConfig

config = AnalysisConfig(
    pdb_file="path/to/structure.pdb",
    chains=["A"],
    output=OutputConfig(output_dir="path/to/results"),
)

result = FlexGeo2App().run(config)
```

## What comes back

The returned `result` contains in-memory data frames:

- `result.raw_df`: normalized Melodia descriptors for each model, chain, and residue
- `result.residue_summary_df`: residue-level curvature, torsion, and `dmax` summaries
- `result.model_summary_df`: per-chain model summaries
- `result.overall_model_summary_df`: model summaries across all selected chains
- `result.outputs`: paths to files written by FleXgeo2

## Default output files

A standard run writes:

- `geometry_descriptors.csv`
- `residue_summary.csv`
- `model_summary_overall.csv`
- `plots/ensemble_overview.png`

The overview plot shows curvature, torsion, and `dmax` along the residue order. For multi-model PDB files, FleXgeo2 summarizes the ensemble mean and variation, and can overlay individual model traces.

For a runnable notebook version of this workflow, see `notebooks/01_quickstart_analysis.ipynb`.
