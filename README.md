# FleXgeo2

This is the first `FleXgeo2` prototype, continuing your earlier software with a new CLI built around [Melodia_py](https://github.com/rwmontalvao/Melodia_py) to compute differential geometry descriptors from PDB files and generate ensemble-aware curvature and torsion outputs.

`FleXgeo2` can now be used both as a CLI and as a Python library. The high-level library entrypoint is `FlexGeo2App`, and advanced users can also import service classes such as `GeometryService`, `DistanceService`, and `ClusteringService`.

## What it does

- reads a PDB file
- computes Melodia descriptors
- writes raw and summarised descriptor tables to CSV
- plots curvature and torsion by residue
- for multi-model PDBs, plots the ensemble mean with a shaded standard deviation band
- overlays individual model traces to help compare conformers
- computes per-residue Euclidean distances in `(curvature, torsion)` space to a reference state
- exports residue-distance matrices and heatmaps
- can cluster conformations residue-by-residue with HDBSCAN in `(curvature, torsion)` space
- writes one clustering scatter plot per residue
- can cluster conformations using a whole residue range as one combined geometric signature
- writes separate outputs for each chain

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Usage

```bash
flexgeo2 path/to/structure.pdb
```

Library usage:

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App

config = AnalysisConfig(pdb_file="ensemble.pdb")
result = FlexGeo2App().run(config)
```

Service-level usage:

```python
from flexgeo2.geometry import GeometryService
from flexgeo2.clustering import ClusteringService

geometry = GeometryService()
raw_df = geometry.load_structure("ensemble.pdb", n_jobs=4)
raw_df = geometry.filter_chains(raw_df, ["A"])
raw_df = geometry.normalize(raw_df)
summary_df = geometry.summarize(raw_df)

clusters = ClusteringService()
assignments_df, cluster_summary_df = clusters.cluster_residues(
    raw_df,
    min_cluster_size=5,
    min_samples=None,
)
```

Lower-level writer/plotter usage:

```python
from flexgeo2.config import OutputConfig
from flexgeo2.outputs import OutputWriter
from flexgeo2.plotting import DistanceHeatmapPlotter

plotter = DistanceHeatmapPlotter()
plotter.plot(distance_long_df, "distance_heatmap.png", title="Reference comparison")

writer = OutputWriter(OutputConfig(output_dir="results", verbose=False))
artifacts = writer.write(
    result,
    max_models_in_plot=12,
    hide_model_traces=False,
)
```

By default, `FleXgeo2` writes a lean set of top-level summary files and plots. Use `--output-verbose` when you want detailed intermediate tables, matrix exports, and duplicated per-chain output folders.

Optional outputs:

```bash
flexgeo2 path/to/ensemble.pdb \
  --output-dir results \
  --chain A \
  --n-jobs -1
```

To write the full detailed output set:

```bash
flexgeo2 path/to/ensemble.pdb --output-verbose
```

To hide individual model overlays:

```bash
flexgeo2 path/to/ensemble.pdb --hide-model-traces
```

To compare the ensemble against a reference model already present in the input:

```bash
flexgeo2 path/to/ensemble.pdb --reference-model 1
```

To compare the ensemble against a state from another PDB:

```bash
flexgeo2 path/to/ensemble.pdb \
  --reference-pdb path/to/reference_state.pdb \
  --reference-pdb-model 1
```

To cluster conformations independently for each residue:

```bash
flexgeo2 path/to/ensemble.pdb --cluster-residues
```

You can tune HDBSCAN if needed:

```bash
flexgeo2 path/to/ensemble.pdb \
  --cluster-residues \
  --cluster-min-size 8 \
  --cluster-min-samples 4
```

To cluster conformations using a biologically interesting residue window:

```bash
flexgeo2 path/to/ensemble.pdb --cluster-residue-range 45-54
```

You can repeat the option to analyze multiple windows:

```bash
flexgeo2 path/to/ensemble.pdb \
  --chain A \
  --cluster-residue-range 45-54 \
  --cluster-residue-range 90-99
```

Default outputs are written to the chosen directory:

- `geometry_descriptors.csv`
- `residue_summary.csv`
- `model_summary_overall.csv`
- `distance_to_reference_summary.csv`
- `residue_cluster_summary.csv`
- `cluster_plots/<chain>_<residue>_clusters.png`
- `residue_range_cluster_summary.csv`
- `range_cluster_plots/<chain>_<start-end>_clusters.png`
- `plots/ensemble_overview.png`
- `plots/distance_to_reference_heatmap.png`

Verbose mode adds:

- `model_summary_by_chain.csv`
- `distance_to_reference_long.csv`
- `distance_matrices/<chain>_distance_matrix.csv`
- `residue_cluster_assignments.csv`
- `residue_range_cluster_assignments.csv`
- `chains/<chain>/geometry_descriptors.csv`
- `chains/<chain>/residue_summary.csv`
- `chains/<chain>/model_summary.csv`
- `chains/<chain>/curvature_torsion.png`
- `chains/<chain>/distance_to_reference_long.csv`
- `chains/<chain>/distance_to_reference_summary.csv`
- `chains/<chain>/distance_to_reference_matrix.csv`
- `chains/<chain>/distance_to_reference_heatmap.png`
- `chains/<chain>/residue_cluster_assignments.csv`
- `chains/<chain>/residue_cluster_summary.csv`
- `chains/<chain>/cluster_plots/<residue>_clusters.png`
- `chains/<chain>/residue_range_cluster_assignments.csv`
- `chains/<chain>/residue_range_cluster_summary.csv`
- `chains/<chain>/range_cluster_plots/<start-end>_clusters.png`

## Notes

- The prototype expects Melodia to return columns including `model`, `chain`, `order`, `name`, `curvature`, and `torsion`.
- Residue positions are plotted using Melodia's `order` column.
- The per-chain model summary includes mean absolute deviation from the ensemble mean, which is useful as a first-pass conformational variability signal.
- Distance matrices use rows for ensemble models and columns for residues, with each cell storing the Euclidean distance to the chosen reference in `(curvature, torsion)` space.
- Residue clustering treats each residue independently and clusters the ensemble conformations using only that residue's curvature and torsion values.
- Residue-range clustering concatenates curvature and torsion values across the selected window into one feature vector per conformation, then runs one HDBSCAN solution for that whole region.
