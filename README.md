
![logo](./flexgeo_logo_wide.png)

[![CI](https://img.shields.io/github/actions/workflow/status/AMarinhoSN/FleXgeo2/ci.yml?branch=main&label=CI)](https://github.com/AMarinhoSN/FleXgeo2/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-%3E%3D3.10-blue)
![Version](https://img.shields.io/badge/version-2.0.0.dev0-orange)
![Code style](https://img.shields.io/badge/code%20style-ruff-46aef7)

---

# FleXgeo2

`FleXgeo2` is a complete refactor of [FleXgeo](https://github.com/AMarinhoSN/FleXgeo), this new version is built around [Melodia_py](https://github.com/rwmontalvao/Melodia_py) to compute differential geometry descriptors from PDB files and generate ensemble-aware curvature and torsion outputs.

`FleXgeo2` can now be used both as a CLI and as a Python library. The high-level library entrypoint is `FlexGeo2App`, and advanced users can also import service classes such as `GeometryService`, `DistanceService`, and `ClusteringService`. Check the [full documentation](https://github.com/AMarinhoSN/FleXgeo2/wiki) for more details.

> **FleXgeo2 2.0.0 is in development** (version `2.0.0.dev0`). If you used the April 2026
> preview tagged `v2.0.0`, note that the output folder layout, some options and model
> numbering have changed since: see [CHANGELOG.md](CHANGELOG.md) for what changed and where
> each file went.

---

## What it does

- reads a PDB file
- computes Melodia descriptors
- writes raw and summarised descriptor tables to CSV
- plots curvature and torsion by residue
- for multi-model PDBs, plots the ensemble mean with a shaded standard deviation band
- overlays individual model traces to help compare conformers
- computes per-residue `dmax`, an outlier-trimmed maximum spread in `(curvature, torsion)` space
- computes per-residue Euclidean distances in `(curvature, torsion)` space to a reference state
- plots those distances as a heatmap, and can export them as matrices
- can cluster conformations residue-by-residue with HDBSCAN in `(curvature, torsion)` space
- maps the cluster of every model at every residue in one figure
- can cluster conformations using a whole residue range as one combined geometric signature
- writes a `README.md` into every output folder, explaining the results and every file

 ---

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

---
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

writer = OutputWriter(OutputConfig(output_dir="results"))
artifacts = writer.write(
    result,
    max_models_in_plot=12,
    hide_model_traces=False,
)
```

By default, `FleXgeo2` writes a lean set of summary tables and plots, organised in one folder per analysis (see the output layout below), and prints a short summary. For example, `flexgeo2 pdb2lj5.pdb --reference-model 1` prints:

```
FleXgeo2 analysed pdb2lj5.pdb: 301 models, 1 chain (A), 76 residues.

Most flexible residues (dmax): A LEU71 (2.308), A GLY75 (1.934), A GLY76 (1.934)
Furthest from input model 1 (mean distance): A LEU71 (1.634), A LEU67 (1.116), A GLY10 (0.764)

Results: results/ (9 files); start with README.md
```

Common options:

```bash
flexgeo2 path/to/ensemble.pdb \
  --output-dir results \
  --chain A \
  --n-jobs -1
```

To also write the distances to the reference as one models x residues matrix per chain:

```bash
flexgeo2 path/to/ensemble.pdb --reference-model 1 --distance-matrices
```

To plot curvature against torsion for chosen residues (one point per model, coloured by cluster with `--cluster-residues`, with the reference marked when one is given):

```bash
flexgeo2 path/to/ensemble.pdb --cluster-residues --plot-residues 45,A:50-52
```

A residue number without a chain applies to every analysed chain that has it.

Figures are PNG by default. For publication, write them as vector PDF or SVG, with text that stays editable in Illustrator or Inkscape:

```bash
flexgeo2 path/to/ensemble.pdb --plot-format pdf
```

To hide individual model overlays:

```bash
flexgeo2 path/to/ensemble.pdb --hide-model-traces
```

To tune `dmax` outlier trimming, set the fraction cutoff for sparse extreme histogram bins:

```bash
flexgeo2 path/to/ensemble.pdb --dmax-outlier-fraction 0.02
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

`--reference-pdb-model` is only valid together with `--reference-pdb`.

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

Outputs are organised with one folder per analysis. Folders for optional analyses are
only created when that analysis runs. By default the output directory contains:

```
results/
├── README.md                    # start here: what ran, key results, guide to every file
├── run.json                     # parameters, package versions, input checksums
├── overview.png                 # per-residue results along the sequence
├── geometry/
│   ├── descriptors.csv          # per model and residue: curvature, torsion, ...
│   ├── residues.csv             # per residue: mean, SD, range and dmax
│   └── models.csv               # per model: deviation from the ensemble mean
├── reference/                   # with --reference-model or --reference-pdb
│   ├── distances.csv            # per model and residue: distance to the reference
│   ├── residues.csv             # per residue: distance statistics
│   └── heatmap.png
├── clusters/                    # with --cluster-residues
│   ├── assignments.csv          # per model and residue: cluster label, probability
│   ├── residues.csv             # per residue: number of clusters, noise fraction
│   └── clusters.png             # cluster of every model at every residue
└── range_clusters/              # with --cluster-residue-range
    ├── assignments.csv          # per model and range: cluster label, PCA coordinates
    ├── ranges.csv               # per range: number of clusters, noise fraction
    └── <chain>_<start-end>.png
```

Figures use the extension of `--plot-format` (`png` by default, `pdf` or `svg`).

All tables are in long ("tidy") format with a `chain` column, so a single chain can be
selected by filtering that column.

Some files are written only when they apply:

- `geometry/models_by_chain.csv`: the per-model summary computed separately per chain,
  when the input has more than one chain
- `reference/matrices/<chain>.csv`: distances as a models x residues matrix, with
  `--distance-matrices`
- `residue_plots/<chain>_<residue number>_<name>.png`: curvature vs torsion of each
  residue chosen with `--plot-residues`
- `overview_<chain>.png`, `reference/heatmap_<chain>.png` and
  `clusters/clusters_<chain>.png`: with more than 4 chains, these figures are written
  one per chain instead of one for all chains, which would be too tall to read

In file names, `<chain>` is the chain ID. When two chain IDs differ only by case (e.g.
`A` and `a` in a large assembly), the one that is not upper case gets a `_lower` suffix
(`overview_a_lower.png`), because macOS and Windows treat `overview_A.png` and
`overview_a.png` as the same file.

FleXgeo2 will not write into an output folder that already contains files, so results
from different runs never mix. To rerun into the same folder, add `--overwrite`
(`OutputConfig(overwrite=True)` in Python): the files of the earlier run are removed
first, and any other files in the folder are kept.

---

## Notes

- The prototype expects Melodia to return columns including `model`, `chain`, `order`, `name`, `curvature`, and `torsion`.
- Model identifiers are the PDB `MODEL` serial numbers (1-based), so `--reference-model 1` selects `MODEL 1`. Files without `MODEL` records are reported as model `1`.
- Residue positions are plotted using Melodia's `order` column, which is the author residue number from the PDB file.
- Chains, reference models and residue ranges are validated before Melodia runs, so input mistakes fail fast with a one-line error.
- `geometry/residues.csv` includes `dmax`. The trimmed extrema used to compute it (`curvature_dmax_min`, `curvature_dmax_max`, `torsion_dmax_min`, `torsion_dmax_max`) and the histogram bin widths are in `AnalysisResult.residue_summary_df` when using FleXgeo2 from Python (see notebook 06).
- `dmax` trims only sparse extreme histogram bins. The default threshold is `0.01`, meaning only extreme bins with less than 1% of a residue's observations are ignored.
- The model summaries (`geometry/models.csv`, and `geometry/models_by_chain.csv` for multi-chain input) include the mean absolute deviation from the ensemble mean, which is useful as a first-pass conformational variability signal.
- Distance matrices use rows for ensemble models and columns for residues, with each cell storing the Euclidean distance to the chosen reference in `(curvature, torsion)` space.
- Residue clustering treats each residue independently and clusters the ensemble conformations using only that residue's curvature and torsion values.
- Residue-range clustering concatenates curvature and torsion values across the selected window into one feature vector per conformation, then runs one HDBSCAN solution for that whole region.

---

## How to cite?

> da Silva Neto AM, Silva SR, Vendruscolo M, Camilloni C, Montalvão RW. A superposition free method for protein conformational ensemble analyses and local clustering based on a differential geometry representation of backbone. Proteins. 2019; 87: 302–312. https://doi.org/10.1002/prot.25652

