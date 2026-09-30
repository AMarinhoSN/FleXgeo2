# Changelog

All notable changes to FleXgeo2 are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

Work towards FleXgeo2 2.0.0, the first release of FleXgeo2, a rewrite of
[FleXgeo](https://github.com/AMarinhoSN/FleXgeo) built around
[Melodia_py](https://github.com/rwmontalvao/Melodia_py). Until then the package version
is `2.0.0.dev0`.

A preview of FleXgeo2 was tagged `v2.0.0` on 2026-04-24 and later withdrawn; 2.0.0 will
be released from this work instead. The entries below are the changes since that
preview. It had a different output layout, an `--output-verbose` option and off-by-one
model numbers, so results and scripts written for it need updating: see
[Upgrading from the April 2026 preview](#upgrading-from-the-april-2026-preview).

### Breaking changes

- **Model numbers now match the PDB file.** The preview used Biopython's 0-based model index,
  so every `model` value was one lower than the PDB `MODEL` record and
  `--reference-model 1` selected `MODEL 2`. Model IDs are now the `MODEL` serial numbers;
  files without `MODEL` records report model `1`.
- **New output layout**: one folder per analysis (`geometry/`, `reference/`,
  `clusters/`, `range_clusters/`), with new file names. The per-chain copies under
  `chains/<chain>/` are gone; every table has a `chain` column instead.
- **`--output-verbose` is removed.** Per-model tables (`clusters/assignments.csv`,
  `reference/distances.csv`, `range_clusters/assignments.csv`) are now written by
  default. Distance matrices need `--distance-matrices`. `geometry/models_by_chain.csv`
  is written whenever the input has more than one chain. Passing `--output-verbose`
  exits with a message naming these replacements.
- **Per-residue cluster scatter plots are removed.** `clusters/clusters.png` shows every
  residue in one figure; notebook 04 shows how to plot a chosen residue from Python.
- **FleXgeo2 refuses to write into an output folder that already has files**, so results
  of different runs cannot mix. Use `--overwrite` to replace an earlier run.
- **Table columns**:
  - Melodia's `id` and `code` columns are dropped. Per-model tables start with
    `chain, model, order, name, residue_label`.
  - `geometry/residues.csv` no longer has the dmax internals (`curvature_dmax_min`,
    `curvature_dmax_max`, `torsion_dmax_min`, `torsion_dmax_max`,
    `curvature_dmax_bin_width`, `torsion_dmax_bin_width`). They are still in
    `AnalysisResult.residue_summary_df`.
  - `clusters/assignments.csv` keeps only the row keys, `curvature`, `torsion`,
    `cluster` and `cluster_probability`; the other descriptors are in
    `geometry/descriptors.csv`.
  - Counts are named alike in every table: `n_conformations` is now `models`, and
    `n_residues` is now `residues`.
  - `model` in range cluster assignments is an integer, as in every other table (it was
    a string).
- **Python API**:
  - `OutputConfig.verbose` is replaced by `OutputConfig.distance_matrices`;
    `OutputConfig.overwrite` is new.
  - `OutputArtifacts`: `chains_dir` and `cluster_plots_dir` are removed; `readme`,
    `run_manifest` and `cluster_map_plot` are new.
  - `OutputWriter` no longer takes `chain_plotter` or `residue_cluster_plotter`, and
    takes a new `cluster_map_plotter`. `ResidueClusterPlotter` is removed.
  - The column changes above apply to the result data frames too (`raw_df`, the
    cluster `assignments_df` and `summary_df` tables).

### Added

- Every output folder has a `README.md` (what ran, key results, and a guide to every file
  and column) and a `run.json` (FleXgeo2 version, parameters, Python and package
  versions, SHA-256 of the input and reference files, list of files written).
- `clusters/clusters.png`: the cluster of every model at every residue, with the number
  of clusters per residue above it.
- `--distance-matrices` and `--overwrite` options.
- A short terminal summary: the input, the headline result of each analysis, and where
  the results are. It replaces the list of every output path.
- Up-front validation: chains, reference models and residue ranges are checked before
  Melodia runs, and mistakes are reported as a one-line `flexgeo2: error: ...` (exit
  status 1, or 130 on Ctrl-C) instead of a traceback.
- The CLI rejects invalid values for `--cluster-min-size`, `--cluster-min-samples`,
  `--max-models-in-plot` and `--n-jobs`.
- `AnalysisResult.config` records the configuration that produced a result.

### Changed

- `overview.png` shows every per-residue result along the sequence, one panel each:
  curvature and torsion (mean +/- SD with model traces), `dmax`, and, when those analyses
  ran, clusters per residue and the mean distance to the reference.
- A cluster label has the same colour in every plot, from an 18-colour palette that
  keeps grey for noise.

### Fixed

- The distance heatmap and distance matrices ordered residues alphabetically by label
  (e.g. VAL17 before THR66); they now follow the sequence.
- `--hide-model-traces` and `--max-models-in-plot` had no effect on the default overview
  plot.
- `--reference-pdb-model` without `--reference-pdb` was silently ignored; it is now an
  error. A `ReferenceConfig` with both `model_id` and `pdb_file` is rejected instead of
  silently using `model_id`.
- Spurious floating-point warnings from the PCA projection on macOS.
- `biopython` is declared as a dependency.

### Upgrading from the April 2026 preview

Output files:

| Preview | Now |
| --- | --- |
| `geometry_descriptors.csv` | `geometry/descriptors.csv` |
| `residue_summary.csv` | `geometry/residues.csv` |
| `model_summary_overall.csv` | `geometry/models.csv` |
| `model_summary_by_chain.csv` (verbose) | `geometry/models_by_chain.csv` (more than one chain) |
| `plots/ensemble_overview.png` | `overview.png` |
| `distance_to_reference_long.csv` (verbose) | `reference/distances.csv` |
| `distance_to_reference_summary.csv` | `reference/residues.csv` |
| `plots/distance_to_reference_heatmap.png` | `reference/heatmap.png` |
| `distance_matrices/<chain>_distance_matrix.csv` (verbose) | `reference/matrices/<chain>.csv` (`--distance-matrices`) |
| `residue_cluster_assignments.csv` (verbose) | `clusters/assignments.csv` |
| `residue_cluster_summary.csv` | `clusters/residues.csv` |
| `cluster_plots/<chain>_<residue>_clusters.png` | removed; see `clusters/clusters.png` |
| `residue_range_cluster_assignments.csv` (verbose) | `range_clusters/assignments.csv` |
| `residue_range_cluster_summary.csv` | `range_clusters/ranges.csv` |
| `range_cluster_plots/<chain>_<start-end>_clusters.png` | `range_clusters/<chain>_<start-end>.png` |
| `chains/<chain>/...` (verbose) | removed; filter the `chain` column |

Command line:

| Preview | Now |
| --- | --- |
| `--output-verbose` | `--distance-matrices` for distance matrices; everything else is written by default or removed (see above) |
| rerunning into the same `--output-dir` | add `--overwrite` |

Model numbers: the preview numbered models 0, 1, 2, ... in file order. When the PDB
`MODEL` records are numbered 1, 2, 3, ... (the usual case), add 1 to a preview model ID to
get the current one; otherwise the current ID is the `MODEL` record itself.
