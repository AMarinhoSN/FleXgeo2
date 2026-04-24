# Understanding Outputs

FleXgeo2 writes a compact default output set for most analyses. Use `--output-verbose` when you want detailed intermediate tables, distance matrices, and duplicated per-chain output folders.

## Default outputs

A standard run writes:

- `geometry_descriptors.csv`
- `residue_summary.csv`
- `model_summary_overall.csv`
- `plots/ensemble_overview.png`

If reference-state comparison is enabled, FleXgeo2 also writes:

- `distance_to_reference_summary.csv`
- `plots/distance_to_reference_heatmap.png`

If per-residue clustering is enabled, FleXgeo2 also writes:

- `residue_cluster_summary.csv`
- `cluster_plots/<chain>_<residue>_clusters.png`

If residue-range clustering is enabled, FleXgeo2 also writes:

- `residue_range_cluster_summary.csv`
- `range_cluster_plots/<chain>_<start-end>_clusters.png`

## Verbose outputs

With `--output-verbose`, FleXgeo2 adds:

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

## Important tables

`geometry_descriptors.csv` contains the normalized per-model, per-chain, per-residue Melodia descriptors. FleXgeo2 expects Melodia-derived columns including `model`, `chain`, `order`, `name`, `curvature`, and `torsion`.

`residue_summary.csv` summarizes each residue across models. It includes curvature and torsion summaries plus `dmax` and the trimmed extrema used to compute it:

- `curvature_dmax_min`
- `curvature_dmax_max`
- `torsion_dmax_min`
- `torsion_dmax_max`
- `dmax`

`model_summary_overall.csv` summarizes each model across the selected chains. It is useful for finding conformations that differ strongly from the ensemble average.

`model_summary_by_chain.csv`, available in verbose mode, gives the same model-level view split by chain.

`distance_to_reference_summary.csv` summarizes each residue's distance to the selected reference state.

`distance_to_reference_long.csv`, available in verbose mode, keeps the model-by-residue distance values used to build summaries and heatmaps.

`residue_cluster_summary.csv` reports per-residue HDBSCAN clustering results, including cluster counts and noise fractions.

`residue_cluster_assignments.csv`, available in verbose mode, lists the cluster label and probability for each conformation at each residue.

`residue_range_cluster_summary.csv` reports HDBSCAN clustering results for each selected residue window.

`residue_range_cluster_assignments.csv`, available in verbose mode, lists range-level cluster labels, probabilities, and PCA plot coordinates for each conformation.

## Important plots

`plots/ensemble_overview.png` shows residue-level geometry summaries along the sequence. Residue positions use Melodia's `order` column.

`plots/distance_to_reference_heatmap.png` shows models on one axis and residues on the other. Brighter regions indicate residues whose local curvature/torsion geometry is farther from the reference state.

`cluster_plots/<chain>_<residue>_clusters.png` shows per-residue conformation clusters in `(curvature, torsion)` space.

`range_cluster_plots/<chain>_<start-end>_clusters.png` shows residue-range clustering in a two-dimensional PCA projection of the combined geometry signature.

## Choosing default or verbose mode

Use the default output set when you want compact summary artifacts for figures, reports, or first-pass analysis.

Use `--output-verbose` when you need per-chain files, full model-by-residue distance values, cluster assignment tables, or CSV matrices for downstream statistical analysis.
