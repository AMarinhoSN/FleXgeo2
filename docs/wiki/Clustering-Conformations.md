# Clustering Conformations

FleXgeo2 uses HDBSCAN to cluster conformations from local backbone geometry. It supports two clustering workflows:

- per-residue clustering, where each residue is clustered independently in `(curvature, torsion)` space
- residue-range clustering, where a selected residue window is converted into one combined geometric signature per conformation

HDBSCAN labels `-1` mean noise or unassigned conformations.

## Per-residue clustering

Run:

```bash
flexgeo2 path/to/ensemble.pdb --cluster-residues
```

Python equivalent:

```python
from flexgeo2 import AnalysisConfig, ClusteringConfig, FlexGeo2App

result = FlexGeo2App().run(
    AnalysisConfig(
        pdb_file="path/to/ensemble.pdb",
        clustering=ClusteringConfig(cluster_residues=True),
    )
)
```

Each residue is clustered using only that residue's curvature and torsion values across models.

## Residue-range clustering

Run:

```bash
flexgeo2 path/to/ensemble.pdb --cluster-residue-range 45-54
```

You can repeat the option:

```bash
flexgeo2 path/to/ensemble.pdb \
  --chain A \
  --cluster-residue-range 45-54 \
  --cluster-residue-range 90-99
```

Python equivalent:

```python
from flexgeo2 import AnalysisConfig, ClusteringConfig, FlexGeo2App

result = FlexGeo2App().run(
    AnalysisConfig(
        pdb_file="path/to/ensemble.pdb",
        chains=["A"],
        clustering=ClusteringConfig(
            cluster_residue_ranges=["45-54", "90-99"],
        ),
    )
)
```

Residue-range clustering concatenates curvature and torsion values across the selected window into one feature vector per conformation, then runs one HDBSCAN solution for the whole region.

## HDBSCAN parameters

Tune cluster size and sample behavior with:

```bash
flexgeo2 path/to/ensemble.pdb \
  --cluster-residues \
  --cluster-min-size 8 \
  --cluster-min-samples 4
```

Python equivalent:

```python
ClusteringConfig(
    cluster_residues=True,
    min_cluster_size=8,
    min_samples=4,
)
```

`min_cluster_size` defaults to `5`. `min_samples` defaults to `None`, which lets HDBSCAN use its own default behavior.

## In-memory results

Per-residue clustering results are available at `result.residue_clustering`.

```python
residue_clusters = result.residue_clustering
assert residue_clusters is not None

residue_clusters.summary_df.sort_values(
    ["n_clusters", "noise_fraction"],
    ascending=[False, True],
).head(15)
```

Residue-range clustering results are available at `result.residue_range_clustering`.

```python
range_clusters = result.residue_range_clustering
assert range_clusters is not None

range_clusters.summary_df
```

## Output files

Default per-residue clustering outputs:

- `residue_cluster_summary.csv`
- `cluster_plots/<chain>_<residue>_clusters.png`

Verbose per-residue clustering output:

- `residue_cluster_assignments.csv`

Default residue-range clustering outputs:

- `residue_range_cluster_summary.csv`
- `range_cluster_plots/<chain>_<start-end>_clusters.png`

Verbose residue-range clustering output:

- `residue_range_cluster_assignments.csv`

## Important columns

`n_clusters` counts non-noise HDBSCAN clusters.

`noise_fraction` is the fraction of conformations labeled `-1`.

`cluster` is the assigned HDBSCAN cluster label for a conformation.

`cluster_probability` is HDBSCAN's confidence-like membership probability.

For residue-range clustering, `pc1` and `pc2` are PCA projection coordinates used for plotting the high-dimensional residue-window signature.

For a runnable version of this workflow, see `notebooks/04_clustering_conformations.ipynb`.
