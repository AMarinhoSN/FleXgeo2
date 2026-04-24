# Reference-State Comparison

FleXgeo2 can compare every model in an ensemble to a reference state. The distance is the Euclidean distance between each model and the reference in each residue's `(curvature, torsion)` space.

Large distances mark residues whose local backbone geometry differs strongly from the selected reference.

## Reference model from the input ensemble

Use `--reference-model` when the reference state is already present in the input PDB.

```bash
flexgeo2 path/to/ensemble.pdb --reference-model 1
```

Python equivalent:

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App, ReferenceConfig

result = FlexGeo2App().run(
    AnalysisConfig(
        pdb_file="path/to/ensemble.pdb",
        reference=ReferenceConfig(model_id="1"),
    )
)
```

## External reference PDB

Use `--reference-pdb` when the reference state is stored in another PDB file.

```bash
flexgeo2 path/to/ensemble.pdb \
  --reference-pdb path/to/reference_state.pdb \
  --reference-pdb-model 1
```

Python equivalent:

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App, ReferenceConfig

result = FlexGeo2App().run(
    AnalysisConfig(
        pdb_file="path/to/ensemble.pdb",
        reference=ReferenceConfig(
            pdb_file="path/to/reference_state.pdb",
            pdb_model_id="1",
        ),
    )
)
```

If `pdb_model_id` is omitted, FleXgeo2 uses the first model found in the reference PDB.

## In-memory results

When reference comparison is enabled, `result.distance_result` contains:

- `long_df`: distance values for each model, chain, and residue
- `summary_df`: residue-level distance summaries
- `reference_label`: human-readable reference state label

Example:

```python
distance_result = result.distance_result
assert distance_result is not None

distance_result.summary_df.sort_values(
    "distance_mean",
    ascending=False,
).head(10)
```

## Output files

Default reference outputs:

- `distance_to_reference_summary.csv`
- `plots/distance_to_reference_heatmap.png`

Verbose reference outputs:

- `distance_to_reference_long.csv`
- `distance_matrices/<chain>_distance_matrix.csv`
- `chains/<chain>/distance_to_reference_long.csv`
- `chains/<chain>/distance_to_reference_summary.csv`
- `chains/<chain>/distance_to_reference_matrix.csv`
- `chains/<chain>/distance_to_reference_heatmap.png`

## How to read the heatmap

The distance heatmap places models on one axis and residues on the other. Brighter regions show larger distance to the reference state in `(curvature, torsion)` space. Localized bright columns can indicate residues with strong reference-state differences; bands across many residues can indicate model-level conformational shifts.

For a runnable version of this workflow, see `notebooks/03_reference_comparison.ipynb`.
