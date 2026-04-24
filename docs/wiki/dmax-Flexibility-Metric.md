# dmax Flexibility Metric

`dmax` is a per-residue flexibility metric. It measures the outlier-trimmed maximum spread of an ensemble in `(curvature, torsion)` space.

## Interpretation

A high `dmax` value marks a residue whose local backbone geometry spans a wide region across the ensemble. This can highlight flexible residues, residues involved in conformational change, or local regions where models disagree.

`dmax` is computed from trimmed extrema rather than raw extrema so that sparse extreme observations have less influence on the final value.

## Output columns

`residue_summary.csv` includes:

- `curvature_dmax_min`
- `curvature_dmax_max`
- `torsion_dmax_min`
- `torsion_dmax_max`
- `dmax`

The four extrema columns show the trimmed curvature and torsion bounds used in the calculation.

## Command-line tuning

The default outlier fraction is `0.01`.

```bash
flexgeo2 path/to/ensemble.pdb --dmax-outlier-fraction 0.01
```

Increase the value to trim more sparse extreme bins:

```bash
flexgeo2 path/to/ensemble.pdb --dmax-outlier-fraction 0.05
```

The fraction must be in `[0, 1)`.

## Ranking flexible residues

After a Python run:

```python
from flexgeo2 import AnalysisConfig, FlexGeo2App, OutputConfig

result = FlexGeo2App().run(
    AnalysisConfig(
        pdb_file="path/to/ensemble.pdb",
        output=OutputConfig(write_files=False),
    )
)

top_dmax = result.residue_summary_df.sort_values("dmax", ascending=False).head(10)
top_dmax[["chain", "order", "residue_label", "models", "dmax"]]
```

These rows identify residues with the widest trimmed curvature/torsion spread.

## Tuning without rerunning Melodia

Once `raw_df` exists, you can recompute residue summaries with different `dmax_outlier_fraction` values without loading the PDB again.

```python
from flexgeo2.geometry import GeometryService

geometry = GeometryService()

summary_1pct = geometry.summarize(result.raw_df, dmax_outlier_fraction=0.01)
summary_5pct = geometry.summarize(result.raw_df, dmax_outlier_fraction=0.05)

comparison = summary_1pct[["chain", "order", "residue_label", "dmax"]].merge(
    summary_5pct[["chain", "order", "residue_label", "dmax"]],
    on=["chain", "order", "residue_label"],
    suffixes=("_1pct", "_5pct"),
)

comparison["delta"] = comparison["dmax_5pct"] - comparison["dmax_1pct"]
comparison.reindex(comparison["delta"].abs().sort_values(ascending=False).index).head(10)
```

Large absolute `delta` values mark residues whose `dmax` is sensitive to the trimming threshold.

For a runnable version of this workflow, see `notebooks/06_dmax_library_api.ipynb`.
