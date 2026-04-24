# Command Line Guide

The `flexgeo2` command computes Melodia differential-geometry descriptors from a PDB file and writes CSV tables and plots.

```bash
flexgeo2 path/to/structure.pdb [options]
```

## Basic input and output

```bash
flexgeo2 path/to/structure.pdb
```

Options:

- `pdb_file`: input PDB file
- `--output-dir PATH`: directory for CSV and plot outputs, default `results`
- `--output-verbose`: write detailed intermediate tables and per-chain output folders

Examples:

```bash
flexgeo2 path/to/structure.pdb --output-dir results
flexgeo2 path/to/structure.pdb --output-verbose
```

## Chain filtering and workers

Use `--chain` to keep specific chains. The option can be repeated.

```bash
flexgeo2 path/to/ensemble.pdb --chain A
flexgeo2 path/to/ensemble.pdb --chain A --chain B
```

Use `--n-jobs` to set the worker count passed to Melodia for multi-model files.

```bash
flexgeo2 path/to/ensemble.pdb --chain A --n-jobs 4
```

## Plot controls

FleXgeo2 overlays individual model traces on ensemble plots by default, up to the configured limit.

```bash
flexgeo2 path/to/ensemble.pdb --max-models-in-plot 20
```

To hide individual model overlays and show only the ensemble mean and standard deviation band:

```bash
flexgeo2 path/to/ensemble.pdb --hide-model-traces
```

## dmax tuning

`dmax` is computed from outlier-trimmed curvature and torsion extrema. The default trimming threshold is `0.01`, meaning sparse extreme bins with less than 1% of a residue's observations are ignored.

```bash
flexgeo2 path/to/ensemble.pdb --dmax-outlier-fraction 0.02
```

The value must be greater than or equal to `0` and less than `1`.

## Reference-state comparison

Use `--reference-model` to compare every model in the input ensemble to one model from the same input PDB.

```bash
flexgeo2 path/to/ensemble.pdb --reference-model 1
```

Use `--reference-pdb` to compare against a state from another PDB file.

```bash
flexgeo2 path/to/ensemble.pdb \
  --reference-pdb path/to/reference_state.pdb \
  --reference-pdb-model 1
```

`--reference-model` and `--reference-pdb` are mutually exclusive. `--reference-pdb-model` requires `--reference-pdb`; if it is omitted, FleXgeo2 uses the first model in the reference PDB.

## Per-residue clustering

Use HDBSCAN to cluster conformations independently for each residue in `(curvature, torsion)` space.

```bash
flexgeo2 path/to/ensemble.pdb --cluster-residues
```

Tune HDBSCAN with:

```bash
flexgeo2 path/to/ensemble.pdb \
  --cluster-residues \
  --cluster-min-size 8 \
  --cluster-min-samples 4
```

## Residue-range clustering

Use one combined geometric signature across a residue window.

```bash
flexgeo2 path/to/ensemble.pdb --cluster-residue-range 45-54
```

Repeat the option to analyze multiple windows.

```bash
flexgeo2 path/to/ensemble.pdb \
  --chain A \
  --cluster-residue-range 45-54 \
  --cluster-residue-range 90-99
```

## Option reference

| Option | Purpose | Default |
| --- | --- | --- |
| `pdb_file` | Input PDB file | Required |
| `--output-dir` | Directory for CSV and plot outputs | `results` |
| `--chain` | Chain ID to keep; can be repeated | All chains |
| `--n-jobs` | Workers passed to Melodia | `1` |
| `--max-models-in-plot` | Maximum individual model traces per chain plot | `12` |
| `--hide-model-traces` | Hide individual model traces | `False` |
| `--dmax-outlier-fraction` | Sparse extreme-bin cutoff for `dmax` | `0.01` |
| `--reference-model` | Input model to use as reference state | None |
| `--reference-pdb` | External PDB to use as reference state | None |
| `--reference-pdb-model` | Model from external reference PDB | First model |
| `--cluster-residues` | Run per-residue HDBSCAN clustering | `False` |
| `--cluster-min-size` | HDBSCAN `min_cluster_size` | `5` |
| `--cluster-min-samples` | HDBSCAN `min_samples` | None |
| `--cluster-residue-range` | Residue window for range clustering; can be repeated | None |
| `--output-verbose` | Write detailed outputs and per-chain folders | `False` |
