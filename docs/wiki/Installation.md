# Installation

FleXgeo2 requires Python `>=3.10`.

## Editable install

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

This installs the `flexgeo2` command-line entry point and makes the package importable from Python.

## Main dependencies

FleXgeo2 uses:

- `melodia-py` for differential-geometry descriptors
- `pandas` and `numpy` for tables and numeric operations
- `matplotlib` for plots
- `hdbscan` for conformation clustering
- `nglview` and `numba` through the analysis stack

## Verify the install

Run:

```bash
flexgeo2 --help
```

You should see the command-line options for input files, output directories, chain filtering, plotting, reference-state comparison, clustering, and verbose output.

## Development extras

For local testing and linting:

```bash
pip install -e ".[dev]"
```

The development extras include `pytest`, `ruff`, and `pre-commit`.
