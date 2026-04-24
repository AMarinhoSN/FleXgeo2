# Notebook Tutorials

The `notebooks/` directory contains executable examples that complement the wiki pages. Use them when you want to inspect data frames directly, adapt plotting code, or build a reproducible analysis workflow.

## Available notebooks

- `notebooks/01_quickstart_analysis.ipynb`: first complete analysis from Python, including returned result objects and default outputs
- `notebooks/02_working_with_results.ipynb`: working with FleXgeo2 result data frames in memory
- `notebooks/03_reference_comparison.ipynb`: comparing an ensemble to a reference model or reference PDB
- `notebooks/04_clustering_conformations.ipynb`: per-residue and residue-range HDBSCAN clustering
- `notebooks/05_cli_to_python_mapping.ipynb`: translating command-line examples into Python config objects
- `notebooks/06_dmax_library_api.ipynb`: computing, inspecting, and retuning `dmax`

## Suggested reading order

Start with `01_quickstart_analysis.ipynb` if you are new to FleXgeo2.

Use `05_cli_to_python_mapping.ipynb` when you already know the command-line workflow and want to move into Python.

Use `03_reference_comparison.ipynb`, `04_clustering_conformations.ipynb`, and `06_dmax_library_api.ipynb` for focused scientific workflows.

Use `02_working_with_results.ipynb` when you want to build custom tables, figures, or downstream analyses from the returned data frames.

## Notes

The notebooks may create temporary output directories while running. When adapting examples for your own project, replace temporary paths with project-specific paths such as `results/`, `analysis/results/`, or another directory under your working folder.
