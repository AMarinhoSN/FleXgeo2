# FleXgeo2 Wiki

FleXgeo2 is a refactor of [FleXgeo](https://github.com/AMarinhoSN/FleXgeo) built around [Melodia_py](https://github.com/rwmontalvao/Melodia_py). It computes protein backbone differential-geometry descriptors from PDB files and helps researchers study ensemble-aware curvature, torsion, local flexibility, reference-state differences, and conformational clusters.

## Pages

- [Installation](Installation)
- [Quickstart Analysis](Quickstart-Analysis)
- [Command Line Guide](Command-Line-Guide)
- [Python Library Guide](Python-Library-Guide)
- [Understanding Outputs](Understanding-Outputs)
- [dmax Flexibility Metric](dmax-Flexibility-Metric)
- [Reference-State Comparison](Reference-State-Comparison)
- [Clustering Conformations](Clustering-Conformations)
- [Notebook Tutorials](Notebook-Tutorials)
- [Citation](Citation)

## When to use FleXgeo2

Use FleXgeo2 when you want to compare local protein backbone geometry across a PDB ensemble without relying on structural superposition alone. The main workflows are:

- finding residues with large curvature/torsion variation across models
- summarizing ensemble flexibility with per-residue `dmax`
- comparing every model to a chosen reference state
- clustering conformations residue-by-residue or across a residue window
- exporting reproducible CSV tables and plots for downstream analysis

For the fastest path through the software, start with [Installation](Installation) and [Quickstart Analysis](Quickstart-Analysis). For deeper analysis examples, see [Notebook Tutorials](Notebook-Tutorials).
