"""Self-describing output folders: README.md (human guide) and run.json (manifest)."""

from __future__ import annotations

import dataclasses
import hashlib
import json
import platform
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from flexgeo2.models import AnalysisResult

PACKAGES = ("FleXgeo2", "melodia-py", "biopython", "numpy", "pandas", "hdbscan", "matplotlib")

_RESIDUE = {
    "chain": "Chain identifier.",
    "order": "Residue number from the PDB file (author numbering).",
    "name": "Residue name.",
    "residue_label": "Residue name and number, e.g. ALA45.",
}
_MODEL_RESIDUE = {
    "chain": "Chain identifier.",
    "model": "PDB MODEL number of the conformation.",
    "order": "Residue number from the PDB file (author numbering).",
    "name": "Residue name.",
    "residue_label": "Residue name and number, e.g. ALA45.",
}

_CURVATURE_TORSION = {
    "curvature": "Frenet-Serret curvature of the C-alpha spline at this residue (1/A).",
    "torsion": "Frenet-Serret torsion of the C-alpha spline at this residue (1/A).",
}

_DESCRIPTORS = {
    **_MODEL_RESIDUE,
    **_CURVATURE_TORSION,
    "arc_length": "C-alpha spline arc length over a 3-residue window (A).",
    "writhing": "Gauss writhing number over a 5-residue window.",
    "phi": "Backbone phi dihedral (degrees; empty for the first residue).",
    "psi": "Backbone psi dihedral (degrees; empty for the last residue).",
}
_MODEL_SUMMARY = {
    "model": "PDB MODEL number of the conformation.",
    "residues": "Number of residues in the model.",
    "curvature_mean": "Mean curvature over the model's residues.",
    "torsion_mean": "Mean torsion over the model's residues.",
    "curvature_std": "Standard deviation of curvature over the model's residues.",
    "torsion_std": "Standard deviation of torsion over the model's residues.",
    "curvature_mean_abs_deviation": (
        "Mean absolute difference from the ensemble-mean curvature, over residues. "
        "High values flag conformations that differ from the rest."
    ),
    "torsion_mean_abs_deviation": (
        "Mean absolute difference from the ensemble-mean torsion, over residues."
    ),
}
_RANGE = {
    "chain": "Chain identifier.",
    "range_start": "First residue number of the range.",
    "range_end": "Last residue number of the range.",
    "range_label": "Range as START-END.",
}


@dataclass(frozen=True)
class OutputFile:
    """One entry of the file guide; ``pattern`` is a glob relative to the output folder."""

    pattern: str
    display: str
    description: str
    columns: dict[str, str] = field(default_factory=dict)


FILE_GUIDE: tuple[OutputFile, ...] = (
    OutputFile("README.md", "README.md", "This guide."),
    OutputFile(
        "run.json",
        "run.json",
        "Parameters, package versions, input checksums and the list of files written.",
    ),
    OutputFile(
        "overview.png",
        "overview.png",
        "Per-residue results along the sequence, one panel each: curvature and torsion "
        "(ensemble mean +/- standard deviation, with individual model traces), dmax, and, "
        "when those analyses ran, clusters per residue and the mean distance to the "
        "reference (+/- standard deviation).",
    ),
    OutputFile(
        "geometry/descriptors.csv",
        "geometry/descriptors.csv",
        "Melodia descriptors for every model and residue.",
        _DESCRIPTORS,
    ),
    OutputFile(
        "geometry/residues.csv",
        "geometry/residues.csv",
        "Per-residue statistics over the ensemble, including dmax.",
        {
            **_RESIDUE,
            "curvature_mean": "Mean curvature over models.",
            "curvature_std": "Standard deviation of curvature over models.",
            "curvature_min": "Minimum curvature over models.",
            "curvature_max": "Maximum curvature over models.",
            "torsion_mean": "Mean torsion over models.",
            "torsion_std": "Standard deviation of torsion over models.",
            "torsion_min": "Minimum torsion over models.",
            "torsion_max": "Maximum torsion over models.",
            "models": "Number of models with this residue.",
            "dmax": (
                "Spread of the residue in (curvature, torsion) space: diagonal of the "
                "curvature x torsion range after trimming sparse extreme histogram bins. "
                "Higher means more flexible. The trimmed ranges are in "
                "AnalysisResult.residue_summary_df when using FleXgeo2 from Python."
            ),
        },
    ),
    OutputFile(
        "geometry/models.csv",
        "geometry/models.csv",
        "Per-model summary over all chains.",
        {k: v for k, v in _MODEL_SUMMARY.items() if not k.endswith("_std")},
    ),
    OutputFile(
        "geometry/models_by_chain.csv",
        "geometry/models_by_chain.csv",
        "Per-model summary computed separately for each chain (written when the input has "
        "more than one chain).",
        {"chain": "Chain identifier.", **_MODEL_SUMMARY},
    ),
    OutputFile(
        "reference/distances.csv",
        "reference/distances.csv",
        "Distance to the reference in (curvature, torsion) space, per model and residue.",
        {
            "model": "PDB MODEL number of the conformation.",
            **_RESIDUE,
            "curvature": "Curvature of this model at this residue (1/A).",
            "torsion": "Torsion of this model at this residue (1/A).",
            "reference_curvature": "Curvature of the reference at this residue.",
            "reference_torsion": "Torsion of the reference at this residue.",
            "distance_to_reference": "Euclidean distance in (curvature, torsion) space.",
            "reference_label": "Which structure was used as the reference.",
        },
    ),
    OutputFile(
        "reference/residues.csv",
        "reference/residues.csv",
        "Per-residue statistics of the distance to the reference over models.",
        {
            **_RESIDUE,
            "distance_mean": "Mean distance to the reference over models.",
            "distance_std": "Standard deviation of the distance over models.",
            "distance_min": "Minimum distance over models.",
            "distance_max": "Maximum distance over models.",
            "models": "Number of models compared at this residue.",
        },
    ),
    OutputFile(
        "reference/heatmap.png",
        "reference/heatmap.png",
        "Distance to the reference for every model (x) and residue (y).",
    ),
    OutputFile(
        "reference/matrices/*.csv",
        "reference/matrices/<chain>.csv",
        "Distances as a models x residues matrix, one file per chain (with --distance-matrices).",
    ),
    OutputFile(
        "clusters/assignments.csv",
        "clusters/assignments.csv",
        "HDBSCAN cluster of every model at every residue, clustered independently per "
        "residue in (curvature, torsion) space.",
        {
            **_MODEL_RESIDUE,
            **_CURVATURE_TORSION,
            "cluster": "Cluster label at this residue; -1 means noise (no cluster).",
            "cluster_probability": "HDBSCAN membership strength (0 for noise).",
        },
    ),
    OutputFile(
        "clusters/residues.csv",
        "clusters/residues.csv",
        "Per-residue clustering summary.",
        {
            **_RESIDUE,
            "models": "Number of models clustered.",
            "n_clusters": "Number of clusters found (noise excluded).",
            "noise_fraction": "Fraction of models labelled as noise.",
        },
    ),
    OutputFile(
        "clusters/clusters.png",
        "clusters/clusters.png",
        "Cluster of every model (rows) at every residue (columns); grey is noise. The bars "
        "above show the number of clusters per residue. Labels are assigned independently "
        "at each residue, so cluster 0 at one residue is unrelated to cluster 0 at another.",
    ),
    OutputFile(
        "range_clusters/assignments.csv",
        "range_clusters/assignments.csv",
        "HDBSCAN cluster of every model for each residue range, using the curvature and "
        "torsion of all residues in the range together.",
        {
            **_RANGE,
            "model": "PDB MODEL number of the conformation.",
            "cluster": "Cluster label for this range; -1 means noise (no cluster).",
            "cluster_probability": "HDBSCAN membership strength (0 for noise).",
            "pc1": "First principal component of the range's features (for plotting).",
            "pc2": "Second principal component of the range's features (for plotting).",
        },
    ),
    OutputFile(
        "range_clusters/ranges.csv",
        "range_clusters/ranges.csv",
        "Per-range clustering summary.",
        {
            **_RANGE,
            "residues": "Number of residues in the range.",
            "models": "Number of models clustered.",
            "n_clusters": "Number of clusters found (noise excluded).",
            "noise_fraction": "Fraction of models labelled as noise.",
        },
    ),
    OutputFile(
        "range_clusters/*.png",
        "range_clusters/<chain>_<start-end>.png",
        "Models projected on the first two principal components, coloured by cluster.",
    ),
)


def written_files(output_dir: Path) -> list[str]:
    """Relative paths of files matched by the guide (plus README.md and run.json)."""
    found = set()
    for entry in FILE_GUIDE:
        found.update(
            path.relative_to(output_dir).as_posix()
            for path in output_dir.glob(entry.pattern)
            if path.is_file()
        )
    return sorted(found | {"README.md", "run.json"})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def package_versions() -> dict[str, str | None]:
    versions = {}
    for package in PACKAGES:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, list | tuple):
        return [_jsonable(item) for item in value]
    return value


def build_manifest(result: AnalysisResult, output_dir: Path, created: datetime) -> dict:
    raw_df = result.raw_df
    config = result.config
    reference = None
    if result.distance_result is not None:
        reference = {"label": result.distance_result.reference_label}
        if config is not None and config.reference is not None and config.reference.pdb_file:
            reference_pdb = Path(config.reference.pdb_file).resolve()
            reference["pdb_file"] = str(reference_pdb)
            reference["sha256"] = sha256(reference_pdb)
    return {
        "flexgeo2_version": package_versions()["FleXgeo2"],
        "created_utc": created.isoformat(timespec="seconds"),
        "input": {
            "pdb_file": str(result.pdb_file),
            "sha256": sha256(result.pdb_file) if Path(result.pdb_file).is_file() else None,
            "models": int(raw_df["model"].nunique()),
            "chains": sorted(str(chain) for chain in raw_df["chain"].unique()),
            "residues": int(raw_df.groupby(["chain", "order"]).ngroups),
        },
        "reference": reference,
        "parameters": _jsonable(dataclasses.asdict(config)) if config is not None else None,
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": package_versions(),
        },
        "outputs": written_files(output_dir),
    }


def _table(headers: list[str], rows: list[list]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]
    return lines


def _residue(row) -> str:
    return f"{row['chain']} {row['residue_label']}" if row["chain"] else row["residue_label"]


def _analyses(result: AnalysisResult) -> list[str]:
    config = result.config
    lines = []
    fraction = config.dmax_outlier_fraction if config is not None else None
    lines.append(
        "- Backbone geometry: curvature, torsion and dmax"
        + (f" (dmax outlier fraction {fraction})." if fraction is not None else ".")
    )
    if result.distance_result is not None:
        lines.append(f"- Distance to reference: {result.distance_result.reference_label}.")
    clustering = config.clustering if config is not None else None
    hdbscan = (
        f" (HDBSCAN, min_cluster_size {clustering.min_cluster_size}, "
        f"min_samples {clustering.min_samples if clustering.min_samples else 'default'})"
        if clustering is not None
        else ""
    )
    if result.residue_clustering is not None:
        lines.append(f"- Per-residue clustering{hdbscan}.")
    if result.residue_range_clustering is not None:
        ranges = ", ".join(result.residue_range_clustering.summary_df["range_label"].unique())
        lines.append(f"- Residue-range clustering of {ranges}{hdbscan}.")
    return lines


def _key_results(result: AnalysisResult, top: int = 5) -> list[str]:
    lines = ["### Most flexible residues (highest dmax)", ""]
    flexible = result.residue_summary_df.nlargest(top, "dmax")
    lines += _table(
        ["Residue", "dmax"],
        [[_residue(row), f"{row['dmax']:.3f}"] for _, row in flexible.iterrows()],
    )

    if result.distance_result is not None:
        furthest = result.distance_result.summary_df.nlargest(top, "distance_mean")
        lines += ["", "### Residues furthest from the reference (mean distance)", ""]
        lines += _table(
            ["Residue", "Mean distance"],
            [[_residue(row), f"{row['distance_mean']:.3f}"] for _, row in furthest.iterrows()],
        )

    if result.residue_clustering is not None:
        summary = result.residue_clustering.summary_df
        split = summary[summary["n_clusters"] >= 2].sort_values(
            ["n_clusters", "noise_fraction"], ascending=[False, True]
        )
        all_noise = int((summary["noise_fraction"] == 1.0).sum())
        lines += ["", "### Per-residue clustering", ""]
        lines.append(f"- {len(split)} of {len(summary)} residues split into two or more clusters.")
        if not split.empty:
            shown = ", ".join(
                f"{_residue(row)} ({int(row['n_clusters'])})"
                for _, row in split.head(10).iterrows()
            )
            more = f", and {len(split) - 10} more" if len(split) > 10 else ""
            lines.append(f"- Residues with most clusters: {shown}{more}.")
        lines.append(f"- {all_noise} residues have every conformation labelled as noise.")
        lines.append("- Map of every model's cluster at every residue: `clusters/clusters.png`.")

    if result.residue_range_clustering is not None:
        ranges = result.residue_range_clustering.summary_df
        lines += ["", "### Residue-range clustering", ""]
        lines += _table(
            ["Chain", "Range", "Clusters", "Noise fraction"],
            [
                [
                    row["chain"],
                    row["range_label"],
                    int(row["n_clusters"]),
                    f"{row['noise_fraction']:.2f}",
                ]
                for _, row in ranges.iterrows()
            ],
        )
        if (ranges["n_clusters"] == 0).any():
            lines += [
                "",
                "Ranges with 0 clusters had every conformation labelled as noise; consider a "
                "smaller `--cluster-min-size`.",
            ]
    return lines


def _file_guide(output_dir: Path) -> list[str]:
    lines = []
    for entry in FILE_GUIDE:
        if entry.pattern not in ("README.md", "run.json") and not any(
            path.is_file() for path in output_dir.glob(entry.pattern)
        ):
            continue
        lines += [f"### `{entry.display}`", "", entry.description, ""]
        if entry.columns:
            lines += _table(
                ["Column", "Description"], [[f"`{k}`", v] for k, v in entry.columns.items()]
            )
            lines.append("")
    return lines


def render_readme(result: AnalysisResult, output_dir: Path, created: datetime) -> str:
    raw_df = result.raw_df
    chains = sorted(str(chain) for chain in raw_df["chain"].unique())
    n_models = int(raw_df["model"].nunique())
    n_residues = int(raw_df.groupby(["chain", "order"]).ngroups)
    version = package_versions()["FleXgeo2"] or "unknown"
    lines = [
        "# FleXgeo2 results",
        "",
        f"Generated {created:%Y-%m-%d %H:%M} UTC by FleXgeo2 {version} from "
        f"`{Path(result.pdb_file).name}`. Parameters, package versions and input checksums "
        "are in `run.json`.",
        "",
        "Start with `overview.png`, then the key results below.",
        "",
        "## Input",
        "",
        f"- File: `{result.pdb_file}`",
        f"- {n_models} models, {len(chains)} chain(s) ({', '.join(chains)}), {n_residues} residues",
        "",
        "## Analyses",
        "",
        *_analyses(result),
        "",
        "## Key results",
        "",
        *_key_results(result),
        "",
        "## Files",
        "",
        *_file_guide(output_dir),
    ]
    return "\n".join(lines).rstrip() + "\n"


def write_report(result: AnalysisResult, output_dir: Path) -> tuple[Path, Path]:
    """Write README.md and run.json into ``output_dir``; call after all other outputs."""
    created = datetime.now(timezone.utc)
    readme = output_dir / "README.md"
    manifest = output_dir / "run.json"
    readme.write_text(render_readme(result, output_dir, created))
    manifest.write_text(json.dumps(build_manifest(result, output_dir, created), indent=2) + "\n")
    return readme, manifest
