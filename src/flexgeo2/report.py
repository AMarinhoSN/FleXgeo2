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

from flexgeo2.clustering import clustering_parameters
from flexgeo2.config import MAX_CHAINS_PER_FIGURE, PLOT_FORMATS
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


_OVERVIEW = (
    "Per-residue results along the sequence, one panel each: curvature and torsion "
    "(ensemble mean +/- standard deviation, with individual model traces), dmax, and, "
    "when those analyses ran, clusters per residue and the mean distance to the "
    "reference (+/- standard deviation)."
)
_HEATMAP = (
    "Distance to the reference for every residue (x) and model (y), laid out like "
    "the cluster map. The colour scale ends at the 99th percentile of each chain's "
    "distances; larger distances share the top colour (arrow on the colour bar)."
)
_CLUSTER_MAP = (
    "Cluster of every model (rows) at every residue (columns); light grey is noise, and "
    'clusters 18 and above share dark grey ("Other clusters"). The bars above show '
    "the number of clusters per residue. Labels are assigned independently at each "
    "residue, so cluster 0 at one residue is unrelated to cluster 0 at another."
)
_PER_CHAIN = (
    f" One file per chain: with more than {MAX_CHAINS_PER_FIGURE} chains, a single figure "
    "for all of them would be too tall to read."
)

# Figure entries use "{ext}" for the plot format; see file_guide().
FILE_GUIDE: tuple[OutputFile, ...] = (
    OutputFile("README.md", "README.md", "This guide."),
    OutputFile(
        "run.json",
        "run.json",
        "Parameters, package versions, input checksums and the list of files written.",
    ),
    OutputFile("overview.{ext}", "overview.{ext}", _OVERVIEW),
    OutputFile("overview_*.{ext}", "overview_<chain>.{ext}", _OVERVIEW + _PER_CHAIN),
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
    OutputFile("reference/heatmap.{ext}", "reference/heatmap.{ext}", _HEATMAP),
    OutputFile(
        "reference/heatmap_*.{ext}", "reference/heatmap_<chain>.{ext}", _HEATMAP + _PER_CHAIN
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
            "noise_fraction": "Fraction of models labelled as noise: outside the dense core "
            "of every cluster; often around half the models, even for a residue with one "
            "broad state.",
        },
    ),
    OutputFile("clusters/clusters.{ext}", "clusters/clusters.{ext}", _CLUSTER_MAP),
    OutputFile(
        "clusters/clusters_*.{ext}", "clusters/clusters_<chain>.{ext}", _CLUSTER_MAP + _PER_CHAIN
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
            "noise_fraction": "Fraction of models labelled as noise: outside the dense core "
            "of every cluster.",
        },
    ),
    OutputFile(
        "residue_plots/*.{ext}",
        "residue_plots/<chain>_<residue number>_<name>.{ext}",
        "Curvature vs torsion of each residue chosen with --plot-residues, one point per "
        "model. Points are coloured by cluster when per-residue clustering ran, and the "
        "reference is marked with a star when one was given. The title gives the residue's "
        "dmax.",
    ),
    OutputFile(
        "range_clusters/*.{ext}",
        "range_clusters/<chain>_<start-end>.{ext}",
        "Models projected on the first two principal components, coloured by cluster.",
    ),
)


def file_guide(plot_format: str = "png") -> tuple[OutputFile, ...]:
    """The file guide with figure names in ``plot_format``."""
    return tuple(
        dataclasses.replace(
            entry,
            pattern=entry.pattern.replace("{ext}", plot_format),
            display=entry.display.replace("{ext}", plot_format),
        )
        for entry in FILE_GUIDE
    )


def written_files(output_dir: Path) -> list[str]:
    """Relative paths of files matched by the guide (plus README.md and run.json).

    Figures are matched in every plot format, so an earlier run's figures are found
    even when it used another format.
    """
    found = set()
    for plot_format in PLOT_FORMATS:
        for entry in file_guide(plot_format):
            found.update(
                path.relative_to(output_dir).as_posix()
                for path in output_dir.glob(entry.pattern)
                if path.is_file()
            )
    return sorted(found | {"README.md", "run.json"})


def _plot_format(result: AnalysisResult) -> str:
    return result.config.output.plot_format if result.config is not None else "png"


def _figure_name(result: AnalysisResult, stem: str) -> str:
    """File name of a figure with one panel per chain: one file, or one per chain."""
    _, chains, _ = _input_counts(result)
    per_chain = "_<chain>" if len(chains) > MAX_CHAINS_PER_FIGURE else ""
    return f"{stem}{per_chain}.{_plot_format(result)}"


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
        # The HDBSCAN settings actually used, with defaults filled in.
        "hdbscan": _hdbscan_settings(result),
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


def _input_counts(result: AnalysisResult) -> tuple[int, list[str], int]:
    raw_df = result.raw_df
    chains = sorted(str(chain) for chain in raw_df["chain"].unique())
    return int(raw_df["model"].nunique()), chains, int(raw_df.groupby(["chain", "order"]).ngroups)


def _most_flexible(result: AnalysisResult, top: int):
    return result.residue_summary_df.nlargest(top, "dmax")


def _furthest_from_reference(result: AnalysisResult, top: int):
    return result.distance_result.summary_df.nlargest(top, "distance_mean")


def _split_residues(result: AnalysisResult):
    summary = result.residue_clustering.summary_df
    return summary[summary["n_clusters"] >= 2].sort_values(
        ["n_clusters", "noise_fraction"], ascending=[False, True]
    )


def _all_noise_ranges(result: AnalysisResult):
    ranges = result.residue_range_clustering.summary_df
    return ranges[ranges["n_clusters"] == 0]


ALL_NOISE_HINT = (
    "Ranges with 0 clusters had every conformation labelled as noise; consider a "
    "smaller `--cluster-min-size`."
)


def _hdbscan_settings(result: AnalysisResult) -> dict | None:
    """HDBSCAN settings used for clustering, or None when no clustering ran."""
    clustering = result.config.clustering if result.config is not None else None
    ran = result.residue_clustering is not None or result.residue_range_clustering is not None
    if clustering is None or not ran:
        return None
    size, samples = clustering_parameters(
        clustering.min_cluster_size, clustering.min_samples, result.raw_df["model"].nunique()
    )
    return {"min_cluster_size": size, "min_samples": samples, "allow_single_cluster": True}


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
    settings = _hdbscan_settings(result)
    hdbscan = (
        f" (HDBSCAN, min_cluster_size {settings['min_cluster_size']}, "
        f"min_samples {settings['min_samples']}, single cluster allowed)"
        if settings is not None
        else ""
    )
    if result.residue_clustering is not None:
        lines.append(f"- Per-residue clustering{hdbscan}.")
    if result.residue_range_clustering is not None:
        ranges = ", ".join(result.residue_range_clustering.summary_df["range_label"].unique())
        lines.append(f"- Residue-range clustering of {ranges}{hdbscan}.")
    if config is not None and config.output.plot_residues:
        lines.append(f"- Residue plots for {', '.join(config.output.plot_residues)}.")
    return lines


def _key_results(result: AnalysisResult, top: int = 5) -> list[str]:
    lines = ["### Most flexible residues (highest dmax)", ""]
    flexible = _most_flexible(result, top)
    lines += _table(
        ["Residue", "dmax"],
        [[_residue(row), f"{row['dmax']:.3f}"] for _, row in flexible.iterrows()],
    )

    if result.distance_result is not None:
        furthest = _furthest_from_reference(result, top)
        lines += ["", "### Residues furthest from the reference (mean distance)", ""]
        lines += _table(
            ["Residue", "Mean distance"],
            [[_residue(row), f"{row['distance_mean']:.3f}"] for _, row in furthest.iterrows()],
        )

    if result.residue_clustering is not None:
        summary = result.residue_clustering.summary_df
        split = _split_residues(result)
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
        lines.append(
            "- Map of every model's cluster at every residue: "
            f"`clusters/{_figure_name(result, 'clusters')}`."
        )

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
        if not _all_noise_ranges(result).empty:
            lines += ["", ALL_NOISE_HINT]
    return lines


def _file_guide(output_dir: Path, plot_format: str) -> list[str]:
    lines = []
    for entry in file_guide(plot_format):
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
    n_models, chains, n_residues = _input_counts(result)
    version = package_versions()["FleXgeo2"] or "unknown"
    lines = [
        "# FleXgeo2 results",
        "",
        f"Generated {created:%Y-%m-%d %H:%M} UTC by FleXgeo2 {version} from "
        f"`{Path(result.pdb_file).name}`. Parameters, package versions and input checksums "
        "are in `run.json`.",
        "",
        f"Start with `{_figure_name(result, 'overview')}`, then the key results below.",
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
        *_file_guide(output_dir, _plot_format(result)),
    ]
    return "\n".join(lines).rstrip() + "\n"


def _display_path(path: Path) -> str:
    """``path`` relative to the working directory when it is inside it."""
    try:
        return str(path.relative_to(Path.cwd()))
    except ValueError:
        return str(path)


def render_terminal_summary(result: AnalysisResult, top: int = 3) -> str:
    """A few lines for the terminal: what was analysed, headline results, where to look."""
    n_models, chains, n_residues = _input_counts(result)
    chain_text = f"{len(chains)} chain{'s' if len(chains) != 1 else ''} ({', '.join(chains)})"
    lines = [
        f"FleXgeo2 analysed {Path(result.pdb_file).name}: {n_models} models, {chain_text}, "
        f"{n_residues} residues.",
        "",
    ]

    def residues(frame, column) -> str:
        return ", ".join(f"{_residue(row)} ({row[column]:.3f})" for _, row in frame.iterrows())

    lines.append(f"Most flexible residues (dmax): {residues(_most_flexible(result, top), 'dmax')}")
    if result.distance_result is not None:
        furthest = _furthest_from_reference(result, top)
        lines.append(
            f"Furthest from {result.distance_result.reference_label} (mean distance): "
            f"{residues(furthest, 'distance_mean')}"
        )
    if result.residue_clustering is not None:
        n_split = len(_split_residues(result))
        n_total = len(result.residue_clustering.summary_df)
        lines.append(
            f"Per-residue clustering: {n_split} of {n_total} residues split into two or more "
            "clusters"
        )
    if result.residue_range_clustering is not None:
        for _, row in result.residue_range_clustering.summary_df.iterrows():
            chain = f"{row['chain']} " if row["chain"] else ""
            n_clusters = int(row["n_clusters"])
            hint = " (try a smaller --cluster-min-size)" if n_clusters == 0 else ""
            lines.append(
                f"Range clustering {chain}{row['range_label']}: {n_clusters} "
                f"cluster{'' if n_clusters == 1 else 's'}, {row['noise_fraction']:.0%} noise{hint}"
            )

    plots_dir = result.outputs.residue_plots_dir if result.outputs is not None else None
    if plots_dir is not None:
        n_plots = sum(1 for path in plots_dir.iterdir() if path.is_file())
        lines.append(f"Residue plots: {n_plots} in {_display_path(plots_dir)}/")

    readme = result.outputs.readme
    if readme is not None:
        output_dir = readme.parent
        n_files = len(written_files(output_dir))
        lines += [
            "",
            f"Results: {_display_path(output_dir)}/ ({n_files} files); start with README.md",
        ]
    return "\n".join(lines)


def write_report(result: AnalysisResult, output_dir: Path) -> tuple[Path, Path]:
    """Write README.md and run.json into ``output_dir``; call after all other outputs."""
    created = datetime.now(timezone.utc)
    readme = output_dir / "README.md"
    manifest = output_dir / "run.json"
    readme.write_text(render_readme(result, output_dir, created))
    manifest.write_text(json.dumps(build_manifest(result, output_dir, created), indent=2) + "\n")
    return readme, manifest
