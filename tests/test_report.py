from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import pandas as pd
import pytest

from flexgeo2 import AnalysisConfig, ClusteringConfig, FlexGeo2App, OutputConfig
from flexgeo2.cli.main import main
from flexgeo2.report import file_guide, render_readme, render_terminal_summary, written_files

pytest.importorskip("melodia_py")

MINI_ENSEMBLE = Path(__file__).parent / "data" / "mini_ensemble.pdb"


def copy_chain_a(path: Path, chains: str) -> Path:
    """Write the three-model fixture with its chain A copied as each of ``chains``."""
    lines: list[str] = []
    chain_a: list[str] = []
    for line in MINI_ENSEMBLE.read_text().splitlines():
        if line.startswith("ENDMDL"):
            lines += [f"{atom[:21]}{chain}{atom[22:]}" for chain in chains for atom in chain_a]
            chain_a = []
        if line.startswith("ATOM"):
            chain_a.append(line)
        lines.append(line)
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture(scope="module")
def two_chain_ensemble(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The three-model fixture with its chain A copied as chain B."""
    return copy_chain_a(tmp_path_factory.mktemp("input") / "two_chains.pdb", "B")


@pytest.fixture(scope="module")
def full_run(tmp_path_factory: pytest.TempPathFactory, two_chain_ensemble: Path) -> Path:
    """Every analysis and every optional table, on a two-chain, three-model ensemble."""
    output_dir = tmp_path_factory.mktemp("report") / "out"
    exit_code = main(
        [
            str(two_chain_ensemble),
            "--output-dir",
            str(output_dir),
            "--reference-pdb",
            str(two_chain_ensemble),
            "--reference-pdb-model",
            "2",
            "--cluster-residues",
            "--cluster-residue-range",
            "2-5",
            "--cluster-min-size",
            "2",
            "--distance-matrices",
            "--plot-residues",
            "3,A:5",
        ]
    )
    assert exit_code == 0
    return output_dir


def files_on_disk(output_dir: Path) -> list[str]:
    return sorted(
        path.relative_to(output_dir).as_posix() for path in output_dir.rglob("*") if path.is_file()
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_every_csv_column_is_described(full_run: Path) -> None:
    guide = {entry.pattern: entry for entry in file_guide()}
    csv_files = [name for name in files_on_disk(full_run) if name.endswith(".csv")]
    assert csv_files, "the run should have written CSV files"

    for name in csv_files:
        if name.startswith("reference/matrices/"):
            continue  # one column per residue; described as a whole
        columns = pd.read_csv(full_run / name, nrows=0).columns.tolist()
        assert name in guide, f"{name} has no entry in the file guide"
        assert set(columns) == set(guide[name].columns), name


def test_guide_and_manifest_cover_exactly_the_files_written(full_run: Path) -> None:
    manifest = json.loads((full_run / "run.json").read_text())

    assert manifest["outputs"] == files_on_disk(full_run)
    assert written_files(full_run) == files_on_disk(full_run)


PER_CHAIN_FIGURES = {
    "overview_<chain>.png",
    "reference/heatmap_<chain>.png",
    "clusters/clusters_<chain>.png",
}


def test_readme_describes_every_output_that_was_written(full_run: Path) -> None:
    readme = (full_run / "README.md").read_text()

    for entry in file_guide("png"):
        if entry.display in PER_CHAIN_FIGURES:
            # Two chains share one figure; see the many-chain test below.
            assert f"`{entry.display}`" not in readme, entry.display
            continue
        assert f"`{entry.display}`" in readme, entry.display
        for column in entry.columns:
            assert f"`{column}`" in readme, (entry.display, column)
    assert "Start with `overview.png`" in readme


def test_more_than_four_chains_get_one_figure_per_chain(tmp_path: Path) -> None:
    pdb = copy_chain_a(tmp_path / "five_chains.pdb", "BCDE")
    output_dir = tmp_path / "out"
    exit_code = main(
        [str(pdb), "--output-dir", str(output_dir), "--reference-model", "2"]
        + ["--cluster-residues", "--cluster-min-size", "2"]
    )
    assert exit_code == 0

    figures = [name for name in files_on_disk(output_dir) if name.endswith(".png")]
    assert figures == sorted(
        f"{stem}_{chain}.png"
        for stem in ("overview", "reference/heatmap", "clusters/clusters")
        for chain in "ABCDE"
    )
    readme = (output_dir / "README.md").read_text()
    assert "Start with `overview_<chain>.png`" in readme
    for display in PER_CHAIN_FIGURES:
        assert f"`{display}`" in readme, display
    for display in ("overview.png", "reference/heatmap.png", "clusters/clusters.png"):
        assert f"`{display}`" not in readme, display
    assert written_files(output_dir) == files_on_disk(output_dir)


def test_readme_skips_outputs_of_analyses_that_did_not_run(tmp_path: Path) -> None:
    output_dir = tmp_path / "out"
    assert main([str(MINI_ENSEMBLE), "--output-dir", str(output_dir)]) == 0

    readme = (output_dir / "README.md").read_text()

    assert "`geometry/residues.csv`" in readme
    for absent in ("reference/", "clusters/", "range_clusters/", "models_by_chain"):
        assert absent not in readme, absent
    assert "Distance to reference" not in readme
    assert json.loads((output_dir / "run.json").read_text())["hdbscan"] is None


def test_readme_gives_the_hdbscan_settings_used(tmp_path: Path) -> None:
    output_dir = tmp_path / "out"
    assert main([str(MINI_ENSEMBLE), "--output-dir", str(output_dir), "--cluster-residues"]) == 0

    readme = (output_dir / "README.md").read_text()

    # Three models: the default size is the floor of 5 and min_samples the default 5.
    assert (
        "Per-residue clustering (HDBSCAN, min_cluster_size 5, min_samples 5, "
        "single cluster allowed)." in readme
    )


def test_manifest_records_inputs_parameters_and_versions(
    full_run: Path, two_chain_ensemble: Path
) -> None:
    manifest = json.loads((full_run / "run.json").read_text())

    assert manifest["input"] == {
        "pdb_file": str(two_chain_ensemble.resolve()),
        "sha256": sha256(two_chain_ensemble),
        "models": 3,
        "chains": ["A", "B"],
        "residues": 20,
    }
    assert manifest["reference"] == {
        "label": "two_chains.pdb model 2",
        "pdb_file": str(two_chain_ensemble.resolve()),
        "sha256": sha256(two_chain_ensemble),
    }
    parameters = manifest["parameters"]
    assert parameters["clustering"] == {
        "cluster_residues": True,
        "cluster_residue_ranges": ["2-5"],
        "min_cluster_size": 2,
        "min_samples": None,
    }
    # min_samples was not given: it defaults to 5, or min_cluster_size if smaller.
    assert manifest["hdbscan"] == {
        "min_cluster_size": 2,
        "min_samples": 2,
        "allow_single_cluster": True,
    }
    assert parameters["reference"]["pdb_model_id"] == "2"
    assert parameters["output"]["distance_matrices"] is True
    assert manifest["flexgeo2_version"] == metadata.version("FleXgeo2")
    assert manifest["environment"]["packages"]["FleXgeo2"] == metadata.version("FleXgeo2")
    assert datetime.fromisoformat(manifest["created_utc"]).tzinfo is not None


def test_readme_key_results_match_the_tables(full_run: Path) -> None:
    readme = (full_run / "README.md").read_text()
    residues = pd.read_csv(full_run / "geometry" / "residues.csv")
    clusters = pd.read_csv(full_run / "clusters" / "residues.csv")

    most_flexible = residues.nlargest(1, "dmax").iloc[0]
    assert f"| A {most_flexible['residue_label']} | {most_flexible['dmax']:.3f} |" in readme
    n_split = int((clusters["n_clusters"] >= 2).sum())
    assert f"{n_split} of {len(clusters)} residues split into two or more clusters" in readme
    assert "- Distance to reference: two_chains.pdb model 2." in readme


def test_readme_suggests_smaller_clusters_when_a_range_is_all_noise(tmp_path: Path) -> None:
    config = AnalysisConfig(
        pdb_file=MINI_ENSEMBLE,
        clustering=ClusteringConfig(cluster_residue_ranges=["2-5"], min_cluster_size=5),
        output=OutputConfig(write_files=False),
    )
    result = FlexGeo2App().run(config)
    # Three models with min_cluster_size 5 cannot form a cluster.
    assert result.residue_range_clustering.summary_df["n_clusters"].tolist() == [0]

    readme = render_readme(result, tmp_path, datetime.now(timezone.utc))

    assert "| A | 2-5 | 0 | 1.00 |" in readme
    assert "consider a smaller `--cluster-min-size`" in readme


def top_residues(table: pd.DataFrame, column: str) -> str:
    rows = table.nlargest(3, column)
    return ", ".join(
        f"{row['chain']} {row['residue_label']} ({row[column]:.3f})" for _, row in rows.iterrows()
    )


def test_terminal_summary_gives_headline_results_not_file_paths(capsys) -> None:
    # The working directory is tmp_path (conftest), so "out" is shown relative to it.
    arguments = [str(MINI_ENSEMBLE), "--output-dir", "out", "--reference-model", "1"]
    arguments += ["--cluster-residues", "--cluster-residue-range", "2-5"]
    assert main([*arguments, "--cluster-min-size", "2"]) == 0

    lines = capsys.readouterr().out.splitlines()
    out = Path("out")
    residues = pd.read_csv(out / "geometry" / "residues.csv")
    distances = pd.read_csv(out / "reference" / "residues.csv")
    clusters = pd.read_csv(out / "clusters" / "residues.csv")
    [range_row] = pd.read_csv(out / "range_clusters" / "ranges.csv").to_dict("records")
    n_split = int((clusters["n_clusters"] >= 2).sum())
    hint = " (try a smaller --cluster-min-size)" if range_row["n_clusters"] == 0 else ""

    assert lines == [
        "FleXgeo2 analysed mini_ensemble.pdb: 3 models, 1 chain (A), 10 residues.",
        "",
        f"Most flexible residues (dmax): {top_residues(residues, 'dmax')}",
        f"Furthest from input model 1 (mean distance): {top_residues(distances, 'distance_mean')}",
        f"Per-residue clustering: {n_split} of 10 residues split into two or more clusters",
        f"Range clustering A 2-5: {range_row['n_clusters']} "
        f"cluster{'' if range_row['n_clusters'] == 1 else 's'}, "
        f"{range_row['noise_fraction']:.0%} noise{hint}",
        "",
        f"Results: out/ ({len(files_on_disk(out))} files); start with README.md",
    ]


def test_terminal_summary_mentions_only_the_analyses_that_ran(capsys) -> None:
    assert main([str(MINI_ENSEMBLE), "--output-dir", "out"]) == 0

    output = capsys.readouterr().out
    assert "Most flexible residues (dmax):" in output
    for absent in ("Furthest from", "Per-residue clustering", "Range clustering"):
        assert absent not in output


def test_terminal_summary_shows_absolute_path_outside_working_dir(
    tmp_path_factory: pytest.TempPathFactory, capsys
) -> None:
    output_dir = tmp_path_factory.mktemp("elsewhere") / "out"
    assert main([str(MINI_ENSEMBLE), "--output-dir", str(output_dir)]) == 0

    last_line = capsys.readouterr().out.splitlines()[-1]
    assert last_line.startswith(f"Results: {output_dir.resolve()}/ (")


def test_terminal_summary_hints_at_all_noise_ranges_and_skips_missing_files() -> None:
    config = AnalysisConfig(
        pdb_file=MINI_ENSEMBLE,
        clustering=ClusteringConfig(cluster_residue_ranges=["2-5"], min_cluster_size=5),
        output=OutputConfig(write_files=False),
    )
    result = FlexGeo2App().run(config)
    # Add a range that did cluster: it must not get the hint.
    ranges = result.residue_range_clustering.summary_df
    clustered = ranges.assign(range_label="6-9", n_clusters=2, noise_fraction=0.25)
    single = ranges.assign(range_label="7-10", n_clusters=1, noise_fraction=0.4)
    result.residue_range_clustering.summary_df = pd.concat([ranges, clustered, single])

    lines = render_terminal_summary(result).splitlines()

    assert "Range clustering A 2-5: 0 clusters, 100% noise (try a smaller --cluster-min-size)" in (
        lines
    )
    assert "Range clustering A 6-9: 2 clusters, 25% noise" in lines
    assert "Range clustering A 7-10: 1 cluster, 40% noise" in lines
    assert not any(line.startswith("Results:") for line in lines)


@pytest.mark.parametrize(("plot_format", "magic"), [("pdf", b"%PDF"), ("svg", b"<?xml")])
def test_vector_plot_formats_are_written_and_described(plot_format: str, magic: bytes) -> None:
    arguments = [str(MINI_ENSEMBLE), "--output-dir", "out", "--reference-model", "1"]
    arguments += ["--cluster-residues", "--cluster-min-size", "2", "--plot-residues", "3"]
    assert main([*arguments, "--plot-format", plot_format]) == 0

    out = Path("out")
    figures = [name for name in files_on_disk(out) if not name.endswith((".csv", ".md", ".json"))]
    assert figures == [
        f"clusters/clusters.{plot_format}",
        f"overview.{plot_format}",
        f"reference/heatmap.{plot_format}",
        f"residue_plots/A_0003_{pd.read_csv(out / 'geometry' / 'residues.csv')['name'][2]}"
        f".{plot_format}",
    ]
    for name in figures:
        assert (out / name).read_bytes().startswith(magic), name

    readme = (out / "README.md").read_text()
    assert f"Start with `overview.{plot_format}`" in readme
    assert f"### `reference/heatmap.{plot_format}`" in readme
    assert ".png" not in readme
    assert (
        json.loads((out / "run.json").read_text())["parameters"]["output"]["plot_format"]
        == plot_format
    )


def test_vector_figures_keep_text_editable() -> None:
    assert main([str(MINI_ENSEMBLE), "--output-dir", "pdf", "--plot-format", "pdf"]) == 0
    assert main([str(MINI_ENSEMBLE), "--output-dir", "svg", "--plot-format", "svg"]) == 0

    pdf = Path("pdf/overview.pdf").read_bytes()
    assert b"/FontFile2" in pdf  # embedded TrueType font
    assert b"/Type3" not in pdf
    svg = Path("svg/overview.svg").read_text()
    assert "<text" in svg  # text kept as text, not converted to paths
