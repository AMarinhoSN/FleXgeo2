"""End-to-end tests that run the real Melodia on a tiny three-model ensemble."""

from __future__ import annotations

import json
from dataclasses import fields, replace
from pathlib import Path

import pandas as pd
import pytest

from flexgeo2 import (
    AnalysisConfig,
    ClusteringConfig,
    FlexGeo2App,
    OutputConfig,
    ReferenceConfig,
    save_figure,
)
from flexgeo2.cli.main import main
from flexgeo2.geometry import GeometryService
from flexgeo2.models import AnalysisResult, OutputArtifacts
from flexgeo2.outputs import OutputDirectoryNotEmptyError

pytest.importorskip("melodia_py")


def test_load_structure_reports_pdb_model_serials(mini_ensemble_pdb: Path) -> None:
    df = GeometryService().load_structure(mini_ensemble_pdb)

    assert sorted(df["model"].unique().tolist()) == [1, 2, 3]
    assert sorted(df["order"].unique().tolist()) == list(range(1, 11))


def test_models_without_model_records_start_at_one(mini_ensemble_pdb: Path, tmp_path: Path) -> None:
    first_model_atoms = []
    for line in mini_ensemble_pdb.read_text().splitlines(keepends=True):
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            first_model_atoms.append(line)
    single_model_pdb = tmp_path / "single.pdb"
    single_model_pdb.write_text("".join(first_model_atoms) + "END\n")

    df = GeometryService().load_structure(single_model_pdb)

    assert df["model"].unique().tolist() == [1]


def test_describe_structure_lists_models_and_residues(mini_ensemble_pdb: Path) -> None:
    service = GeometryService()

    info = service.describe_structure(service.parse_structure(mini_ensemble_pdb))

    assert info.model_ids == [1, 2, 3]
    assert info.residues_by_chain == {"A": set(range(1, 11))}


def test_cli_reference_model_uses_pdb_model_numbering(
    mini_ensemble_pdb: Path, tmp_path: Path
) -> None:
    output_dir = tmp_path / "out"

    exit_code = main(
        [
            str(mini_ensemble_pdb),
            "--output-dir",
            str(output_dir),
            "--reference-model",
            "1",
            "--distance-matrices",
        ]
    )

    assert exit_code == 0
    long_df = pd.read_csv(output_dir / "reference" / "distances.csv")
    model_one = long_df[long_df["model"] == 1]
    assert model_one["distance_to_reference"].tolist() == pytest.approx([0.0] * 10)
    assert (long_df[long_df["model"] == 2]["distance_to_reference"] > 0).any()

    matrix = pd.read_csv(output_dir / "reference" / "matrices" / "A.csv", index_col=0)
    orders = [int("".join(ch for ch in label if ch.isdigit())) for label in matrix.columns]
    assert orders == list(range(1, 11))
    assert matrix.index.tolist() == [1, 2, 3]

    assert (output_dir / "overview.png").is_file()
    assert (output_dir / "reference" / "heatmap.png").is_file()


def test_cli_reference_pdb_compares_against_external_model(
    mini_ensemble_pdb: Path, tmp_path: Path
) -> None:
    output_dir = tmp_path / "out"

    exit_code = main(
        [
            str(mini_ensemble_pdb),
            "--output-dir",
            str(output_dir),
            "--reference-pdb",
            str(mini_ensemble_pdb),
            "--reference-pdb-model",
            "2",
        ]
    )

    assert exit_code == 0
    long_df = pd.read_csv(output_dir / "reference" / "distances.csv")
    assert set(long_df["reference_label"]) == {"mini_ensemble.pdb model 2"}
    distances_by_model = long_df.groupby("model")["distance_to_reference"].max()
    assert distances_by_model.loc[2] == pytest.approx(0.0)
    assert distances_by_model.loc[1] > 0.0
    assert distances_by_model.loc[3] > 0.0


def test_cli_runs_both_clustering_modes(mini_ensemble_pdb: Path, tmp_path: Path) -> None:
    output_dir = tmp_path / "out"

    exit_code = main(
        [
            str(mini_ensemble_pdb),
            "--output-dir",
            str(output_dir),
            "--cluster-residues",
            "--cluster-residue-range",
            "2-5",
            "--cluster-min-size",
            "2",
        ]
    )

    assert exit_code == 0
    residue_summary = pd.read_csv(output_dir / "clusters" / "residues.csv")
    assert residue_summary["order"].tolist() == list(range(1, 11))
    assert set(residue_summary["models"]) == {3}
    assert residue_summary["noise_fraction"].between(0.0, 1.0).all()
    assert (output_dir / "clusters" / "clusters.png").is_file()

    # Per-model answers are written by default.
    assignments = pd.read_csv(output_dir / "clusters" / "assignments.csv")
    assert len(assignments) == 3 * 10
    assert set(assignments["model"]) == {1, 2, 3}
    keys = ["chain", "model", "order", "name", "residue_label"]
    assert assignments.columns.tolist() == [
        *keys,
        "curvature",
        "torsion",
        "cluster",
        "cluster_probability",
    ]
    descriptors = pd.read_csv(output_dir / "geometry" / "descriptors.csv", nrows=0)
    assert descriptors.columns.tolist() == [
        *keys,
        "curvature",
        "torsion",
        "arc_length",
        "writhing",
        "phi",
        "psi",
    ]

    range_summary = pd.read_csv(output_dir / "range_clusters" / "ranges.csv")
    assert range_summary[["chain", "range_label", "residues", "models"]].to_dict("records") == [
        {"chain": "A", "range_label": "2-5", "residues": 4, "models": 3}
    ]
    assert (output_dir / "range_clusters" / "A_2-5.png").is_file()

    range_assignments = pd.read_csv(output_dir / "range_clusters" / "assignments.csv")
    assert range_assignments["model"].tolist() == [1, 2, 3]
    assert set(range_assignments["range_label"]) == {"2-5"}


def test_cli_plots_chosen_residues(mini_ensemble_pdb: Path, tmp_path: Path, capsys) -> None:
    output_dir = tmp_path / "out"

    exit_code = main(
        [
            str(mini_ensemble_pdb),
            "--output-dir",
            str(output_dir),
            "--cluster-residues",
            "--cluster-min-size",
            "2",
            "--reference-model",
            "1",
            "--plot-residues",
            "3,A:5-6",
        ]
    )

    assert exit_code == 0
    plots = sorted(path.name for path in (output_dir / "residue_plots").iterdir())
    residues = pd.read_csv(output_dir / "geometry" / "residues.csv").set_index("order")
    assert plots == [f"A_{order:04d}_{residues.loc[order, 'name']}.png" for order in (3, 5, 6)]
    assert "Residue plots: 3 in out/residue_plots/" in capsys.readouterr().out


def test_cli_rejects_residue_plots_outside_the_structure_before_melodia(
    monkeypatch: pytest.MonkeyPatch, mini_ensemble_pdb: Path, tmp_path: Path, capsys
) -> None:
    def fail_compute_geometry(*args, **kwargs):
        raise AssertionError("Melodia should not run when the input is invalid.")

    monkeypatch.setattr(GeometryService, "compute_geometry", fail_compute_geometry)

    exit_code = main(
        [str(mini_ensemble_pdb), "--output-dir", str(tmp_path / "out"), "--plot-residues", "9-12"]
    )

    assert exit_code == 1
    assert "residue(s) 11, 12 not found in chain(s) A" in capsys.readouterr().err


def saving_config(pdb_file: Path, output_dir: Path | None = None) -> AnalysisConfig:
    return AnalysisConfig(
        pdb_file=pdb_file,
        reference=ReferenceConfig(model_id="1"),
        clustering=ClusteringConfig(cluster_residues=True, min_cluster_size=2),
        output=OutputConfig(output_dir=output_dir, plot_format="svg", plot_residues=["3"]),
    )


def output_files(output_dir: Path) -> list[str]:
    return sorted(
        p.relative_to(output_dir).as_posix() for p in output_dir.rglob("*") if p.is_file()
    )


def test_library_run_writes_nothing_by_default(mini_ensemble_pdb: Path, tmp_path: Path) -> None:
    # The autouse fixture runs each test in tmp_path, the old default folder's parent.
    result = FlexGeo2App().run(AnalysisConfig(pdb_file=mini_ensemble_pdb))

    assert result.outputs == OutputArtifacts()
    assert list(tmp_path.iterdir()) == []


def test_save_writes_what_a_run_with_an_output_dir_writes(
    mini_ensemble_pdb: Path, tmp_path: Path
) -> None:
    written = tmp_path / "written"
    FlexGeo2App().run(saving_config(mini_ensemble_pdb, written))
    result = FlexGeo2App().run(saving_config(mini_ensemble_pdb))
    saved = tmp_path / "saved"

    artifacts = result.save(saved)

    assert output_files(saved) == output_files(written)
    assert "residue_plots/A_0003_ILE.svg" in output_files(saved)  # the run's output settings
    assert result.outputs is artifacts
    assert artifacts.overview_plot == saved.resolve() / "overview.svg"
    manifest = json.loads((saved / "run.json").read_text())
    assert manifest["parameters"]["output"]["output_dir"] == str(saved)
    assert result.config.output.output_dir == saved


def test_save_refuses_a_folder_with_files_unless_overwriting(
    mini_ensemble_pdb: Path, tmp_path: Path
) -> None:
    result = FlexGeo2App().run(saving_config(mini_ensemble_pdb))
    folder = tmp_path / "out"
    folder.mkdir()
    (folder / "notes.txt").write_text("keep me")

    with pytest.raises(OutputDirectoryNotEmptyError):
        result.save(folder)
    assert result.config.output.output_dir is None  # unchanged by the refused save

    result.save(folder, overwrite=True)
    result.save(folder, overwrite=True)  # replaces its own earlier outputs

    assert (folder / "notes.txt").read_text() == "keep me"
    assert (folder / "README.md").is_file()


def test_run_and_save_leave_the_callers_matplotlib_settings_alone(
    mini_ensemble_pdb: Path, tmp_path: Path
) -> None:
    import matplotlib.pyplot as plt

    with plt.rc_context({"axes.titleweight": "normal", "pdf.fonttype": 3}):
        user_settings = dict(plt.rcParams)
        result = FlexGeo2App().run(saving_config(mini_ensemble_pdb, tmp_path / "run"))
        assert dict(plt.rcParams) == user_settings
        result.save(tmp_path / "saved")
        assert dict(plt.rcParams) == user_settings


@pytest.fixture(scope="module")
def full_run(mini_ensemble_pdb: Path, tmp_path_factory: pytest.TempPathFactory):
    """A run with every analysis and PNG figures, and the folder it wrote."""
    output_dir = tmp_path_factory.mktemp("full_run")
    config = AnalysisConfig(
        pdb_file=mini_ensemble_pdb,
        reference=ReferenceConfig(model_id="1"),
        clustering=ClusteringConfig(
            cluster_residues=True, cluster_residue_ranges=["2-4"], min_cluster_size=2
        ),
        output=OutputConfig(output_dir=output_dir, plot_residues=["3"]),
        max_models_in_plot=2,  # of 3, so the overview shows the run's settings are used
    )
    return FlexGeo2App().run(config), output_dir


@pytest.mark.parametrize(
    ("draw", "written", "dpi"),
    [
        (lambda result: result.plot_overview(), "overview.png", 300),
        (lambda result: result.plot_distance_heatmap(), "reference/heatmap.png", 300),
        (lambda result: result.plot_cluster_map(), "clusters/clusters.png", 300),
        (lambda result: result.plot_residue(3), "residue_plots/A_0003_ILE.png", 250),
        (lambda result: result.plot_residue("A:3"), "residue_plots/A_0003_ILE.png", 250),
        (lambda result: result.plot_residue_range("2-4"), "range_clusters/A_2-4.png", 250),
        (lambda result: result.plot_overview(chains="A"), "overview.png", 300),
    ],
)
def test_plot_methods_draw_the_figures_a_run_writes(
    full_run, tmp_path: Path, draw, written: str, dpi: int
) -> None:
    import matplotlib.pyplot as plt

    result, output_dir = full_run
    figure = draw(result)

    save_figure(figure, tmp_path / "figure.png", dpi=dpi)
    plt.close(figure)

    assert (tmp_path / "figure.png").read_bytes() == (output_dir / written).read_bytes()


@pytest.mark.parametrize(
    ("draw", "chains"),
    [
        (lambda result: result.plot_overview(), ["A", "B"]),
        (lambda result: result.plot_overview(chains="B"), ["B"]),
        (lambda result: result.plot_distance_heatmap(chains=["B"]), ["B"]),
        (lambda result: result.plot_cluster_map(chains=["A"]), ["A"]),
    ],
)
def test_plot_methods_draw_the_chosen_chains(full_run, draw, chains: list[str]) -> None:
    import matplotlib.pyplot as plt

    result, _ = full_run
    two_chains = with_copy_as_chain_b(result)

    figure = draw(two_chains)

    titles = " ".join(axis.get_title() for axis in figure.axes)
    plt.close(figure)
    assert [chain for chain in ("A", "B") if f"Chain {chain}" in titles] == chains


def with_copy_as_chain_b(result: AnalysisResult) -> AnalysisResult:
    """The result with every table's chain A rows repeated as chain B."""

    def doubled(df: pd.DataFrame) -> pd.DataFrame:
        return pd.concat([df, df.assign(chain="B")], ignore_index=True)

    def doubled_tables(analysis):
        return replace(
            analysis,
            **{
                field.name: doubled(getattr(analysis, field.name))
                for field in fields(analysis)
                if field.name.endswith("_df")
            },
        )

    copy = doubled_tables(result)
    copy.distance_result = doubled_tables(result.distance_result)
    copy.residue_clustering = doubled_tables(result.residue_clustering)
    return copy


@pytest.mark.parametrize(
    ("draw", "message"),
    [
        (lambda result: result.plot_overview(chains=["B"]), r"Chain\(s\) B not in this result"),
        (lambda result: result.plot_residue("2-3"), "matches 2 residues"),
        (lambda result: result.plot_residue("99"), "99 not found"),
        (lambda result: result.plot_residue_range("3-5"), "is not a clustered range.*A:2-4"),
        (lambda result: result.plot_residue_range("2-4,3-5"), "Give one residue range"),
    ],
)
def test_plot_methods_explain_bad_selections(full_run, draw, message: str) -> None:
    result, _ = full_run

    with pytest.raises(ValueError, match=message):
        draw(result)


@pytest.mark.parametrize(
    ("draw", "message"),
    [
        (lambda result: result.plot_distance_heatmap(), "needs a run with a reference"),
        (lambda result: result.plot_cluster_map(), "needs a run with per-residue clustering"),
        (lambda result: result.plot_residue_range("2-4"), "needs a run with range clustering"),
    ],
)
def test_plot_methods_name_the_analysis_they_need(
    mini_ensemble_pdb: Path, draw, message: str
) -> None:
    result = FlexGeo2App().run(AnalysisConfig(pdb_file=mini_ensemble_pdb))

    with pytest.raises(ValueError, match=message):
        draw(result)


def test_save_figure_embeds_editable_fonts_whatever_the_session_uses(
    full_run, tmp_path: Path
) -> None:
    import matplotlib.pyplot as plt

    result, _ = full_run
    figure = result.plot_residue(3)

    with plt.rc_context({"pdf.fonttype": 3}):
        save_figure(figure, tmp_path / "residue.pdf")
        figure.savefig(tmp_path / "session.pdf")

    assert plt.fignum_exists(figure.number)  # left open for the caller
    plt.close(figure)
    assert b"/Type3" not in (tmp_path / "residue.pdf").read_bytes()
    assert b"/Type3" in (tmp_path / "session.pdf").read_bytes()
