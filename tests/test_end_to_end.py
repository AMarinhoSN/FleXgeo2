"""End-to-end tests that run the real Melodia on a tiny three-model ensemble."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from flexgeo2.cli.main import main
from flexgeo2.geometry import GeometryService

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
