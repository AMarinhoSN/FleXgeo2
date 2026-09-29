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
            "--output-verbose",
        ]
    )

    assert exit_code == 0
    long_df = pd.read_csv(output_dir / "distance_to_reference_long.csv")
    model_one = long_df[long_df["model"] == 1]
    assert model_one["distance_to_reference"].tolist() == pytest.approx([0.0] * 10)
    assert (long_df[long_df["model"] == 2]["distance_to_reference"] > 0).any()

    matrix = pd.read_csv(output_dir / "distance_matrices" / "A_distance_matrix.csv", index_col=0)
    orders = [int("".join(ch for ch in label if ch.isdigit())) for label in matrix.columns]
    assert orders == list(range(1, 11))
    assert matrix.index.tolist() == [1, 2, 3]

    assert (output_dir / "plots" / "ensemble_overview.png").is_file()
    assert (output_dir / "plots" / "distance_to_reference_heatmap.png").is_file()
