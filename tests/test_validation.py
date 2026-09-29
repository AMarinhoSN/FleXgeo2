from __future__ import annotations

from pathlib import Path

import pytest

from flexgeo2.config import AnalysisConfig, ClusteringConfig, ReferenceConfig
from flexgeo2.geometry import StructureInfo
from flexgeo2.validation import validate_against_structure, validate_config


@pytest.fixture
def pdb_file(tmp_path: Path) -> Path:
    path = tmp_path / "ensemble.pdb"
    path.write_text("HEADER test\n")
    return path


@pytest.fixture
def structure_info() -> StructureInfo:
    return StructureInfo(
        model_ids=[1, 2, 3],
        residues_by_chain={"A": set(range(1, 11)), "B": set(range(5, 8))},
    )


def test_validate_config_accepts_defaults(pdb_file: Path) -> None:
    validate_config(AnalysisConfig(pdb_file=pdb_file))


def test_validate_config_rejects_missing_input(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="Input PDB file not found"):
        validate_config(AnalysisConfig(pdb_file=tmp_path / "missing.pdb"))


def test_validate_config_rejects_missing_reference_pdb(pdb_file: Path, tmp_path: Path) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        reference=ReferenceConfig(pdb_file=tmp_path / "missing_ref.pdb"),
    )

    with pytest.raises(FileNotFoundError, match="Reference PDB file not found"):
        validate_config(config)


def test_validate_config_rejects_reference_pdb_model_without_pdb(pdb_file: Path) -> None:
    config = AnalysisConfig(pdb_file=pdb_file, reference=ReferenceConfig(pdb_model_id="1"))

    with pytest.raises(ValueError, match="--reference-pdb-model requires --reference-pdb"):
        validate_config(config)


def test_validate_config_rejects_both_reference_sources(pdb_file: Path) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        reference=ReferenceConfig(model_id="1", pdb_file=pdb_file),
    )

    with pytest.raises(ValueError, match="not both"):
        validate_config(config)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"n_jobs": 0}, "n_jobs"),
        ({"max_models_in_plot": -1}, "max_models_in_plot"),
        ({"dmax_outlier_fraction": 1.0}, "dmax_outlier_fraction"),
        (
            {"clustering": ClusteringConfig(cluster_residues=True, min_cluster_size=1)},
            "min_cluster_size",
        ),
        (
            {"clustering": ClusteringConfig(cluster_residues=True, min_samples=0)},
            "min_samples",
        ),
        (
            {"clustering": ClusteringConfig(cluster_residue_ranges=["10-5"])},
            "End must be greater",
        ),
    ],
)
def test_validate_config_rejects_invalid_options(
    pdb_file: Path, overrides: dict, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        validate_config(AnalysisConfig(pdb_file=pdb_file, **overrides))


def test_validate_against_structure_accepts_valid_selection(
    pdb_file: Path, structure_info: StructureInfo
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        chains=["A"],
        reference=ReferenceConfig(model_id="2"),
        clustering=ClusteringConfig(cluster_residue_ranges=["3-6"]),
    )

    validate_against_structure(config, structure_info)


def test_validate_against_structure_rejects_unknown_chain(
    pdb_file: Path, structure_info: StructureInfo
) -> None:
    config = AnalysisConfig(pdb_file=pdb_file, chains=["A", "Z"])

    with pytest.raises(ValueError, match=r"Chain\(s\) not found in the input: Z"):
        validate_against_structure(config, structure_info)


def test_validate_against_structure_rejects_unknown_reference_model(
    pdb_file: Path, structure_info: StructureInfo
) -> None:
    config = AnalysisConfig(pdb_file=pdb_file, reference=ReferenceConfig(model_id="0"))

    with pytest.raises(ValueError, match="Reference model '0' was not found"):
        validate_against_structure(config, structure_info)


def test_validate_against_structure_checks_reference_pdb_models(
    pdb_file: Path, structure_info: StructureInfo
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        reference=ReferenceConfig(pdb_file=pdb_file, pdb_model_id="7"),
    )
    reference_info = StructureInfo(model_ids=[1], residues_by_chain={"A": {1}})

    with pytest.raises(ValueError, match="not found in the reference PDB"):
        validate_against_structure(config, structure_info, reference_info)


def test_validate_against_structure_rejects_partial_range(
    pdb_file: Path, structure_info: StructureInfo
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        clustering=ClusteringConfig(cluster_residue_ranges=["8-12"]),
    )

    with pytest.raises(ValueError, match=r"8-12 on chain 'A' is incomplete.*11, 12"):
        validate_against_structure(config, structure_info)


def test_validate_against_structure_rejects_range_outside_selected_chains(
    pdb_file: Path, structure_info: StructureInfo
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        chains=["A"],
        clustering=ClusteringConfig(cluster_residue_ranges=["50-60"]),
    )

    with pytest.raises(ValueError, match="does not match any residues"):
        validate_against_structure(config, structure_info)
