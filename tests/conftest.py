from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(autouse=True)
def _run_in_tmp_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Run each test from its own temporary directory.

    Relative default paths (e.g. OutputConfig's ``results``) then land in ``tmp_path``
    instead of the directory pytest was launched from.
    """
    monkeypatch.chdir(tmp_path)


@pytest.fixture
def mini_ensemble_pdb() -> Path:
    """Three NMR models (PDB MODEL 1-3) of residues 1-10 from PDB 2LJ5, chain A."""
    return DATA_DIR / "mini_ensemble.pdb"


@pytest.fixture
def raw_geometry_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model": "2",
                "chain": "B",
                "order": 1,
                "name": "GLY",
                "curvature": 2.0,
                "torsion": 2.5,
            },
            {
                "model": "1",
                "chain": "A",
                "order": 2,
                "name": "GLY",
                "curvature": 0.4,
                "torsion": 0.5,
            },
            {
                "model": "1",
                "chain": "A",
                "order": 1,
                "name": "ALA",
                "curvature": 0.1,
                "torsion": 0.2,
            },
            {
                "model": "2",
                "chain": "A",
                "order": 1,
                "name": "ALA",
                "curvature": 0.2,
                "torsion": 0.3,
            },
            {
                "model": "2",
                "chain": "A",
                "order": 2,
                "name": "GLY",
                "curvature": 0.6,
                "torsion": 0.9,
            },
            {
                "model": "1",
                "chain": "B",
                "order": 1,
                "name": "GLY",
                "curvature": 1.5,
                "torsion": 2.0,
            },
        ]
    )


@pytest.fixture
def normalized_geometry_df(raw_geometry_df: pd.DataFrame) -> pd.DataFrame:
    df = raw_geometry_df.copy()
    df["residue_label"] = [
        f"{name}{int(order)}" for order, name in zip(df["order"], df["name"], strict=False)
    ]
    return df.sort_values(["chain", "model", "order"]).reset_index(drop=True)
