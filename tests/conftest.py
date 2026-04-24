from __future__ import annotations

import pandas as pd
import pytest


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
        f"{name}{int(order)}" for order, name in zip(df["order"], df["name"])
    ]
    return df.sort_values(["chain", "model", "order"]).reset_index(drop=True)
