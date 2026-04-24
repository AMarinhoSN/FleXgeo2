from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from flexgeo2.clustering import ClusteringService


def test_parse_residue_range_accepts_start_and_end() -> None:
    assert ClusteringService.parse_residue_range(" 45-54 ") == (45, 54)


def test_parse_residue_range_rejects_reversed_range() -> None:
    with pytest.raises(ValueError, match="End must be greater than or equal to start"):
        ClusteringService.parse_residue_range("54-45")


def test_parse_residue_range_rejects_non_integer_range() -> None:
    with pytest.raises(ValueError, match="Start and end must be integers"):
        ClusteringService.parse_residue_range("A-45")


def test_compute_pca_projection_returns_two_columns() -> None:
    matrix = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]])

    projection = ClusteringService.compute_pca_projection(matrix)

    assert projection.shape == (3, 2)
    assert projection[:, 1].tolist() == pytest.approx([0.0, 0.0, 0.0])


def test_cluster_residues_marks_small_groups_as_noise(
    monkeypatch: pytest.MonkeyPatch,
    normalized_geometry_df: pd.DataFrame,
) -> None:
    monkeypatch.setitem(sys.modules, "hdbscan", SimpleNamespace())

    assignments_df, summary_df = ClusteringService().cluster_residues(
        raw_df=normalized_geometry_df,
        min_cluster_size=5,
        min_samples=None,
    )

    assert set(assignments_df["cluster"]) == {-1}
    assert set(assignments_df["cluster_probability"]) == {0.0}
    assert set(summary_df["n_clusters"]) == {0}
    assert set(summary_df["noise_fraction"]) == {1.0}


def test_cluster_residue_ranges_rejects_incomplete_ranges(
    monkeypatch: pytest.MonkeyPatch,
    normalized_geometry_df: pd.DataFrame,
) -> None:
    monkeypatch.setitem(sys.modules, "hdbscan", SimpleNamespace())

    with pytest.raises(ValueError, match="is incomplete"):
        ClusteringService().cluster_residue_ranges(
            raw_df=normalized_geometry_df,
            range_texts=["1-2"],
            min_cluster_size=5,
            min_samples=None,
        )


def test_cluster_residue_ranges_marks_small_groups_as_noise(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(sys.modules, "hdbscan", SimpleNamespace())
    df = pd.DataFrame(
        [
            {
                "model": "1",
                "chain": "A",
                "order": 1,
                "name": "ALA",
                "residue_label": "ALA1",
                "curvature": 0.1,
                "torsion": 0.2,
            },
            {
                "model": "2",
                "chain": "A",
                "order": 1,
                "name": "ALA",
                "residue_label": "ALA1",
                "curvature": 0.2,
                "torsion": 0.3,
            },
        ]
    )

    assignments_df, summary_df = ClusteringService().cluster_residue_ranges(
        raw_df=df,
        range_texts=["1-1"],
        min_cluster_size=5,
        min_samples=None,
    )

    assert assignments_df["cluster"].tolist() == [-1, -1]
    assert assignments_df["cluster_probability"].tolist() == [0.0, 0.0]
    assert assignments_df["range_label"].tolist() == ["1-1", "1-1"]
    assert summary_df.loc[0, "n_clusters"] == 0
    assert summary_df.loc[0, "noise_fraction"] == 1.0
