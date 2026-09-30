from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from flexgeo2.clustering import ClusteringService, clustering_parameters


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


@pytest.mark.filterwarnings("error")
def test_compute_pca_projection_emits_no_floating_point_warnings() -> None:
    # A 40 x 30 matrix made matmul raise spurious divide/overflow/invalid warnings with
    # Apple Accelerate BLAS (numpy 2.2). On other BLAS builds this passes trivially.
    matrix = np.random.default_rng(0).normal(size=(40, 30))

    projection = ClusteringService.compute_pca_projection(matrix)

    centered = matrix - matrix.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    reference = np.column_stack(
        [np.einsum("ij,j->i", centered, vh[0]), np.einsum("ij,j->i", centered, vh[1])]
    )
    assert projection == pytest.approx(reference)


@pytest.mark.parametrize("bad_value", [np.nan, np.inf])
def test_compute_pca_projection_rejects_non_finite_input(bad_value: float) -> None:
    matrix = np.ones((5, 4))
    matrix[2, 1] = bad_value

    with pytest.raises(ValueError, match="requires finite"):
        ClusteringService.compute_pca_projection(matrix)


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


def _blob_frame(
    centers_by_order: dict[int, list[tuple[float, float]]],
    points_per_blob: int = 20,
    spread: float = 0.01,
    seed: int = 0,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Build chain A geometry where each model sits in one blob per residue.

    Model ``i`` belongs to blob ``i // points_per_blob`` at every residue, so the
    returned mapping gives the expected group of each model.
    """
    rng = np.random.default_rng(seed)
    n_blobs = max(len(centers) for centers in centers_by_order.values())
    n_models = n_blobs * points_per_blob
    rows = []
    for order, centers in centers_by_order.items():
        for model_index in range(n_models):
            blob = (model_index // points_per_blob) % len(centers)
            curvature, torsion = centers[blob] + rng.normal(0.0, spread, size=2)
            rows.append(
                {
                    "model": model_index + 1,
                    "chain": "A",
                    "order": order,
                    "name": "ALA",
                    "residue_label": f"ALA{order}",
                    "curvature": float(curvature),
                    "torsion": float(torsion),
                }
            )
    expected_groups = {
        model_index + 1: model_index // points_per_blob for model_index in range(n_models)
    }
    return pd.DataFrame(rows), expected_groups


def _assert_labels_match_groups(labels: pd.Series, groups: pd.Series) -> None:
    """Each true group maps to exactly one cluster label, and labels are not shared."""
    group_to_labels = labels.groupby(groups).unique()
    assert all(len(group_labels) == 1 for group_labels in group_to_labels)
    assigned = [int(group_labels[0]) for group_labels in group_to_labels]
    assert -1 not in assigned
    assert len(set(assigned)) == len(assigned)


def test_cluster_residues_recovers_separated_blobs_with_hdbscan() -> None:
    raw_df, _ = _blob_frame(
        {
            1: [(0.2, -0.5), (0.8, 0.5)],
            2: [(0.1, 0.1), (0.5, -0.5), (0.9, 0.4)],
        }
    )

    assignments_df, summary_df = ClusteringService().cluster_residues(
        raw_df=raw_df,
        min_cluster_size=5,
        min_samples=None,
    )

    summary = summary_df.set_index("order")
    assert summary.loc[1, "n_clusters"] == 2
    assert summary.loc[2, "n_clusters"] == 3
    assert summary["noise_fraction"].tolist() == [0.0, 0.0]
    assert summary["models"].tolist() == [60, 60]

    assert len(assignments_df) == len(raw_df)
    for order, n_blobs in ((1, 2), (2, 3)):
        residue_df = assignments_df[assignments_df["order"] == order]
        true_group = (residue_df["model"] - 1) // 20 % n_blobs
        _assert_labels_match_groups(residue_df["cluster"], true_group)

    probabilities = assignments_df["cluster_probability"]
    assert probabilities.between(0.0, 1.0).all()
    assert (probabilities > 0.0).all()


def test_cluster_residues_keeps_only_row_keys_features_and_labels() -> None:
    raw_df, _ = _blob_frame({1: [(0.2, -0.5), (0.8, 0.5)]})
    raw_df = raw_df.assign(arc_length=3.8, writhing=0.0, phi=-60.0, psi=-45.0)

    for min_cluster_size in (5, 1000):  # clustered, and too few models to cluster
        assignments_df, _ = ClusteringService().cluster_residues(
            raw_df=raw_df, min_cluster_size=min_cluster_size, min_samples=None
        )
        assert assignments_df.columns.tolist() == [
            "chain",
            "model",
            "order",
            "name",
            "residue_label",
            "curvature",
            "torsion",
            "cluster",
            "cluster_probability",
        ]


def test_cluster_residues_passes_hdbscan_parameters(monkeypatch: pytest.MonkeyPatch) -> None:
    import hdbscan

    created = []
    real_hdbscan = hdbscan.HDBSCAN

    def recording_hdbscan(**kwargs):
        created.append(kwargs)
        return real_hdbscan(**kwargs)

    monkeypatch.setattr(hdbscan, "HDBSCAN", recording_hdbscan)
    raw_df, _ = _blob_frame({1: [(0.2, -0.5), (0.8, 0.5)]})

    ClusteringService().cluster_residues(raw_df=raw_df, min_cluster_size=7, min_samples=3)

    assert created == [{"min_cluster_size": 7, "min_samples": 3, "allow_single_cluster": True}]


@pytest.mark.parametrize(
    "centers",
    [
        # Groups differ only in curvature.
        {1: [(0.2, 0.3), (0.8, 0.3)]},
        # Groups differ only in torsion.
        {1: [(0.3, 0.2), (0.3, 0.8)]},
    ],
    ids=["curvature-only", "torsion-only"],
)
def test_cluster_residues_uses_curvature_and_torsion(centers: dict) -> None:
    raw_df, expected_groups = _blob_frame(centers)

    assignments_df, summary_df = ClusteringService().cluster_residues(
        raw_df=raw_df,
        min_cluster_size=5,
        min_samples=None,
    )

    assert summary_df.iloc[0]["n_clusters"] == 2
    _assert_labels_match_groups(
        assignments_df["cluster"], assignments_df["model"].map(expected_groups)
    )


def test_cluster_residue_ranges_recovers_window_signatures_with_hdbscan() -> None:
    raw_df, expected_groups = _blob_frame(
        {
            10: [(0.2, -0.5), (0.8, 0.5)],
            11: [(0.4, 0.1), (0.4, 0.9)],
            12: [(0.6, -0.2), (0.1, -0.2)],
            13: [(0.3, 0.3), (0.3, 0.3)],
        }
    )
    # Unequal group sizes (20 vs 10) so misaligned PCA rows cannot go unnoticed.
    raw_df = raw_df[raw_df["model"] <= 30]

    assignments_df, summary_df = ClusteringService().cluster_residue_ranges(
        raw_df=raw_df,
        range_texts=["10-12"],
        min_cluster_size=5,
        min_samples=None,
    )

    assert pd.api.types.is_integer_dtype(assignments_df["model"])
    summary = summary_df.iloc[0]
    assert summary["range_label"] == "10-12"
    assert summary["residues"] == 3
    assert summary["models"] == 30
    assert summary["n_clusters"] == 2
    assert summary["noise_fraction"] == 0.0

    true_group = assignments_df["model"].map(expected_groups)
    _assert_labels_match_groups(assignments_df["cluster"], true_group)

    # The groups differ strongly across the window, so PC1 alone separates them.
    pc1_by_group = assignments_df.groupby(true_group)["pc1"].agg(["min", "max"])
    assert pc1_by_group.loc[0, "max"] < pc1_by_group.loc[1, "min"] or (
        pc1_by_group.loc[1, "max"] < pc1_by_group.loc[0, "min"]
    )


@pytest.mark.parametrize(
    "centers",
    [
        # Groups differ only in curvature.
        {1: [(0.2, 0.3), (0.8, 0.3)], 2: [(0.4, -0.1), (0.9, -0.1)]},
        # Groups differ only in torsion.
        {1: [(0.3, 0.2), (0.3, 0.8)], 2: [(-0.1, 0.4), (-0.1, 0.9)]},
    ],
    ids=["curvature-only", "torsion-only"],
)
def test_cluster_residue_ranges_uses_curvature_and_torsion(centers: dict) -> None:
    raw_df, expected_groups = _blob_frame(centers)

    assignments_df, summary_df = ClusteringService().cluster_residue_ranges(
        raw_df=raw_df,
        range_texts=["1-2"],
        min_cluster_size=5,
        min_samples=None,
    )

    assert summary_df.iloc[0]["n_clusters"] == 2
    _assert_labels_match_groups(
        assignments_df["cluster"], assignments_df["model"].map(expected_groups)
    )


def test_cluster_residue_ranges_clusters_each_window_independently() -> None:
    raw_df, _ = _blob_frame(
        {
            1: [(0.2, -0.5), (0.8, 0.5), (0.5, 0.0)],
            2: [(0.2, -0.5), (0.8, 0.5), (0.5, 0.0)],
            3: [(0.3, 0.3), (0.9, 0.9), (0.3, 0.3)],
            4: [(0.3, 0.3), (0.9, 0.9), (0.3, 0.3)],
        }
    )

    _, summary_df = ClusteringService().cluster_residue_ranges(
        raw_df=raw_df,
        range_texts=["1-2", "3-4"],
        min_cluster_size=5,
        min_samples=None,
    )

    assert summary_df.set_index("range_label")["n_clusters"].to_dict() == {"1-2": 3, "3-4": 2}


@pytest.mark.parametrize(
    ("min_cluster_size", "min_samples", "n_models", "expected"),
    [
        (None, None, 20, (5, 5)),
        (None, None, 100, (5, 5)),
        (None, None, 301, (15, 5)),
        (None, None, 1000, (50, 5)),
        (8, None, 1000, (8, 5)),
        (3, None, 1000, (3, 3)),
        (None, 2, 1000, (50, 2)),
        (20, 12, 30, (20, 12)),
    ],
)
def test_clustering_parameters_fill_in_the_defaults(
    min_cluster_size, min_samples, n_models: int, expected: tuple[int, int]
) -> None:
    # Smallest cluster 5% of the models (at least 5); min_samples 5 (or the size if smaller).
    assert clustering_parameters(min_cluster_size, min_samples, n_models) == expected


def test_default_parameters_scale_with_the_number_of_models(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import hdbscan

    created = []
    real_hdbscan = hdbscan.HDBSCAN

    def recording_hdbscan(**kwargs):
        created.append(kwargs)
        return real_hdbscan(**kwargs)

    monkeypatch.setattr(hdbscan, "HDBSCAN", recording_hdbscan)
    raw_df, _ = _blob_frame({1: [(0.2, -0.5)], 2: [(0.4, 0.1)]}, points_per_blob=200)
    service = ClusteringService()

    service.cluster_residues(raw_df=raw_df, min_cluster_size=None, min_samples=None)
    service.cluster_residue_ranges(
        raw_df=raw_df, range_texts=["1-2"], min_cluster_size=None, min_samples=None
    )

    expected = {"min_cluster_size": 10, "min_samples": 5, "allow_single_cluster": True}
    assert created == [expected] * 3  # two residues, one range
