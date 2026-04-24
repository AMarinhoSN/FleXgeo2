from __future__ import annotations

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from flexgeo2.geometry import GeometryService


def test_filter_chains_keeps_requested_chains(normalized_geometry_df: pd.DataFrame) -> None:
    result = GeometryService().filter_chains(normalized_geometry_df, ["B"])

    assert result["chain"].tolist() == ["B", "B"]
    assert result["model"].tolist() == ["1", "2"]


def test_filter_chains_rejects_empty_selection(normalized_geometry_df: pd.DataFrame) -> None:
    with pytest.raises(ValueError, match="No residues found"):
        GeometryService().filter_chains(normalized_geometry_df, ["Z"])


def test_normalize_adds_residue_labels_and_sorts(raw_geometry_df: pd.DataFrame) -> None:
    result = GeometryService().normalize(raw_geometry_df)

    assert result[["chain", "model", "order"]].to_dict("records") == [
        {"chain": "A", "model": "1", "order": 1},
        {"chain": "A", "model": "1", "order": 2},
        {"chain": "A", "model": "2", "order": 1},
        {"chain": "A", "model": "2", "order": 2},
        {"chain": "B", "model": "1", "order": 1},
        {"chain": "B", "model": "2", "order": 1},
    ]
    assert result["residue_label"].tolist() == ["ALA1", "GLY2", "ALA1", "GLY2", "GLY1", "GLY1"]
    assert str(result["order"].dtype) == "int64"


def test_summarize_groups_residues(normalized_geometry_df: pd.DataFrame) -> None:
    result = GeometryService().summarize(normalized_geometry_df)

    expected = pd.DataFrame(
        [
            {
                "chain": "A",
                "order": 1,
                "name": "ALA",
                "residue_label": "ALA1",
                "curvature_mean": 0.15,
                "curvature_std": 0.0707106781,
                "curvature_min": 0.1,
                "curvature_max": 0.2,
                "torsion_mean": 0.25,
                "torsion_std": 0.0707106781,
                "torsion_min": 0.2,
                "torsion_max": 0.3,
                "models": 2,
                "curvature_dmax_min": 0.1,
                "curvature_dmax_max": 0.2,
                "torsion_dmax_min": 0.2,
                "torsion_dmax_max": 0.3,
                "curvature_dmax_bin_width": 0.05,
                "torsion_dmax_bin_width": 0.05,
                "dmax": 0.1414213562,
            },
            {
                "chain": "A",
                "order": 2,
                "name": "GLY",
                "residue_label": "GLY2",
                "curvature_mean": 0.5,
                "curvature_std": 0.1414213562,
                "curvature_min": 0.4,
                "curvature_max": 0.6,
                "torsion_mean": 0.7,
                "torsion_std": 0.2828427125,
                "torsion_min": 0.5,
                "torsion_max": 0.9,
                "models": 2,
                "curvature_dmax_min": 0.4,
                "curvature_dmax_max": 0.6,
                "torsion_dmax_min": 0.5,
                "torsion_dmax_max": 0.9,
                "curvature_dmax_bin_width": 0.1,
                "torsion_dmax_bin_width": 0.2,
                "dmax": 0.4472135955,
            },
            {
                "chain": "B",
                "order": 1,
                "name": "GLY",
                "residue_label": "GLY1",
                "curvature_mean": 1.75,
                "curvature_std": 0.3535533906,
                "curvature_min": 1.5,
                "curvature_max": 2.0,
                "torsion_mean": 2.25,
                "torsion_std": 0.3535533906,
                "torsion_min": 2.0,
                "torsion_max": 2.5,
                "models": 2,
                "curvature_dmax_min": 1.5,
                "curvature_dmax_max": 2.0,
                "torsion_dmax_min": 2.0,
                "torsion_dmax_max": 2.5,
                "curvature_dmax_bin_width": 0.25,
                "torsion_dmax_bin_width": 0.25,
                "dmax": 0.7071067812,
            },
        ]
    )
    assert_frame_equal(result, expected, check_exact=False, rtol=1e-9, atol=1e-12)


def test_summarize_dmax_defaults_to_raw_range_when_no_extreme_bins_are_sparse(
    normalized_geometry_df: pd.DataFrame,
) -> None:
    result = GeometryService().summarize(normalized_geometry_df)

    assert result["dmax"].tolist() == pytest.approx(
        [
            ((0.2 - 0.1) ** 2 + (0.3 - 0.2) ** 2) ** 0.5,
            ((0.6 - 0.4) ** 2 + (0.9 - 0.5) ** 2) ** 0.5,
            ((2.0 - 1.5) ** 2 + (2.5 - 2.0) ** 2) ** 0.5,
        ]
    )


def test_summarize_dmax_trims_sparse_extreme_bins() -> None:
    df = pd.DataFrame(
        {
            "model": [str(index) for index in range(12)],
            "chain": ["A"] * 12,
            "order": [1] * 12,
            "name": ["ALA"] * 12,
            "residue_label": ["ALA1"] * 12,
            "curvature": [-100, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 100],
            "torsion": [-50, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50],
        }
    )

    result = GeometryService().summarize(df, dmax_outlier_fraction=0.1)
    row = result.iloc[0]

    assert row["curvature_dmax_min"] == 0.0
    assert row["curvature_dmax_max"] == 9.0
    assert row["torsion_dmax_min"] == 10.0
    assert row["torsion_dmax_max"] == 19.0
    assert row["dmax"] == pytest.approx(((9.0 - 0.0) ** 2 + (19.0 - 10.0) ** 2) ** 0.5)


def test_summarize_dmax_zero_threshold_disables_trimming() -> None:
    df = pd.DataFrame(
        {
            "model": [str(index) for index in range(12)],
            "chain": ["A"] * 12,
            "order": [1] * 12,
            "name": ["ALA"] * 12,
            "residue_label": ["ALA1"] * 12,
            "curvature": [-100, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 100],
            "torsion": [-50, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 50],
        }
    )

    result = GeometryService().summarize(df, dmax_outlier_fraction=0.0)
    row = result.iloc[0]

    assert row["curvature_dmax_min"] == -100.0
    assert row["curvature_dmax_max"] == 100.0
    assert row["torsion_dmax_min"] == -50.0
    assert row["torsion_dmax_max"] == 50.0
    assert row["dmax"] == pytest.approx(((100.0 + 100.0) ** 2 + (50.0 + 50.0) ** 2) ** 0.5)


def test_summarize_dmax_constant_residue_is_zero() -> None:
    df = pd.DataFrame(
        [
            {
                "model": "1",
                "chain": "A",
                "order": 1,
                "name": "ALA",
                "residue_label": "ALA1",
                "curvature": 1.5,
                "torsion": 2.5,
            },
            {
                "model": "2",
                "chain": "A",
                "order": 1,
                "name": "ALA",
                "residue_label": "ALA1",
                "curvature": 1.5,
                "torsion": 2.5,
            },
        ]
    )

    result = GeometryService().summarize(df)

    assert result.loc[0, "curvature_dmax_bin_width"] == 0.0
    assert result.loc[0, "torsion_dmax_bin_width"] == 0.0
    assert result.loc[0, "dmax"] == 0.0


def test_build_model_summary_computes_chain_and_overall_values(
    normalized_geometry_df: pd.DataFrame,
) -> None:
    service = GeometryService()
    residue_summary_df = service.summarize(normalized_geometry_df)

    model_summary_df, overall_summary_df = service.build_model_summary(
        normalized_geometry_df, residue_summary_df
    )

    expected_model_summary = pd.DataFrame(
        [
            {
                "chain": "A",
                "model": "1",
                "residues": 2,
                "curvature_mean": 0.25,
                "torsion_mean": 0.35,
                "curvature_std": 0.2121320344,
                "torsion_std": 0.2121320344,
                "curvature_mean_abs_deviation": 0.075,
                "torsion_mean_abs_deviation": 0.125,
            },
            {
                "chain": "A",
                "model": "2",
                "residues": 2,
                "curvature_mean": 0.4,
                "torsion_mean": 0.6,
                "curvature_std": 0.2828427125,
                "torsion_std": 0.4242640687,
                "curvature_mean_abs_deviation": 0.075,
                "torsion_mean_abs_deviation": 0.125,
            },
            {
                "chain": "B",
                "model": "1",
                "residues": 1,
                "curvature_mean": 1.5,
                "torsion_mean": 2.0,
                "curvature_std": 0.0,
                "torsion_std": 0.0,
                "curvature_mean_abs_deviation": 0.25,
                "torsion_mean_abs_deviation": 0.25,
            },
            {
                "chain": "B",
                "model": "2",
                "residues": 1,
                "curvature_mean": 2.0,
                "torsion_mean": 2.5,
                "curvature_std": 0.0,
                "torsion_std": 0.0,
                "curvature_mean_abs_deviation": 0.25,
                "torsion_mean_abs_deviation": 0.25,
            },
        ]
    )
    assert_frame_equal(
        model_summary_df,
        expected_model_summary,
        check_exact=False,
        rtol=1e-9,
        atol=1e-12,
    )

    expected_overall = pd.DataFrame(
        [
            {
                "model": "1",
                "residues": 3,
                "curvature_mean": 0.6666666667,
                "torsion_mean": 0.9,
                "curvature_mean_abs_deviation": 0.1333333333,
                "torsion_mean_abs_deviation": 0.1666666667,
            },
            {
                "model": "2",
                "residues": 3,
                "curvature_mean": 0.9333333333,
                "torsion_mean": 1.2333333333,
                "curvature_mean_abs_deviation": 0.1333333333,
                "torsion_mean_abs_deviation": 0.1666666667,
            },
        ]
    )
    assert_frame_equal(
        overall_summary_df,
        expected_overall,
        check_exact=False,
        rtol=1e-9,
        atol=1e-12,
    )


def test_summarize_single_model_std_is_zero() -> None:
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
            }
        ]
    )

    result = GeometryService().summarize(df)

    assert_frame_equal(
        result[["curvature_std", "torsion_std"]],
        pd.DataFrame({"curvature_std": [0.0], "torsion_std": [0.0]}),
    )
    assert result.loc[0, "dmax"] == 0.0
