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
            },
        ]
    )
    assert_frame_equal(result, expected, check_exact=False, rtol=1e-9, atol=1e-12)


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
