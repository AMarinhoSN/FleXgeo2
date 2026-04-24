from __future__ import annotations

import pandas as pd
import pytest

from flexgeo2.distances import DistanceService


def test_select_reference_rows_defaults_to_first_model(
    normalized_geometry_df: pd.DataFrame,
) -> None:
    selected, label = DistanceService().select_reference_rows(normalized_geometry_df, None)

    assert label == "1"
    assert selected["model"].drop_duplicates().tolist() == ["1"]


def test_select_reference_rows_rejects_missing_model(normalized_geometry_df: pd.DataFrame) -> None:
    with pytest.raises(ValueError, match="Reference model '9' was not found"):
        DistanceService().select_reference_rows(normalized_geometry_df, "9")


def test_distance_compute_uses_chain_order_residue_overlap(
    normalized_geometry_df: pd.DataFrame,
) -> None:
    service = DistanceService()
    reference_df, _ = service.select_reference_rows(normalized_geometry_df, "1")

    long_df, summary_df = service.compute(
        raw_df=normalized_geometry_df,
        reference_df=reference_df,
        reference_label="input model 1",
    )

    model_two_a2 = long_df[
        (long_df["chain"] == "A") & (long_df["model"] == "2") & (long_df["order"] == 2)
    ].iloc[0]
    assert model_two_a2["distance_to_reference"] == pytest.approx(
        ((0.6 - 0.4) ** 2 + (0.9 - 0.5) ** 2) ** 0.5
    )
    assert set(long_df["reference_label"]) == {"input model 1"}

    a1_summary = summary_df[(summary_df["chain"] == "A") & (summary_df["order"] == 1)].iloc[0]
    assert a1_summary["distance_min"] == 0.0
    assert a1_summary["distance_max"] == pytest.approx(((0.2 - 0.1) ** 2 + (0.3 - 0.2) ** 2) ** 0.5)
    assert a1_summary["models"] == 2


def test_distance_compute_rejects_no_overlap(normalized_geometry_df: pd.DataFrame) -> None:
    reference_df = normalized_geometry_df[normalized_geometry_df["model"] == "1"].copy()
    reference_df["chain"] = "Z"

    with pytest.raises(ValueError, match="No overlapping residues"):
        DistanceService().compute(
            raw_df=normalized_geometry_df,
            reference_df=reference_df,
            reference_label="external",
        )
