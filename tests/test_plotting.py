from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from flexgeo2.geometry import GeometryService
from flexgeo2.plotting import DistanceHeatmapPlotter, OverviewPlotter


@pytest.fixture
def distance_long_df() -> pd.DataFrame:
    # Residue labels deliberately sort differently from sequence order.
    residues = [(1, "VAL1"), (2, "ALA2"), (10, "MET10")]
    return pd.DataFrame(
        [
            {
                "chain": "A",
                "model": model,
                "order": order,
                "residue_label": label,
                "distance_to_reference": float(model * order),
            }
            for model in (1, 2)
            for order, label in residues
        ]
    )


def test_heatmap_residue_axis_follows_sequence_order(distance_long_df: pd.DataFrame) -> None:
    fig = DistanceHeatmapPlotter().render(distance_long_df, "test")
    try:
        axis = fig.axes[0]
        labels = [tick.get_text() for tick in axis.get_yticklabels()]
        assert labels == ["VAL1", "ALA2", "MET10"]
        assert not any(line.get_visible() for line in axis.get_xgridlines())
    finally:
        plt.close(fig)


def test_heatmap_plot_writes_file(distance_long_df: pd.DataFrame, tmp_path) -> None:
    output = tmp_path / "heatmap.png"

    DistanceHeatmapPlotter().plot(distance_long_df, output, "test")

    assert output.is_file()


@pytest.mark.parametrize(
    ("show_model_traces", "max_models", "expected_traces"),
    [(True, 12, 2), (True, 1, 1), (False, 12, 0)],
)
def test_overview_model_traces_follow_flags(
    normalized_geometry_df: pd.DataFrame,
    show_model_traces: bool,
    max_models: int,
    expected_traces: int,
) -> None:
    summary_df = GeometryService().summarize(normalized_geometry_df)
    fig = OverviewPlotter().render(
        summary_df,
        raw_df=normalized_geometry_df,
        show_model_traces=show_model_traces,
        max_models_in_plot=max_models,
    )
    try:
        curvature_axis = fig.axes[0]
        # One line for the ensemble mean plus one per model trace.
        assert len(curvature_axis.get_lines()) == 1 + expected_traces
    finally:
        plt.close(fig)
