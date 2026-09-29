from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.figure import Figure

from flexgeo2.geometry import GeometryService
from flexgeo2.plotting import (
    ChainGeometryPlotter,
    DistanceHeatmapPlotter,
    OverviewPlotter,
    ResidueClusterPlotter,
    ResidueRangeClusterPlotter,
)


@pytest.fixture(autouse=True)
def _no_leaked_figures():
    yield
    open_figures = plt.get_fignums()
    plt.close("all")
    assert open_figures == [], "plotter left matplotlib figures open"


@pytest.fixture
def saved_figures(monkeypatch: pytest.MonkeyPatch) -> list[Figure]:
    """Record every figure a plotter saves; figures stay inspectable after plt.close."""
    figures: list[Figure] = []
    original_savefig = Figure.savefig

    def recording_savefig(self, *args, **kwargs):
        figures.append(self)
        return original_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", recording_savefig)
    return figures


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


def residue_cluster_df(clusters: list[int], chain: str = "A") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chain": chain,
            "model": range(1, len(clusters) + 1),
            "order": 1,
            "name": "ALA",
            "residue_label": "ALA1",
            "curvature": [0.1 * index for index in range(len(clusters))],
            "torsion": [0.05 * index for index in range(len(clusters))],
            "cluster": clusters,
            "cluster_probability": 1.0,
        }
    )


def range_cluster_df(clusters: list[int], chain: str = "A") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chain": chain,
            "range_start": 1,
            "range_end": 2,
            "range_label": "1-2",
            "model": [str(model) for model in range(1, len(clusters) + 1)],
            "cluster": clusters,
            "cluster_probability": 1.0,
            "pc1": [0.1 * index for index in range(len(clusters))],
            "pc2": [-0.1 * index for index in range(len(clusters))],
        }
    )


def chain_frames(normalized_geometry_df: pd.DataFrame, chain: str = "A"):
    summary_df = GeometryService().summarize(normalized_geometry_df)
    return (
        normalized_geometry_df[normalized_geometry_df["chain"] == chain],
        summary_df[summary_df["chain"] == chain],
    )


def two_chain_distances(distance_long_df: pd.DataFrame) -> pd.DataFrame:
    return pd.concat([distance_long_df, distance_long_df.assign(chain="B")], ignore_index=True)


def legend_labels(figure: Figure) -> list[str]:
    return [text.get_text() for text in figure.axes[0].get_legend().get_texts()]


@pytest.mark.parametrize(
    "plotter",
    ["chain", "overview", "heatmap", "residue_cluster", "range_cluster"],
)
def test_every_plotter_writes_a_readable_png(
    plotter: str,
    normalized_geometry_df: pd.DataFrame,
    distance_long_df: pd.DataFrame,
    tmp_path,
) -> None:
    output = tmp_path / f"{plotter}.png"
    chain_raw_df, chain_summary_df = chain_frames(normalized_geometry_df)
    calls = {
        "chain": lambda: ChainGeometryPlotter().plot(
            chain_raw_df, chain_summary_df, output, show_model_traces=True, max_models_in_plot=12
        ),
        "overview": lambda: OverviewPlotter().plot(
            GeometryService().summarize(normalized_geometry_df),
            output,
            raw_df=normalized_geometry_df,
        ),
        "heatmap": lambda: DistanceHeatmapPlotter().plot(distance_long_df, output, "title"),
        "residue_cluster": lambda: ResidueClusterPlotter().plot(
            residue_cluster_df([-1, 0, 0, 1]), output
        ),
        "range_cluster": lambda: ResidueRangeClusterPlotter().plot(
            range_cluster_df([-1, 0, 0, 1]), output
        ),
    }

    calls[plotter]()

    assert output.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    height, width, _ = plt.imread(output).shape
    assert height > 100 and width > 100


@pytest.mark.parametrize(
    ("show_model_traces", "max_models", "expected_traces"),
    [(True, 12, 2), (True, 1, 1), (True, 0, 0), (False, 12, 0)],
)
def test_chain_plot_titles_and_model_traces(
    normalized_geometry_df: pd.DataFrame,
    saved_figures: list[Figure],
    tmp_path,
    show_model_traces: bool,
    max_models: int,
    expected_traces: int,
) -> None:
    chain_raw_df, chain_summary_df = chain_frames(normalized_geometry_df)

    ChainGeometryPlotter().plot(
        chain_raw_df,
        chain_summary_df,
        tmp_path / "chain.png",
        show_model_traces=show_model_traces,
        max_models_in_plot=max_models,
    )

    [figure] = saved_figures
    curvature_axis, torsion_axis = figure.axes
    assert curvature_axis.get_title() == "Chain A: Curvature"
    assert torsion_axis.get_title() == "Chain A: Torsion"
    for axis in (curvature_axis, torsion_axis):
        # One line for the ensemble mean plus one per model trace.
        assert len(axis.get_lines()) == 1 + expected_traces


@pytest.mark.parametrize(
    ("plotter", "frame"),
    [(ResidueClusterPlotter, residue_cluster_df), (ResidueRangeClusterPlotter, range_cluster_df)],
    ids=["residue", "range"],
)
def test_cluster_plots_label_noise_and_clusters(
    plotter, frame, saved_figures: list[Figure], tmp_path
) -> None:
    plotter().plot(frame([1, -1, 0, 0, 1, -1]), tmp_path / "clusters.png")

    [figure] = saved_figures
    assert legend_labels(figure) == ["Noise", "Cluster 0", "Cluster 1"]
    noise_collection = figure.axes[0].collections[0]
    assert noise_collection.get_facecolor()[0][:3] == pytest.approx((0x9E / 255,) * 3)


@pytest.mark.parametrize(
    ("plotter", "frame", "title"),
    [
        (ResidueClusterPlotter, residue_cluster_df, "Chain: ALA1"),
        (ResidueRangeClusterPlotter, range_cluster_df, "Chain: residues 1-2"),
    ],
    ids=["residue", "range"],
)
def test_cluster_plots_handle_blank_chain_and_many_clusters(
    plotter, frame, title: str, saved_figures: list[Figure], tmp_path
) -> None:
    # More clusters than colours in the tab10 palette.
    plotter().plot(frame(list(range(12)), chain=""), tmp_path / "clusters.png")

    [figure] = saved_figures
    assert figure.axes[0].get_title() == title
    assert legend_labels(figure) == [f"Cluster {index}" for index in range(12)]


def test_overview_has_one_row_per_chain(normalized_geometry_df: pd.DataFrame) -> None:
    summary_df = GeometryService().summarize(normalized_geometry_df)

    fig = OverviewPlotter().render(summary_df, raw_df=normalized_geometry_df)
    try:
        assert [axis.get_title() for axis in fig.axes] == [
            "Chain A: Curvature",
            "Chain A: Torsion",
            "Chain B: Curvature",
            "Chain B: Torsion",
        ]
    finally:
        plt.close(fig)


def test_heatmap_has_one_panel_per_chain(distance_long_df: pd.DataFrame) -> None:
    fig = DistanceHeatmapPlotter().render(two_chain_distances(distance_long_df), "title")
    try:
        panels = [axis for axis in fig.axes if axis.get_images()]
        assert [axis.get_title() for axis in panels] == [
            "Chain A: Distance to reference",
            "Chain B: Distance to reference",
        ]
    finally:
        plt.close(fig)


def test_heatmap_rejects_all_missing_distances(distance_long_df: pd.DataFrame) -> None:
    missing = distance_long_df.assign(distance_to_reference=float("nan"))

    with pytest.raises(ValueError, match="empty after alignment"):
        DistanceHeatmapPlotter().render(missing, "title")
