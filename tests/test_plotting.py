from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure

from flexgeo2.geometry import GeometryService
from flexgeo2.plotting import (
    NOISE_COLOR,
    ChainGeometryPlotter,
    ClusterMapPlotter,
    DistanceHeatmapPlotter,
    OverviewPlotter,
    ResidueRangeClusterPlotter,
    cluster_color,
    cluster_palette,
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
    ["chain", "overview", "heatmap", "range_cluster", "cluster_map"],
)
def test_every_plotter_writes_a_readable_png(
    plotter: str,
    normalized_geometry_df: pd.DataFrame,
    distance_long_df: pd.DataFrame,
    tmp_path,
) -> None:
    output = tmp_path / f"{plotter}.png"
    chain_raw_df, chain_summary_df = chain_frames(normalized_geometry_df)
    map_summary, map_assignments = cluster_map_frames({1: [0, 0, 1], 2: [-1, 0, 0]})
    overview_summary = GeometryService().summarize(normalized_geometry_df)
    calls = {
        "chain": lambda: ChainGeometryPlotter().plot(
            chain_raw_df, chain_summary_df, output, show_model_traces=True, max_models_in_plot=12
        ),
        "overview": lambda: OverviewPlotter().plot(
            overview_summary,
            output,
            raw_df=normalized_geometry_df,
            cluster_summary_df=overview_extras(overview_summary)[0],
            distance_summary_df=overview_extras(overview_summary)[1],
        ),
        "heatmap": lambda: DistanceHeatmapPlotter().plot(distance_long_df, output, "title"),
        "range_cluster": lambda: ResidueRangeClusterPlotter().plot(
            range_cluster_df([-1, 0, 0, 1]), output
        ),
        "cluster_map": lambda: ClusterMapPlotter().plot(map_summary, output, map_assignments),
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


def test_range_cluster_plot_labels_noise_and_clusters(
    saved_figures: list[Figure], tmp_path
) -> None:
    ResidueRangeClusterPlotter().plot(range_cluster_df([1, -1, 0, 0, 1, -1]), tmp_path / "r.png")

    [figure] = saved_figures
    assert legend_labels(figure) == ["Noise", "Cluster 0", "Cluster 1"]
    noise_collection = figure.axes[0].collections[0]
    assert noise_collection.get_facecolor()[0][:3] == pytest.approx((0x9E / 255,) * 3)


def test_range_cluster_plot_handles_blank_chain_and_many_clusters(
    saved_figures: list[Figure], tmp_path
) -> None:
    # More clusters than colours in the tab10 palette.
    ResidueRangeClusterPlotter().plot(
        range_cluster_df(list(range(12)), chain=""), tmp_path / "r.png"
    )

    [figure] = saved_figures
    assert figure.axes[0].get_title() == "Chain: residues 1-2"
    assert legend_labels(figure) == [f"Cluster {index}" for index in range(12)]


def overview_extras(summary_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cluster and distance summaries for every residue, in reverse sequence order."""
    residues = summary_df[["chain", "order"]].iloc[::-1].reset_index(drop=True)
    clusters = residues.assign(n_clusters=range(1, len(residues) + 1))
    distances = residues.assign(
        distance_mean=[0.5 * (index + 1) for index in range(len(residues))],
        distance_std=0.1,
    )
    return clusters, distances


def panel_labels(figure: Figure) -> list[str]:
    return [axis.get_ylabel() for axis in figure.axes]


@pytest.mark.parametrize(
    ("with_clusters", "with_distances", "expected"),
    [
        (False, False, ["Curvature (1/Å)", "Torsion (1/Å)", "dmax"]),
        (True, False, ["Curvature (1/Å)", "Torsion (1/Å)", "dmax", "Clusters"]),
        (False, True, ["Curvature (1/Å)", "Torsion (1/Å)", "dmax", "Distance to\nreference"]),
        (
            True,
            True,
            ["Curvature (1/Å)", "Torsion (1/Å)", "dmax", "Clusters", "Distance to\nreference"],
        ),
    ],
)
def test_overview_adds_panels_for_the_analyses_that_ran(
    normalized_geometry_df: pd.DataFrame,
    with_clusters: bool,
    with_distances: bool,
    expected: list[str],
) -> None:
    summary_df = GeometryService().summarize(normalized_geometry_df)
    clusters, distances = overview_extras(summary_df)

    fig = OverviewPlotter().render(
        summary_df,
        cluster_summary_df=clusters if with_clusters else None,
        distance_summary_df=distances if with_distances else None,
    )
    try:
        # Two chains, each with the same stack of panels.
        assert panel_labels(fig) == expected * 2
        titles = [axis.get_title() for axis in fig.axes]
        assert titles == ["Chain A"] + [""] * (len(expected) - 1) + ["Chain B"] + [""] * (
            len(expected) - 1
        )
        # Only the bottom panel of each chain shows residue labels.
        bottoms = {len(expected) - 1, 2 * len(expected) - 1}
        for index, axis in enumerate(fig.axes):
            assert axis.xaxis.get_tick_params()["labelbottom"] is (index in bottoms), index
    finally:
        plt.close(fig)


def test_overview_panels_line_up_by_residue(normalized_geometry_df: pd.DataFrame) -> None:
    summary_df = GeometryService().summarize(normalized_geometry_df)
    clusters, distances = overview_extras(summary_df)
    # Chain A residue 2 was not compared with the reference.
    distances = distances.drop(
        distances[(distances["chain"] == "A") & (distances["order"] == 2)].index
    )

    fig = OverviewPlotter().render(
        summary_df, cluster_summary_df=clusters, distance_summary_df=distances
    )
    try:
        _, _, dmax_axis, clusters_axis, distance_axis = fig.axes[:5]
        chain_a = summary_df[summary_df["chain"] == "A"].sort_values("order")
        assert [bar.get_x() + bar.get_width() / 2 for bar in dmax_axis.patches] == [1, 2]
        assert dmax_axis.get_xlim() == pytest.approx((0.4, 2.6))
        assert [bar.get_height() for bar in dmax_axis.patches] == pytest.approx(
            chain_a["dmax"].tolist()
        )
        # overview_extras numbers residues in reverse: B1 -> 1, A2 -> 2, A1 -> 3.
        assert [bar.get_height() for bar in clusters_axis.patches] == [3, 2]
        [mean_line] = distance_axis.get_lines()
        assert mean_line.get_xdata().tolist() == [1, 2]
        assert mean_line.get_ydata()[0] == pytest.approx(1.5)
        assert np.isnan(mean_line.get_ydata()[1])
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


def cluster_map_frames(
    labels_by_order: dict[int, list[int]], chain: str = "A"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summary and assignments for residues {order: [label of model 1, 2, ...]}."""
    names = {1: "VAL1", 2: "ALA2", 10: "MET10"}
    assignments = pd.DataFrame(
        [
            {
                "chain": chain,
                "model": model,
                "order": order,
                "residue_label": names[order],
                "cluster": label,
            }
            for order, labels in labels_by_order.items()
            for model, label in enumerate(labels, start=1)
        ]
    )
    summary = pd.DataFrame(
        [
            {
                "chain": chain,
                "order": order,
                "residue_label": names[order],
                "n_clusters": len({label for label in labels if label >= 0}),
            }
            for order, labels in labels_by_order.items()
        ]
    )
    return summary, assignments


def test_cluster_map_colours_each_cell_by_its_label() -> None:
    # Given out of sequence order; residue labels also sort differently from order.
    labels = {10: [2, 2, -1], 1: [0, 0, 1], 2: [-1, 0, 0]}
    summary, assignments = cluster_map_frames(labels)

    fig = ClusterMapPlotter().render(summary, assignments)
    try:
        strip_axis, map_axis = fig.axes
        image = map_axis.get_images()[0].get_array()
        for column, order in enumerate([1, 2, 10]):
            for row, label in enumerate(labels[order]):
                assert tuple(image[row, column]) == to_rgba(cluster_color(label)), (order, row)

        assert [bar.get_height() for bar in strip_axis.patches] == [2, 1, 1]
        assert [tick.get_text() for tick in map_axis.get_xticklabels()] == [
            "VAL1",
            "ALA2",
            "MET10",
        ]
        assert [tick.get_text() for tick in map_axis.get_yticklabels()] == ["1", "2", "3"]
        assert [text.get_text() for text in fig.legends[0].get_texts()] == [
            "Noise",
            "Cluster 0",
            "Cluster 1",
            "Cluster 2",
        ]
    finally:
        plt.close(fig)


def test_cluster_map_marks_missing_residues_and_draws_each_chain() -> None:
    summary_a, assignments_a = cluster_map_frames({1: [0, 0, 1], 2: [0, 0, 0]})
    summary_b, assignments_b = cluster_map_frames({1: [0, 1, 1]}, chain="B")
    # Model 3 lacks residue 2 of chain A.
    assignments_a = assignments_a.drop(
        assignments_a[(assignments_a["order"] == 2) & (assignments_a["model"] == 3)].index
    )
    summary = pd.concat([summary_a, summary_b], ignore_index=True)
    assignments = pd.concat([assignments_a, assignments_b], ignore_index=True)

    fig = ClusterMapPlotter().render(summary, assignments)
    try:
        titles = [axis.get_title() for axis in fig.axes]
        assert titles == ["Chain A: clusters per residue", "", "Chain B: clusters per residue", ""]
        image_a = fig.axes[1].get_images()[0].get_array()
        assert tuple(image_a[2, 1]) == to_rgba("white")
        assert fig.legends[0].get_texts()[-1].get_text() == "Not present"
    finally:
        plt.close(fig)


def test_cluster_colours_depend_only_on_the_label(saved_figures: list[Figure], tmp_path) -> None:
    ResidueRangeClusterPlotter().plot(range_cluster_df([0, 0, 1]), tmp_path / "no_noise.png")
    ResidueRangeClusterPlotter().plot(range_cluster_df([-1, 0, 1]), tmp_path / "noise.png")

    def colours(figure: Figure) -> dict[str, tuple]:
        axis = figure.axes[0]
        return {
            text.get_text(): tuple(collection.get_facecolor()[0])
            for text, collection in zip(
                axis.get_legend().get_texts(), axis.collections, strict=True
            )
        }

    without_noise, with_noise = (colours(figure) for figure in saved_figures)
    assert with_noise["Cluster 0"] == without_noise["Cluster 0"]
    assert with_noise["Cluster 1"] == without_noise["Cluster 1"]


def test_cluster_palette_is_distinct_and_avoids_the_noise_grey() -> None:
    palette = [to_rgba(color) for color in cluster_palette()]

    assert len(set(palette)) == len(palette) >= 18
    for color in palette:
        red, green, blue, _ = color
        assert not (np.isclose(red, green) and np.isclose(green, blue)), color
    assert cluster_color(0) != cluster_color(10)
    assert cluster_color(-1) == NOISE_COLOR
