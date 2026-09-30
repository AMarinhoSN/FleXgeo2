from __future__ import annotations

from dataclasses import fields
from pathlib import Path

import pandas as pd
import pytest

from flexgeo2.config import OutputConfig
from flexgeo2.distances import DistanceService
from flexgeo2.geometry import DMAX_DETAIL_COLUMNS, GeometryService
from flexgeo2.models import (
    AnalysisResult,
    DistanceResult,
    OutputArtifacts,
    ResidueClusteringResult,
    ResidueRangeClusteringResult,
)
from flexgeo2.outputs import OutputDirectoryNotEmptyError, OutputWriter


class RecordingPlotter:
    """Stand-in plotter that records calls and creates an empty output file."""

    def __init__(self) -> None:
        self.calls: list[dict] = []

    def plot(self, *args, **kwargs) -> None:
        output_path = Path(kwargs["output_path"] if "output_path" in kwargs else args[1])
        output_path.touch()
        self.calls.append({"args": args, "kwargs": kwargs, "output_path": output_path})


@pytest.fixture
def plotters() -> dict[str, RecordingPlotter]:
    return {
        name: RecordingPlotter()
        for name in (
            "overview_plotter",
            "distance_plotter",
            "residue_range_cluster_plotter",
            "cluster_map_plotter",
            "residue_plotter",
        )
    }


def make_writer(
    output_dir: Path | None,
    plotters: dict[str, RecordingPlotter],
    distance_matrices: bool = False,
    write_files: bool = True,
    overwrite: bool = False,
    plot_residues: list[str] | None = None,
) -> OutputWriter:
    config = OutputConfig(
        output_dir=output_dir,
        distance_matrices=distance_matrices,
        write_files=write_files,
        overwrite=overwrite,
        plot_residues=plot_residues or [],
    )
    return OutputWriter(config, **plotters)


@pytest.fixture
def base_result(normalized_geometry_df: pd.DataFrame) -> AnalysisResult:
    geometry = GeometryService()
    residue_summary_df = geometry.summarize(normalized_geometry_df)
    model_summary_df, overall_model_summary_df = geometry.build_model_summary(
        normalized_geometry_df, residue_summary_df
    )
    return AnalysisResult(
        pdb_file=Path("ensemble.pdb"),
        raw_df=normalized_geometry_df,
        residue_summary_df=residue_summary_df,
        model_summary_df=model_summary_df,
        overall_model_summary_df=overall_model_summary_df,
    )


@pytest.fixture
def full_result(base_result: AnalysisResult) -> AnalysisResult:
    raw_df = base_result.raw_df
    distances = DistanceService()
    reference_df, _ = distances.select_reference_rows(raw_df, "1")
    long_df, summary_df = distances.compute(raw_df, reference_df, "input model 1")
    base_result.distance_result = DistanceResult(
        long_df=long_df, summary_df=summary_df, reference_label="input model 1"
    )

    assignments_df = raw_df.assign(cluster=0, cluster_probability=1.0)
    residue_summary = (
        raw_df.groupby(["chain", "order", "name", "residue_label"])
        .size()
        .reset_index(name="models")
        .assign(n_clusters=1, noise_fraction=0.0)
    )
    base_result.residue_clustering = ResidueClusteringResult(
        assignments_df=assignments_df, summary_df=residue_summary
    )

    # Range clustering only on chain A, so chain B gets no range outputs.
    range_assignments = pd.DataFrame(
        {
            "chain": "A",
            "range_start": 1,
            "range_end": 2,
            "range_label": "1-2",
            "model": [1, 2],
            "cluster": [0, 0],
            "cluster_probability": [1.0, 1.0],
            "pc1": [-0.5, 0.5],
            "pc2": [0.0, 0.0],
        }
    )
    range_summary = pd.DataFrame(
        [
            {
                "chain": "A",
                "range_start": 1,
                "range_end": 2,
                "range_label": "1-2",
                "residues": 2,
                "models": 2,
                "n_clusters": 1,
                "noise_fraction": 0.0,
            }
        ]
    )
    base_result.residue_range_clustering = ResidueRangeClusteringResult(
        assignments_df=range_assignments, summary_df=range_summary
    )
    return base_result


def written_files(output_dir: Path) -> set[str]:
    return {
        path.relative_to(output_dir).as_posix() for path in output_dir.rglob("*") if path.is_file()
    }


def artifact_paths(artifacts: OutputArtifacts) -> dict[str, Path]:
    return {
        field.name: getattr(artifacts, field.name)
        for field in fields(artifacts)
        if getattr(artifacts, field.name) is not None
    }


# The fixture has two chains, so the per-chain model summary is written too.
BASE_FILES = {
    "README.md",
    "run.json",
    "overview.png",
    "geometry/descriptors.csv",
    "geometry/residues.csv",
    "geometry/models.csv",
    "geometry/models_by_chain.csv",
}

FULL_DEFAULT_FILES = BASE_FILES | {
    "reference/distances.csv",
    "reference/residues.csv",
    "reference/heatmap.png",
    "clusters/assignments.csv",
    "clusters/residues.csv",
    "clusters/clusters.png",
    "range_clusters/assignments.csv",
    "range_clusters/ranges.csv",
    "range_clusters/A_1-2.png",
}

FULL_WITH_MATRICES_FILES = FULL_DEFAULT_FILES | {
    "reference/matrices/A.csv",
    "reference/matrices/B.csv",
}


def test_write_files_disabled_writes_nothing(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    output_dir = tmp_path / "out"

    artifacts = make_writer(output_dir, plotters, write_files=False).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert artifacts == OutputArtifacts()
    assert not output_dir.exists()
    assert all(not plotter.calls for plotter in plotters.values())


def test_missing_output_dir_is_rejected(plotters: dict, base_result: AnalysisResult) -> None:
    with pytest.raises(ValueError, match="output_dir must be set"):
        make_writer(None, plotters).write(
            base_result, max_models_in_plot=12, hide_model_traces=False
        )


def test_default_mode_writes_only_core_outputs(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters).write(
        base_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert written_files(tmp_path) == BASE_FILES
    assert set(artifact_paths(artifacts)) == {
        "raw_csv",
        "residue_summary_csv",
        "overall_model_summary_csv",
        "model_summary_csv",
        "overview_plot",
        "readme",
        "run_manifest",
    }


def test_default_mode_with_all_analyses(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert written_files(tmp_path) == FULL_DEFAULT_FILES
    assert artifacts.distance_matrix_dir is None


def test_distance_matrices_option_adds_one_matrix_per_chain(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters, distance_matrices=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert written_files(tmp_path) == FULL_WITH_MATRICES_FILES
    assert artifacts.distance_matrix_dir == tmp_path.resolve() / "reference" / "matrices"


@pytest.mark.parametrize("distance_matrices", [False, True])
def test_reported_artifacts_exist_on_disk(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult, distance_matrices: bool
) -> None:
    artifacts = make_writer(tmp_path, plotters, distance_matrices=distance_matrices).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    for name, path in artifact_paths(artifacts).items():
        assert path.exists(), f"{name} points to missing path {path}"
        assert path.is_relative_to(tmp_path.resolve())


def test_csv_outputs_round_trip(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters, distance_matrices=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    expected_tables = {
        artifacts.raw_csv: full_result.raw_df,
        artifacts.residue_summary_csv: full_result.residue_summary_df.drop(
            columns=list(DMAX_DETAIL_COLUMNS)
        ),
        artifacts.model_summary_csv: full_result.model_summary_df,
        artifacts.overall_model_summary_csv: full_result.overall_model_summary_df,
        artifacts.distance_long_csv: full_result.distance_result.long_df,
        artifacts.distance_summary_csv: full_result.distance_result.summary_df,
        artifacts.cluster_assignments_csv: full_result.residue_clustering.assignments_df,
        artifacts.cluster_summary_csv: full_result.residue_clustering.summary_df,
        artifacts.range_cluster_assignments_csv: (
            full_result.residue_range_clustering.assignments_df
        ),
        artifacts.range_cluster_summary_csv: full_result.residue_range_clustering.summary_df,
    }
    for path, expected in expected_tables.items():
        written = pd.read_csv(path, dtype={"model": str, "chain": str})
        if "model" in expected.columns:
            expected = expected.astype({"model": str})
        pd.testing.assert_frame_equal(
            written,
            expected.reset_index(drop=True),
            check_dtype=False,
            obj=path.name,
        )


def test_residue_table_leaves_out_dmax_details(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters).write(
        base_result, max_models_in_plot=12, hide_model_traces=False
    )

    columns = pd.read_csv(artifacts.residue_summary_csv, nrows=0).columns.tolist()
    assert "dmax" in columns
    assert not set(DMAX_DETAIL_COLUMNS) & set(columns)
    # The library result keeps them, e.g. to draw the trimmed ranges.
    assert set(DMAX_DETAIL_COLUMNS) <= set(base_result.residue_summary_df.columns)


def test_distance_matrices_are_split_by_chain(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters, distance_matrices=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    matrix_a = pd.read_csv(artifacts.distance_matrix_dir / "A.csv", index_col="model")
    matrix_b = pd.read_csv(artifacts.distance_matrix_dir / "B.csv", index_col="model")
    assert matrix_a.columns.tolist() == ["ALA1", "GLY2"]
    assert matrix_b.columns.tolist() == ["GLY1"]
    assert matrix_a.index.tolist() == matrix_b.index.tolist() == [1, 2]


@pytest.mark.parametrize(
    ("hide_model_traces", "max_models"),
    [(False, 7), (True, 3)],
)
def test_trace_options_reach_overview_plot(
    tmp_path: Path,
    plotters: dict,
    base_result: AnalysisResult,
    hide_model_traces: bool,
    max_models: int,
) -> None:
    make_writer(tmp_path, plotters).write(
        base_result, max_models_in_plot=max_models, hide_model_traces=hide_model_traces
    )

    [overview_call] = plotters["overview_plotter"].calls
    assert overview_call["kwargs"]["show_model_traces"] is (not hide_model_traces)
    assert overview_call["kwargs"]["max_models_in_plot"] == max_models
    pd.testing.assert_frame_equal(overview_call["kwargs"]["raw_df"], base_result.raw_df)


def test_distance_heatmap_title_names_reference(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    make_writer(tmp_path, plotters, distance_matrices=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    titles = [call["args"][2] for call in plotters["distance_plotter"].calls]
    assert titles == ["Distance to reference: input model 1"]


def test_blank_chain_ids_use_unassigned_in_file_names(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    frames = [
        full_result.raw_df,
        full_result.residue_summary_df,
        full_result.model_summary_df,
        full_result.distance_result.long_df,
        full_result.distance_result.summary_df,
        full_result.residue_clustering.assignments_df,
        full_result.residue_clustering.summary_df,
    ]
    for frame in frames:
        frame.loc[frame["chain"] == "B", "chain"] = ""

    make_writer(tmp_path, plotters, distance_matrices=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert (tmp_path / "reference" / "matrices" / "unassigned.csv").is_file()


def test_non_empty_output_dir_is_refused(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    (tmp_path / "notes.txt").write_text("mine")

    with pytest.raises(OutputDirectoryNotEmptyError, match="overwrite=True"):
        make_writer(tmp_path, plotters).write(
            base_result, max_models_in_plot=12, hide_model_traces=False
        )

    assert written_files(tmp_path) == {"notes.txt"}


def test_hidden_files_do_not_make_output_dir_non_empty(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    (tmp_path / ".DS_Store").write_text("")

    make_writer(tmp_path, plotters).write(
        base_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert written_files(tmp_path) == BASE_FILES | {".DS_Store"}


def test_overwrite_removes_stale_outputs_and_keeps_other_files(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    make_writer(tmp_path, plotters, distance_matrices=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )
    (tmp_path / "notes.txt").write_text("mine")
    (tmp_path / "clusters" / "picked.txt").write_text("mine too")
    full_result.distance_result = None
    full_result.residue_range_clustering = None

    make_writer(tmp_path, plotters, overwrite=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    # Reference and range-clustering outputs are gone; residue clusters are rewritten.
    assert written_files(tmp_path) == {
        name
        for name in FULL_DEFAULT_FILES
        if not name.startswith(("reference/", "range_clusters/"))
    } | {"notes.txt", "clusters/picked.txt"}
    assert not (tmp_path / "reference").exists()
    assert not (tmp_path / "range_clusters").exists()
    assert (tmp_path / "notes.txt").read_text() == "mine"


def test_overwrite_refuses_a_readme_flexgeo2_did_not_write(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    (tmp_path / "README.md").write_text("# My project")

    with pytest.raises(ValueError, match="README.md that FleXgeo2 did not write"):
        make_writer(tmp_path, plotters, overwrite=True).write(
            base_result, max_models_in_plot=12, hide_model_traces=False
        )

    assert (tmp_path / "README.md").read_text() == "# My project"


def test_output_path_that_is_a_file_is_rejected(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    output_path = tmp_path / "results"
    output_path.write_text("")

    with pytest.raises(ValueError, match="exists and is not a folder"):
        make_writer(output_path, plotters).write(
            base_result, max_models_in_plot=12, hide_model_traces=False
        )


def test_overview_gets_the_summaries_of_the_analyses_that_ran(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    make_writer(tmp_path, plotters).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    [overview_call] = plotters["overview_plotter"].calls
    kwargs = overview_call["kwargs"]
    assert kwargs["cluster_summary_df"] is full_result.residue_clustering.summary_df
    assert kwargs["distance_summary_df"] is full_result.distance_result.summary_df


def test_overview_without_optional_analyses_gets_no_summaries(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    make_writer(tmp_path, plotters).write(
        base_result, max_models_in_plot=12, hide_model_traces=False
    )

    [overview_call] = plotters["overview_plotter"].calls
    assert overview_call["kwargs"]["cluster_summary_df"] is None
    assert overview_call["kwargs"]["distance_summary_df"] is None


def test_single_chain_input_has_no_per_chain_model_summary(
    tmp_path: Path, plotters: dict, normalized_geometry_df: pd.DataFrame
) -> None:
    chain_a = normalized_geometry_df[normalized_geometry_df["chain"] == "A"]
    geometry = GeometryService()
    residue_summary_df = geometry.summarize(chain_a)
    model_summary_df, overall_model_summary_df = geometry.build_model_summary(
        chain_a, residue_summary_df
    )
    result = AnalysisResult(
        pdb_file=Path("ensemble.pdb"),
        raw_df=chain_a,
        residue_summary_df=residue_summary_df,
        model_summary_df=model_summary_df,
        overall_model_summary_df=overall_model_summary_df,
    )

    artifacts = make_writer(tmp_path, plotters).write(
        result, max_models_in_plot=12, hide_model_traces=False
    )

    assert artifacts.model_summary_csv is None
    assert written_files(tmp_path) == BASE_FILES - {"geometry/models_by_chain.csv"}


def test_chosen_residues_are_plotted_with_clusters_reference_and_dmax(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    # Residue 1 exists in both chains; A:2 only in chain A.
    artifacts = make_writer(tmp_path, plotters, plot_residues=["1", "A:2"]).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert artifacts.residue_plots_dir == tmp_path.resolve() / "residue_plots"
    assert written_files(tmp_path) == FULL_DEFAULT_FILES | {
        "residue_plots/A_0001_ALA.png",
        "residue_plots/A_0002_GLY.png",
        "residue_plots/B_0001_GLY.png",
    }
    calls = plotters["residue_plotter"].calls
    assert len(calls) == 3
    summary = full_result.residue_summary_df.set_index(["chain", "order"])
    distances = full_result.distance_result.long_df
    for call in calls:
        points = call["args"][0]
        chain, order = points["chain"].iloc[0], points["order"].iloc[0]
        assert set(points["chain"]) == {chain} and set(points["order"]) == {order}
        assert len(points) == 2  # both models
        assert "cluster" in points.columns
        assert call["kwargs"]["dmax"] == summary.loc[(chain, order), "dmax"]
        reference = distances[(distances["chain"] == chain) & (distances["order"] == order)]
        assert call["kwargs"]["reference"] == (
            reference["reference_curvature"].iloc[0],
            reference["reference_torsion"].iloc[0],
        )


def test_residue_plots_without_clustering_or_reference(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    make_writer(tmp_path, plotters, plot_residues=["B:1"]).write(
        base_result, max_models_in_plot=12, hide_model_traces=False
    )

    [call] = plotters["residue_plotter"].calls
    assert call["output_path"].name == "B_0001_GLY.png"
    assert "cluster" not in call["args"][0].columns
    assert call["kwargs"]["reference"] is None


def test_no_residue_plots_unless_requested(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert artifacts.residue_plots_dir is None
    assert plotters["residue_plotter"].calls == []
