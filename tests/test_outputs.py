from __future__ import annotations

from dataclasses import fields
from pathlib import Path

import pandas as pd
import pytest

from flexgeo2.config import OutputConfig
from flexgeo2.distances import DistanceService
from flexgeo2.geometry import GeometryService
from flexgeo2.models import (
    AnalysisResult,
    DistanceResult,
    OutputArtifacts,
    ResidueClusteringResult,
    ResidueRangeClusteringResult,
)
from flexgeo2.outputs import OutputWriter


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
            "chain_plotter",
            "distance_plotter",
            "residue_cluster_plotter",
            "residue_range_cluster_plotter",
        )
    }


def make_writer(
    output_dir: Path | None,
    plotters: dict[str, RecordingPlotter],
    verbose: bool = False,
    write_files: bool = True,
) -> OutputWriter:
    config = OutputConfig(output_dir=output_dir, verbose=verbose, write_files=write_files)
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
        .reset_index(name="n_conformations")
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
            "model": ["1", "2"],
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
                "n_conformations": 2,
                "n_residues": 2,
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


BASE_FILES = {
    "geometry_descriptors.csv",
    "residue_summary.csv",
    "model_summary_overall.csv",
    "plots/ensemble_overview.png",
}

FULL_DEFAULT_FILES = BASE_FILES | {
    "distance_to_reference_summary.csv",
    "plots/distance_to_reference_heatmap.png",
    "residue_cluster_summary.csv",
    "cluster_plots/A_ALA1_clusters.png",
    "cluster_plots/A_GLY2_clusters.png",
    "cluster_plots/B_GLY1_clusters.png",
    "residue_range_cluster_summary.csv",
    "range_cluster_plots/A_1-2_clusters.png",
}

VERBOSE_CHAIN_A_FILES = {
    "chains/A/geometry_descriptors.csv",
    "chains/A/residue_summary.csv",
    "chains/A/model_summary.csv",
    "chains/A/curvature_torsion.png",
    "chains/A/distance_to_reference_long.csv",
    "chains/A/distance_to_reference_summary.csv",
    "chains/A/distance_to_reference_matrix.csv",
    "chains/A/distance_to_reference_heatmap.png",
    "chains/A/residue_cluster_assignments.csv",
    "chains/A/residue_cluster_summary.csv",
    "chains/A/cluster_plots/ALA1_clusters.png",
    "chains/A/cluster_plots/GLY2_clusters.png",
    "chains/A/residue_range_cluster_assignments.csv",
    "chains/A/residue_range_cluster_summary.csv",
    "chains/A/range_cluster_plots/1-2_clusters.png",
}

VERBOSE_CHAIN_B_FILES = {
    "chains/B/geometry_descriptors.csv",
    "chains/B/residue_summary.csv",
    "chains/B/model_summary.csv",
    "chains/B/curvature_torsion.png",
    "chains/B/distance_to_reference_long.csv",
    "chains/B/distance_to_reference_summary.csv",
    "chains/B/distance_to_reference_matrix.csv",
    "chains/B/distance_to_reference_heatmap.png",
    "chains/B/residue_cluster_assignments.csv",
    "chains/B/residue_cluster_summary.csv",
    "chains/B/cluster_plots/GLY1_clusters.png",
}

FULL_VERBOSE_FILES = (
    FULL_DEFAULT_FILES
    | {
        "model_summary_by_chain.csv",
        "distance_to_reference_long.csv",
        "distance_matrices/A_distance_matrix.csv",
        "distance_matrices/B_distance_matrix.csv",
        "residue_cluster_assignments.csv",
        "residue_range_cluster_assignments.csv",
    }
    | VERBOSE_CHAIN_A_FILES
    | VERBOSE_CHAIN_B_FILES
)


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
        "overview_plot",
    }


def test_default_mode_with_all_analyses(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert written_files(tmp_path) == FULL_DEFAULT_FILES
    assert artifacts.chains_dir is None
    assert artifacts.model_summary_csv is None
    assert artifacts.distance_long_csv is None
    assert artifacts.distance_matrix_dir is None
    assert artifacts.cluster_assignments_csv is None
    assert artifacts.range_cluster_assignments_csv is None
    assert not plotters["chain_plotter"].calls


def test_verbose_mode_with_all_analyses(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters, verbose=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert written_files(tmp_path) == FULL_VERBOSE_FILES
    assert artifacts.chains_dir == tmp_path.resolve() / "chains"


@pytest.mark.parametrize("verbose", [False, True])
def test_reported_artifacts_exist_on_disk(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult, verbose: bool
) -> None:
    artifacts = make_writer(tmp_path, plotters, verbose=verbose).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    for name, path in artifact_paths(artifacts).items():
        assert path.exists(), f"{name} points to missing path {path}"
        assert path.is_relative_to(tmp_path.resolve())


def test_csv_outputs_round_trip(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    artifacts = make_writer(tmp_path, plotters, verbose=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    expected_tables = {
        artifacts.raw_csv: full_result.raw_df,
        artifacts.residue_summary_csv: full_result.residue_summary_df,
        artifacts.model_summary_csv: full_result.model_summary_df,
        artifacts.overall_model_summary_csv: full_result.overall_model_summary_df,
        artifacts.distance_long_csv: full_result.distance_result.long_df,
        artifacts.distance_summary_csv: full_result.distance_result.summary_df,
        artifacts.cluster_summary_csv: full_result.residue_clustering.summary_df,
    }
    for path, expected in expected_tables.items():
        written = pd.read_csv(path, dtype={"model": str, "chain": str})
        pd.testing.assert_frame_equal(
            written,
            expected.reset_index(drop=True),
            check_dtype=False,
            obj=path.name,
        )


def test_chain_outputs_contain_only_that_chain(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    make_writer(tmp_path, plotters, verbose=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    for chain in ("A", "B"):
        for csv_path in (tmp_path / "chains" / chain).glob("*.csv"):
            if csv_path.name == "distance_to_reference_matrix.csv":
                continue
            chains = pd.read_csv(csv_path, dtype={"chain": str})["chain"].unique().tolist()
            assert chains == [chain], csv_path.name


@pytest.mark.parametrize(
    ("hide_model_traces", "max_models"),
    [(False, 7), (True, 3)],
)
def test_trace_options_reach_overview_and_chain_plots(
    tmp_path: Path,
    plotters: dict,
    base_result: AnalysisResult,
    hide_model_traces: bool,
    max_models: int,
) -> None:
    make_writer(tmp_path, plotters, verbose=True).write(
        base_result, max_models_in_plot=max_models, hide_model_traces=hide_model_traces
    )

    [overview_call] = plotters["overview_plotter"].calls
    assert overview_call["kwargs"]["show_model_traces"] is (not hide_model_traces)
    assert overview_call["kwargs"]["max_models_in_plot"] == max_models
    pd.testing.assert_frame_equal(overview_call["kwargs"]["raw_df"], base_result.raw_df)

    chain_calls = plotters["chain_plotter"].calls
    assert len(chain_calls) == 2
    for call in chain_calls:
        assert call["kwargs"]["show_model_traces"] is (not hide_model_traces)
        assert call["kwargs"]["max_models_in_plot"] == max_models


def test_distance_heatmap_title_names_reference(
    tmp_path: Path, plotters: dict, full_result: AnalysisResult
) -> None:
    make_writer(tmp_path, plotters, verbose=True).write(
        full_result, max_models_in_plot=12, hide_model_traces=False
    )

    titles = [call["args"][2] for call in plotters["distance_plotter"].calls]
    assert titles == [
        "Distance to reference: input model 1",
        "Chain A: distance to input model 1",
        "Chain B: distance to input model 1",
    ]


def test_blank_chain_ids_use_unassigned_folder(
    tmp_path: Path, plotters: dict, base_result: AnalysisResult
) -> None:
    for frame in (
        base_result.raw_df,
        base_result.residue_summary_df,
        base_result.model_summary_df,
    ):
        frame.loc[frame["chain"] == "B", "chain"] = ""

    make_writer(tmp_path, plotters, verbose=True).write(
        base_result, max_models_in_plot=12, hide_model_traces=False
    )

    assert {path.name for path in (tmp_path / "chains").iterdir()} == {"A", "unassigned"}
    assert (tmp_path / "chains" / "unassigned" / "geometry_descriptors.csv").is_file()
