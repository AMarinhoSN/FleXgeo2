from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from flexgeo2.config import AnalysisConfig, ClusteringConfig, OutputConfig, ReferenceConfig
from flexgeo2.geometry import StructureInfo
from flexgeo2.models import OutputArtifacts
from flexgeo2.pipeline import FlexGeo2App


class FakeGeometryService:
    def __init__(self, raw_df: pd.DataFrame) -> None:
        self.raw_df = raw_df
        self.calls: list[str] = []

    def ensure_dependencies(self) -> None:
        self.calls.append("ensure_dependencies")

    def parse_structure(self, pdb_file: Path) -> str:
        self.calls.append(f"parse_structure:{Path(pdb_file).name}")
        return Path(pdb_file).name

    def describe_structure(self, structure: str) -> StructureInfo:
        self.calls.append(f"describe_structure:{structure}")
        return StructureInfo(
            model_ids=sorted(int(model) for model in self.raw_df["model"].unique()),
            residues_by_chain={
                chain: set(chain_df["order"]) for chain, chain_df in self.raw_df.groupby("chain")
            },
        )

    def compute_geometry(self, structure: str, n_jobs: int = 1) -> pd.DataFrame:
        self.calls.append(f"compute_geometry:{structure}:{n_jobs}")
        return self.raw_df.copy()

    def filter_chains(self, df: pd.DataFrame, chains: list[str] | None) -> pd.DataFrame:
        self.calls.append(f"filter_chains:{chains}")
        return df[df["chain"].isin(chains)].copy() if chains else df.copy()

    def normalize(self, df: pd.DataFrame) -> pd.DataFrame:
        self.calls.append("normalize")
        return df.copy()

    def summarize(
        self,
        df: pd.DataFrame,
        dmax_outlier_fraction: float = 0.01,
    ) -> pd.DataFrame:
        self.calls.append(f"summarize:{dmax_outlier_fraction}")
        return pd.DataFrame(
            [
                {
                    "chain": "A",
                    "order": 1,
                    "name": "ALA",
                    "residue_label": "ALA1",
                    "curvature_mean": df["curvature"].mean(),
                    "curvature_std": 0.0,
                    "curvature_min": df["curvature"].min(),
                    "curvature_max": df["curvature"].max(),
                    "torsion_mean": df["torsion"].mean(),
                    "torsion_std": 0.0,
                    "torsion_min": df["torsion"].min(),
                    "torsion_max": df["torsion"].max(),
                    "models": df["model"].nunique(),
                    "curvature_dmax_min": df["curvature"].min(),
                    "curvature_dmax_max": df["curvature"].max(),
                    "torsion_dmax_min": df["torsion"].min(),
                    "torsion_dmax_max": df["torsion"].max(),
                    "curvature_dmax_bin_width": 0.0,
                    "torsion_dmax_bin_width": 0.0,
                    "dmax": 0.0,
                }
            ]
        )

    def build_model_summary(
        self, raw_df: pd.DataFrame, residue_summary_df: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        self.calls.append("build_model_summary")
        return (
            pd.DataFrame(
                [
                    {
                        "chain": "A",
                        "model": "1",
                        "residues": 1,
                        "curvature_mean": 0.1,
                        "torsion_mean": 0.2,
                        "curvature_std": 0.0,
                        "torsion_std": 0.0,
                        "curvature_mean_abs_deviation": 0.0,
                        "torsion_mean_abs_deviation": 0.0,
                    }
                ]
            ),
            pd.DataFrame(
                [
                    {
                        "model": "1",
                        "residues": 1,
                        "curvature_mean": 0.1,
                        "torsion_mean": 0.2,
                        "curvature_mean_abs_deviation": 0.0,
                        "torsion_mean_abs_deviation": 0.0,
                    }
                ]
            ),
        )


class FakeDistanceService:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def select_reference_rows(
        self, df: pd.DataFrame, model_id: str | None
    ) -> tuple[pd.DataFrame, str]:
        self.calls.append(f"select_reference_rows:{model_id}")
        if model_id is None:
            model_id = df["model"].iloc[0]
        return df[df["model"] == model_id].copy(), str(model_id)

    def compute(
        self, raw_df: pd.DataFrame, reference_df: pd.DataFrame, reference_label: str
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        self.calls.append(f"compute:{reference_label}")
        long_df = raw_df.copy()
        long_df["reference_curvature"] = reference_df["curvature"].iloc[0]
        long_df["reference_torsion"] = reference_df["torsion"].iloc[0]
        long_df["distance_to_reference"] = 0.0
        long_df["reference_label"] = reference_label
        return long_df, pd.DataFrame({"reference_label": [reference_label]})


class FakeClusteringService:
    def __init__(self) -> None:
        self.calls: list[tuple] = []

    def cluster_residues(
        self, raw_df: pd.DataFrame, min_cluster_size: int, min_samples: int | None
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        self.calls.append(("cluster_residues", raw_df, min_cluster_size, min_samples))
        return (
            pd.DataFrame({"source": ["residue_assignments"]}),
            pd.DataFrame({"source": ["residue_summary"]}),
        )

    def cluster_residue_ranges(
        self,
        raw_df: pd.DataFrame,
        range_texts: list[str],
        min_cluster_size: int,
        min_samples: int | None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        self.calls.append(
            ("cluster_residue_ranges", raw_df, range_texts, min_cluster_size, min_samples)
        )
        return (
            pd.DataFrame({"source": ["range_assignments"]}),
            pd.DataFrame({"source": ["range_summary"]}),
        )


class FakeOutputWriter:
    instances: list[FakeOutputWriter] = []

    def __init__(self, config: OutputConfig) -> None:
        self.config = config
        self.write_calls: list[tuple[int, bool]] = []
        FakeOutputWriter.instances.append(self)

    def write(
        self,
        result,
        max_models_in_plot: int,
        hide_model_traces: bool,
    ) -> OutputArtifacts:
        self.write_calls.append((max_models_in_plot, hide_model_traces))
        return OutputArtifacts(raw_csv=Path("raw.csv"))


@pytest.fixture(autouse=True)
def _isolate_pipeline(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("flexgeo2.pipeline.PlotStyle.apply", lambda: None)
    FakeOutputWriter.instances = []


@pytest.fixture
def pdb_file(tmp_path: Path) -> Path:
    path = tmp_path / "ensemble.pdb"
    path.write_text("HEADER test\n")
    return path


@pytest.fixture
def reference_pdb(tmp_path: Path) -> Path:
    path = tmp_path / "reference.pdb"
    path.write_text("HEADER reference\n")
    return path


@pytest.fixture
def geometry(normalized_geometry_df: pd.DataFrame) -> FakeGeometryService:
    return FakeGeometryService(normalized_geometry_df)


@pytest.fixture
def distances() -> FakeDistanceService:
    return FakeDistanceService()


@pytest.fixture
def clustering() -> FakeClusteringService:
    return FakeClusteringService()


@pytest.fixture
def app(
    geometry: FakeGeometryService,
    distances: FakeDistanceService,
    clustering: FakeClusteringService,
) -> FlexGeo2App:
    return FlexGeo2App(
        geometry_service=geometry,
        distance_service=distances,
        clustering_service=clustering,
        output_writer_cls=FakeOutputWriter,
    )


def test_app_run_with_reference_model_wires_distance_result(
    app: FlexGeo2App,
    geometry: FakeGeometryService,
    distances: FakeDistanceService,
    pdb_file: Path,
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        chains=["A"],
        n_jobs=2,
        max_models_in_plot=4,
        hide_model_traces=True,
        dmax_outlier_fraction=0.05,
        reference=ReferenceConfig(model_id="1"),
        output=OutputConfig(write_files=False),
    )

    result = app.run(config)

    assert result.pdb_file == pdb_file.resolve()
    assert result.distance_result is not None
    assert result.distance_result.reference_label == "input model 1"
    assert geometry.calls == [
        "ensure_dependencies",
        "parse_structure:ensemble.pdb",
        "describe_structure:ensemble.pdb",
        "compute_geometry:ensemble.pdb:2",
        "filter_chains:['A']",
        "normalize",
        "summarize:0.05",
        "build_model_summary",
    ]
    assert distances.calls == ["select_reference_rows:1", "compute:input model 1"]
    assert FakeOutputWriter.instances[0].config.write_files is False
    assert FakeOutputWriter.instances[0].write_calls == [(4, True)]
    assert result.outputs == OutputArtifacts(raw_csv=Path("raw.csv"))


def test_app_run_without_optional_analyses(
    app: FlexGeo2App,
    distances: FakeDistanceService,
    clustering: FakeClusteringService,
    pdb_file: Path,
) -> None:
    config = AnalysisConfig(pdb_file=pdb_file, output=OutputConfig(write_files=False))

    result = app.run(config)

    assert result.distance_result is None
    assert result.residue_clustering is None
    assert result.residue_range_clustering is None
    assert distances.calls == []
    assert clustering.calls == []
    assert result.raw_df["chain"].unique().tolist() == ["A", "B"]


def test_app_run_with_reference_pdb_loads_and_filters_reference(
    app: FlexGeo2App,
    geometry: FakeGeometryService,
    distances: FakeDistanceService,
    pdb_file: Path,
    reference_pdb: Path,
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        chains=["A"],
        n_jobs=3,
        reference=ReferenceConfig(pdb_file=reference_pdb, pdb_model_id="2"),
        output=OutputConfig(write_files=False),
    )

    result = app.run(config)

    assert geometry.calls == [
        "ensure_dependencies",
        "parse_structure:ensemble.pdb",
        "parse_structure:reference.pdb",
        "describe_structure:ensemble.pdb",
        "describe_structure:reference.pdb",
        "compute_geometry:ensemble.pdb:3",
        "filter_chains:['A']",
        "normalize",
        "summarize:0.01",
        "build_model_summary",
        "compute_geometry:reference.pdb:3",
        "filter_chains:['A']",
        "normalize",
    ]
    assert distances.calls == ["select_reference_rows:2", "compute:reference.pdb model 2"]
    assert result.distance_result.reference_label == "reference.pdb model 2"


def test_app_run_with_reference_pdb_defaults_to_first_model(
    app: FlexGeo2App,
    distances: FakeDistanceService,
    pdb_file: Path,
    reference_pdb: Path,
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        reference=ReferenceConfig(pdb_file=reference_pdb),
        output=OutputConfig(write_files=False),
    )

    result = app.run(config)

    assert distances.calls == ["select_reference_rows:None", "compute:reference.pdb model 1"]
    assert result.distance_result.reference_label == "reference.pdb model 1"


def test_app_run_wires_both_clustering_modes(
    app: FlexGeo2App,
    clustering: FakeClusteringService,
    pdb_file: Path,
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        chains=["A"],
        clustering=ClusteringConfig(
            cluster_residues=True,
            cluster_residue_ranges=["1-2"],
            min_cluster_size=3,
            min_samples=2,
        ),
        output=OutputConfig(write_files=False),
    )

    result = app.run(config)

    [residue_call, range_call] = clustering.calls
    assert residue_call[0] == "cluster_residues"
    assert residue_call[2:] == (3, 2)
    assert range_call[0] == "cluster_residue_ranges"
    assert range_call[2:] == (["1-2"], 3, 2)
    for call in (residue_call, range_call):
        pd.testing.assert_frame_equal(call[1], result.raw_df)
        assert call[1]["chain"].unique().tolist() == ["A"]

    assert result.residue_clustering.assignments_df["source"].tolist() == ["residue_assignments"]
    assert result.residue_clustering.summary_df["source"].tolist() == ["residue_summary"]
    assert result.residue_range_clustering.assignments_df["source"].tolist() == [
        "range_assignments"
    ]
    assert result.residue_range_clustering.summary_df["source"].tolist() == ["range_summary"]


def test_app_run_passes_output_config_and_result_to_writer(
    app: FlexGeo2App, pdb_file: Path, tmp_path: Path
) -> None:
    output = OutputConfig(output_dir=tmp_path / "out", overwrite=True)
    config = AnalysisConfig(pdb_file=pdb_file, output=output)

    result = app.run(config)

    [writer] = FakeOutputWriter.instances
    assert writer.config is output
    assert writer.write_calls == [(12, False)]
    assert result.outputs == OutputArtifacts(raw_csv=Path("raw.csv"))


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"reference": ReferenceConfig(model_id="9")}, "Reference model '9' was not found"),
        ({"chains": ["Z"]}, r"Chain\(s\) not found"),
        (
            {"clustering": ClusteringConfig(cluster_residue_ranges=["1-5"])},
            "1-5 on chain 'A' is incomplete",
        ),
    ],
)
def test_app_run_validates_before_computing_geometry(
    app: FlexGeo2App,
    geometry: FakeGeometryService,
    pdb_file: Path,
    overrides: dict,
    message: str,
) -> None:
    config = AnalysisConfig(pdb_file=pdb_file, output=OutputConfig(write_files=False), **overrides)

    with pytest.raises(ValueError, match=message):
        app.run(config)

    assert not any(call.startswith("compute_geometry") for call in geometry.calls)
    assert FakeOutputWriter.instances == []


def test_app_run_validates_reference_pdb_model_before_computing_geometry(
    app: FlexGeo2App,
    geometry: FakeGeometryService,
    pdb_file: Path,
    reference_pdb: Path,
) -> None:
    config = AnalysisConfig(
        pdb_file=pdb_file,
        reference=ReferenceConfig(pdb_file=reference_pdb, pdb_model_id="9"),
        output=OutputConfig(write_files=False),
    )

    with pytest.raises(ValueError, match="not found in the reference PDB"):
        app.run(config)

    assert not any(call.startswith("compute_geometry") for call in geometry.calls)
