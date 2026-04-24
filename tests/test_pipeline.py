from __future__ import annotations

from pathlib import Path

import pandas as pd

from flexgeo2.config import AnalysisConfig, OutputConfig, ReferenceConfig
from flexgeo2.models import OutputArtifacts
from flexgeo2.pipeline import FlexGeo2App


class FakeGeometryService:
    def __init__(self, raw_df: pd.DataFrame) -> None:
        self.raw_df = raw_df
        self.calls: list[str] = []

    def ensure_dependencies(self) -> None:
        self.calls.append("ensure_dependencies")

    def load_structure(self, pdb_file: Path, n_jobs: int = 1) -> pd.DataFrame:
        self.calls.append(f"load_structure:{Path(pdb_file).name}:{n_jobs}")
        return self.raw_df.copy()

    def filter_chains(self, df: pd.DataFrame, chains: list[str] | None) -> pd.DataFrame:
        self.calls.append(f"filter_chains:{chains}")
        return df[df["chain"].isin(chains)].copy() if chains else df.copy()

    def normalize(self, df: pd.DataFrame) -> pd.DataFrame:
        self.calls.append("normalize")
        return df.copy()

    def summarize(self, df: pd.DataFrame) -> pd.DataFrame:
        self.calls.append("summarize")
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
    pass


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


def test_app_run_with_reference_model_wires_distance_result(
    monkeypatch,
    normalized_geometry_df: pd.DataFrame,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr("flexgeo2.pipeline.PlotStyle.apply", lambda: None)
    FakeOutputWriter.instances = []
    pdb_file = tmp_path / "ensemble.pdb"
    pdb_file.write_text("HEADER test\n")
    geometry = FakeGeometryService(normalized_geometry_df)
    distances = FakeDistanceService()
    app = FlexGeo2App(
        geometry_service=geometry,
        distance_service=distances,
        clustering_service=FakeClusteringService(),
        output_writer_cls=FakeOutputWriter,
    )
    config = AnalysisConfig(
        pdb_file=pdb_file,
        chains=["A"],
        n_jobs=2,
        max_models_in_plot=4,
        hide_model_traces=True,
        reference=ReferenceConfig(model_id="1"),
        output=OutputConfig(write_files=False),
    )

    result = app.run(config)

    assert result.pdb_file == pdb_file.resolve()
    assert result.distance_result is not None
    assert result.distance_result.reference_label == "input model 1"
    assert geometry.calls == [
        "ensure_dependencies",
        "load_structure:ensemble.pdb:2",
        "filter_chains:['A']",
        "normalize",
        "summarize",
        "build_model_summary",
    ]
    assert distances.calls == ["select_reference_rows:1", "compute:input model 1"]
    assert FakeOutputWriter.instances[0].config.write_files is False
    assert FakeOutputWriter.instances[0].write_calls == [(4, True)]
    assert result.outputs == OutputArtifacts(raw_csv=Path("raw.csv"))
