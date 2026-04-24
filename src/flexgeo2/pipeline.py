from __future__ import annotations

from pathlib import Path

from flexgeo2.clustering import ClusteringService
from flexgeo2.config import AnalysisConfig
from flexgeo2.distances import DistanceService
from flexgeo2.geometry import GeometryService
from flexgeo2.models import (
    AnalysisResult,
    DistanceResult,
    ResidueClusteringResult,
    ResidueRangeClusteringResult,
)
from flexgeo2.outputs import OutputWriter
from flexgeo2.plotting import PlotStyle


class FlexGeo2App:
    """Library-facing application facade for complete FleXgeo2 analyses."""

    def __init__(
        self,
        geometry_service: GeometryService | None = None,
        distance_service: DistanceService | None = None,
        clustering_service: ClusteringService | None = None,
        output_writer_cls: type[OutputWriter] = OutputWriter,
    ) -> None:
        self.geometry = geometry_service or GeometryService()
        self.distances = distance_service or DistanceService()
        self.clustering = clustering_service or ClusteringService()
        self.output_writer_cls = output_writer_cls

    def run(self, config: AnalysisConfig) -> AnalysisResult:
        self.geometry.ensure_dependencies()
        PlotStyle.apply()

        if config.reference and config.reference.pdb_model_id and not config.reference.pdb_file:
            raise ValueError("--reference-pdb-model requires --reference-pdb.")

        pdb_file = Path(config.pdb_file).resolve()
        if not pdb_file.is_file():
            raise FileNotFoundError(f"Input PDB file not found: {pdb_file}")

        raw_df = self.geometry.load_structure(pdb_file, n_jobs=config.n_jobs)
        raw_df = self.geometry.filter_chains(raw_df, config.chains)
        raw_df = self.geometry.normalize(raw_df)
        residue_summary_df = self.geometry.summarize(
            raw_df,
            dmax_outlier_fraction=config.dmax_outlier_fraction,
        )
        model_summary_df, overall_model_summary_df = self.geometry.build_model_summary(
            raw_df, residue_summary_df
        )

        distance_result = self._build_distance_result(config, raw_df)
        residue_clustering = self._build_residue_clustering(config, raw_df)
        residue_range_clustering = self._build_residue_range_clustering(config, raw_df)

        result = AnalysisResult(
            pdb_file=pdb_file,
            raw_df=raw_df,
            residue_summary_df=residue_summary_df,
            model_summary_df=model_summary_df,
            overall_model_summary_df=overall_model_summary_df,
            distance_result=distance_result,
            residue_clustering=residue_clustering,
            residue_range_clustering=residue_range_clustering,
        )

        writer = self.output_writer_cls(config.output)
        result.outputs = writer.write(
            result,
            max_models_in_plot=config.max_models_in_plot,
            hide_model_traces=config.hide_model_traces,
        )
        return result

    def _build_distance_result(self, config: AnalysisConfig, raw_df):
        if config.reference is None:
            return None

        if config.reference.model_id:
            reference_rows, reference_model_label = self.distances.select_reference_rows(
                raw_df, config.reference.model_id
            )
            reference_label = f"input model {reference_model_label}"
        elif config.reference.pdb_file:
            reference_pdb = Path(config.reference.pdb_file).resolve()
            if not reference_pdb.is_file():
                raise FileNotFoundError(f"Reference PDB file not found: {reference_pdb}")
            reference_df = self.geometry.load_structure(reference_pdb, n_jobs=config.n_jobs)
            reference_df = self.geometry.filter_chains(reference_df, config.chains)
            reference_df = self.geometry.normalize(reference_df)
            reference_rows, reference_model_label = self.distances.select_reference_rows(
                reference_df, config.reference.pdb_model_id
            )
            reference_label = f"{reference_pdb.name} model {reference_model_label}"
        else:
            return None

        distance_long_df, distance_summary_df = self.distances.compute(
            raw_df=raw_df,
            reference_df=reference_rows,
            reference_label=reference_label,
        )
        return DistanceResult(
            long_df=distance_long_df,
            summary_df=distance_summary_df,
            reference_label=reference_label,
        )

    def _build_residue_clustering(self, config: AnalysisConfig, raw_df):
        if not config.clustering.cluster_residues:
            return None
        assignments_df, summary_df = self.clustering.cluster_residues(
            raw_df=raw_df,
            min_cluster_size=config.clustering.min_cluster_size,
            min_samples=config.clustering.min_samples,
        )
        return ResidueClusteringResult(assignments_df=assignments_df, summary_df=summary_df)

    def _build_residue_range_clustering(self, config: AnalysisConfig, raw_df):
        if not config.clustering.cluster_residue_ranges:
            return None
        assignments_df, summary_df = self.clustering.cluster_residue_ranges(
            raw_df=raw_df,
            range_texts=config.clustering.cluster_residue_ranges,
            min_cluster_size=config.clustering.min_cluster_size,
            min_samples=config.clustering.min_samples,
        )
        return ResidueRangeClusteringResult(
            assignments_df=assignments_df,
            summary_df=summary_df,
        )
