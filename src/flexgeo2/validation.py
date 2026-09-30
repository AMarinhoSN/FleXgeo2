"""Up-front validation so user mistakes surface before Melodia runs."""

from __future__ import annotations

from pathlib import Path

from flexgeo2.clustering import ClusteringService
from flexgeo2.config import AnalysisConfig
from flexgeo2.geometry import GeometryService, StructureInfo
from flexgeo2.outputs import check_output_dir
from flexgeo2.selection import parse_residue_selections, select_residues


def validate_config(config: AnalysisConfig) -> None:
    """Check options that do not depend on the structure contents."""
    pdb_file = Path(config.pdb_file).resolve()
    if not pdb_file.is_file():
        raise FileNotFoundError(f"Input PDB file not found: {pdb_file}")

    if config.n_jobs == 0:
        raise ValueError("n_jobs must be a positive integer or -1 (all CPUs), not 0.")
    if config.max_models_in_plot < 0:
        raise ValueError("max_models_in_plot must be zero or a positive integer.")
    GeometryService._validate_dmax_outlier_fraction(config.dmax_outlier_fraction)

    reference = config.reference
    if reference is not None:
        if reference.model_id is not None and reference.pdb_file is not None:
            raise ValueError(
                "Choose either a reference model from the input (--reference-model) or an "
                "external reference PDB (--reference-pdb), not both."
            )
        if reference.pdb_model_id is not None and reference.pdb_file is None:
            raise ValueError("--reference-pdb-model requires --reference-pdb.")
        if reference.pdb_file is not None:
            reference_pdb = Path(reference.pdb_file).resolve()
            if not reference_pdb.is_file():
                raise FileNotFoundError(f"Reference PDB file not found: {reference_pdb}")

    clustering = config.clustering
    if clustering.cluster_residues or clustering.cluster_residue_ranges:
        if clustering.min_cluster_size < 2:
            raise ValueError("min_cluster_size must be at least 2.")
        if clustering.min_samples is not None and clustering.min_samples < 1:
            raise ValueError("min_samples must be a positive integer.")
    for range_text in clustering.cluster_residue_ranges:
        ClusteringService.parse_residue_range(range_text)

    parse_residue_selections(config.output.plot_residues)

    if config.output.distance_matrices and reference is None:
        raise ValueError(
            "Distance matrices need a reference (--reference-model or --reference-pdb)."
        )

    check_output_dir(config.output)


def _check_model(model_id: str, info: StructureInfo, source: str) -> None:
    available = [str(model) for model in info.model_ids]
    if str(model_id) not in available:
        raise ValueError(
            f"Reference model '{model_id}' was not found in {source}. "
            f"Available models: {_format_ids(available)}"
        )


def _format_ids(ids: list[str], limit: int = 20) -> str:
    if len(ids) <= limit:
        return ", ".join(ids)
    return f"{', '.join(ids[:limit])}, ... ({len(ids)} total)"


def validate_against_structure(
    config: AnalysisConfig,
    info: StructureInfo,
    reference_info: StructureInfo | None = None,
) -> None:
    """Check chains, reference models and residue selections against parsed structures."""
    available_chains = sorted(info.residues_by_chain)
    if config.chains:
        missing = [chain for chain in config.chains if chain not in info.residues_by_chain]
        if missing:
            raise ValueError(
                f"Chain(s) not found in the input: {', '.join(missing)}. "
                f"Available chains: {', '.join(available_chains)}"
            )
        selected_chains = list(config.chains)
    else:
        selected_chains = available_chains

    reference = config.reference
    if reference is not None:
        if reference.model_id is not None:
            _check_model(reference.model_id, info, "the input ensemble")
        if reference_info is not None and reference.pdb_model_id is not None:
            _check_model(reference.pdb_model_id, reference_info, "the reference PDB")

    for range_text in config.clustering.cluster_residue_ranges:
        start, end = ClusteringService.parse_residue_range(range_text)
        expected = set(range(start, end + 1))
        complete_somewhere = False
        for chain in selected_chains:
            present = expected & info.residues_by_chain.get(chain, set())
            if not present:
                continue
            if present != expected:
                missing_residues = sorted(expected - present)
                raise ValueError(
                    f"Residue range {start}-{end} on chain '{chain}' is incomplete; "
                    f"missing residue(s): {_format_ids([str(r) for r in missing_residues])}."
                )
            complete_somewhere = True
        if not complete_somewhere:
            raise ValueError(
                f"Residue range {start}-{end} does not match any residues in the selected "
                f"chain(s): {', '.join(selected_chains)}."
            )

    select_residues(
        parse_residue_selections(config.output.plot_residues),
        {chain: info.residues_by_chain.get(chain, set()) for chain in selected_chains},
    )
