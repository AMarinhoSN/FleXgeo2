from __future__ import annotations

import argparse
import math
from pathlib import Path

from flexgeo2.config import AnalysisConfig, ClusteringConfig, OutputConfig, ReferenceConfig
from flexgeo2.pipeline import FlexGeo2App


def fraction_in_unit_interval(value: str) -> float:
    try:
        fraction = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number in [0, 1).") from exc
    if not math.isfinite(fraction) or fraction < 0.0 or fraction >= 1.0:
        raise argparse.ArgumentTypeError("must be a number in [0, 1).")
    return fraction


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute Melodia differential geometry descriptors from a PDB file "
            "and generate ensemble-aware curvature and torsion outputs."
        )
    )
    parser.add_argument("pdb_file", type=Path, help="Input PDB file.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for CSV and plot outputs. Default: results",
    )
    parser.add_argument(
        "--chain",
        action="append",
        dest="chains",
        help="Chain ID to keep. Can be repeated, e.g. --chain A --chain B",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of workers passed to Melodia for multi-model files.",
    )
    parser.add_argument(
        "--max-models-in-plot",
        type=int,
        default=12,
        help="Maximum number of individual model traces to overlay per chain plot.",
    )
    parser.add_argument(
        "--hide-model-traces",
        action="store_true",
        help="Plot only the ensemble mean and standard deviation band.",
    )
    parser.add_argument(
        "--dmax-outlier-fraction",
        type=fraction_in_unit_interval,
        default=0.01,
        help=(
            "Extreme histogram bins with less than this fraction of a residue's observations "
            "are ignored when computing dmax. Default: 0.01"
        ),
    )
    reference_group = parser.add_mutually_exclusive_group()
    reference_group.add_argument(
        "--reference-model",
        help=(
            "Model identifier from the input ensemble to use as the reference state "
            "for curvature/torsion distance calculations."
        ),
    )
    reference_group.add_argument(
        "--reference-pdb",
        type=Path,
        help="External PDB file providing the reference state for distance calculations.",
    )
    parser.add_argument(
        "--reference-pdb-model",
        help=(
            "Model identifier to use from --reference-pdb. Defaults to the first model "
            "found in that file."
        ),
    )
    parser.add_argument(
        "--cluster-residues",
        action="store_true",
        help="Run HDBSCAN independently for each residue in curvature/torsion space.",
    )
    parser.add_argument(
        "--cluster-min-size",
        type=int,
        default=5,
        help="Minimum cluster size passed to HDBSCAN. Default: 5",
    )
    parser.add_argument(
        "--cluster-min-samples",
        type=int,
        default=None,
        help="Optional min_samples value passed to HDBSCAN.",
    )
    parser.add_argument(
        "--cluster-residue-range",
        action="append",
        dest="cluster_residue_ranges",
        help=(
            "Residue range for one combined HDBSCAN solution, formatted as START-END "
            "(for example 45-54). Can be repeated."
        ),
    )
    parser.add_argument(
        "--output-verbose",
        action="store_true",
        help="Write detailed intermediate tables and per-chain output folders.",
    )
    return parser


def build_config(args: argparse.Namespace) -> AnalysisConfig:
    reference = None
    if args.reference_model or args.reference_pdb:
        reference = ReferenceConfig(
            model_id=args.reference_model,
            pdb_file=args.reference_pdb,
            pdb_model_id=args.reference_pdb_model,
        )

    clustering = ClusteringConfig(
        cluster_residues=args.cluster_residues,
        cluster_residue_ranges=args.cluster_residue_ranges or [],
        min_cluster_size=args.cluster_min_size,
        min_samples=args.cluster_min_samples,
    )

    output = OutputConfig(
        output_dir=args.output_dir,
        verbose=args.output_verbose,
        write_files=True,
    )

    return AnalysisConfig(
        pdb_file=args.pdb_file,
        chains=args.chains,
        n_jobs=args.n_jobs,
        max_models_in_plot=args.max_models_in_plot,
        hide_model_traces=args.hide_model_traces,
        dmax_outlier_fraction=args.dmax_outlier_fraction,
        reference=reference,
        clustering=clustering,
        output=output,
    )


def print_run_summary(result) -> None:
    model_count = int(result.raw_df["model"].nunique())
    chain_count = int(result.raw_df["chain"].nunique())
    residue_count = int(result.raw_df.groupby(["chain", "order"]).ngroups)
    outputs = result.outputs

    print(f"Processed: {result.pdb_file}")
    print(f"Models: {model_count}")
    print(f"Chains: {chain_count}")
    print(f"Residues: {residue_count}")
    print(f"Raw descriptors: {outputs.raw_csv}")
    print(f"Residue summary: {outputs.residue_summary_csv}")
    print(f"Overall model summary: {outputs.overall_model_summary_csv}")
    print(f"Overview plot: {outputs.overview_plot}")
    if outputs.model_summary_csv is not None:
        print(f"Model summary by chain: {outputs.model_summary_csv}")
    if outputs.distance_summary_csv is not None:
        print(f"Distance summary: {outputs.distance_summary_csv}")
        print(f"Distance heatmap: {outputs.distance_heatmap}")
    if outputs.distance_long_csv is not None:
        print(f"Distance details: {outputs.distance_long_csv}")
        print(f"Distance matrices by chain: {outputs.distance_matrix_dir}")
    if outputs.cluster_summary_csv is not None:
        print(f"Residue cluster summary: {outputs.cluster_summary_csv}")
        print(f"Residue cluster plots: {outputs.cluster_plots_dir}")
    if outputs.cluster_assignments_csv is not None:
        print(f"Residue cluster assignments: {outputs.cluster_assignments_csv}")
    if outputs.range_cluster_summary_csv is not None:
        print(f"Residue-range cluster summary: {outputs.range_cluster_summary_csv}")
        print(f"Residue-range cluster plots: {outputs.range_cluster_plots_dir}")
    if outputs.range_cluster_assignments_csv is not None:
        print(f"Residue-range cluster assignments: {outputs.range_cluster_assignments_csv}")
    if outputs.chains_dir is not None:
        print(f"Per-chain outputs: {outputs.chains_dir}")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    app = FlexGeo2App()
    result = app.run(build_config(args))
    print_run_summary(result)


if __name__ == "__main__":
    main()
