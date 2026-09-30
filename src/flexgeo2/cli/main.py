from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from flexgeo2.config import (
    PLOT_FORMATS,
    AnalysisConfig,
    ClusteringConfig,
    OutputConfig,
    ReferenceConfig,
)
from flexgeo2.outputs import OutputDirectoryNotEmptyError
from flexgeo2.pipeline import FlexGeo2App
from flexgeo2.report import render_terminal_summary


def fraction_in_unit_interval(value: str) -> float:
    try:
        fraction = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number in [0, 1).") from exc
    if not math.isfinite(fraction) or fraction < 0.0 or fraction >= 1.0:
        raise argparse.ArgumentTypeError("must be a number in [0, 1).")
    return fraction


def int_at_least(minimum: int):
    def parse(value: str) -> int:
        try:
            number = int(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"must be an integer >= {minimum}.") from exc
        if number < minimum:
            raise argparse.ArgumentTypeError(f"must be an integer >= {minimum}.")
        return number

    return parse


def n_jobs_value(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer or -1.") from exc
    if number == 0 or number < -1:
        raise argparse.ArgumentTypeError("must be a positive integer or -1.")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="flexgeo2",
        description=(
            "Compute Melodia differential geometry descriptors from a PDB file "
            "and generate ensemble-aware curvature and torsion outputs."
        ),
    )
    parser.add_argument("pdb_file", type=Path, help="Input PDB file.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for CSV and plot outputs. Default: results",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Replace the outputs of an earlier run in --output-dir. Without it, FleXgeo2 "
            "refuses to write into a folder that is not empty. Other files are kept."
        ),
    )
    parser.add_argument(
        "--chain",
        action="append",
        dest="chains",
        help="Chain ID to keep. Can be repeated, e.g. --chain A --chain B",
    )
    parser.add_argument(
        "--n-jobs",
        type=n_jobs_value,
        default=1,
        help="Number of workers passed to Melodia for multi-model files (-1 = all CPUs).",
    )
    parser.add_argument(
        "--max-models-in-plot",
        type=int_at_least(0),
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
            "PDB MODEL number from the input ensemble to use as the reference state "
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
            "PDB MODEL number to use from --reference-pdb. Defaults to the first model "
            "found in that file. Requires --reference-pdb."
        ),
    )
    parser.add_argument(
        "--cluster-residues",
        action="store_true",
        help="Run HDBSCAN independently for each residue in curvature/torsion space.",
    )
    parser.add_argument(
        "--cluster-min-size",
        type=int_at_least(2),
        default=None,
        help="Smallest cluster HDBSCAN reports, in models. Default: 5%% of the models, at least 5.",
    )
    parser.add_argument(
        "--cluster-min-samples",
        type=int_at_least(1),
        default=None,
        help="HDBSCAN min_samples: larger values label more models as noise. Default: 5, "
        "or --cluster-min-size if smaller.",
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
        "--plot-residues",
        action="append",
        dest="plot_residues",
        metavar="SELECTION",
        help=(
            "Plot curvature vs torsion for chosen residues, one point per model, e.g. 45, "
            "45-50, A:45 or A:45-50 (comma-separated; can be repeated). Points are coloured "
            "by cluster with --cluster-residues, and the reference is marked when given."
        ),
    )
    parser.add_argument(
        "--plot-format",
        choices=PLOT_FORMATS,
        default="png",
        help=(
            "File format of every figure. pdf and svg are vector formats with editable "
            "text, for publication. Default: png"
        ),
    )
    parser.add_argument(
        "--distance-matrices",
        action="store_true",
        help=(
            "Also write the distances to the reference as one models x residues CSV per "
            "chain (reference/matrices/). Requires --reference-model or --reference-pdb."
        ),
    )
    # Removed option, kept hidden so old command lines get a pointer to its replacements.
    parser.add_argument("--output-verbose", action="store_true", help=argparse.SUPPRESS)
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
        distance_matrices=args.distance_matrices,
        plot_residues=args.plot_residues or [],
        plot_format=args.plot_format,
        write_files=True,
        overwrite=args.overwrite,
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
    print(render_terminal_summary(result))


def parse_args(parser: argparse.ArgumentParser, argv: list[str] | None = None):
    args = parser.parse_args(argv)
    if args.output_verbose:
        parser.error(
            "--output-verbose was removed. Distance matrices: --distance-matrices. "
            "geometry/models_by_chain.csv is now written whenever there is more than one "
            "chain. Per-residue cluster plots are replaced by clusters/clusters.png and "
            "--plot-residues for chosen residues."
        )
    if args.reference_pdb_model is not None and args.reference_pdb is None:
        parser.error("--reference-pdb-model requires --reference-pdb.")
    if args.distance_matrices and args.reference_model is None and args.reference_pdb is None:
        parser.error("--distance-matrices requires --reference-model or --reference-pdb.")
    return args


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parse_args(parser, argv)
    try:
        result = FlexGeo2App().run(build_config(args))
    except OutputDirectoryNotEmptyError as exc:
        print(
            f"{parser.prog}: error: output folder {exc.output_dir} is not empty. "
            "Use --overwrite to replace the outputs of an earlier run, or choose another "
            "--output-dir.",
            file=sys.stderr,
        )
        return 1
    except (FileNotFoundError, ValueError) as exc:
        print(f"{parser.prog}: error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print(f"{parser.prog}: interrupted", file=sys.stderr)
        return 130
    print_run_summary(result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
