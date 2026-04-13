from __future__ import annotations

import argparse
import math
from pathlib import Path

REQUIRED_COLUMNS = {"model", "chain", "order", "name", "curvature", "torsion"}


def parse_args() -> argparse.Namespace:
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
    return parser.parse_args()


def ensure_dependencies() -> None:
    try:
        import hdbscan  # noqa: F401
        import matplotlib  # noqa: F401
        import melodia_py  # noqa: F401
        import pandas  # noqa: F401
    except ModuleNotFoundError as exc:
        missing = exc.name or "required dependency"
        raise SystemExit(
            f"Missing dependency: {missing}. Install the project dependencies first, "
            "for example with: pip install -e ."
        ) from exc


def compute_geometry(pdb_file: Path, n_jobs: int):
    import melodia_py as mel

    df = mel.geometry_from_structure_file(str(pdb_file), n_jobs=n_jobs)
    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        missing_csv = ", ".join(sorted(missing))
        raise ValueError(
            f"Melodia output is missing required columns: {missing_csv}. "
            f"Available columns: {', '.join(df.columns)}"
        )
    return df.copy()


def filter_chains(df, chains: list[str] | None):
    if not chains:
        return df
    filtered = df[df["chain"].isin(chains)].copy()
    if filtered.empty:
        raise ValueError(f"No residues found for requested chain(s): {', '.join(chains)}")
    return filtered


def normalise_geometry_table(df):
    normalised = df.copy()
    normalised["order"] = normalised["order"].astype(int)
    normalised["residue_label"] = [
        f"{name}{int(order)}" for order, name in zip(normalised["order"], normalised["name"])
    ]
    return normalised.sort_values(["chain", "model", "order"]).reset_index(drop=True)


def summarise_geometry(df):
    grouped = (
        df.groupby(["chain", "order", "name", "residue_label"], dropna=False)
        .agg(
            curvature_mean=("curvature", "mean"),
            curvature_std=("curvature", "std"),
            curvature_min=("curvature", "min"),
            curvature_max=("curvature", "max"),
            torsion_mean=("torsion", "mean"),
            torsion_std=("torsion", "std"),
            torsion_min=("torsion", "min"),
            torsion_max=("torsion", "max"),
            models=("model", "nunique"),
        )
        .reset_index()
        .sort_values(["chain", "order"])
    )
    return grouped.fillna({"curvature_std": 0.0, "torsion_std": 0.0})


def build_model_summary(raw_df, residue_summary_df):
    merged = raw_df.merge(
        residue_summary_df[
            [
                "chain",
                "order",
                "curvature_mean",
                "torsion_mean",
            ]
        ],
        on=["chain", "order"],
        how="left",
    )
    merged["curvature_abs_deviation"] = (
        merged["curvature"] - merged["curvature_mean"]
    ).abs()
    merged["torsion_abs_deviation"] = (merged["torsion"] - merged["torsion_mean"]).abs()

    model_summary = (
        merged.groupby(["chain", "model"], dropna=False)
        .agg(
            residues=("order", "count"),
            curvature_mean=("curvature", "mean"),
            torsion_mean=("torsion", "mean"),
            curvature_std=("curvature", "std"),
            torsion_std=("torsion", "std"),
            curvature_mean_abs_deviation=("curvature_abs_deviation", "mean"),
            torsion_mean_abs_deviation=("torsion_abs_deviation", "mean"),
        )
        .reset_index()
        .sort_values(["chain", "model"])
    )

    overall_summary = (
        merged.groupby("model", dropna=False)
        .agg(
            residues=("order", "count"),
            curvature_mean=("curvature", "mean"),
            torsion_mean=("torsion", "mean"),
            curvature_mean_abs_deviation=("curvature_abs_deviation", "mean"),
            torsion_mean_abs_deviation=("torsion_abs_deviation", "mean"),
        )
        .reset_index()
        .sort_values("model")
    )

    model_summary = model_summary.fillna({"curvature_std": 0.0, "torsion_std": 0.0})
    overall_summary = overall_summary.fillna(0.0)
    return model_summary, overall_summary


def select_reference_rows(df, model_id):
    available_models = [str(model) for model in df["model"].drop_duplicates().tolist()]
    if model_id is None:
        first_model = df["model"].drop_duplicates().iloc[0]
        return df[df["model"] == first_model].copy(), str(first_model)

    model_mask = df["model"].astype(str) == str(model_id)
    selected = df[model_mask].copy()
    if selected.empty:
        raise ValueError(
            f"Reference model '{model_id}' was not found. Available models: "
            f"{', '.join(available_models)}"
        )
    return selected, str(model_id)


def compute_reference_distance_tables(raw_df, reference_df, reference_label):
    comparison = raw_df.merge(
        reference_df[
            [
                "chain",
                "order",
                "residue_label",
                "curvature",
                "torsion",
            ]
        ].rename(
            columns={
                "curvature": "reference_curvature",
                "torsion": "reference_torsion",
            }
        ),
        on=["chain", "order", "residue_label"],
        how="inner",
    )

    if comparison.empty:
        raise ValueError(
            "No overlapping residues were found between the ensemble and the reference state."
        )

    comparison["distance_to_reference"] = (
        (
            comparison["curvature"] - comparison["reference_curvature"]
        ) ** 2
        + (comparison["torsion"] - comparison["reference_torsion"]) ** 2
    ) ** 0.5
    comparison["reference_label"] = reference_label

    distance_long = comparison[
        [
            "chain",
            "model",
            "order",
            "name",
            "residue_label",
            "curvature",
            "torsion",
            "reference_curvature",
            "reference_torsion",
            "distance_to_reference",
            "reference_label",
        ]
    ].sort_values(["chain", "model", "order"])

    distance_summary = (
        distance_long.groupby(["chain", "order", "name", "residue_label"], dropna=False)
        .agg(
            distance_mean=("distance_to_reference", "mean"),
            distance_std=("distance_to_reference", "std"),
            distance_min=("distance_to_reference", "min"),
            distance_max=("distance_to_reference", "max"),
            models=("model", "nunique"),
        )
        .reset_index()
        .sort_values(["chain", "order"])
        .fillna({"distance_std": 0.0})
    )

    return distance_long, distance_summary


def cluster_residues(raw_df, min_cluster_size: int, min_samples: int | None):
    import hdbscan
    import numpy as np
    import pandas as pd

    cluster_frames = []
    summary_rows = []

    for (chain, order, name, residue_label), residue_df in raw_df.groupby(
        ["chain", "order", "name", "residue_label"], dropna=False
    ):
        residue_points = residue_df[["curvature", "torsion"]].copy()
        if len(residue_points) < max(2, min_cluster_size):
            residue_result = residue_df.copy()
            residue_result["cluster"] = -1
            residue_result["cluster_probability"] = 0.0
            cluster_frames.append(residue_result)
            summary_rows.append(
                {
                    "chain": chain,
                    "order": order,
                    "name": name,
                    "residue_label": residue_label,
                    "n_conformations": len(residue_df),
                    "n_clusters": 0,
                    "noise_fraction": 1.0,
                }
            )
            continue

        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
        )
        feature_matrix = np.array(residue_points.to_numpy(), dtype=float, order="C", copy=True)
        labels = clusterer.fit_predict(feature_matrix)
        probabilities = getattr(clusterer, "probabilities_", None)

        residue_result = residue_df.copy()
        residue_result["cluster"] = labels
        residue_result["cluster_probability"] = (
            probabilities if probabilities is not None else 0.0
        )
        cluster_frames.append(residue_result)

        non_noise_clusters = sorted({int(label) for label in labels if int(label) >= 0})
        noise_fraction = float((residue_result["cluster"] == -1).mean())
        summary_rows.append(
            {
                "chain": chain,
                "order": order,
                "name": name,
                "residue_label": residue_label,
                "n_conformations": len(residue_df),
                "n_clusters": len(non_noise_clusters),
                "noise_fraction": noise_fraction,
            }
        )

    clustered_df = (
        pd.concat(cluster_frames, ignore_index=True)
        .sort_values(["chain", "order", "model"])
        .reset_index(drop=True)
    )
    cluster_summary_df = pd.DataFrame(summary_rows).sort_values(["chain", "order"])
    return clustered_df, cluster_summary_df


def parse_residue_range(range_text: str) -> tuple[int, int]:
    cleaned = range_text.strip()
    if "-" not in cleaned:
        raise ValueError(
            f"Invalid residue range '{range_text}'. Expected format START-END, e.g. 45-54."
        )
    start_text, end_text = cleaned.split("-", 1)
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError as exc:
        raise ValueError(
            f"Invalid residue range '{range_text}'. Start and end must be integers."
        ) from exc
    if end < start:
        raise ValueError(
            f"Invalid residue range '{range_text}'. End must be greater than or equal to start."
        )
    return start, end


def compute_pca_projection(matrix):
    import numpy as np

    centered = matrix - matrix.mean(axis=0, keepdims=True)
    if centered.shape[0] < 2:
        return np.column_stack([centered[:, 0], np.zeros(centered.shape[0])])

    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    if vh.shape[0] == 1:
        pc1 = centered @ vh[0]
        pc2 = np.zeros(centered.shape[0])
    else:
        pc1 = centered @ vh[0]
        pc2 = centered @ vh[1]
    return np.column_stack([pc1, pc2])


def cluster_residue_ranges(
    raw_df,
    range_texts: list[str],
    min_cluster_size: int,
    min_samples: int | None,
):
    import hdbscan
    import numpy as np
    import pandas as pd

    assignment_frames = []
    summary_rows = []

    parsed_ranges = [parse_residue_range(range_text) for range_text in range_texts]

    for chain, chain_df in raw_df.groupby("chain", dropna=False):
        for start_order, end_order in parsed_ranges:
            range_df = chain_df[
                (chain_df["order"] >= start_order) & (chain_df["order"] <= end_order)
            ].copy()
            if range_df.empty:
                continue

            residue_orders = sorted(range_df["order"].drop_duplicates().tolist())
            expected_orders = list(range(start_order, end_order + 1))
            if residue_orders != expected_orders:
                raise ValueError(
                    f"Residue range {start_order}-{end_order} on chain '{chain}' is incomplete "
                    "after filtering. Make sure all residues in the range are present."
                )

            feature_table = (
                range_df.pivot_table(
                    index="model",
                    columns="order",
                    values=["curvature", "torsion"],
                    aggfunc="first",
                )
                .sort_index(axis=0)
                .sort_index(axis=1)
            )
            if feature_table.isna().any().any():
                raise ValueError(
                    f"Residue range {start_order}-{end_order} on chain '{chain}' contains "
                    "missing curvature/torsion values for at least one conformation."
                )

            feature_matrix = np.array(feature_table.to_numpy(), dtype=float, order="C", copy=True)
            range_label = f"{start_order}-{end_order}"
            if len(feature_table) < max(2, min_cluster_size):
                labels = [-1] * len(feature_table)
                probabilities = [0.0] * len(feature_table)
            else:
                clusterer = hdbscan.HDBSCAN(
                    min_cluster_size=min_cluster_size,
                    min_samples=min_samples,
                )
                labels = clusterer.fit_predict(feature_matrix)
                probabilities = getattr(clusterer, "probabilities_", [0.0] * len(feature_table))

            projection = compute_pca_projection(feature_matrix)
            assignment_df = pd.DataFrame(
                {
                    "chain": chain,
                    "range_start": start_order,
                    "range_end": end_order,
                    "range_label": range_label,
                    "model": feature_table.index.astype(str),
                    "cluster": labels,
                    "cluster_probability": probabilities,
                    "pc1": projection[:, 0],
                    "pc2": projection[:, 1],
                }
            )
            assignment_frames.append(assignment_df)

            non_noise_clusters = sorted(
                {int(label) for label in assignment_df["cluster"].tolist() if int(label) >= 0}
            )
            summary_rows.append(
                {
                    "chain": chain,
                    "range_start": start_order,
                    "range_end": end_order,
                    "range_label": range_label,
                    "n_conformations": len(feature_table),
                    "n_residues": len(expected_orders),
                    "n_clusters": len(non_noise_clusters),
                    "noise_fraction": float((assignment_df["cluster"] == -1).mean()),
                }
            )

    if not assignment_frames:
        raise ValueError("No valid residue-range clustering solutions could be generated.")

    assignments_df = pd.concat(assignment_frames, ignore_index=True).sort_values(
        ["chain", "range_start", "model"]
    )
    summary_df = pd.DataFrame(summary_rows).sort_values(["chain", "range_start"])
    return assignments_df, summary_df


def sanitize_chain_id(chain) -> str:
    if chain is None or chain == "":
        return "unassigned"
    return str(chain).replace("/", "_")


def configure_plot_style() -> None:
    import matplotlib.pyplot as plt

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#c7c7c7",
            "axes.titleweight": "bold",
            "axes.labelsize": 11,
            "axes.titlesize": 13,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.frameon": False,
        }
    )


def apply_residue_ticks(axis, x_values, residue_labels) -> None:
    if len(residue_labels) == 0:
        return
    tick_step = max(1, math.ceil(len(residue_labels) / 18))
    axis.set_xticks(x_values[::tick_step])
    axis.set_xticklabels(residue_labels[::tick_step], rotation=45, ha="right")


def plot_chain_geometry(
    chain_raw_df,
    chain_summary_df,
    output_path: Path,
    show_model_traces: bool,
    max_models_in_plot: int,
) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), constrained_layout=True)

    x_values = chain_summary_df["order"].to_numpy()
    residue_labels = chain_summary_df["residue_label"].tolist()
    chain_id = chain_summary_df["chain"].iloc[0] if not chain_summary_df.empty else None
    title_suffix = f"Chain {chain_id}" if chain_id not in (None, "") else "Chain"

    if show_model_traces:
        model_ids = list(chain_raw_df["model"].drop_duplicates())[:max_models_in_plot]
        for model_id in model_ids:
            model_df = chain_raw_df[chain_raw_df["model"] == model_id]
            axes[0].plot(
                model_df["order"],
                model_df["curvature"],
                color="#4c78a8",
                alpha=0.22,
                linewidth=1,
            )
            axes[1].plot(
                model_df["order"],
                model_df["torsion"],
                color="#e45756",
                alpha=0.22,
                linewidth=1,
            )

    axes[0].fill_between(
        x_values,
        chain_summary_df["curvature_mean"] - chain_summary_df["curvature_std"],
        chain_summary_df["curvature_mean"] + chain_summary_df["curvature_std"],
        color="#4c78a8",
        alpha=0.18,
        label="Ensemble SD",
    )
    axes[0].plot(
        x_values,
        chain_summary_df["curvature_mean"],
        color="#1f4e79",
        linewidth=2.5,
        label="Ensemble mean",
    )
    axes[0].set_title(f"{title_suffix}: Curvature")
    axes[0].set_ylabel("Curvature")
    axes[0].legend(loc="upper right")

    axes[1].fill_between(
        x_values,
        chain_summary_df["torsion_mean"] - chain_summary_df["torsion_std"],
        chain_summary_df["torsion_mean"] + chain_summary_df["torsion_std"],
        color="#e45756",
        alpha=0.18,
        label="Ensemble SD",
    )
    axes[1].plot(
        x_values,
        chain_summary_df["torsion_mean"],
        color="#b22222",
        linewidth=2.5,
        label="Ensemble mean",
    )
    axes[1].set_title(f"{title_suffix}: Torsion")
    axes[1].set_ylabel("Torsion")
    axes[1].set_xlabel("Residue")
    axes[1].legend(loc="upper right")

    for axis in axes:
        axis.grid(alpha=0.3)
        apply_residue_ticks(axis, x_values, residue_labels)

    fig.suptitle("Backbone Differential Geometry", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_overview(summary_df, output_path: Path) -> None:
    import matplotlib.pyplot as plt

    chains = list(summary_df["chain"].drop_duplicates())
    fig, axes = plt.subplots(
        nrows=len(chains),
        ncols=2,
        figsize=(14, max(4, 3.8 * len(chains))),
        constrained_layout=True,
    )

    if len(chains) == 1:
        axes = [axes]

    for row_axes, chain in zip(axes, chains):
        chain_df = summary_df[summary_df["chain"] == chain]
        x_values = chain_df["order"].to_numpy()
        residue_labels = chain_df["residue_label"].tolist()
        title_suffix = f"Chain {chain}" if chain not in (None, "") else "Chain"

        row_axes[0].plot(x_values, chain_df["curvature_mean"], color="#1f4e79", linewidth=2)
        row_axes[0].fill_between(
            x_values,
            chain_df["curvature_mean"] - chain_df["curvature_std"],
            chain_df["curvature_mean"] + chain_df["curvature_std"],
            color="#4c78a8",
            alpha=0.18,
        )
        row_axes[0].set_title(f"{title_suffix}: Curvature")
        row_axes[0].set_ylabel("Curvature")

        row_axes[1].plot(x_values, chain_df["torsion_mean"], color="#b22222", linewidth=2)
        row_axes[1].fill_between(
            x_values,
            chain_df["torsion_mean"] - chain_df["torsion_std"],
            chain_df["torsion_mean"] + chain_df["torsion_std"],
            color="#e45756",
            alpha=0.18,
        )
        row_axes[1].set_title(f"{title_suffix}: Torsion")
        row_axes[1].set_ylabel("Torsion")

        for axis in row_axes:
            axis.set_xlabel("Residue")
            axis.grid(alpha=0.3)
            apply_residue_ticks(axis, x_values, residue_labels)

    fig.suptitle("Ensemble Overview", fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_distance_heatmap(distance_long_df, output_path: Path, title: str) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd

    chains = list(distance_long_df["chain"].drop_duplicates())
    fig, axes = plt.subplots(
        nrows=len(chains),
        ncols=1,
        figsize=(14, max(4, 3.8 * len(chains))),
        constrained_layout=True,
    )

    if len(chains) == 1:
        axes = [axes]

    for axis, chain in zip(axes, chains):
        chain_df = distance_long_df[distance_long_df["chain"] == chain].copy()
        matrix = (
            chain_df.pivot(index="model", columns="residue_label", values="distance_to_reference")
            .sort_index()
        )
        matrix = matrix.apply(pd.to_numeric, errors="coerce")
        if matrix.empty or matrix.isna().all().all():
            raise ValueError(
                f"Distance heatmap for chain '{chain}' is empty after alignment. "
                "This usually means the chosen reference does not overlap with the "
                "ensemble on chain/residue numbering."
            )

        image = axis.imshow(
            matrix.to_numpy().T,
            aspect="auto",
            cmap="magma",
            interpolation="nearest",
            origin="lower",
        )
        axis.set_title(
            f"Chain {chain}: Distance to reference"
            if chain not in (None, "")
            else "Chain: Distance to reference"
        )
        axis.set_xlabel("Conformation")
        axis.set_ylabel("Residue")

        model_labels = [str(model) for model in matrix.index]
        if model_labels:
            model_tick_step = max(1, math.ceil(len(model_labels) / 20))
            model_tick_positions = list(range(0, len(model_labels), model_tick_step))
            axis.set_xticks(model_tick_positions)
            axis.set_xticklabels(
                [model_labels[index] for index in model_tick_positions],
                rotation=45,
                ha="right",
            )

        residue_labels = list(matrix.columns)
        if residue_labels:
            residue_tick_step = max(1, math.ceil(len(residue_labels) / 25))
            residue_tick_positions = list(range(0, len(residue_labels), residue_tick_step))
            axis.set_yticks(residue_tick_positions)
            axis.set_yticklabels([residue_labels[index] for index in residue_tick_positions])

        fig.colorbar(image, ax=axis, fraction=0.024, pad=0.02, label="Euclidean distance")

    fig.suptitle(title, fontsize=16, fontweight="bold")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_distance_matrix_csv(distance_long_df, output_path: Path) -> None:
    matrix_df = (
        distance_long_df.pivot(index="model", columns="residue_label", values="distance_to_reference")
        .sort_index()
    )
    matrix_df.to_csv(output_path)


def plot_residue_cluster_solution(residue_cluster_df, output_path: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    residue_label = residue_cluster_df["residue_label"].iloc[0]
    chain = residue_cluster_df["chain"].iloc[0]
    title_prefix = f"Chain {chain}" if chain not in (None, "") else "Chain"

    unique_clusters = sorted(residue_cluster_df["cluster"].drop_duplicates())
    cmap = plt.get_cmap("tab10")

    for idx, cluster_label in enumerate(unique_clusters):
        cluster_points = residue_cluster_df[residue_cluster_df["cluster"] == cluster_label]
        if int(cluster_label) == -1:
            color = "#9e9e9e"
            legend_label = "Noise"
        else:
            color = cmap(idx % 10)
            legend_label = f"Cluster {int(cluster_label)}"
        ax.scatter(
            cluster_points["curvature"],
            cluster_points["torsion"],
            s=42,
            alpha=0.85,
            c=[color],
            label=legend_label,
            edgecolors="none",
        )

    ax.set_title(f"{title_prefix}: {residue_label}")
    ax.set_xlabel("Curvature")
    ax.set_ylabel("Torsion")
    ax.legend(loc="best")
    ax.grid(alpha=0.3)
    fig.savefig(output_path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def plot_residue_range_cluster_solution(range_cluster_df, output_path: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    chain = range_cluster_df["chain"].iloc[0]
    range_label = range_cluster_df["range_label"].iloc[0]
    title_prefix = f"Chain {chain}" if chain not in (None, "") else "Chain"

    unique_clusters = sorted(range_cluster_df["cluster"].drop_duplicates())
    cmap = plt.get_cmap("tab10")

    for idx, cluster_label in enumerate(unique_clusters):
        cluster_points = range_cluster_df[range_cluster_df["cluster"] == cluster_label]
        if int(cluster_label) == -1:
            color = "#9e9e9e"
            legend_label = "Noise"
        else:
            color = cmap(idx % 10)
            legend_label = f"Cluster {int(cluster_label)}"
        ax.scatter(
            cluster_points["pc1"],
            cluster_points["pc2"],
            s=48,
            alpha=0.88,
            c=[color],
            label=legend_label,
            edgecolors="none",
        )

    for _, row in range_cluster_df.iterrows():
        ax.text(row["pc1"], row["pc2"], str(row["model"]), fontsize=7, alpha=0.7)

    ax.set_title(f"{title_prefix}: residues {range_label}")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(loc="best")
    ax.grid(alpha=0.3)
    fig.savefig(output_path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def write_outputs(
    raw_df,
    residue_summary_df,
    model_summary_df,
    overall_model_summary_df,
    distance_long_df,
    distance_summary_df,
    residue_clusters_df,
    residue_cluster_summary_df,
    residue_range_clusters_df,
    residue_range_cluster_summary_df,
    output_dir: Path,
    show_model_traces: bool,
    max_models_in_plot: int,
    output_verbose: bool,
) -> dict[str, Path]:
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    chains_dir = output_dir / "chains" if output_verbose else None
    if chains_dir is not None:
        chains_dir.mkdir(parents=True, exist_ok=True)

    raw_csv = output_dir / "geometry_descriptors.csv"
    residue_summary_csv = output_dir / "residue_summary.csv"
    model_summary_csv = output_dir / "model_summary_by_chain.csv"
    overall_model_summary_csv = output_dir / "model_summary_overall.csv"
    overview_plot = plots_dir / "ensemble_overview.png"
    distance_long_csv = output_dir / "distance_to_reference_long.csv"
    distance_summary_csv = output_dir / "distance_to_reference_summary.csv"
    distance_heatmap = plots_dir / "distance_to_reference_heatmap.png"
    distance_matrix_dir = output_dir / "distance_matrices"
    cluster_summary_csv = output_dir / "residue_cluster_summary.csv"
    cluster_assignments_csv = output_dir / "residue_cluster_assignments.csv"
    cluster_plots_dir = output_dir / "cluster_plots"
    range_cluster_summary_csv = output_dir / "residue_range_cluster_summary.csv"
    range_cluster_assignments_csv = output_dir / "residue_range_cluster_assignments.csv"
    range_cluster_plots_dir = output_dir / "range_cluster_plots"

    raw_df.to_csv(raw_csv, index=False)
    residue_summary_df.to_csv(residue_summary_csv, index=False)
    overall_model_summary_df.to_csv(overall_model_summary_csv, index=False)
    plot_overview(residue_summary_df, overview_plot)
    if output_verbose:
        model_summary_df.to_csv(model_summary_csv, index=False)
    if distance_long_df is not None and distance_summary_df is not None:
        distance_summary_df.to_csv(distance_summary_csv, index=False)
        reference_label = distance_long_df["reference_label"].iloc[0]
        plot_distance_heatmap(
            distance_long_df=distance_long_df,
            output_path=distance_heatmap,
            title=f"Distance to reference: {reference_label}",
        )
        if output_verbose:
            distance_matrix_dir.mkdir(parents=True, exist_ok=True)
            distance_long_df.to_csv(distance_long_csv, index=False)
            for chain in distance_long_df["chain"].drop_duplicates():
                chain_key = sanitize_chain_id(chain)
                chain_distance_long_df = distance_long_df[
                    distance_long_df["chain"] == chain
                ].copy()
                write_distance_matrix_csv(
                    chain_distance_long_df,
                    distance_matrix_dir / f"{chain_key}_distance_matrix.csv",
                )
    if residue_clusters_df is not None and residue_cluster_summary_df is not None:
        cluster_plots_dir.mkdir(parents=True, exist_ok=True)
        residue_cluster_summary_df.to_csv(cluster_summary_csv, index=False)
        if output_verbose:
            residue_clusters_df.to_csv(cluster_assignments_csv, index=False)
        for (chain, residue_label), residue_cluster_df in residue_clusters_df.groupby(
            ["chain", "residue_label"], dropna=False
        ):
            chain_key = sanitize_chain_id(chain)
            plot_residue_cluster_solution(
                residue_cluster_df=residue_cluster_df,
                output_path=cluster_plots_dir / f"{chain_key}_{residue_label}_clusters.png",
            )
    if residue_range_clusters_df is not None and residue_range_cluster_summary_df is not None:
        range_cluster_plots_dir.mkdir(parents=True, exist_ok=True)
        residue_range_cluster_summary_df.to_csv(range_cluster_summary_csv, index=False)
        if output_verbose:
            residue_range_clusters_df.to_csv(range_cluster_assignments_csv, index=False)
        for (chain, range_label), range_cluster_df in residue_range_clusters_df.groupby(
            ["chain", "range_label"], dropna=False
        ):
            chain_key = sanitize_chain_id(chain)
            safe_range_label = str(range_label).replace("/", "_")
            plot_residue_range_cluster_solution(
                range_cluster_df=range_cluster_df,
                output_path=range_cluster_plots_dir
                / f"{chain_key}_{safe_range_label}_clusters.png",
            )

    if output_verbose and chains_dir is not None:
        for chain in residue_summary_df["chain"].drop_duplicates():
            chain_key = sanitize_chain_id(chain)
            chain_output_dir = chains_dir / chain_key
            chain_output_dir.mkdir(parents=True, exist_ok=True)

            chain_raw_df = raw_df[raw_df["chain"] == chain].copy()
            chain_summary_df = residue_summary_df[residue_summary_df["chain"] == chain].copy()
            chain_model_summary_df = model_summary_df[model_summary_df["chain"] == chain].copy()

            chain_raw_df.to_csv(chain_output_dir / "geometry_descriptors.csv", index=False)
            chain_summary_df.to_csv(chain_output_dir / "residue_summary.csv", index=False)
            chain_model_summary_df.to_csv(chain_output_dir / "model_summary.csv", index=False)
            plot_chain_geometry(
                chain_raw_df=chain_raw_df,
                chain_summary_df=chain_summary_df,
                output_path=chain_output_dir / "curvature_torsion.png",
                show_model_traces=show_model_traces,
                max_models_in_plot=max_models_in_plot,
            )
            if distance_long_df is not None and distance_summary_df is not None:
                chain_distance_long_df = distance_long_df[distance_long_df["chain"] == chain].copy()
                chain_distance_summary_df = distance_summary_df[
                    distance_summary_df["chain"] == chain
                ].copy()
                if not chain_distance_long_df.empty:
                    chain_distance_long_df.to_csv(
                        chain_output_dir / "distance_to_reference_long.csv", index=False
                    )
                    chain_distance_summary_df.to_csv(
                        chain_output_dir / "distance_to_reference_summary.csv", index=False
                    )
                    write_distance_matrix_csv(
                        chain_distance_long_df,
                        chain_output_dir / "distance_to_reference_matrix.csv",
                    )
                    chain_reference_label = chain_distance_long_df["reference_label"].iloc[0]
                    plot_distance_heatmap(
                        distance_long_df=chain_distance_long_df,
                        output_path=chain_output_dir / "distance_to_reference_heatmap.png",
                        title=f"Chain {chain}: distance to {chain_reference_label}",
                    )
            if residue_clusters_df is not None and residue_cluster_summary_df is not None:
                chain_cluster_df = residue_clusters_df[
                    residue_clusters_df["chain"] == chain
                ].copy()
                chain_cluster_summary_df = residue_cluster_summary_df[
                    residue_cluster_summary_df["chain"] == chain
                ].copy()
                if not chain_cluster_df.empty:
                    chain_cluster_plot_dir = chain_output_dir / "cluster_plots"
                    chain_cluster_plot_dir.mkdir(parents=True, exist_ok=True)
                    chain_cluster_df.to_csv(
                        chain_output_dir / "residue_cluster_assignments.csv", index=False
                    )
                    chain_cluster_summary_df.to_csv(
                        chain_output_dir / "residue_cluster_summary.csv", index=False
                    )
                    for residue_label, residue_cluster_df in chain_cluster_df.groupby(
                        "residue_label", dropna=False
                    ):
                        plot_residue_cluster_solution(
                            residue_cluster_df=residue_cluster_df,
                            output_path=chain_cluster_plot_dir / f"{residue_label}_clusters.png",
                        )
            if (
                residue_range_clusters_df is not None
                and residue_range_cluster_summary_df is not None
            ):
                chain_range_cluster_df = residue_range_clusters_df[
                    residue_range_clusters_df["chain"] == chain
                ].copy()
                chain_range_cluster_summary_df = residue_range_cluster_summary_df[
                    residue_range_cluster_summary_df["chain"] == chain
                ].copy()
                if not chain_range_cluster_df.empty:
                    chain_range_cluster_plot_dir = chain_output_dir / "range_cluster_plots"
                    chain_range_cluster_plot_dir.mkdir(parents=True, exist_ok=True)
                    chain_range_cluster_df.to_csv(
                        chain_output_dir / "residue_range_cluster_assignments.csv", index=False
                    )
                    chain_range_cluster_summary_df.to_csv(
                        chain_output_dir / "residue_range_cluster_summary.csv", index=False
                    )
                    for range_label, range_cluster_df in chain_range_cluster_df.groupby(
                        "range_label", dropna=False
                    ):
                        safe_range_label = str(range_label).replace("/", "_")
                        plot_residue_range_cluster_solution(
                            range_cluster_df=range_cluster_df,
                            output_path=chain_range_cluster_plot_dir
                            / f"{safe_range_label}_clusters.png",
                        )

    return {
        "raw_csv": raw_csv,
        "residue_summary_csv": residue_summary_csv,
        "model_summary_csv": model_summary_csv if output_verbose else None,
        "overall_model_summary_csv": overall_model_summary_csv,
        "overview_plot": overview_plot,
        "chains_dir": chains_dir,
        "distance_long_csv": (
            distance_long_csv if distance_long_df is not None and output_verbose else None
        ),
        "distance_summary_csv": (
            distance_summary_csv if distance_summary_df is not None else None
        ),
        "distance_heatmap": distance_heatmap if distance_long_df is not None else None,
        "distance_matrix_dir": (
            distance_matrix_dir if distance_long_df is not None and output_verbose else None
        ),
        "cluster_assignments_csv": (
            cluster_assignments_csv
            if residue_clusters_df is not None and output_verbose
            else None
        ),
        "cluster_summary_csv": (
            cluster_summary_csv if residue_cluster_summary_df is not None else None
        ),
        "cluster_plots_dir": cluster_plots_dir if residue_clusters_df is not None else None,
        "range_cluster_assignments_csv": (
            range_cluster_assignments_csv
            if residue_range_clusters_df is not None and output_verbose
            else None
        ),
        "range_cluster_summary_csv": (
            range_cluster_summary_csv
            if residue_range_cluster_summary_df is not None
            else None
        ),
        "range_cluster_plots_dir": (
            range_cluster_plots_dir if residue_range_clusters_df is not None else None
        ),
    }


def print_run_summary(pdb_file: Path, raw_df, outputs: dict[str, Path]) -> None:
    model_count = int(raw_df["model"].nunique())
    chain_count = int(raw_df["chain"].nunique())
    residue_count = int(raw_df.groupby(["chain", "order"]).ngroups)

    print(f"Processed: {pdb_file}")
    print(f"Models: {model_count}")
    print(f"Chains: {chain_count}")
    print(f"Residues: {residue_count}")
    print(f"Raw descriptors: {outputs['raw_csv']}")
    print(f"Residue summary: {outputs['residue_summary_csv']}")
    print(f"Overall model summary: {outputs['overall_model_summary_csv']}")
    print(f"Overview plot: {outputs['overview_plot']}")
    if outputs["model_summary_csv"] is not None:
        print(f"Model summary by chain: {outputs['model_summary_csv']}")
    if outputs["distance_summary_csv"] is not None:
        print(f"Distance summary: {outputs['distance_summary_csv']}")
        print(f"Distance heatmap: {outputs['distance_heatmap']}")
    if outputs["distance_long_csv"] is not None:
        print(f"Distance details: {outputs['distance_long_csv']}")
        print(f"Distance matrices by chain: {outputs['distance_matrix_dir']}")
    if outputs["cluster_summary_csv"] is not None:
        print(f"Residue cluster summary: {outputs['cluster_summary_csv']}")
        print(f"Residue cluster plots: {outputs['cluster_plots_dir']}")
    if outputs["cluster_assignments_csv"] is not None:
        print(f"Residue cluster assignments: {outputs['cluster_assignments_csv']}")
    if outputs["range_cluster_summary_csv"] is not None:
        print(f"Residue-range cluster summary: {outputs['range_cluster_summary_csv']}")
        print(f"Residue-range cluster plots: {outputs['range_cluster_plots_dir']}")
    if outputs["range_cluster_assignments_csv"] is not None:
        print(f"Residue-range cluster assignments: {outputs['range_cluster_assignments_csv']}")
    if outputs["chains_dir"] is not None:
        print(f"Per-chain outputs: {outputs['chains_dir']}")


def main() -> None:
    args = parse_args()
    ensure_dependencies()
    configure_plot_style()

    if args.reference_pdb_model and not args.reference_pdb:
        raise ValueError("--reference-pdb-model requires --reference-pdb.")

    pdb_file = args.pdb_file.resolve()
    output_dir = args.output_dir.resolve()

    if not pdb_file.is_file():
        raise FileNotFoundError(f"Input PDB file not found: {pdb_file}")

    output_dir.mkdir(parents=True, exist_ok=True)

    raw_df = compute_geometry(pdb_file, n_jobs=args.n_jobs)
    raw_df = filter_chains(raw_df, args.chains)
    raw_df = normalise_geometry_table(raw_df)
    residue_summary_df = summarise_geometry(raw_df)
    model_summary_df, overall_model_summary_df = build_model_summary(
        raw_df, residue_summary_df
    )
    distance_long_df = None
    distance_summary_df = None
    residue_clusters_df = None
    residue_cluster_summary_df = None
    residue_range_clusters_df = None
    residue_range_cluster_summary_df = None

    if args.reference_model or args.reference_pdb:
        if args.reference_model:
            reference_rows, reference_model_label = select_reference_rows(
                raw_df, args.reference_model
            )
            reference_label = f"input model {reference_model_label}"
        else:
            reference_pdb = args.reference_pdb.resolve()
            if not reference_pdb.is_file():
                raise FileNotFoundError(f"Reference PDB file not found: {reference_pdb}")
            reference_df = compute_geometry(reference_pdb, n_jobs=args.n_jobs)
            reference_df = filter_chains(reference_df, args.chains)
            reference_df = normalise_geometry_table(reference_df)
            reference_rows, reference_model_label = select_reference_rows(
                reference_df, args.reference_pdb_model
            )
            reference_label = f"{reference_pdb.name} model {reference_model_label}"

        distance_long_df, distance_summary_df = compute_reference_distance_tables(
            raw_df=raw_df,
            reference_df=reference_rows,
            reference_label=reference_label,
        )

    if args.cluster_residues:
        residue_clusters_df, residue_cluster_summary_df = cluster_residues(
            raw_df=raw_df,
            min_cluster_size=args.cluster_min_size,
            min_samples=args.cluster_min_samples,
        )
    if args.cluster_residue_ranges:
        residue_range_clusters_df, residue_range_cluster_summary_df = cluster_residue_ranges(
            raw_df=raw_df,
            range_texts=args.cluster_residue_ranges,
            min_cluster_size=args.cluster_min_size,
            min_samples=args.cluster_min_samples,
        )

    outputs = write_outputs(
        raw_df=raw_df,
        residue_summary_df=residue_summary_df,
        model_summary_df=model_summary_df,
        overall_model_summary_df=overall_model_summary_df,
        distance_long_df=distance_long_df,
        distance_summary_df=distance_summary_df,
        residue_clusters_df=residue_clusters_df,
        residue_cluster_summary_df=residue_cluster_summary_df,
        residue_range_clusters_df=residue_range_clusters_df,
        residue_range_cluster_summary_df=residue_range_cluster_summary_df,
        output_dir=output_dir,
        show_model_traces=not args.hide_model_traces,
        max_models_in_plot=args.max_models_in_plot,
        output_verbose=args.output_verbose,
    )
    print_run_summary(pdb_file, raw_df, outputs)


if __name__ == "__main__":
    main()
