from __future__ import annotations

import math
from pathlib import Path


def sanitize_chain_id(chain) -> str:
    if chain is None or chain == "":
        return "unassigned"
    return str(chain).replace("/", "_")


class PlotStyle:
    """Global matplotlib styling for FleXgeo2 plots."""

    @staticmethod
    def apply() -> None:
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


class BasePlotter:
    @staticmethod
    def apply_residue_ticks(axis, x_values, residue_labels) -> None:
        if len(residue_labels) == 0:
            return
        tick_step = max(1, math.ceil(len(residue_labels) / 18))
        axis.set_xticks(x_values[::tick_step])
        axis.set_xticklabels(residue_labels[::tick_step], rotation=45, ha="right")


class ChainGeometryPlotter(BasePlotter):
    def plot(
        self,
        chain_raw_df,
        chain_summary_df,
        output_path: str | Path,
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
            self.apply_residue_ticks(axis, x_values, residue_labels)

        fig.suptitle("Backbone Differential Geometry", fontsize=16, fontweight="bold")
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


class OverviewPlotter(BasePlotter):
    def plot(self, summary_df, output_path: str | Path) -> None:
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

        for row_axes, chain in zip(axes, chains, strict=False):
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
                self.apply_residue_ticks(axis, x_values, residue_labels)

        fig.suptitle("Ensemble Overview", fontsize=16, fontweight="bold")
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


class DistanceHeatmapPlotter:
    def plot(self, distance_long_df, output_path: str | Path, title: str) -> None:
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

        for axis, chain in zip(axes, chains, strict=False):
            chain_df = distance_long_df[distance_long_df["chain"] == chain].copy()
            matrix = chain_df.pivot(
                index="model", columns="residue_label", values="distance_to_reference"
            ).sort_index()
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


class ResidueClusterPlotter:
    def plot(self, residue_cluster_df, output_path: str | Path) -> None:
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


class ResidueRangeClusterPlotter:
    def plot(self, range_cluster_df, output_path: str | Path) -> None:
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
