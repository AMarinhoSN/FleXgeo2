from __future__ import annotations

from pathlib import Path

from flexgeo2.distances import DistanceService


def sanitize_chain_id(chain) -> str:
    if chain is None or chain == "":
        return "unassigned"
    return str(chain).replace("/", "_")


NOISE_COLOR = "#9e9e9e"
# Clusters past the end of the palette share this colour and one "Other clusters" entry
# in legends: dark grey, distinct from the palette, the noise grey and the black
# reference star.
OTHER_CLUSTERS_COLOR = "#404040"


def cluster_palette() -> list[tuple[float, float, float]]:
    """Cluster colours: tab10, then tab20's light shades. Greys are left out for noise."""
    import matplotlib.pyplot as plt

    tab10 = [color for index, color in enumerate(plt.get_cmap("tab10").colors) if index != 7]
    light = [
        color
        for index, color in enumerate(plt.get_cmap("tab20").colors)
        if index % 2 == 1 and index != 15
    ]
    return tab10 + light


def cluster_color(label):
    """Colour for a cluster label, the same in every plot; -1 (noise) is grey."""
    if int(label) < 0:
        return NOISE_COLOR
    palette = cluster_palette()
    return palette[int(label)] if int(label) < len(palette) else OTHER_CLUSTERS_COLOR


def cluster_legend_label(label) -> str:
    if int(label) < 0:
        return "Noise"
    n_colors = len(cluster_palette())
    return f"Cluster {int(label)}" if int(label) < n_colors else f"Other clusters ({n_colors}+)"


def cluster_legend_groups(labels) -> dict[str, list[int]]:
    """``{legend entry: cluster labels}`` in label order; one entry for all "other" clusters."""
    groups: dict[str, list[int]] = {}
    for label in sorted({int(label) for label in labels}):
        groups.setdefault(cluster_legend_label(label), []).append(label)
    return groups


def label_rows(axis, values) -> None:
    """Label an image axis that has one row (or column) per value, e.g. model IDs.

    Labels go on round values (multiples of 1, 2, 5, 10, 20, 50, ...), as many as fit
    the axis once the figure is laid out (at most 10), so they do not overlap however
    many rows there are. Values need not be contiguous.
    """
    from matplotlib.ticker import FuncFormatter, Locator

    values = [int(value) for value in values]

    class RoundValueLocator(Locator):
        def __call__(self):
            if not values:
                return []
            # As many labels as fit (matplotlib's estimate), at most 10 like its other axes.
            room = max(1, min(self.axis.get_tick_space(), 10))
            step = _round_step(max(values) - min(values), room)
            ticks = [row for row, value in enumerate(values) if value % step == 0]
            by_row = list(range(0, len(values), _round_step(len(values) - 1, room)))
            # Irregular values may have few round ones; then space the labels by row.
            return ticks if len(ticks) >= len(by_row) / 2 else by_row

    def format_row(position, _) -> str:
        row = round(position)
        return str(values[row]) if 0 <= row < len(values) else ""

    axis.set_major_locator(RoundValueLocator())
    axis.set_major_formatter(FuncFormatter(format_row))


# Residue bars keep a gap between them up to this many residues per chain. On longer
# chains the gaps shrink below ~2 px at 300 dpi and alias into stripes, so bars touch.
MAX_RESIDUES_WITH_BAR_GAPS = 400


def residue_bar_width(n_residues: int) -> float:
    """Width of per-residue bars (1.0 = touching) for a chain of ``n_residues``."""
    return 0.8 if n_residues <= MAX_RESIDUES_WITH_BAR_GAPS else 1.0


def model_axis_height(n_models: int) -> float:
    """Height in inches of a models x residues image, shared by the heatmap and cluster map."""
    return min(6.0, max(2.0, 0.12 * n_models))


def _round_step(span: int, max_ticks: int) -> int:
    """Smallest step of 1, 2, 5, 10, 20, 50, ... giving at most ``max_ticks`` ticks."""
    scale = 1
    while True:
        for base in (1, 2, 5):
            step = base * scale
            if span // step + 1 <= max_ticks:
                return step
        scale *= 10


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
                # Vector figures stay editable: TrueType fonts in PDF (journals often
                # reject matplotlib's default Type 3), real text rather than paths in SVG.
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
            }
        )


class BasePlotter:
    @staticmethod
    def apply_residue_ticks(axis) -> None:
        """Ticks on round residue numbers for an x axis in residue-number coordinates."""
        from matplotlib.ticker import MaxNLocator

        axis.xaxis.set_major_locator(MaxNLocator(nbins="auto", integer=True, steps=[1, 2, 5, 10]))

    @staticmethod
    def plot_model_traces(curvature_axis, torsion_axis, chain_raw_df, max_models_in_plot) -> None:
        model_ids = list(chain_raw_df["model"].drop_duplicates())[:max_models_in_plot]
        for model_id in model_ids:
            model_df = chain_raw_df[chain_raw_df["model"] == model_id]
            curvature_axis.plot(
                model_df["order"],
                model_df["curvature"],
                color="#4c78a8",
                alpha=0.22,
                linewidth=1,
            )
            torsion_axis.plot(
                model_df["order"],
                model_df["torsion"],
                color="#e45756",
                alpha=0.22,
                linewidth=1,
            )


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
        chain_id = chain_summary_df["chain"].iloc[0] if not chain_summary_df.empty else None
        title_suffix = f"Chain {chain_id}" if chain_id not in (None, "") else "Chain"

        if show_model_traces:
            self.plot_model_traces(axes[0], axes[1], chain_raw_df, max_models_in_plot)

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
            self.apply_residue_ticks(axis)

        fig.suptitle("Backbone Differential Geometry", fontsize=16, fontweight="bold")
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


class OverviewPlotter(BasePlotter):
    """Per-residue results along the sequence, one panel per quantity, per chain.

    Curvature, torsion and dmax are always shown; clusters per residue and the distance
    to the reference are added when those analyses ran.
    """

    PANEL_HEIGHTS = {
        "curvature": 2.4,
        "torsion": 2.4,
        "dmax": 1.4,
        "clusters": 1.2,
        "distance": 1.6,
    }

    def plot(
        self,
        summary_df,
        output_path: str | Path,
        raw_df=None,
        show_model_traces: bool = True,
        max_models_in_plot: int = 12,
        cluster_summary_df=None,
        distance_summary_df=None,
    ) -> None:
        import matplotlib.pyplot as plt

        fig = self.render(
            summary_df,
            raw_df=raw_df,
            show_model_traces=show_model_traces,
            max_models_in_plot=max_models_in_plot,
            cluster_summary_df=cluster_summary_df,
            distance_summary_df=distance_summary_df,
        )
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def render(
        self,
        summary_df,
        raw_df=None,
        show_model_traces: bool = True,
        max_models_in_plot: int = 12,
        cluster_summary_df=None,
        distance_summary_df=None,
    ):
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator

        panels = ["curvature", "torsion", "dmax"]
        if cluster_summary_df is not None:
            panels.append("clusters")
        if distance_summary_df is not None:
            panels.append("distance")
        heights = [self.PANEL_HEIGHTS[panel] for panel in panels]

        chains = list(summary_df["chain"].drop_duplicates())
        fig = plt.figure(
            figsize=(14, 0.8 + len(chains) * (sum(heights) + 0.6)), constrained_layout=True
        )
        grid = fig.add_gridspec(
            nrows=len(chains) * len(panels), ncols=1, height_ratios=heights * len(chains)
        )

        for chain_index, chain in enumerate(chains):
            chain_df = summary_df[summary_df["chain"] == chain].sort_values("order")
            x_values = chain_df["order"].to_numpy()
            axes = {}
            for panel_index, panel in enumerate(panels):
                shared = axes.get("curvature")
                axes[panel] = fig.add_subplot(
                    grid[chain_index * len(panels) + panel_index], sharex=shared
                )

            if show_model_traces and raw_df is not None:
                self.plot_model_traces(
                    axes["curvature"],
                    axes["torsion"],
                    raw_df[raw_df["chain"] == chain],
                    max_models_in_plot,
                )
            for name, line_color, band_color in (
                ("curvature", "#1f4e79", "#4c78a8"),
                ("torsion", "#b22222", "#e45756"),
            ):
                self._mean_and_band(
                    axes[name],
                    x_values,
                    chain_df[f"{name}_mean"],
                    chain_df[f"{name}_std"],
                    line_color,
                    band_color,
                )
            axes["curvature"].set_ylabel("Curvature (1/Å)")
            axes["torsion"].set_ylabel("Torsion (1/Å)")

            bar_width = residue_bar_width(len(x_values))
            axes["dmax"].bar(x_values, chain_df["dmax"], width=bar_width, color="#6a3d9a")
            axes["dmax"].set_ylabel("dmax")

            if cluster_summary_df is not None:
                clusters = self._align(chain_df, cluster_summary_df, chain, ["n_clusters"])
                axes["clusters"].bar(
                    x_values, clusters["n_clusters"].fillna(0), width=bar_width, color="#4c4c4c"
                )
                axes["clusters"].set_ylabel("Clusters")
                axes["clusters"].yaxis.set_major_locator(MaxNLocator(integer=True))

            if distance_summary_df is not None:
                distances = self._align(
                    chain_df, distance_summary_df, chain, ["distance_mean", "distance_std"]
                )
                self._mean_and_band(
                    axes["distance"],
                    x_values,
                    distances["distance_mean"],
                    distances["distance_std"],
                    "#b35806",
                    "#fdb863",
                )
                axes["distance"].set_ylabel("Distance to\nreference")

            title = f"Chain {chain}" if chain not in (None, "") else "Chain"
            axes["curvature"].set_title(title)
            for panel in panels:
                axis = axes[panel]
                axis.grid(alpha=0.3)
                if panel in ("dmax", "clusters"):
                    axis.grid(axis="x", visible=False)
                axis.tick_params(labelbottom=False)
            if len(x_values):
                # Bars would otherwise pad the shared axis with wide empty margins.
                axes["curvature"].set_xlim(x_values.min() - 0.6, x_values.max() + 0.6)
            bottom = axes[panels[-1]]
            bottom.tick_params(labelbottom=True)
            bottom.set_xlabel("Residue")
            self.apply_residue_ticks(bottom)

        fig.suptitle("Ensemble Overview", fontsize=16, fontweight="bold")
        return fig

    @staticmethod
    def _mean_and_band(axis, x_values, mean, std, line_color, band_color) -> None:
        axis.plot(x_values, mean, color=line_color, linewidth=2)
        axis.fill_between(x_values, mean - std, mean + std, color=band_color, alpha=0.18)

    @staticmethod
    def _align(chain_df, other_df, chain, columns):
        """Values of ``columns`` for each residue of ``chain_df``; NaN where absent."""
        chain_rows = other_df[other_df["chain"] == chain][["order", *columns]]
        return chain_df[["order"]].merge(chain_rows, on="order", how="left")


class DistanceHeatmapPlotter:
    # The colour scale ends at this percentile of each chain's distances. A few very
    # large distances (e.g. torsion spikes) would otherwise stretch it until the rest
    # of the map is black; the larger values share the top colour.
    COLOR_PERCENTILE = 99

    def plot(self, distance_long_df, output_path: str | Path, title: str) -> None:
        import matplotlib.pyplot as plt

        fig = self.render(distance_long_df, title)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def render(self, distance_long_df, title: str):
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        # Build and validate every chain's matrix before creating the figure, so an
        # invalid input does not leave an open figure behind.
        matrices = {}
        for chain in distance_long_df["chain"].drop_duplicates():
            chain_df = distance_long_df[distance_long_df["chain"] == chain]
            matrix = DistanceService.to_matrix(chain_df).apply(pd.to_numeric, errors="coerce")
            if matrix.empty or matrix.isna().all().all():
                raise ValueError(
                    f"Distance heatmap for chain '{chain}' is empty after alignment. "
                    "This usually means the chosen reference does not overlap with the "
                    "ensemble on chain/residue numbering."
                )
            matrices[chain] = matrix

        # Laid out like the cluster map: residues along x, models down y (first at top).
        heights = [model_axis_height(len(matrix.index)) for matrix in matrices.values()]
        fig, axes = plt.subplots(
            nrows=len(matrices),
            ncols=1,
            figsize=(14, 0.4 + sum(height + 1.0 for height in heights)),
            height_ratios=heights,
            constrained_layout=True,
            squeeze=False,
        )

        for axis, (chain, matrix) in zip(axes[:, 0], matrices.items(), strict=True):
            values = matrix.to_numpy(dtype=float)
            # "higher" picks an observed distance, so small matrices are not capped.
            vmax = float(np.nanpercentile(values, self.COLOR_PERCENTILE, method="higher"))
            image = axis.imshow(
                values,
                aspect="auto",
                cmap="magma",
                interpolation="nearest",
                vmin=0,
                # All distances zero: keep a 0-1 scale (matplotlib would pick -0.1 to 0.1).
                vmax=vmax if vmax > 0 else 1.0,
            )
            axis.set_title(
                f"Chain {chain}: Distance to reference"
                if chain not in (None, "")
                else "Chain: Distance to reference"
            )
            axis.grid(False)
            axis.set_xlabel("Residue")
            axis.set_ylabel("Model")

            # to_matrix puts residues in sequence order, one column per residue number.
            chain_orders = distance_long_df.loc[distance_long_df["chain"] == chain, "order"]
            label_rows(axis.xaxis, sorted(chain_orders.unique()))
            label_rows(axis.yaxis, matrix.index)

            fig.colorbar(
                image,
                ax=axis,
                fraction=0.024,
                pad=0.02,
                label="Euclidean distance",
                # An arrow at the top of the colour bar marks a capped scale.
                extend="max" if np.nanmax(values) > vmax else "neither",
            )

        fig.suptitle(title, fontsize=16, fontweight="bold")
        return fig


class ClusterMapPlotter(BasePlotter):
    """Cluster label of every model at every residue, with clusters per residue above."""

    MISSING_COLOR = "white"

    def plot(self, summary_df, output_path: str | Path, assignments_df) -> None:
        import matplotlib.pyplot as plt

        fig = self.render(summary_df, assignments_df)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def render(self, summary_df, assignments_df):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.colors import to_rgba
        from matplotlib.patches import Patch
        from matplotlib.ticker import MaxNLocator

        chains = list(summary_df["chain"].drop_duplicates())
        models = list(assignments_df["model"].drop_duplicates())
        map_height = model_axis_height(len(models))
        fig = plt.figure(figsize=(14, len(chains) * (map_height + 1.6)), constrained_layout=True)
        grid = fig.add_gridspec(
            nrows=2 * len(chains), ncols=1, height_ratios=[1.2, map_height] * len(chains)
        )

        labels_shown: set[int] = set()
        any_missing = False
        for index, chain in enumerate(chains):
            chain_summary = summary_df[summary_df["chain"] == chain].sort_values("order")
            chain_assignments = assignments_df[assignments_df["chain"] == chain]
            # Models x residues, in sequence order; NaN where a model lacks the residue.
            labels = (
                chain_assignments.pivot(index="model", columns="order", values="cluster")
                .reindex(index=models, columns=chain_summary["order"])
                .to_numpy(dtype=float)
            )
            image = np.empty(labels.shape + (4,))
            image[:] = to_rgba(self.MISSING_COLOR)
            for label in np.unique(labels[~np.isnan(labels)]):
                image[labels == label] = to_rgba(cluster_color(label))
                labels_shown.add(int(label))
            any_missing = any_missing or bool(np.isnan(labels).any())

            positions = np.arange(len(chain_summary))
            strip_axis = fig.add_subplot(grid[2 * index])
            map_axis = fig.add_subplot(grid[2 * index + 1], sharex=strip_axis)

            strip_axis.bar(
                positions,
                chain_summary["n_clusters"],
                width=residue_bar_width(len(positions)),
                color="#4c4c4c",
            )
            strip_axis.set_ylabel("Clusters")
            strip_axis.yaxis.set_major_locator(MaxNLocator(integer=True))
            strip_axis.tick_params(labelbottom=False)
            strip_axis.grid(axis="x", visible=False)
            title_suffix = f"Chain {chain}" if chain not in (None, "") else "Chain"
            strip_axis.set_title(f"{title_suffix}: clusters per residue")

            map_axis.imshow(image, aspect="auto", interpolation="nearest")
            map_axis.grid(False)
            map_axis.set_xlabel("Residue")
            map_axis.set_ylabel("Model")
            label_rows(map_axis.xaxis, chain_summary["order"])
            label_rows(map_axis.yaxis, models)

        handles = [
            Patch(facecolor=cluster_color(members[0]), label=entry)
            for entry, members in cluster_legend_groups(labels_shown).items()
        ]
        if any_missing:
            handles.append(
                Patch(facecolor=self.MISSING_COLOR, edgecolor="#c7c7c7", label="Not present")
            )
        fig.legend(handles=handles, loc="outside right upper")
        fig.suptitle(
            "Per-residue clusters (labels are assigned independently at each residue)",
            fontsize=14,
            fontweight="bold",
        )
        return fig


class ResiduePlotter:
    """Curvature vs torsion of one residue, one point per model.

    Points are coloured by cluster when ``residue_df`` has a ``cluster`` column, and the
    reference state is marked when ``reference`` (curvature, torsion) is given.
    """

    def plot(self, residue_df, output_path: str | Path, dmax=None, reference=None) -> None:
        import matplotlib.pyplot as plt

        fig = self.render(residue_df, dmax=dmax, reference=reference)
        fig.savefig(output_path, dpi=250, bbox_inches="tight")
        plt.close(fig)

    def render(self, residue_df, dmax=None, reference=None):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
        if "cluster" in residue_df.columns:
            for entry, members in cluster_legend_groups(residue_df["cluster"]).items():
                points = residue_df[residue_df["cluster"].isin(members)]
                ax.scatter(
                    points["curvature"],
                    points["torsion"],
                    s=36,
                    alpha=0.8,
                    c=[cluster_color(members[0])],
                    label=entry,
                    edgecolors="none",
                )
        else:
            ax.scatter(
                residue_df["curvature"],
                residue_df["torsion"],
                s=36,
                alpha=0.8,
                c="#4c78a8",
                label="Models",
                edgecolors="none",
            )
        if reference is not None:
            ax.scatter(
                [reference[0]],
                [reference[1]],
                marker="*",
                s=260,
                c="black",
                label="Reference",
                zorder=3,
            )

        chain = residue_df["chain"].iloc[0]
        title = residue_df["residue_label"].iloc[0]
        if chain not in (None, ""):
            title = f"Chain {chain}: {title}"
        if dmax is not None:
            title = f"{title} (dmax {dmax:.3f})"
        ax.set_title(title)
        ax.set_xlabel("Curvature (1/Å)")
        ax.set_ylabel("Torsion (1/Å)")
        if len(ax.get_legend_handles_labels()[0]) > 1:
            ax.legend(loc="best")
        ax.grid(alpha=0.3)
        return fig


class ResidueRangeClusterPlotter:
    def plot(self, range_cluster_df, output_path: str | Path) -> None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
        chain = range_cluster_df["chain"].iloc[0]
        range_label = range_cluster_df["range_label"].iloc[0]
        title_prefix = f"Chain {chain}" if chain not in (None, "") else "Chain"

        for entry, members in cluster_legend_groups(range_cluster_df["cluster"]).items():
            cluster_points = range_cluster_df[range_cluster_df["cluster"].isin(members)]
            ax.scatter(
                cluster_points["pc1"],
                cluster_points["pc2"],
                s=48,
                alpha=0.88,
                c=[cluster_color(members[0])],
                label=entry,
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
