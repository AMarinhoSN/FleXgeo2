from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

# Columns that identify a row of per-model, per-residue tables, in output order.
DESCRIPTOR_KEYS = ("chain", "model", "order", "name", "residue_label")

# residue_summary_df columns that explain how dmax was computed; kept for the library,
# left out of geometry/residues.csv.
DMAX_DETAIL_COLUMNS = (
    "curvature_dmax_min",
    "curvature_dmax_max",
    "torsion_dmax_min",
    "torsion_dmax_max",
    "curvature_dmax_bin_width",
    "torsion_dmax_bin_width",
)


@dataclass(slots=True)
class StructureInfo:
    """Lightweight description of a parsed structure, used for early validation."""

    model_ids: list[int]
    residues_by_chain: dict[str, set[int]]


class GeometryService:
    """Load and summarize Melodia backbone geometry."""

    required_columns = {"model", "chain", "order", "name", "curvature", "torsion"}

    @staticmethod
    def parse_structure(pdb_file: str | Path):
        """Parse a PDB file with Biopython, without computing any geometry."""
        from Bio.PDB import PDBParser

        path = Path(pdb_file)
        return PDBParser(QUIET=True).get_structure(path.stem, str(path))

    @staticmethod
    def model_serial_map(structure) -> dict[int, int]:
        """Map Biopython's 0-based model index to the PDB MODEL serial number.

        Files without MODEL records have a serial of 0; those fall back to index + 1 so
        model identifiers are always 1-based, matching PDB conventions.
        """
        return {
            int(model.id): int(model.serial_num) if model.serial_num else int(model.id) + 1
            for model in structure
        }

    @classmethod
    def describe_structure(cls, structure) -> StructureInfo:
        """Summarise model identifiers and first-model residue numbers per chain."""
        serial_map = cls.model_serial_map(structure)
        residues_by_chain: dict[str, set[int]] = {}
        first_model = next(iter(structure), None)
        if first_model is not None:
            for chain in first_model:
                residues = {int(res.id[1]) for res in chain.get_residues() if res.id[0] == " "}
                if residues:
                    residues_by_chain[str(chain.id)] = residues
        return StructureInfo(
            model_ids=[serial_map[key] for key in sorted(serial_map)],
            residues_by_chain=residues_by_chain,
        )

    def compute_geometry(self, structure, n_jobs: int = 1):
        """Run Melodia on a parsed structure and report PDB MODEL serials as model IDs."""
        import melodia_py as mel

        df = mel.geometry_from_structure(structure, n_jobs=n_jobs)
        missing = self.required_columns.difference(df.columns)
        if missing:
            missing_csv = ", ".join(sorted(missing))
            raise ValueError(
                f"Melodia output is missing required columns: {missing_csv}. "
                f"Available columns: {', '.join(df.columns)}"
            )
        df = df.copy()
        df["model"] = df["model"].astype(int).map(self.model_serial_map(structure))
        return df

    def load_structure(self, pdb_file: str | Path, n_jobs: int = 1):
        return self.compute_geometry(self.parse_structure(pdb_file), n_jobs=n_jobs)

    def filter_chains(self, df, chains: list[str] | None):
        if not chains:
            return df
        filtered = df[df["chain"].isin(chains)].copy()
        if filtered.empty:
            raise ValueError(f"No residues found for requested chain(s): {', '.join(chains)}")
        return filtered

    def normalize(self, df):
        # Melodia's row index ("id") and file-name code ("code") carry no information.
        normalised = df.drop(columns=["id", "code"], errors="ignore")
        normalised["order"] = normalised["order"].astype(int)
        normalised["residue_label"] = [
            f"{name}{int(order)}"
            for order, name in zip(normalised["order"], normalised["name"], strict=False)
        ]
        keys = [column for column in DESCRIPTOR_KEYS if column in normalised.columns]
        normalised = normalised[keys + [c for c in normalised.columns if c not in keys]]
        return normalised.sort_values(["chain", "model", "order"]).reset_index(drop=True)

    @staticmethod
    def _validate_dmax_outlier_fraction(dmax_outlier_fraction: float) -> float:
        fraction = float(dmax_outlier_fraction)
        if not math.isfinite(fraction) or fraction < 0.0 or fraction >= 1.0:
            raise ValueError("dmax_outlier_fraction must be a finite fraction in [0, 1).")
        return fraction

    @staticmethod
    def _optimise_histogram_bin_width(values) -> tuple[float, int]:
        import numpy as np

        finite_values = np.asarray(values, dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]
        if finite_values.size <= 1:
            return 0.0, 1

        value_min = float(finite_values.min())
        value_max = float(finite_values.max())
        value_range = value_max - value_min
        if value_range <= 0.0:
            return 0.0, 1

        max_bins = min(200, max(2, finite_values.size))
        candidate_bins = np.arange(2, max_bins + 1)
        bin_widths = value_range / candidate_bins
        costs = []
        for n_bins, bin_width in zip(candidate_bins, bin_widths, strict=True):
            counts, _ = np.histogram(
                finite_values,
                bins=int(n_bins),
                range=(value_min, value_max),
            )
            mean_count = float(counts.mean())
            biased_variance = float(((counts - mean_count) ** 2).sum() / n_bins)
            costs.append((2.0 * mean_count - biased_variance) / (bin_width**2))

        best_index = int(np.argmin(np.asarray(costs, dtype=float)))
        return float(bin_widths[best_index]), int(candidate_bins[best_index])

    @classmethod
    def _compute_dmax_axis_extremes(
        cls,
        values,
        dmax_outlier_fraction: float,
    ) -> tuple[float, float, float]:
        import numpy as np

        finite_values = np.asarray(values, dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]
        if finite_values.size == 0:
            return math.nan, math.nan, 0.0
        if finite_values.size == 1:
            value = float(finite_values[0])
            return value, value, 0.0

        value_min = float(finite_values.min())
        value_max = float(finite_values.max())
        if value_max == value_min:
            return value_min, value_max, 0.0

        bin_width, n_bins = cls._optimise_histogram_bin_width(finite_values)
        counts, edges = np.histogram(
            finite_values,
            bins=n_bins,
            range=(value_min, value_max),
        )

        first_bin = 0
        last_bin = len(counts) - 1
        total_count = float(finite_values.size)
        while (
            first_bin < last_bin
            and (float(counts[first_bin]) / total_count) < dmax_outlier_fraction
        ):
            first_bin += 1
        while (
            last_bin > first_bin and (float(counts[last_bin]) / total_count) < dmax_outlier_fraction
        ):
            last_bin -= 1

        retained_min_edge = edges[first_bin]
        retained_max_edge = edges[last_bin + 1]
        retained_values = finite_values[
            (finite_values >= retained_min_edge) & (finite_values <= retained_max_edge)
        ]
        if retained_values.size == 0:
            retained_values = finite_values

        return float(retained_values.min()), float(retained_values.max()), bin_width

    def summarize(self, df, dmax_outlier_fraction: float = 0.01):
        import pandas as pd

        dmax_outlier_fraction = self._validate_dmax_outlier_fraction(dmax_outlier_fraction)
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
        grouped = grouped.fillna({"curvature_std": 0.0, "torsion_std": 0.0})

        dmax_rows = []
        for (chain, order, name, residue_label), residue_df in df.groupby(
            ["chain", "order", "name", "residue_label"], dropna=False
        ):
            curvature_min, curvature_max, curvature_bin_width = self._compute_dmax_axis_extremes(
                residue_df["curvature"],
                dmax_outlier_fraction=dmax_outlier_fraction,
            )
            torsion_min, torsion_max, torsion_bin_width = self._compute_dmax_axis_extremes(
                residue_df["torsion"],
                dmax_outlier_fraction=dmax_outlier_fraction,
            )
            dmax_rows.append(
                {
                    "chain": chain,
                    "order": order,
                    "name": name,
                    "residue_label": residue_label,
                    "curvature_dmax_min": curvature_min,
                    "curvature_dmax_max": curvature_max,
                    "torsion_dmax_min": torsion_min,
                    "torsion_dmax_max": torsion_max,
                    "curvature_dmax_bin_width": curvature_bin_width,
                    "torsion_dmax_bin_width": torsion_bin_width,
                    "dmax": math.sqrt(
                        (curvature_min - curvature_max) ** 2 + (torsion_min - torsion_max) ** 2
                    ),
                }
            )

        dmax_df = pd.DataFrame(dmax_rows)
        return grouped.merge(
            dmax_df,
            on=["chain", "order", "name", "residue_label"],
            how="left",
        )

    def build_model_summary(self, raw_df, residue_summary_df):
        merged = raw_df.merge(
            residue_summary_df[["chain", "order", "curvature_mean", "torsion_mean"]],
            on=["chain", "order"],
            how="left",
        )
        merged["curvature_abs_deviation"] = (merged["curvature"] - merged["curvature_mean"]).abs()
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
