from __future__ import annotations

from pathlib import Path


class GeometryService:
    """Load and summarize Melodia backbone geometry."""

    required_columns = {"model", "chain", "order", "name", "curvature", "torsion"}

    @staticmethod
    def ensure_dependencies() -> None:
        try:
            import matplotlib  # noqa: F401
            import melodia_py  # noqa: F401
            import pandas  # noqa: F401
        except ModuleNotFoundError as exc:
            missing = exc.name or "required dependency"
            raise SystemExit(
                f"Missing dependency: {missing}. Install the project dependencies first, "
                "for example with: pip install -e ."
            ) from exc

    def load_structure(self, pdb_file: str | Path, n_jobs: int = 1):
        import melodia_py as mel

        df = mel.geometry_from_structure_file(str(pdb_file), n_jobs=n_jobs)
        missing = self.required_columns.difference(df.columns)
        if missing:
            missing_csv = ", ".join(sorted(missing))
            raise ValueError(
                f"Melodia output is missing required columns: {missing_csv}. "
                f"Available columns: {', '.join(df.columns)}"
            )
        return df.copy()

    def filter_chains(self, df, chains: list[str] | None):
        if not chains:
            return df
        filtered = df[df["chain"].isin(chains)].copy()
        if filtered.empty:
            raise ValueError(f"No residues found for requested chain(s): {', '.join(chains)}")
        return filtered

    def normalize(self, df):
        normalised = df.copy()
        normalised["order"] = normalised["order"].astype(int)
        normalised["residue_label"] = [
            f"{name}{int(order)}" for order, name in zip(normalised["order"], normalised["name"])
        ]
        return normalised.sort_values(["chain", "model", "order"]).reset_index(drop=True)

    def summarize(self, df):
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

    def build_model_summary(self, raw_df, residue_summary_df):
        merged = raw_df.merge(
            residue_summary_df[["chain", "order", "curvature_mean", "torsion_mean"]],
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
