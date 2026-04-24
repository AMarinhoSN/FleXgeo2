from __future__ import annotations


class DistanceService:
    """Reference-state selection and curvature/torsion distance analysis."""

    def select_reference_rows(self, df, model_id: str | None):
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

    def compute(self, raw_df, reference_df, reference_label: str):
        comparison = raw_df.merge(
            reference_df[
                ["chain", "order", "residue_label", "curvature", "torsion"]
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
            (comparison["curvature"] - comparison["reference_curvature"]) ** 2
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
