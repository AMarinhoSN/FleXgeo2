from __future__ import annotations


class ClusteringService:
    """HDBSCAN clustering helpers for residues and residue ranges."""

    @staticmethod
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

    @staticmethod
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

    def cluster_residues(self, raw_df, min_cluster_size: int, min_samples: int | None):
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

    def cluster_residue_ranges(
        self,
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
        parsed_ranges = [self.parse_residue_range(range_text) for range_text in range_texts]

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

                feature_matrix = np.array(
                    feature_table.to_numpy(), dtype=float, order="C", copy=True
                )
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

                projection = self.compute_pca_projection(feature_matrix)
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
