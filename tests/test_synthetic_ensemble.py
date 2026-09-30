"""Synthetic switch ensembles: known local switches, analysed end to end.

The fixtures in ``tests/data`` (made by ``tests/data/generate_synthetic_ensembles.py``)
are ideal poly-alanine alpha helices in which one residue switches conformation between
blocks of models. Each switch swings the rest of the chain as a rigid body. FleXgeo2
should separate the states at the switch and, being superposition free, see no change
far from it. State definitions are read from the REMARK 999 header of each file.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

from flexgeo2 import AnalysisConfig, ClusteringConfig, FlexGeo2App
from flexgeo2.geometry import GeometryService

pytest.importorskip("melodia_py")

DATA_DIR = Path(__file__).parent / "data"
MIN_CLUSTER_SIZE = 10

# Expected clusters per residue, keyed by offset from the switch residue. Measured with
# Melodia on the noise-free states and checked for robustness over 50 noise seeds:
# - two states (alpha/beta): the switched residue (~0.12) and the next one (~0.43)
#   separate cleanly;
# - three states (alpha/beta/alpha-L): only the switched residue separates all three
#   (0.12 / 1.40 / 1.50); one residue downstream beta and alpha-L are too close.
EXPECTED_CLUSTERS = {
    "two_state_switch.pdb": {0: 2, 1: 2},
    "three_state_switch.pdb": {0: 3},
}


@dataclass(frozen=True)
class SwitchState:
    first_model: int
    last_model: int
    switch_phi_psi: tuple[float, float]


@dataclass(frozen=True)
class SwitchEnsemble:
    path: Path
    n_residues: int
    switch_residue: int
    noise_deg: float
    helix_phi_psi: tuple[float, float]
    states: tuple[SwitchState, ...]

    @classmethod
    def from_pdb(cls, path: Path) -> SwitchEnsemble:
        meta = {}
        for line in path.read_text().splitlines():
            if line.startswith("REMARK 999 "):
                key, value = line[len("REMARK 999 ") :].split(" ", 1)
                meta[key] = value

        def pair(text: str) -> tuple[float, float]:
            phi, psi = (float(value) for value in text.split())
            return phi, psi

        helix = pair(meta["state_a_phi_psi"])
        states = []
        for letter in "abc":
            if f"state_{letter}_models" not in meta:
                break
            first, last = (int(value) for value in meta[f"state_{letter}_models"].split("-"))
            switch = pair(meta[f"state_{letter}_switch_phi_psi"]) if letter != "a" else helix
            states.append(SwitchState(first, last, switch))
        return cls(
            path=path,
            n_residues=int(meta["n_residues"]),
            switch_residue=int(meta["switch_residue"]),
            noise_deg=float(meta["dihedral_noise_deg"]),
            helix_phi_psi=helix,
            states=tuple(states),
        )

    def state_of(self, model: int) -> int:
        for index, state in enumerate(self.states):
            if state.first_model <= model <= state.last_model:
                return index
        raise KeyError(model)

    @property
    def switch_windows(self) -> list[str]:
        """Residue ranges containing the switch: tight, wider, and the whole chain."""
        switch = self.switch_residue
        return [f"{switch}-{switch + 1}", f"{switch - 2}-{switch + 3}", f"1-{self.n_residues}"]

    @property
    def far_windows(self) -> list[str]:
        """Residue ranges covering the distant residues on either side of the switch."""
        switch = self.switch_residue
        return [f"1-{switch - 3}", f"{switch + 4}-{self.n_residues}"]

    @property
    def distant_residues(self) -> tuple[int, ...]:
        """Residues at least three positions before or four after the switch."""
        return tuple(range(1, self.switch_residue - 2)) + tuple(
            range(self.switch_residue + 4, self.n_residues + 1)
        )


@pytest.fixture(scope="module", params=sorted(EXPECTED_CLUSTERS))
def analysed(request):
    ensemble = SwitchEnsemble.from_pdb(DATA_DIR / request.param)
    config = AnalysisConfig(
        pdb_file=ensemble.path,
        clustering=ClusteringConfig(
            cluster_residues=True,
            cluster_residue_ranges=ensemble.switch_windows + ensemble.far_windows,
            min_cluster_size=MIN_CLUSTER_SIZE,
        ),
    )
    return ensemble, FlexGeo2App().run(config)


def test_fixture_matches_its_header(analysed) -> None:
    ensemble, result = analysed
    raw_df = result.raw_df
    n_models = ensemble.states[-1].last_model
    assert sorted(raw_df["model"].unique()) == list(range(1, n_models + 1))
    assert sorted(raw_df["order"].unique()) == list(range(1, ensemble.n_residues + 1))

    # Melodia leaves phi undefined for the first residue and psi for the last.
    inner = raw_df[raw_df["order"].between(2, ensemble.n_residues - 1)]
    states = inner["model"].map(ensemble.state_of)
    at_switch = inner["order"] == ensemble.switch_residue
    for axis, column in enumerate(("phi", "psi")):
        switch_values = states.map(lambda s, axis=axis: ensemble.states[s].switch_phi_psi[axis])
        expected = np.where(at_switch, switch_values, ensemble.helix_phi_psi[axis])
        assert np.abs(inner[column] - expected).max() < 5 * ensemble.noise_deg, column


def test_each_switch_swings_the_rest_of_the_chain(analysed) -> None:
    ensemble, _ = analysed
    structure = GeometryService.parse_structure(ensemble.path)
    ca = np.array([[residue["CA"].get_coord() for residue in model["A"]] for model in structure])
    states = np.array([ensemble.state_of(serial) for serial in range(1, len(ca) + 1)])
    reference = ca[states == 0].mean(axis=0)

    # Mean CA displacement from the all-helix state. Residues before the switch stay put
    # (noise alone shifts them by < 0.5 A). After it, the largest displacement measured
    # 20.3 A for beta and 11.1 A for alpha-L, against ~2 A within-state noise spread.
    # The residue that moves most varies: an alpha-L kink returns the chain end close
    # to its helix position.
    before = slice(0, ensemble.switch_residue - 1)
    after = slice(ensemble.switch_residue, None)
    for state in range(1, len(ensemble.states)):
        shift = np.linalg.norm(ca[states == state].mean(axis=0) - reference, axis=1)
        assert shift[before].max() < 0.5, f"state {state}"
        assert shift[after].max() > 10.0, (
            f"state {state} moves the chain by only {shift[after].max():.1f} A"
        )


def test_clustering_recovers_every_state_at_the_switch(analysed) -> None:
    ensemble, result = analysed
    summary = result.residue_clustering.summary_df.set_index("order")
    assignments = result.residue_clustering.assignments_df
    expected = EXPECTED_CLUSTERS[ensemble.path.name]

    for offset, n_states in expected.items():
        residue = ensemble.switch_residue + offset
        assert summary.loc[residue, "n_clusters"] == n_states, f"residue {residue}"
        assert summary.loc[residue, "noise_fraction"] == 0.0, f"residue {residue}"

        residue_df = assignments[assignments["order"] == residue]
        labels_by_state = residue_df.groupby(residue_df["model"].map(ensemble.state_of))[
            "cluster"
        ].unique()
        assert all(len(labels) == 1 for labels in labels_by_state), f"residue {residue}"
        assert len({labels[0] for labels in labels_by_state}) == n_states, f"residue {residue}"


def test_descriptors_far_from_the_switch_are_state_independent(analysed) -> None:
    ensemble, result = analysed
    raw_df = result.raw_df
    means = (
        raw_df.assign(state=raw_df["model"].map(ensemble.state_of))
        .groupby(["order", "state"])[["curvature", "torsion"]]
        .mean()
    )

    def largest_state_shift(residue: int) -> float:
        return max(
            float(np.hypot(*(means.loc[(residue, a)] - means.loc[(residue, b)])))
            for a, b in combinations(range(len(ensemble.states)), 2)
        )

    switch_shift = max(largest_state_shift(ensemble.switch_residue + offset) for offset in (0, 1))
    distant_shift = max(largest_state_shift(residue) for residue in ensemble.distant_residues)

    assert switch_shift > 0.3
    assert distant_shift < 0.05
    assert distant_shift < switch_shift / 10


# Range clustering, measured over 50 noise seeds for both fixtures: every window that
# contains the switch (up to the whole chain) recovered all states exactly. Windows away
# from it were pure noise until single clusters were allowed; they are now one cluster.


def range_assignments(result, window: str):
    assignments = result.residue_range_clustering.assignments_df
    window_df = assignments[assignments["range_label"] == window]
    # Range assignments store model IDs as strings (see TODO.md).
    return window_df.assign(model=window_df["model"].astype(int))


def test_range_clustering_recovers_every_state_for_windows_covering_the_switch(
    analysed,
) -> None:
    ensemble, result = analysed
    summary = result.residue_range_clustering.summary_df.set_index("range_label")
    n_states = len(ensemble.states)

    for window in ensemble.switch_windows:
        assert summary.loc[window, "n_clusters"] == n_states, window
        assert summary.loc[window, "noise_fraction"] == 0.0, window

        window_df = range_assignments(result, window)
        labels_by_state = window_df.groupby(window_df["model"].map(ensemble.state_of))[
            "cluster"
        ].unique()
        assert all(len(labels) == 1 for labels in labels_by_state), window
        assert len({labels[0] for labels in labels_by_state}) == n_states, window


def test_range_clustering_finds_one_state_far_from_the_switch(analysed) -> None:
    ensemble, result = analysed
    summary = result.residue_range_clustering.summary_df.set_index("range_label")

    # One state. HDBSCAN keeps only its dense core and labels the rest noise (65-90% here).
    for window in ensemble.far_windows:
        assert summary.loc[window, "n_clusters"] == 1, window


def test_single_state_residues_form_one_cluster(analysed) -> None:
    ensemble, result = analysed
    summary = result.residue_clustering.summary_df.set_index("order")

    distant = summary.loc[list(ensemble.distant_residues), "n_clusters"]
    assert (distant == 1).all(), distant[distant != 1].to_dict()


def test_range_pca_coordinates_are_finite(analysed) -> None:
    ensemble, result = analysed

    for window in ensemble.switch_windows + ensemble.far_windows:
        window_df = range_assignments(result, window)
        assert np.isfinite(window_df[["pc1", "pc2"]].to_numpy()).all(), window
