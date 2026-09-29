"""Synthetic two-state ensemble: a known local switch, analysed end to end.

The fixture (``tests/data/two_state_switch.pdb``, made by
``tests/data/generate_synthetic_ensembles.py``) is an ideal poly-alanine alpha helix in
which one residue flips to beta in the second half of the models. The flip swings the
rest of the chain as a rigid body. FleXgeo2 should find two clusters at the switch and,
being superposition free, see no change far from it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from flexgeo2 import AnalysisConfig, ClusteringConfig, FlexGeo2App, OutputConfig
from flexgeo2.geometry import GeometryService

pytest.importorskip("melodia_py")

FIXTURE = Path(__file__).parent / "data" / "two_state_switch.pdb"
MIN_CLUSTER_SIZE = 10


def read_remarks(pdb_file: Path) -> dict[str, str]:
    remarks = {}
    for line in pdb_file.read_text().splitlines():
        if line.startswith("REMARK 999 "):
            key, value = line[len("REMARK 999 ") :].split(" ", 1)
            remarks[key] = value
    return remarks


META = read_remarks(FIXTURE)
N_RESIDUES = int(META["n_residues"])
SWITCH_RESIDUE = int(META["switch_residue"])
MODELS_PER_STATE = int(META["models_per_state"])
NOISE_DEG = float(META["dihedral_noise_deg"])
HELIX_PHI_PSI = tuple(float(value) for value in META["state_a_phi_psi"].split())
BETA_PHI_PSI = tuple(float(value) for value in META["state_b_switch_phi_psi"].split())

# Measured with Melodia on the noise-free states: the switch moves curvature/torsion at
# the switched residue (~0.12) and the next one (~0.43). The two residues after that
# respond weakly (~0.02-0.03, spline smoothing) and are left out of both groups;
# everything else stays within ~0.006.
RESPONDING_RESIDUES = (SWITCH_RESIDUE, SWITCH_RESIDUE + 1)
DISTANT_RESIDUES = tuple(range(1, SWITCH_RESIDUE - 2)) + tuple(
    range(SWITCH_RESIDUE + 4, N_RESIDUES + 1)
)


def model_state(model: int) -> int:
    """State of a PDB MODEL number: 0 (all helix) or 1 (switch residue in beta)."""
    return 0 if model <= MODELS_PER_STATE else 1


@pytest.fixture(scope="module")
def result():
    config = AnalysisConfig(
        pdb_file=FIXTURE,
        clustering=ClusteringConfig(cluster_residues=True, min_cluster_size=MIN_CLUSTER_SIZE),
        output=OutputConfig(write_files=False),
    )
    return FlexGeo2App().run(config)


def test_fixture_matches_its_header(result) -> None:
    raw_df = result.raw_df
    assert sorted(raw_df["model"].unique()) == list(range(1, 2 * MODELS_PER_STATE + 1))
    assert sorted(raw_df["order"].unique()) == list(range(1, N_RESIDUES + 1))

    # Melodia leaves phi undefined for the first residue and psi for the last.
    inner = raw_df[raw_df["order"].between(2, N_RESIDUES - 1)]
    is_switch = (inner["order"] == SWITCH_RESIDUE) & (inner["model"].map(model_state) == 1)
    expected_phi = np.where(is_switch, BETA_PHI_PSI[0], HELIX_PHI_PSI[0])
    expected_psi = np.where(is_switch, BETA_PHI_PSI[1], HELIX_PHI_PSI[1])
    assert np.abs(inner["phi"] - expected_phi).max() < 5 * NOISE_DEG
    assert np.abs(inner["psi"] - expected_psi).max() < 5 * NOISE_DEG


def test_switch_swings_the_rest_of_the_chain() -> None:
    structure = GeometryService.parse_structure(FIXTURE)
    ca = np.array([[residue["CA"].get_coord() for residue in model["A"]] for model in structure])
    states = np.array([model_state(serial) for serial in range(1, len(ca) + 1)])

    mean_shift = np.linalg.norm(ca[states == 0].mean(axis=0) - ca[states == 1].mean(axis=0), axis=1)

    assert mean_shift[0] == pytest.approx(0.0, abs=1e-3)
    assert mean_shift[-1] > 15.0


@pytest.mark.parametrize("residue", RESPONDING_RESIDUES)
def test_clustering_recovers_both_states_at_the_switch(result, residue: int) -> None:
    summary = result.residue_clustering.summary_df.set_index("order")
    assignments = result.residue_clustering.assignments_df

    assert summary.loc[residue, "n_clusters"] == 2
    assert summary.loc[residue, "noise_fraction"] == 0.0

    residue_df = assignments[assignments["order"] == residue]
    labels_by_state = residue_df.groupby(residue_df["model"].map(model_state))["cluster"].unique()
    assert all(len(labels) == 1 for labels in labels_by_state)
    assert labels_by_state[0][0] != labels_by_state[1][0]


def test_descriptors_far_from_the_switch_are_state_independent(result) -> None:
    raw_df = result.raw_df
    means = (
        raw_df.assign(state=raw_df["model"].map(model_state))
        .groupby(["order", "state"])[["curvature", "torsion"]]
        .mean()
    )

    def state_shift(residue: int) -> float:
        return float(np.hypot(*(means.loc[(residue, 0)] - means.loc[(residue, 1)])))

    switch_shift = max(state_shift(residue) for residue in RESPONDING_RESIDUES)
    distant_shift = max(state_shift(residue) for residue in DISTANT_RESIDUES)

    assert switch_shift > 0.3
    assert distant_shift < 0.05
    assert distant_shift < switch_shift / 10
