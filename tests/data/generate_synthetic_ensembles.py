"""Generate synthetic backbone ensembles used as test fixtures.

Run manually from the repository root; pytest does not execute this script:

    python tests/data/generate_synthetic_ensembles.py

Each ensemble is an ideal poly-alanine chain built from backbone dihedrals with the
NeRF method (N, CA and C atoms only, which is all Melodia needs). Parameters are
written to REMARK 999 lines in the PDB header, and the tests read them back.

two_state_switch.pdb
    A 15-residue alpha helix in which residue 8 flips to beta in the second half of
    the models. The flip swings the C-terminal half as a rigid body, so FleXgeo2
    should detect two states at the switch and no change far from it.

three_state_switch.pdb
    The same helix with residue 8 in one of three states: alpha (all helix), beta or
    left-handed helix (alpha-L). Polyproline II was tried as the third state but sits
    too close to alpha in curvature/torsion at the switch for robust clustering.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent

# Ideal backbone geometry (Engh & Huber).
BOND_N_CA, BOND_CA_C, BOND_C_N = 1.458, 1.525, 1.329
ANGLE_N_CA_C, ANGLE_CA_C_N, ANGLE_C_N_CA = 111.2, 116.2, 121.7
OMEGA = 180.0

HELIX_PHI_PSI = (-60.0, -45.0)
BETA_PHI_PSI = (-120.0, 130.0)
ALPHA_L_PHI_PSI = (57.0, 47.0)

TWO_STATE = {
    "n_residues": 15,
    "switch_residue": 8,
    "models_per_state": 20,
    "dihedral_noise_deg": 2.0,
    "seed": 0,
}
THREE_STATE = dict(TWO_STATE)


def place_atom(a, b, c, bond: float, angle: float, torsion: float) -> np.ndarray:
    """Place atom D from A, B, C given |CD|, angle BCD and dihedral ABCD (NeRF)."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    normal = np.cross(b - a, bc)
    normal /= np.linalg.norm(normal)
    frame = np.column_stack([bc, np.cross(normal, bc), normal])
    local = bond * np.array(
        [-np.cos(angle), np.sin(angle) * np.cos(torsion), np.sin(angle) * np.sin(torsion)]
    )
    return c + frame @ local


def build_backbone(phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
    """Return an (n_residues, 3, 3) array of N, CA, C coordinates."""
    n = np.zeros(3)
    ca = np.array([BOND_N_CA, 0.0, 0.0])
    theta = np.radians(180.0 - ANGLE_N_CA_C)
    c = ca + BOND_CA_C * np.array([np.cos(theta), np.sin(theta), 0.0])
    residues = [(n, ca, c)]
    for index in range(len(phi) - 1):
        n, ca, c = residues[-1]
        next_n = place_atom(n, ca, c, BOND_C_N, ANGLE_CA_C_N, psi[index])
        next_ca = place_atom(ca, c, next_n, BOND_N_CA, ANGLE_C_N_CA, OMEGA)
        next_c = place_atom(c, next_n, next_ca, BOND_CA_C, ANGLE_N_CA_C, phi[index + 1])
        residues.append((next_n, next_ca, next_c))
    return np.array(residues)


def dihedral(a, b, c, d) -> float:
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def check_backbone(backbone: np.ndarray, phi=None, psi=None) -> None:
    """Assert ideal bonds, planar peptides and (optionally) the requested phi/psi."""
    n, ca, c = backbone[:, 0], backbone[:, 1], backbone[:, 2]
    assert np.allclose(np.linalg.norm(ca - n, axis=1), BOND_N_CA)
    assert np.allclose(np.linalg.norm(c - ca, axis=1), BOND_CA_C)
    assert np.allclose(np.linalg.norm(n[1:] - c[:-1], axis=1), BOND_C_N)
    for i in range(len(backbone) - 1):
        assert np.isclose(abs(dihedral(ca[i], c[i], n[i + 1], ca[i + 1])), OMEGA)
        if psi is not None:
            assert np.isclose(dihedral(n[i], ca[i], c[i], n[i + 1]), psi[i])
        if phi is not None:
            assert np.isclose(dihedral(c[i], n[i + 1], ca[i + 1], c[i + 1]), phi[i + 1])


def switch_dihedrals(switch_phi_psi, params: dict) -> tuple[np.ndarray, np.ndarray]:
    """All-helix dihedrals, with the switch residue set to ``switch_phi_psi`` if given."""
    phi = np.full(params["n_residues"], HELIX_PHI_PSI[0])
    psi = np.full(params["n_residues"], HELIX_PHI_PSI[1])
    if switch_phi_psi is not None:
        phi[params["switch_residue"] - 1], psi[params["switch_residue"] - 1] = switch_phi_psi
    return phi, psi


def pdb_lines(models: list[np.ndarray], remarks: dict) -> list[str]:
    lines = [f"REMARK 999 {key} {value}" for key, value in remarks.items()]
    serial = 1
    for model_number, backbone in enumerate(models, start=1):
        lines.append(f"MODEL     {model_number:>4}")
        for residue_number, atoms in enumerate(backbone, start=1):
            for name, (x, y, z) in zip(("N", "CA", "C"), atoms, strict=True):
                lines.append(
                    f"ATOM  {serial:>5}  {name:<3} ALA A{residue_number:>4}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]:>2}"
                )
                serial += 1
        lines.append("ENDMDL")
    lines.append("END")
    return lines


def switch_ensemble(switches: list, params: dict) -> list[np.ndarray]:
    """Models for each state in turn; ``switches`` holds each state's switch phi/psi."""
    for switch_phi_psi in switches:
        ideal_phi, ideal_psi = switch_dihedrals(switch_phi_psi, params)
        check_backbone(build_backbone(ideal_phi, ideal_psi), ideal_phi, ideal_psi)

    rng = np.random.default_rng(params["seed"])
    noise = params["dihedral_noise_deg"]
    models = []
    for switch_phi_psi in switches:
        for _ in range(params["models_per_state"]):
            phi, psi = switch_dihedrals(switch_phi_psi, params)
            phi = phi + rng.normal(0.0, noise, params["n_residues"])
            psi = psi + rng.normal(0.0, noise, params["n_residues"])
            backbone = build_backbone(phi, psi)
            check_backbone(backbone, phi, psi)
            models.append(backbone)
    return models


def state_remarks(switches: list, params: dict) -> dict:
    """REMARK entries: state_a is all helix; later states give their switch phi/psi."""
    per_state = params["models_per_state"]
    remarks = dict(params)
    for index, letter in enumerate("abc"[: len(switches)]):
        remarks[f"state_{letter}_models"] = f"{index * per_state + 1}-{(index + 1) * per_state}"
    remarks["state_a_phi_psi"] = "{} {}".format(*HELIX_PHI_PSI)
    for letter, switch_phi_psi in zip("bc", switches[1:], strict=False):
        remarks[f"state_{letter}_switch_phi_psi"] = "{} {}".format(*switch_phi_psi)
    return remarks


def generate_two_state(params: dict = TWO_STATE) -> list[str]:
    switches = [None, BETA_PHI_PSI]
    return pdb_lines(switch_ensemble(switches, params), state_remarks(switches, params))


def generate_three_state(params: dict = THREE_STATE) -> list[str]:
    switches = [None, BETA_PHI_PSI, ALPHA_L_PHI_PSI]
    return pdb_lines(switch_ensemble(switches, params), state_remarks(switches, params))


def main() -> None:
    for name, generate in (
        ("two_state_switch.pdb", generate_two_state),
        ("three_state_switch.pdb", generate_three_state),
    ):
        output = DATA_DIR / name
        output.write_text("\n".join(generate()) + "\n")
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
