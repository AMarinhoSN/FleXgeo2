from __future__ import annotations

import pytest

from flexgeo2.selection import ResidueSelection, parse_residue_selections, select_residues


def test_parse_accepts_numbers_ranges_chains_and_lists() -> None:
    selections = parse_residue_selections(["45", "A:50-52, B:-3", "-5--2,7"])

    assert selections == [
        ResidueSelection(None, 45, 45),
        ResidueSelection("A", 50, 52),
        ResidueSelection("B", -3, -3),
        ResidueSelection(None, -5, -2),
        ResidueSelection(None, 7, 7),
    ]
    assert [selection.text for selection in selections] == ["45", "A:50-52", "B:-3", "-5--2", "7"]


@pytest.mark.parametrize("text", ["A", "A:", "45-", "4 5", "A:B:3", "45-50-55", "1.5"])
def test_parse_rejects_malformed_selections(text: str) -> None:
    with pytest.raises(ValueError, match="Invalid residue selection"):
        parse_residue_selections([text])


def test_parse_rejects_reversed_ranges() -> None:
    with pytest.raises(ValueError, match="ends before it starts"):
        parse_residue_selections(["A:52-50"])


RESIDUES = {"A": {1, 2, 3, 4}, "B": {3, 4, 5}}


def test_select_without_chain_uses_every_chain_that_has_the_residue() -> None:
    chosen = select_residues(parse_residue_selections(["3-5", "A:1", "A:1"]), RESIDUES)

    assert chosen == [("A", 1), ("A", 3), ("A", 4), ("B", 3), ("B", 4), ("B", 5)]


def test_select_with_chain_keeps_only_that_chain() -> None:
    # Residue 3 exists in both chains.
    assert select_residues(parse_residue_selections(["B:3"]), RESIDUES) == [("B", 3)]


def test_select_rejects_unknown_chain() -> None:
    with pytest.raises(ValueError, match=r"chain 'C' is not among the analysed chains \(A, B\)"):
        select_residues(parse_residue_selections(["C:3"]), RESIDUES)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("9", r"'9': residue\(s\) 9 not found in chain\(s\) A, B"),
        ("B:1-3", r"'B:1-3': residue\(s\) 1, 2 not found in chain\(s\) B"),
    ],
)
def test_select_names_missing_residues(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        select_residues(parse_residue_selections([text]), RESIDUES)
