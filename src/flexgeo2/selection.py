"""Residue selections such as ``45``, ``45-50``, ``A:45`` or ``A:45-50``."""

from __future__ import annotations

import re
from dataclasses import dataclass, field

# Optional "CHAIN:" prefix, then a residue number or START-END; numbers may be negative.
_SELECTION = re.compile(r"^(?:(?P<chain>[^:\s]+):)?(?P<start>-?\d+)(?:-(?P<end>-?\d+))?$")


@dataclass(frozen=True, slots=True)
class ResidueSelection:
    """Residues ``start`` to ``end`` (inclusive) of ``chain``, or of every chain if None."""

    chain: str | None
    start: int
    end: int
    text: str = field(default="", compare=False)

    @property
    def orders(self) -> range:
        return range(self.start, self.end + 1)


def parse_residue_selections(texts: list[str]) -> list[ResidueSelection]:
    """Parse selections; each text may hold several, separated by commas."""
    selections = []
    for text in texts:
        for item in (part.strip() for part in text.split(",")):
            if not item:
                continue
            match = _SELECTION.match(item)
            if match is None:
                raise ValueError(
                    f"Invalid residue selection '{item}'. Use a residue number or range, "
                    "optionally with a chain: 45, 45-50, A:45 or A:45-50."
                )
            start = int(match["start"])
            end = int(match["end"]) if match["end"] is not None else start
            if end < start:
                raise ValueError(
                    f"Invalid residue selection '{item}': the range ends before it starts."
                )
            selections.append(ResidueSelection(match["chain"], start, end, text=item))
    return selections


def select_residues(
    selections: list[ResidueSelection], residues_by_chain: dict[str, set[int]]
) -> list[tuple[str, int]]:
    """``(chain, residue number)`` pairs chosen by ``selections``, sorted.

    ``residues_by_chain`` holds the residue numbers of each analysed chain. A selection
    without a chain applies to every chain. Every residue a selection names must exist
    in at least one of its chains; otherwise a ``ValueError`` names the missing ones.
    """
    chosen: set[tuple[str, int]] = set()
    for selection in selections:
        if selection.chain is None:
            chains = sorted(residues_by_chain)
        elif selection.chain in residues_by_chain:
            chains = [selection.chain]
        else:
            raise ValueError(
                f"Residue selection '{selection.text}': chain '{selection.chain}' is not "
                f"among the analysed chains ({', '.join(sorted(residues_by_chain))})."
            )
        found = {
            (chain, order)
            for chain in chains
            for order in selection.orders
            if order in residues_by_chain[chain]
        }
        missing = sorted(set(selection.orders) - {order for _, order in found})
        if missing:
            shown = ", ".join(str(order) for order in missing[:20])
            more = f", ... ({len(missing)} total)" if len(missing) > 20 else ""
            raise ValueError(
                f"Residue selection '{selection.text}': residue(s) {shown}{more} not found "
                f"in chain(s) {', '.join(chains)}."
            )
        chosen |= found
    return sorted(chosen)
