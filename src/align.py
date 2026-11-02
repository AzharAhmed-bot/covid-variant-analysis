"""Project spike protein sequences onto Wuhan-Hu-1 reference coordinates."""
from __future__ import annotations

from Bio import Align

AMBIGUOUS = set("XBZJUO*")


def make_aligner() -> Align.PairwiseAligner:
    """Global protein aligner with gap costs that favour few, long gaps."""
    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    for new_name, old_name in (("end_insertion_score", "target_end_gap_score"),
                               ("end_deletion_score", "query_end_gap_score")):
        try:
            setattr(aligner, new_name, 0.0)
        except AttributeError:
            setattr(aligner, old_name, 0.0)
    return aligner


def project_to_reference(
    query: str,
    reference: str,
    aligner: Align.PairwiseAligner | None = None,
) -> tuple[str, int]:
    """Align `query` to `reference` and return it in reference coordinates."""
    aligner = aligner or make_aligner()
    query = query.upper().replace("*", "")

    alignment = aligner.align(reference, query)[0]
    ref_idx, qry_idx = alignment.aligned

    projected = ["-"] * len(reference)
    reference_positions_matched = 0
    for (r_start, r_end), (q_start, q_end) in zip(ref_idx, qry_idx):
        block = query[q_start:q_end]
        for offset, residue in enumerate(block):
            projected[r_start + offset] = "X" if residue in AMBIGUOUS else residue
        reference_positions_matched += r_end - r_start

    n_inserted = len(query) - reference_positions_matched
    return "".join(projected), max(n_inserted, 0)


def mutations_vs_reference(
    projected: str,
    reference: str,
    skip_ambiguous: bool = True,
) -> list[tuple[str, int, str]]:
    """List (ref_aa, 1-based position, alt_aa) differences from the reference."""
    out = []
    for i, (ref_aa, alt_aa) in enumerate(zip(reference, projected), start=1):
        if alt_aa == ref_aa:
            continue
        if skip_ambiguous and alt_aa == "X":
            continue
        out.append((ref_aa, i, alt_aa))
    return out


def format_mutation(ref_aa: str, pos: int, alt_aa: str) -> str:
    """Render a mutation in the conventional notation, e.g. N501Y or H69del."""
    return f"{ref_aa}{pos}del" if alt_aa == "-" else f"{ref_aa}{pos}{alt_aa}"


INDEL_AMBIGUOUS_REGIONS = [(69, 70), (141, 146), (156, 158), (211, 212), (241, 243)]


def is_indel_ambiguous(pos: int) -> bool:
    """True if a deletion call at this position is alignment-ambiguous."""
    return any(start <= pos <= end for start, end in INDEL_AMBIGUOUS_REGIONS)


def normalise_gaps(projected: str, reference: str) -> str:
    """Left-shift gap runs through score-equivalent placements."""
    chars = list(projected)
    n = len(chars)
    i = 0
    while i < n:
        if chars[i] != "-":
            i += 1
            continue
        j = i
        while j + 1 < n and chars[j + 1] == "-":
            j += 1
        start, end = i, j
        while (start > 0 and chars[start - 1] != "-"
               and reference[start - 1] == reference[end]):
            chars[end] = chars[start - 1]
            chars[start - 1] = "-"
            start -= 1
            end -= 1
        i = j + 1
    return "".join(chars)


def deletion_runs(projected: str) -> list[tuple[int, int]]:
    """Maximal runs of deleted positions as (start, end), 1-based inclusive."""
    runs, start = [], None
    for idx, ch in enumerate(projected, start=1):
        if ch == "-" and start is None:
            start = idx
        elif ch != "-" and start is not None:
            runs.append((start, idx - 1))
            start = None
    if start is not None:
        runs.append((start, len(projected)))
    return runs
