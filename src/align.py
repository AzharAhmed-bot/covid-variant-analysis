"""Project spike protein sequences onto Wuhan-Hu-1 reference coordinates.

Every analysed sequence is represented as a fixed-length string of
EXPECTED_SPIKE_LENGTH (1273) characters, one per reference residue:

    an amino-acid letter  the residue aligned to that reference position
    '-'                   deletion relative to the reference
    'X'                   ambiguous or unresolved residue

Insertions relative to the reference are recorded separately rather than
widening the coordinate system, so that position 501 always means S:501 in the
numbering used by the literature. Insertions are rare in spike and are counted
but not used as model features; this is stated as a limitation.
"""
from __future__ import annotations

from Bio import Align

AMBIGUOUS = set("XBZJUO*")


def make_aligner() -> Align.PairwiseAligner:
    """Global protein aligner with gap costs that favour few, long gaps.

    BLOSUM62 with an open/extend penalty chosen so that the known spike
    deletions (for example the 69-70 and 142-144 deletions) come out as single
    contiguous gaps rather than scattered single-residue gaps.
    """
    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    # Do not penalise gaps at the sequence ends: many records are 5'/3' truncated.
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
    """Align `query` to `reference` and return it in reference coordinates.

    Returns (projected, n_inserted_residues) where `projected` has the same
    length as `reference`.
    """
    aligner = aligner or make_aligner()
    query = query.upper().replace("*", "")

    alignment = aligner.align(reference, query)[0]
    ref_idx, qry_idx = alignment.aligned

    projected = ["-"] * len(reference)
    matched_ref = 0
    for (r_start, r_end), (q_start, q_end) in zip(ref_idx, qry_idx):
        block = query[q_start:q_end]
        for offset, residue in enumerate(block):
            projected[r_start + offset] = "X" if residue in AMBIGUOUS else residue
        matched_ref += r_end - r_start

    n_inserted = len(query) - matched_ref
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


# Regions where protein-level alignment cannot reliably place deletions.
# A contiguous deletion run next to a substituted residue admits several
# placements with different BLOSUM scores, so the highest-scoring alignment is
# not necessarily the one curated annotation reports. In BA.1, residues 142-144
# are deleted and 145 is substituted (Y145D); BLOSUM62 scores D against G142
# better than against Y145, so the aligner prefers the non-curated placement.
# Nucleotide-level codon-aware alignment (as used by Nextclade) avoids this.
# Per-position deletion calls in these windows are therefore treated as
# uncertain, and deletions are also summarised per region.
INDEL_AMBIGUOUS_REGIONS = [(69, 70), (141, 146), (156, 158), (211, 212), (241, 243)]


def is_indel_ambiguous(pos: int) -> bool:
    """True if a deletion call at this position is alignment-ambiguous."""
    return any(start <= pos <= end for start, end in INDEL_AMBIGUOUS_REGIONS)


def normalise_gaps(projected: str, reference: str) -> str:
    """Left-shift gap runs through score-equivalent placements.

    A gap run moves one position left when the residue displaced from its left
    edge lands on an identical reference residue, which leaves the alignment
    score unchanged. Applying this consistently makes deletion placement
    deterministic across sequences wherever several alignments tie.
    """
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
