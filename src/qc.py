"""Stage 3: load, project, quality-filter and deduplicate the spike dataset."""
from __future__ import annotations

import collections
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

import align
import spike_reference as ref

MIN_SPIKE_LENGTH_AA = 1200
MAX_AMBIGUOUS_FRACTION = 0.02


def read_fasta(path: str | Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(chunks)
                header, chunks = line[1:].strip().split("|")[0], []
            else:
                chunks.append(line.strip())
    if header:
        seqs[header] = "".join(chunks)
    return seqs


def parse_collection_date(value: str) -> pd.Timestamp | None:
    """Collection dates appear as YYYY, YYYY-MM or YYYY-MM-DD."""
    if not value or not isinstance(value, str):
        return None
    text = value.strip()
    for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
        try:
            return pd.to_datetime(text, format=fmt)
        except (ValueError, TypeError):
            continue
    return None


@dataclass
class QCReport:
    """Record counts at each filtering step, for the Stage 3 before/after table."""
    steps: list[tuple[str, int]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    pre_dedup: "pd.DataFrame | None" = None

    def add(self, label: str, n: int) -> None:
        self.steps.append((label, n))

    def note(self, text: str) -> None:
        self.notes.append(text)

    def to_frame(self) -> pd.DataFrame:
        rows = []
        prev = None
        for label, n in self.steps:
            rows.append({"step": label, "records": n,
                         "removed": "" if prev is None else prev - n})
            prev = n
        return pd.DataFrame(rows)


def build_dataset(
    metadata_csv: str | Path,
    fasta_path: str | Path,
    reference: str,
    min_spike_length_aa: int = MIN_SPIKE_LENGTH_AA,
    max_ambiguous_fraction: float = MAX_AMBIGUOUS_FRACTION,
    verbose: bool = True,
) -> tuple[pd.DataFrame, QCReport]:
    """Return the analysis-ready table plus a record of what was filtered."""
    report = QCReport()
    meta = pd.read_csv(metadata_csv)
    seqs = read_fasta(fasta_path)

    report.add("Downloaded metadata records", len(meta))
    meta = meta.drop_duplicates(subset="accession", keep="first")
    report.add("Unique accessions", len(meta))

    meta["raw_spike"] = meta["accession"].map(seqs)
    meta = meta[meta["raw_spike"].notna()].copy()
    report.add("With a spike sequence", len(meta))

    meta["raw_length"] = meta["raw_spike"].str.len()
    meta = meta[meta["raw_length"] >= min_spike_length_aa].copy()
    report.add(f"Spike length >= {min_spike_length_aa} aa", len(meta))

    meta["collection_date"] = meta["collection_date"].map(parse_collection_date)
    meta = meta[meta["collection_date"].notna()].copy()
    report.add("With a parsable collection date", len(meta))

    aligner = align.make_aligner()
    projected, inserted = [], []
    for n, seq in enumerate(meta["raw_spike"], start=1):
        p, ins = align.project_to_reference(seq, reference, aligner)
        projected.append(align.normalise_gaps(p, reference))
        inserted.append(ins)
        if verbose and n % 250 == 0:
            print(f"    projected {n}/{len(meta)}", flush=True)
    meta["spike"] = projected
    meta["n_inserted"] = inserted

    meta["n_ambiguous"] = meta["spike"].str.count("X")
    meta["ambiguous_fraction"] = meta["n_ambiguous"] / len(reference)
    meta = meta[meta["ambiguous_fraction"] <= max_ambiguous_fraction].copy()
    report.add(f"Ambiguity <= {max_ambiguous_fraction:.0%} of positions", len(meta))

    meta["mutation_profile"] = [
        ";".join(align.format_mutation(*m)
                 for m in align.mutations_vs_reference(s, reference))
        for s in meta["spike"]
    ]

    groups = meta.groupby("mutation_profile")
    dup_groups = [g for _, g in groups if len(g) > 1]
    n_dup_records = sum(len(g) for g in dup_groups)
    cross_region = [g for g in dup_groups if g["region"].nunique() > 1]
    cross_variant = [g for g in dup_groups if g["variant"].nunique() > 1]
    report.note(
        f"{len(dup_groups)} distinct sequences occur more than once, covering "
        f"{n_dup_records} records.")
    report.note(
        f"{len(cross_region)} of those duplicate groups span more than one region. "
        f"Where one spike sequence occurs in several regions, no model reading only "
        f"the sequence can place all of them correctly, so this sets a ceiling on "
        f"geographic accuracy.")
    report.note(
        f"{len(cross_variant)} duplicate groups span more than one variant label.")

    report.pre_dedup = meta.copy()
    meta = meta.drop_duplicates(subset="mutation_profile", keep="first").copy()
    report.add("Unique spike mutation profiles", len(meta))

    meta["n_mutations"] = [
        len(align.mutations_vs_reference(s, reference)) for s in meta["spike"]
    ]
    meta["n_rbd_mutations"] = [
        sum(1 for _, pos, _ in align.mutations_vs_reference(s, reference)
            if ref.in_rbd(pos))
        for s in meta["spike"]
    ]
    meta["year_month"] = meta["collection_date"].dt.to_period("M").astype(str)
    return meta.reset_index(drop=True), report


def duplicate_region_ceiling(meta_before_dedup: pd.DataFrame) -> float:
    """Best accuracy any sequence-only region classifier could reach."""
    total = len(meta_before_dedup)
    if not total:
        return float("nan")
    best = sum(
        collections.Counter(g["region"]).most_common(1)[0][1]
        for _, g in meta_before_dedup.groupby("mutation_profile")
    )
    return best / total
