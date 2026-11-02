"""Reference layer for the SARS-CoV-2 spike project.

Provides:
  * the Wuhan-Hu-1 (MN908947) spike protein sequence, translated from the genome
  * spike domain coordinates used for Stage 6 annotation
  * variant-defining spike mutations parsed from CoVariants

All coordinates are 1-based amino-acid positions in the mature S protein
numbering used throughout the SARS-CoV-2 literature (S starts at M1).
"""
from __future__ import annotations

import json
from pathlib import Path

from Bio.Seq import Seq

# Spike CDS in MN908947 / NC_045512.2, 1-based inclusive.
SPIKE_CDS_START = 21563
SPIKE_CDS_END = 25384
EXPECTED_SPIKE_LENGTH = 1273

# Domain boundaries, 1-based inclusive, in S protein numbering.
# Non-overlapping structural partition used for per-position annotation.
SPIKE_DOMAINS = {
    "Signal peptide": (1, 13),
    "NTD": (14, 305),
    "NTD-RBD linker": (306, 318),
    "RBD": (319, 541),
    "SD1": (542, 590),
    "SD2": (591, 680),
    "S1/S2 furin site": (681, 685),
    "S2 N-terminal": (686, 815),
    "Fusion peptide": (816, 833),
    "FP-HR1 linker": (834, 911),
    "HR1": (912, 984),
    "Central helix": (985, 1035),
    "CD-proximal": (1036, 1075),
    "Connector domain": (1076, 1141),
    "HR2-proximal": (1142, 1162),
    "HR2": (1163, 1213),
    "Transmembrane": (1214, 1237),
    "Cytoplasmic tail": (1238, 1273),
}

# Receptor-binding motif: the RBD sub-region that makes direct ACE2 contact.
RBM = (437, 508)
RBD = SPIKE_DOMAINS["RBD"]
FURIN_SITE = SPIKE_DOMAINS["S1/S2 furin site"]

# Nextstrain clade -> WHO variant label, for the groups this project analyses.
CLADE_TO_VARIANT = {
    "20I (Alpha, V1)": "Alpha",
    "21A (Delta)": "Delta",
    "21I (Delta)": "Delta",
    "21J (Delta)": "Delta",
    "21K (Omicron)": "Omicron",
    "21L (Omicron)": "Omicron",
    "22B (Omicron)": "Omicron",
}

# Pango lineage prefixes -> WHO variant label. Order matters: longest first.
PANGO_TO_VARIANT = [
    ("B.1.1.7", "Alpha"),
    ("Q.", "Alpha"),          # Q.* are Alpha sublineages
    ("B.1.617.2", "Delta"),
    ("AY.", "Delta"),         # AY.* are Delta sublineages
    ("B.1.1.529", "Omicron"),
    ("BA.", "Omicron"),
]


def load_reference_spike(genome_fasta: str | Path) -> str:
    """Translate the spike CDS out of the Wuhan-Hu-1 genome.

    Returns the 1273-residue protein with the terminal stop removed.
    """
    seq_lines = []
    with open(genome_fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                seq_lines.append(line.strip())
    genome = "".join(seq_lines).upper()

    cds = genome[SPIKE_CDS_START - 1 : SPIKE_CDS_END]
    protein = str(Seq(cds).translate())

    if not protein.endswith("*"):
        raise ValueError("Spike CDS did not end in a stop codon; check coordinates.")
    protein = protein[:-1]

    if "*" in protein:
        raise ValueError("Internal stop codon in translated spike; check coordinates.")
    if len(protein) != EXPECTED_SPIKE_LENGTH:
        raise ValueError(
            f"Spike length {len(protein)} != expected {EXPECTED_SPIKE_LENGTH}"
        )
    return protein


def annotate_position(pos: int) -> str:
    """Return the spike domain containing a 1-based residue position."""
    for name, (start, end) in SPIKE_DOMAINS.items():
        if start <= pos <= end:
            return name
    return "Outside annotated range"


def in_rbd(pos: int) -> bool:
    return RBD[0] <= pos <= RBD[1]


def in_rbm(pos: int) -> bool:
    return RBM[0] <= pos <= RBM[1]


def load_covariants_mutations(clusters_json: str | Path) -> dict[str, list[dict]]:
    """Parse CoVariants clusters.json into {clade: [spike mutations]}.

    Each mutation is {"ref": str, "pos": int, "alt": str}. 'alt' of "-" means
    a deletion. Only gene == "S" entries are kept.
    """
    with open(clusters_json) as handle:
        clusters = json.load(handle)["clusters"]

    out: dict[str, list[dict]] = {}
    for cluster in clusters:
        name = cluster.get("display_name")
        muts = (cluster.get("mutations") or {}).get("nonsynonymous") or []
        spike = [
            {"ref": m["left"], "pos": int(m["pos"]), "alt": m["right"]}
            for m in muts
            if m.get("gene") == "S" and m.get("pos") is not None
        ]
        if name:
            out[name] = sorted(spike, key=lambda m: m["pos"])
    return out


def defining_mutations_by_variant(
    clusters_json: str | Path,
) -> dict[str, dict[str, list[dict]]]:
    """Variant-defining spike mutations for the clades this project analyses."""
    all_clades = load_covariants_mutations(clusters_json)
    grouped: dict[str, dict[str, list[dict]]] = {}
    for clade, variant in CLADE_TO_VARIANT.items():
        if clade in all_clades:
            grouped.setdefault(variant, {})[clade] = all_clades[clade]
    return grouped


def variant_from_pango(lineage: str | None) -> str | None:
    """Map a Pango lineage string to Alpha / Delta / Omicron, else None."""
    if not lineage:
        return None
    text = str(lineage).strip()
    if not text or text.lower() in {"nan", "none", "unclassifiable"}:
        return None
    for prefix, variant in PANGO_TO_VARIANT:
        if text == prefix.rstrip(".") or text.startswith(prefix):
            return variant
    return None
