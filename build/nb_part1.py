"""Builds the project notebook. Part 1: cell helpers and Stages 0-2."""
import json
from pathlib import Path

CELLS = []

def md(text):
    CELLS.append({"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")})

def code(text):
    CELLS.append({"cell_type": "code", "metadata": {}, "source": text.strip("\n"),
                  "execution_count": None, "outputs": []})

def module_cell(name):
    """Embed a project module so the notebook is self-contained in Colab."""
    body = Path(f"src/{name}").read_text()
    code(f"%%writefile {name}\n{body}")

md("""
# SARS-CoV-2 Spike Protein Variation Across Alpha, Delta and Omicron

**Computing Bioinformatics — Group Project**

This notebook compares amino-acid variation in the SARS-CoV-2 spike protein
between the Alpha, Delta and Omicron variants, asks how much of the
variant-distinguishing signal lies outside the receptor-binding domain, and
tests whether an isolate's geographic region can be predicted from its spike
sequence.

It runs top to bottom in Google Colab. Set `DOWNLOAD_FRESH = True` in Stage 2
to rebuild the dataset from NCBI instead of using the copy stored with the
project.
""")

md("""
## Stage 0. Setup

The project code is written to files here so the notebook is self-contained:
it does not depend on a checkout being present. Each module is also kept in
`src/` in the repository.
""")

code("""
import sys, subprocess
def pip_install(*pkgs):
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", *pkgs], check=False)

# Colab already has most of these; this is a no-op when they are present.
pip_install("biopython", "scikit-learn", "statsmodels", "seaborn")

import warnings
warnings.filterwarnings("ignore")
print("Python", sys.version.split()[0])
""")

for mod in ["spike_reference.py", "align.py", "qc.py", "analysis.py", "ml.py"]:
    module_cell(mod)

code("""
import collections, json, textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, ".")
import spike_reference as sref
import align, qc, analysis, ml

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 50)
sns.set_theme(style="whitegrid", context="notebook")
PALETTE = {"Alpha": "#4C72B0", "Delta": "#DD8452", "Omicron": "#55A868"}
Path("figures").mkdir(exist_ok=True)
print("modules loaded")
""")

md("""
## Stage 1. Question and objectives

### Research question

How does amino-acid variation in the SARS-CoV-2 spike protein differ between
the Alpha, Delta and Omicron variants relative to the Wuhan-Hu-1 reference,
how much of the variant-distinguishing signal lies outside the
receptor-binding domain (RBD), and can an isolate's geographic region be
predicted from its spike sequence alone?

### Why the question is framed this way

An earlier version of this project asked only whether the variant could be
predicted from the spike sequence, and whether Omicron carries more RBD
mutations than Alpha and Delta. Both parts are answerable without collecting
any data, because Pango lineages are *defined* by their spike mutation
profiles. Counting the defining mutations in the CoVariants catalogue gives:

| Clade | Variant | Defining spike mutations | Of which in the RBD |
|---|---|---|---|
| 20I | Alpha | 10 | 1 |
| 21A / 21I / 21J | Delta | 10–11 | 2 |
| 21K (BA.1) | Omicron | 36 | 15 |
| 21L (BA.2) | Omicron | 31 | 16 |
| 22B (BA.5) | Omicron | 34 | 17 |

A classifier asked to predict variant from spike is therefore being asked to
rediscover its own labelling rule, and a hypothesis that Omicron has more RBD
mutations is arithmetic from the table above rather than an empirical finding.

The question was restructured around parts whose answers are *not* fixed in
advance, while keeping the variant comparison as descriptive work.

### Objectives

1. Retrieve spike protein sequences for Alpha, Delta and Omicron isolates from
   NCBI, with Pango lineage, collection date and geographic region, and filter
   them to a set of distinct, good-quality sequences.
2. Compare the frequency of non-reference residues at each of the 1273 spike
   positions between the three variants, using contingency tests with
   Benjamini-Hochberg correction, and compare RBD mutation burden per isolate
   using Kruskal-Wallis with pairwise follow-up tests.
3. Map the significant positions onto spike domains (NTD, RBD, receptor-binding
   motif, furin cleavage site, HR1, HR2) and compare them with the curated
   variant-defining mutations and the published literature.
4. Train classifiers to predict geographic region from spike sequence within a
   single variant, against a majority-class baseline, and measure how much
   variant signal survives when the RBD is removed.

### Hypotheses

**Mutation burden.** H0: the number of RBD mutations per isolate does not
differ between Alpha, Delta and Omicron. H1: it does differ. This is tested
for completeness, but the answer is effectively determined by the lineage
definitions, so it is reported as a check that the pipeline reproduces known
biology rather than as a discovery.

**Geographic prediction.** H0: spike sequence carries no information about the
geographic region an isolate came from, so a classifier performs no better
than a majority-class baseline. H1: it performs better. This hypothesis is
genuinely open, and a failure to reject H0 would be a valid result.

**Localisation of variant signal.** H0: removing the RBD does not reduce the
accuracy of variant classification, meaning the signal is not localised to the
RBD. This is tested by ablation.
""")

md("""
## Stage 2. Data acquisition

### Sources

| Input | Source | Why |
|---|---|---|
| Wuhan-Hu-1 genome (GenBank MN908947) | Nextstrain `ncov`, `defaults/reference_seq.fasta` | Reference for all coordinates |
| Variant-defining mutations | CoVariants (`hodcroftlab/covariants`), `clusters.json` | Independent annotation for Stage 6 |
| Isolate metadata | NCBI Datasets v2 virus API, `dataset_report` | Pango lineage, collection date, geographic region |
| Spike protein sequences | NCBI E-utilities `efetch`, `rettype=fasta_cds_aa` | Annotated CDS translations |

Two notes on how the sequences were obtained, because the obvious routes do
not work. The NCBI Datasets download endpoint ignores
`include_annotation_type` in its per-accession form and returns only a
metadata report, while its taxon form offers no way to restrict to chosen
accessions and returns roughly 641 MB for a single lineage. E-utilities
`efetch` with `rettype=fasta_cds_aa` returns the annotated CDS translations
for specific accessions, from which the record labelled
`[protein=surface glycoprotein]` is the spike protein.

The reference spike is **translated from the genome** at CDS 21563..25384
rather than downloaded as a separate protein record, so the sequence and its
numbering cannot disagree.
""")

code("""
REFERENCE_DIR = Path("data/reference")
RAW_DIR = Path("data/raw")
REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)

REFERENCE_FILES = {
    "wuhan_hu_1_genome.fasta":
        "https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/reference_seq.fasta",
    "covariants_clusters.json":
        "https://raw.githubusercontent.com/hodcroftlab/covariants/master/web/public/data/clusters.json",
}

import urllib.request
for name, url in REFERENCE_FILES.items():
    target = REFERENCE_DIR / name
    if not target.exists():
        print(f"downloading {name}")
        urllib.request.urlretrieve(url, target)
    print(f"{name}: {target.stat().st_size:,} bytes")
""")

code("""
REFERENCE_SPIKE = sref.load_reference_spike(REFERENCE_DIR / "wuhan_hu_1_genome.fasta")
print(f"reference spike: {len(REFERENCE_SPIKE)} aa")
print(f"  N-terminus: {REFERENCE_SPIKE[:24]}")
print(f"  C-terminus: {REFERENCE_SPIKE[-24:]}")

# The numbering must match the literature, so check residues whose reference
# state is well known. If any of these is wrong, every position below is wrong.
KNOWN = {339: "G", 417: "K", 452: "L", 478: "T", 484: "E", 501: "N",
         570: "A", 614: "D", 681: "P", 716: "T", 982: "S", 1118: "D"}
bad = {p: (aa, REFERENCE_SPIKE[p - 1]) for p, aa in KNOWN.items()
       if REFERENCE_SPIKE[p - 1] != aa}
print(f"\\nknown reference residues checked: {len(KNOWN)}, mismatches: {len(bad)}")
assert not bad, bad
""")

code("""
DEFINING = sref.defining_mutations_by_variant(REFERENCE_DIR / "covariants_clusters.json")
rows = []
for variant, clades in DEFINING.items():
    for clade, muts in clades.items():
        rows.append({
            "clade": clade, "variant": variant,
            "spike_mutations": len(muts),
            "rbd_mutations": sum(1 for m in muts if sref.in_rbd(m["pos"])),
            # Does our translated reference agree with CoVariants at every
            # position this clade annotates?
            "reference_mismatches": sum(
                1 for m in muts
                if len(m["ref"]) == 1 and m["pos"] <= len(REFERENCE_SPIKE)
                and REFERENCE_SPIKE[m["pos"] - 1] != m["ref"]),
        })
defining_table = pd.DataFrame(rows).sort_values(["variant", "clade"])
display(defining_table)
assert defining_table["reference_mismatches"].sum() == 0
print("Every reference residue annotated by CoVariants for these clades matches "
      "the translated reference.")
""")

Path("build/cells_part1.json").write_text(json.dumps(CELLS))
print(f"part 1: {len(CELLS)} cells")
