"""Part 2: dataset construction (Stage 2 cont.), Stage 3 QC, Stage 4 EDA."""
import json
from pathlib import Path

CELLS = []
def md(t): CELLS.append({"cell_type": "markdown", "metadata": {}, "source": t.strip("\n")})
def code(t): CELLS.append({"cell_type": "code", "metadata": {}, "source": t.strip("\n"),
                           "execution_count": None, "outputs": []})
def module_cell(n):
    code(f"%%writefile {n}\n" + Path(f"src/{n}").read_text())

module_cell("build_dataset.py")
module_cell("topup.py")

md("""
### Building the isolate dataset

`DOWNLOAD_FRESH` controls where the data comes from. With the default `False`,
the notebook uses the dataset stored with the project, which is what produced
the results reported here. Set it to `True` to rebuild from NCBI; that takes
roughly an hour because the sample is deepened to find enough *distinct* spike
sequences, for the reason explained in Stage 3.

Sampling is stratified by geographic region. Records come back clustered by
submitting laboratory, so the first pages of any lineage are dominated by a
single country, and an unstratified sample would be almost entirely one region.
""")

code("""
DOWNLOAD_FRESH = False
TARGET_PER_VARIANT = 500

metadata_path = RAW_DIR / "metadata.csv"
fasta_path = RAW_DIR / "spike_proteins.fasta"

if DOWNLOAD_FRESH or not (metadata_path.exists() and fasta_path.exists()):
    print("Building the dataset from NCBI. This takes a while.")
    subprocess.run([sys.executable, "-I", "build_dataset.py", str(TARGET_PER_VARIANT)],
                   check=True)
    # Deepen Alpha and Delta: the geographic task is limited by the number of
    # distinct spike sequences, not by the number of records.
    subprocess.run([sys.executable, "-I", "topup.py", "2500", "120",
                    "B.1.1.7", "B.1.617.2", "AY.4"], check=False)
else:
    print("Using the dataset stored with the project.")

print(f"metadata: {metadata_path.stat().st_size:,} bytes")
print(f"sequences: {fasta_path.stat().st_size:,} bytes")
""")

md("""
### What every variable means

| Column | Type | Meaning |
|---|---|---|
| `accession` | string | NCBI nucleotide accession, unique per isolate |
| `lineage` | string | Pango lineage as assigned by NCBI |
| `variant` | string | WHO label (Alpha, Delta, Omicron) derived from the lineage |
| `country` | string | Country of collection |
| `region` | string | Continent-level region, used as the geographic label |
| `collection_date` | date | When the sample was taken, not when it was deposited |
| `genome_length` | integer | Length of the deposited genome in bases |
| `completeness` | string | NCBI's `COMPLETE` / `PARTIAL` flag for the genome |
| `raw_spike` | string | Spike protein exactly as deposited, variable length |
| `spike` | string | Spike projected onto the 1273 reference positions |
| `n_inserted` | integer | Residues inserted relative to the reference |
| `n_ambiguous` | integer | Unresolved (X) positions after projection |
| `mutation_profile` | string | Semicolon-separated mutations, the deduplication key |
| `n_mutations` | integer | Non-reference residues, excluding X |
| `n_rbd_mutations` | integer | Of those, how many fall in the RBD (319–541) |
""")

code("""
raw_metadata = pd.read_csv(metadata_path)
print(f"{len(raw_metadata):,} records downloaded")
display(raw_metadata.head(3))
print("\\nRecords by variant and region, before any filtering:")
display(raw_metadata.groupby(["variant", "region"]).size().unstack(fill_value=0))
""")

md("""
## Stage 3. Preprocessing and quality control

Four decisions are made here, each for a stated reason.

**Length.** Records whose deposited spike is shorter than 1200 residues are
dropped. These are truncated genes rather than short proteins, and their
absent positions would otherwise be scored as deletions.

**Genome completeness is deliberately *not* used as a filter.** Many
full-length Alpha genomes of about 29,758 bases are flagged `PARTIAL` because
of ambiguous bases at the ends, while the spike CDS is complete and annotated.
Requiring `COMPLETE` would silently remove most Alpha isolates and bias the
comparison. Quality is judged on the spike protein itself.

**Ambiguity.** Positions NCBI could not resolve arrive as `X`. Records with
more than 2% ambiguous positions are dropped, and the `X` positions that
remain are excluded from mutation calls rather than counted as changes.

**Deduplication, on the mutation profile rather than the sequence string.**
Many isolates of one variant differ nowhere in spike, and leaving duplicates
in would place copies of the same sequence in both the training and test folds.
The key matters: because `X` positions are excluded from mutation calls, two
records differing only in *where* their unresolved positions fall have
identical biology but different strings. Keying on the string counts them as
distinct and manufactures diversity that grows as the ambiguity threshold is
loosened. Keying on the mutation profile is invariant to that. On this dataset
the string key gave 512 "unique" sequences where the profile key gives 401, so
the artefact was real and material.
""")

code("""
dataset, report = qc.build_dataset(metadata_path, fasta_path, REFERENCE_SPIKE,
                                   verbose=False)
print("Records remaining after each step:\\n")
display(report.to_frame())
print("\\nWhat the duplicates imply:")
for note in report.notes:
    print(textwrap.fill(f"- {note}", 95, subsequent_indent="  "))
""")

code("""
print(f"{len(dataset):,} distinct spike sequences enter the analysis\\n")
display(dataset.groupby(["variant", "region"]).size().unstack(fill_value=0))

print("Variant totals:", dict(dataset["variant"].value_counts()))
print("Collection dates:", dataset["collection_date"].min().date(),
      "to", dataset["collection_date"].max().date())
""")

code("""
# The ceiling that duplicate sequences place on geographic prediction. Where one
# spike sequence occurs in several regions, no model reading only the sequence
# can place all of those records correctly. `report.pre_dedup` is the filtered
# table before duplicates were collapsed, so no re-projection is needed.
pre_dedup = report.pre_dedup
ceiling_all = qc.duplicate_region_ceiling(pre_dedup)

alpha_pre = pre_dedup[(pre_dedup["variant"] == "Alpha")
                      & (pre_dedup["region"].isin(["Europe", "North America"]))]
ceiling_alpha = qc.duplicate_region_ceiling(alpha_pre)

print(f"Records before deduplication: {len(pre_dedup):,}")
print(f"Best achievable region accuracy from sequence alone, all variants: {ceiling_all:.3f}")
print(f"Same, restricted to Alpha in Europe and North America:            {ceiling_alpha:.3f}")
print()
print("This is an upper bound, not a model score: identical sequences collected in")
print("different regions cannot all be assigned correctly by any model that sees")
print("only the sequence.")
""")

md("""
## Stage 4. Exploratory data analysis
""")

code("""
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
order = ["Alpha", "Delta", "Omicron"]

sns.boxplot(data=dataset, x="variant", y="n_mutations", order=order,
            palette=PALETTE, ax=axes[0], showfliers=False)
sns.stripplot(data=dataset, x="variant", y="n_mutations", order=order,
              color="0.25", size=2, alpha=0.4, ax=axes[0])
axes[0].set(title="Spike mutations per isolate", xlabel="", ylabel="mutations vs Wuhan-Hu-1")

sns.boxplot(data=dataset, x="variant", y="n_rbd_mutations", order=order,
            palette=PALETTE, ax=axes[1], showfliers=False)
sns.stripplot(data=dataset, x="variant", y="n_rbd_mutations", order=order,
              color="0.25", size=2, alpha=0.4, ax=axes[1])
axes[1].set(title="Receptor-binding domain mutations per isolate", xlabel="",
            ylabel="RBD mutations (positions 319-541)")
fig.suptitle("Figure 1. Mutation burden separates Omicron from Alpha and Delta", y=1.02)
fig.tight_layout(); fig.savefig("figures/fig1_burden.png", dpi=150, bbox_inches="tight")
plt.show()

display(dataset.groupby("variant")[["n_mutations", "n_rbd_mutations"]]
        .agg(["median", "mean", "min", "max"]).round(2))
""")

code("""
# Per-position mutation frequency along the protein, by variant.
freq_rows = []
for variant in order:
    sub = dataset[dataset["variant"] == variant]
    if sub.empty:
        continue
    arr = np.zeros(len(REFERENCE_SPIKE))
    for seq in sub["spike"]:
        for i, (c, r) in enumerate(zip(seq, REFERENCE_SPIKE)):
            if c != r and c != "X":
                arr[i] += 1
    freq_rows.append(pd.Series(arr / len(sub), name=variant))
freq_df = pd.concat(freq_rows, axis=1)
freq_df.index = np.arange(1, len(REFERENCE_SPIKE) + 1)

fig, axes = plt.subplots(len(freq_df.columns), 1, figsize=(13, 7), sharex=True)
for ax, variant in zip(np.atleast_1d(axes), freq_df.columns):
    ax.fill_between(freq_df.index, freq_df[variant], color=PALETTE[variant], step="mid")
    ax.set_ylim(0, 1.05); ax.set_ylabel(variant)
    ax.axvspan(*sref.RBD, color="red", alpha=0.08)
    ax.axvspan(*sref.SPIKE_DOMAINS["NTD"], color="blue", alpha=0.06)
np.atleast_1d(axes)[-1].set_xlabel("spike residue position")
np.atleast_1d(axes)[0].set_title(
    "Figure 2. Mutation frequency per position. Blue shading: NTD. Red shading: RBD.")
fig.tight_layout(); fig.savefig("figures/fig2_positions.png", dpi=150, bbox_inches="tight")
plt.show()
""")

code("""
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

counts = dataset.groupby(["variant", "region"]).size().unstack(fill_value=0)
counts.plot(kind="bar", stacked=True, ax=axes[0], colormap="tab20", width=0.7)
axes[0].set(title="Sequences by variant and region", xlabel="", ylabel="distinct sequences")
axes[0].tick_params(axis="x", rotation=0)
axes[0].legend(fontsize=7, title="region")

timeline = (dataset.assign(month=dataset["collection_date"].dt.to_period("M").dt.to_timestamp())
            .groupby(["month", "variant"]).size().unstack(fill_value=0))
timeline.plot(ax=axes[1], color=[PALETTE.get(c, "grey") for c in timeline.columns], lw=2)
axes[1].set(title="Collection dates by variant", xlabel="", ylabel="sequences")

fig.suptitle("Figure 3. Sampling is uneven across regions and variants", y=1.02)
fig.tight_layout(); fig.savefig("figures/fig3_sampling.png", dpi=150, bbox_inches="tight")
plt.show()

print("Region and variant are confounded in GenBank. Omicron records are")
print("overwhelmingly North American and Delta records are concentrated in Europe")
print("and North America, so a model free to use the variant could infer the region")
print("from it. This is why the geographic task in Stage 7 holds the variant fixed.")
""")

Path("build/cells_part2.json").write_text(json.dumps(CELLS))
print(f"part 2: {len(CELLS)} cells")
