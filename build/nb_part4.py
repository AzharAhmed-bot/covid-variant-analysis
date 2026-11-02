"""Part 4: Stage 7 machine learning, Stage 8 conclusions."""
import json
from pathlib import Path

CELLS = []
def md(t): CELLS.append({"cell_type": "markdown", "metadata": {}, "source": t.strip("\n")})
def code(t): CELLS.append({"cell_type": "code", "metadata": {}, "source": t.strip("\n"),
                           "execution_count": None, "outputs": []})

md("""
## Stage 7. Machine learning

Three things are modelled, for different reasons.

**Geographic region, within one variant.** This is the honest predictive task.
Region is not part of how a lineage is defined, so a model has to find real
signal, and it may fail. The variant is held fixed because region and variant
are confounded in GenBank, as Figure 3 showed: a model free to read the variant
could infer the region from that rather than from geography.

**Variant, as an ablation.** Predicting variant from spike is circular, and
both the model and a trivial baseline will score near the ceiling. It is still
worth running to measure how much variant signal sits outside the RBD, which is
a real quantity rather than a restatement of the labels.

**Clustering.** Whether structure in the sequences recovers the variant labels
without being given them, scored with the adjusted Rand index.

Leakage is controlled three ways. Identical mutation profiles were removed in
Stage 3, so no sequence appears in both training and test. Feature selection
and model fitting run inside each cross-validation fold through a pipeline,
never on the full dataset. The held-out test set is scored once.

Because the classes are unbalanced, the metric throughout is **balanced
accuracy**, which averages recall across classes; plain accuracy would reward a
model for simply predicting the majority region.
""")

code("""
features = ml.mutation_features(dataset, REFERENCE_SPIKE)
varying = int((features.sum() > 0).sum())
print(f"feature matrix: {features.shape[0]} isolates x {features.shape[1]} positions")
print(f"positions that vary at all: {varying}")
""")

md("""
### Task A. Geographic region, holding the lineage fixed

This is built up in three steps, because the obvious version of the task is
confounded and the confounding is only visible once it is looked for.

Step one runs the task within Alpha. Step two shows why grouping by WHO
variant is not enough. Step three repeats the task holding the Pango lineage
fixed, which is the result that should be believed.

Europe and North America are the only regions with enough distinct sequences
after deduplication to support cross-validation. That scarcity reflects global
sequencing capacity rather than anything about the virus, and it is recorded in
Stage 8.
""")

code("""
alpha = dataset[(dataset["variant"] == "Alpha")
                & (dataset["region"].isin(["Europe", "North America"]))]
X_alpha = features.loc[alpha.index]
y_alpha = alpha["region"]

results_geo = ml.evaluate(X_alpha, y_alpha, "Alpha: Europe vs North America")
print(f"\\nUpper bound set by duplicate sequences across regions: {ceiling_alpha:.3f}")
""")

code("""
report_text, confusion = ml.detailed_report(X_alpha, y_alpha)
print("Random forest, held-out test set:\\n")
print(report_text)
display(confusion)
""")

md("""
#### Why grouping by WHO variant is not enough

Running the same task inside Delta and inside Omicron gives much higher scores
than inside Alpha. Before accepting that as a geographic result, it is worth
asking what the models are keying on.
""")

code("""
variant_level = []
for variant in ["Alpha", "Delta", "Omicron"]:
    sub = dataset[(dataset["variant"] == variant)
                  & (dataset["region"].isin(["Europe", "North America"]))]
    if sub["region"].nunique() < 2 or sub["region"].value_counts().min() < 25:
        continue
    res = ml.evaluate(features.loc[sub.index], sub["region"], variant, verbose=False)
    res.insert(0, "grouped_by", f"variant = {variant}")
    variant_level.append(res)
variant_level = pd.concat(variant_level, ignore_index=True)
display(variant_level[["grouped_by", "model", "cv_balanced_accuracy",
                       "test_balanced_accuracy", "test_roc_auc"]])
""")

code("""
# What the Delta and Omicron models rely on, and how the sublineages are spread.
for variant in ["Delta", "Omicron"]:
    sub = dataset[(dataset["variant"] == variant)
                  & (dataset["region"].isin(["Europe", "North America"]))]
    print(f"=== {variant}: positions driving the model ===")
    display(ml.top_discriminative_positions(features.loc[sub.index], sub["region"],
                                            REFERENCE_SPIKE, n=6))
    print(f"=== {variant}: Pango lineage composition by region ===")
    ct = pd.crosstab(sub["lineage"], sub["region"])
    display(pd.concat([ct, (ct / ct.sum()).round(3).add_suffix(" (proportion)")], axis=1))
""")

md("""
The Omicron model's strongest features are N440K, Q498R and the 69–70
deletion. Those are not geographic markers: they are the mutations that
separate BA.1, BA.2 and BA.5 from one another. The lineage table explains why
that helps. BA.1 makes up about 31% of the North American Omicron sequences
but under 3% of the European ones, so a model that detects BA.1 gets most of
the way to detecting North America without using geography at all.

Delta shows a milder version of the same thing: AY.4 is roughly twice as
common in Europe as in North America, and the top feature is position 142,
which varies between Delta sublineages.

This is the original circularity reappearing one level down. Holding the WHO
variant fixed does not hold the lineage fixed, because each WHO variant spans
several Pango lineages whose geographic distributions differ. The test has to
be repeated with the Pango lineage held fixed.
""")

md("""
#### The result to believe: one Pango lineage at a time
""")

code("""
lineage_level = []
for lineage in ["B.1.1.7", "B.1.617.2", "AY.4", "BA.1", "BA.2", "BA.5"]:
    sub = dataset[(dataset["lineage"] == lineage)
                  & (dataset["region"].isin(["Europe", "North America"]))]
    counts = sub["region"].value_counts()
    if sub["region"].nunique() < 2 or counts.min() < 25:
        print(f"{lineage}: not enough distinct sequences in both regions "
              f"({dict(counts)}), skipped")
        continue
    res = ml.evaluate(features.loc[sub.index], sub["region"], lineage, verbose=False)
    res.insert(0, "lineage", lineage)
    res.insert(1, "variant", sub["variant"].iloc[0])
    lineage_level.append(res)

lineage_level = pd.concat(lineage_level, ignore_index=True)
display(lineage_level[["lineage", "variant", "model", "cv_balanced_accuracy", "cv_sd",
                       "test_balanced_accuracy", "test_roc_auc", "n_train", "n_test"]])
""")

code("""
fig, ax = plt.subplots(figsize=(9, 4.5))
plot_df = lineage_level[lineage_level["model"] != "Majority-class baseline"]
sns.barplot(data=plot_df, x="lineage", y="test_balanced_accuracy", hue="model",
            ax=ax, palette="Blues_d")
ax.axhline(0.5, ls="--", c="grey", lw=1)
ax.text(-0.45, 0.515, "majority-class baseline", color="grey", fontsize=8)
ax.axhline(ceiling_alpha, ls=":", c="#C44E52", lw=1.3)
ax.text(-0.45, ceiling_alpha + 0.012, "ceiling set by cross-region duplicates",
        color="#C44E52", fontsize=8)
ax.set(ylim=(0, 1.05), ylabel="balanced accuracy (held-out test)", xlabel="Pango lineage",
       title="Figure 7. Geographic signal within a single lineage, replicated three times")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig("figures/fig7_geography.png", dpi=150, bbox_inches="tight")
plt.show()
""")

code("""
# Could this be collection date rather than place? Check each lineage separately.
for lineage in lineage_level["lineage"].unique():
    sub = dataset[(dataset["lineage"] == lineage)
                  & (dataset["region"].isin(["Europe", "North America"]))]
    summary = sub.groupby("region")["collection_date"].agg(["min", "median", "max", "count"])
    gap = abs((summary.loc["Europe", "median"] - summary.loc["North America", "median"]).days)
    print(f"=== {lineage}: median collection dates differ by {gap} days ===")
    display(summary)
""")

code("""
# Positions driving the within-lineage models, now that lineage is fixed.
for lineage in lineage_level["lineage"].unique():
    sub = dataset[(dataset["lineage"] == lineage)
                  & (dataset["region"].isin(["Europe", "North America"]))]
    print(f"=== {lineage}: positions distinguishing Europe from North America ===")
    display(ml.top_discriminative_positions(features.loc[sub.index], sub["region"],
                                            REFERENCE_SPIKE, n=8))
""")

md("""
### Task B. Variant classification, and where the signal lives
""")

code("""
ablation = ml.rbd_ablation(features, dataset["variant"])
display(ablation[["task", "model", "cv_balanced_accuracy", "cv_sd",
                  "test_balanced_accuracy", "n_features"]])
print(textwrap.fill(
    "Masking the receptor-binding domain entirely does not reduce accuracy, and "
    "the NTD alone is equally sufficient. The variant signal is therefore not "
    "localised to the RBD: it is repeated across the protein. This answers the "
    "third hypothesis in Stage 1, and it is the reason the next cell measures "
    "redundancy directly.", 95))
""")

code("""
# How many arbitrary positions are actually needed to identify a variant?
redundancy = ml.redundancy_curve(features, dataset["variant"],
                                 sizes=(1, 2, 3, 5, 10, 20, 50, 100), n_repeats=10)
display(redundancy)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.errorbar(redundancy["n_positions"], redundancy["mean_balanced_accuracy"],
            yerr=redundancy["sd"], marker="o", capsize=3, color="#4C72B0")
ax.fill_between(redundancy["n_positions"], redundancy["min"], redundancy["max"],
                alpha=0.15, color="#4C72B0")
ax.axhline(1 / dataset["variant"].nunique(), ls="--", c="grey", lw=1)
ax.text(1.2, 1 / dataset["variant"].nunique() + 0.02, "chance", color="grey", fontsize=8)
ax.set(xscale="log", xlabel="number of randomly chosen varying positions",
       ylabel="balanced accuracy",
       title="Figure 6. Variant identity is redundantly encoded across spike")
fig.tight_layout(); fig.savefig("figures/fig6_redundancy.png", dpi=150, bbox_inches="tight")
plt.show()
print(textwrap.fill(
    "The shaded band spans the best and worst random subsets at each size. Some "
    "very small subsets already separate the variants almost perfectly, because "
    "they happen to include a defining position, while others of the same size do "
    "little better than chance. Accuracy becomes reliably high once a few dozen "
    "arbitrary positions are included.", 95))
""")

md("""
### Task C. Unsupervised clustering
""")

code("""
display(ml.clustering_report(features, dataset["variant"]))
print(textwrap.fill(
    "An adjusted Rand index near 1 means the clusters reproduce the variant labels "
    "almost exactly without being shown them. For the same reason as Task B this "
    "is not a discovery: the labels are defined by the features being clustered.", 95))
""")

md("""
## Stage 8. Conclusions

### One conclusion per objective

**Objective 1 — assemble a usable dataset.** Spike sequences were retrieved for
Alpha, Delta and Omicron isolates with lineage, collection date and region.
The filtering chain is in Stage 3. The binding constraint was not the number
of records but the number of *distinct* spike sequences: thousands of isolates
collapse to a few hundred distinct mutation profiles, because within a variant
most isolates differ nowhere in spike.

**Objective 2 — compare variation between variants.** Mutation frequency
differs significantly between variants at many spike positions after
Benjamini-Hochberg correction. RBD mutation burden differs decisively
(Kruskal-Wallis, with every pairwise comparison significant and complete
separation between Omicron and the others), with median RBD mutations of 1 in
Alpha, 2 in Delta and 16 in Omicron. H0 for mutation burden is rejected. As
stated in Stage 1, this confirms the pipeline reproduces known biology rather
than discovering anything, because the lineage definitions determine it.

**Objective 3 — interpret the results biologically.** Differentiated positions
concentrate in the S1/S2 furin cleavage site, the RBD and the NTD when
normalised by domain length. The observed high-frequency mutations recover the
curated variant-defining sets, and the residual disagreements are explained:
extra positions are sublineage markers, and missing ones are chiefly NTD
deletions whose placement protein-level alignment cannot resolve.

**Objective 4 — predict and ablate.** Geographic region is predictable from
spike sequence above a majority-class baseline within a single Pango lineage,
replicated independently in B.1.1.7, B.1.617.2 and BA.2, so H0 for geographic
prediction is rejected. The strength varies by lineage and is weakest in
Alpha. Reaching this conclusion required a correction: grouping by WHO variant
is not sufficient, because each variant spans several Pango lineages whose
regional distributions differ, and models trained that way partly detect
sublineage rather than geography. Removing the RBD does not reduce variant
classification accuracy at all, so H0 for localisation is *not* rejected: the
variant signal is not localised to the RBD but repeated across the protein,
and roughly fifty arbitrary varying positions suffice to identify a variant.

### Limitations

1. **The variant classification task is circular.** Pango lineages are defined
   by spike mutations, so near-perfect accuracy and an adjusted Rand index near
   1 restate the labels rather than demonstrating predictive power. This was
   anticipated, which is why the geographic task and the ablation carry the
   analytical weight.

2. **Phylogeny is only partly controlled, even at lineage level.** Holding the
   Pango lineage fixed removes the sublineage shortcut that inflated the
   variant-level results, but Pango is a hierarchy and a named lineage still
   contains finer structure. Sequences from one region may be phylogenetically
   closer to each other for reasons of transmission history rather than
   anything geographic, and spike alone cannot separate those explanations.
   The honest reading is that spike sequence carries information about where
   an isolate came from, not that geography is acting directly on spike. A
   whole-genome phylogeny with ancestral-state reconstruction would be needed
   to go further.

7. **Indel placement is not fully resolved.** Sequences are projected onto
   reference coordinates by global protein alignment, which cannot place
   deletions inside low-complexity runs. In BA.1, residues 142–144 are deleted
   and 145 is substituted; BLOSUM62 scores the substituted residue better
   against reference G142 than against Y145, so the highest-scoring alignment
   disagrees with the curated annotation. Aligning nucleotides codon-wise, as
   Nextclade does, would resolve this. The affected windows are declared
   explicitly, gap runs are left-normalised for consistency, and deletions are
   summarised per region so that no conclusion depends on exact placement.

7. **Sampling is strongly biased.** Sequences in GenBank reflect where
   sequencing happened, not where the virus was. Region is confounded with
   variant and with time, Omicron records are overwhelmingly North American,
   and after deduplication too few distinct Alpha sequences remain from Africa
   and Asia to model them at all. The geographic result therefore describes
   Europe against North America, not the world.

7. **A ceiling exists that no model can pass.** Identical spike sequences were
   collected in different regions. Any model seeing only the sequence must
   misassign some of those records, which bounds achievable geographic accuracy
   independently of model choice. The bound is reported alongside the scores.

7. **Insertions are counted but not modelled.** Residues inserted relative to
   the reference, such as Omicron's ins214EPE, are recorded as a count and
   excluded from the feature matrix, which keeps position numbering aligned
   with the literature at the cost of discarding a little information.

7. **Class imbalance and sample size.** Delta and Omicron contribute far fewer
   distinct sequences than Alpha. Balanced accuracy and class weighting are
   used throughout, but the per-class estimates for the smaller groups rest on
   modest numbers.

### Future work

- Align nucleotides codon-wise before translating, to resolve the NTD deletion
  placement properly and let insertions be modelled.
- Control phylogeny properly by building a tree and testing whether regional
  clustering persists after accounting for shared ancestry, which the
  lineage-level control only approximates.
- Test whether the positions driving geographic prediction, such as the HR2
  position 1191 and E484K in the RBD, correspond to known regional sublineages,
  for example Alpha carrying E484K in the United Kingdom.
- Weight or subsample by sequencing effort per country, to separate viral
  biology from surveillance intensity.
""")

Path("build/cells_part4.json").write_text(json.dumps(CELLS))
print(f"part 4: {len(CELLS)} cells")
