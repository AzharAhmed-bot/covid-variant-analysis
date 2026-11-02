"""Part 3: Stage 5 core analysis, Stage 6 interpretation."""
import json
from pathlib import Path

CELLS = []
def md(t): CELLS.append({"cell_type": "markdown", "metadata": {}, "source": t.strip("\n")})
def code(t): CELLS.append({"cell_type": "code", "metadata": {}, "source": t.strip("\n"),
                           "execution_count": None, "outputs": []})

md("""
## Stage 5. Core analysis

Two tests are run.

**Per position.** For each of the 1273 reference positions, does the proportion
of isolates carrying a non-reference residue differ between the three
variants? Each position gives a 3 x 2 contingency table, tested with a
chi-square test of independence, or Fisher's exact test where expected counts
are too small for the chi-square approximation to be valid. Testing over a
thousand positions inflates false positives, so p-values are corrected with
Benjamini-Hochberg and reported as q-values.

**Mutation burden.** Does the number of RBD mutations per isolate differ
between variants? These are small integer counts and far from normal, so the
omnibus test is Kruskal-Wallis, followed by pairwise Mann-Whitney U tests with
Benjamini-Hochberg correction across the three comparisons.

Significance and effect size are reported separately throughout. With hundreds
of sequences a difference of a few percent can reach significance while
meaning little, so every q-value is accompanied by the frequency difference it
corresponds to.
""")

code("""
mutations = analysis.mutation_matrix(dataset, REFERENCE_SPIKE)
print(f"mutation matrix: {mutations.shape[0]} isolates x {mutations.shape[1]} positions")

position_results = analysis.per_position_tests(mutations, dataset["variant"],
                                               REFERENCE_SPIKE)
significant = position_results[position_results["q_value"] < 0.05]
print(f"positions with any variation: {len(position_results)}")
print(f"positions significant at q < 0.05: {len(significant)}")
""")

code("""
show_cols = ["position", "ref_aa", "domain", "in_rbd",
             "freq_Alpha", "freq_Delta", "freq_Omicron",
             "max_freq_difference", "test", "q_value"]
show_cols = [c for c in show_cols if c in position_results.columns]
print("Twenty most strongly differentiated positions:")
display(position_results[show_cols].head(20).style.format({
    c: "{:.3f}" for c in show_cols if c.startswith("freq_") or c == "max_freq_difference"
}).format({"q_value": "{:.2e}"}))
""")

code("""
# Volcano-style view: effect size against significance.
fig, ax = plt.subplots(figsize=(9, 5.5))
pr = position_results.copy()
pr["neglog_q"] = -np.log10(pr["q_value"].clip(lower=1e-300))
colors = np.where(pr["in_rbd"], "#C44E52", "#4C72B0")
ax.scatter(pr["max_freq_difference"], pr["neglog_q"], c=colors, s=18, alpha=0.75)
ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=1)
ax.text(0.01, -np.log10(0.05) + 2, "q = 0.05", color="grey", fontsize=8)
for _, row in pr.head(12).iterrows():
    ax.annotate(f"{row['ref_aa']}{int(row['position'])}",
                (row["max_freq_difference"], row["neglog_q"]),
                fontsize=7, xytext=(3, 3), textcoords="offset points")
ax.set(xlabel="largest difference in mutation frequency between variants",
       ylabel="-log10(q)",
       title="Figure 4. Effect size against significance. Red: RBD positions.")
fig.tight_layout(); fig.savefig("figures/fig4_volcano.png", dpi=150, bbox_inches="tight")
plt.show()
""")

code("""
burden = analysis.burden_comparison(dataset, value="n_rbd_mutations")
print("Median RBD mutations per isolate:", burden["medians"])
print(f"\\nKruskal-Wallis: H = {burden['kruskal_H']:.1f}, p = {burden['kruskal_p']:.3g}")
print("\\nPairwise Mann-Whitney with Benjamini-Hochberg correction:")
display(burden["pairwise"][["comparison", "n_a", "n_b", "median_a", "median_b",
                            "rank_biserial", "p_value", "q_value"]])
print(textwrap.fill(
    "Rank-biserial correlation is the effect size matching Mann-Whitney. A value "
    "of -1 means complete separation: every isolate of the second group carries "
    "more RBD mutations than every isolate of the first.", 95))
""")

code("""
burden_total = analysis.burden_comparison(dataset, value="n_mutations")
print("Whole-spike mutation burden, for comparison:")
print("Medians:", burden_total["medians"])
print(f"Kruskal-Wallis: H = {burden_total['kruskal_H']:.1f}, "
      f"p = {burden_total['kruskal_p']:.3g}")
display(burden_total["pairwise"][["comparison", "median_a", "median_b",
                                  "rank_biserial", "q_value"]])
""")

md("""
### What Stage 5 shows

The RBD burden hypothesis is rejected decisively, and the direction matches
the lineage definitions: Omicron carries far more RBD mutations than Alpha or
Delta, with complete separation between the groups. As set out in Stage 1 this
is a check that the pipeline reproduces known biology rather than a discovery,
because the lineage definitions already determine it.

The per-position results are the informative part. They show *where* the
variants differ, including positions that are not variant-defining, and the
effect sizes distinguish positions that are nearly fixed within a variant from
those that merely drift.
""")

md("""
## Stage 6. Biological interpretation

Three forms of interpretation are used: mapping significant positions onto
spike domains, comparing the positions found against the independently curated
variant-defining mutations, and relating the strongest signals to published
findings on spike function.
""")

code("""
domains = analysis.domain_summary(position_results)
display(domains)

fig, ax = plt.subplots(figsize=(9, 5))
d = domains[domains["positions_tested"] > 0].sort_values("significant_per_100aa")
ax.barh(d["domain"], d["significant_per_100aa"], color="#4C72B0")
ax.set(xlabel="significant positions per 100 residues",
       title="Figure 5. Differentiated positions concentrate in the furin site, RBD and NTD")
fig.tight_layout(); fig.savefig("figures/fig5_domains.png", dpi=150, bbox_inches="tight")
plt.show()
""")

code("""
# Do the positions found empirically match the curated variant-defining sets?
for variant in ["Alpha", "Delta", "Omicron"]:
    if variant not in set(dataset["variant"]):
        continue
    observed = analysis.top_mutations_by_variant(dataset, REFERENCE_SPIKE, variant, n=12)
    curated_positions = {m["pos"] for clade in DEFINING.get(variant, {}).values()
                         for m in clade}
    observed["curated_defining"] = observed["position"].isin(curated_positions)
    print(f"\\n=== {variant}: most frequent spike mutations "
          f"(n = {(dataset['variant'] == variant).sum()}) ===")
    display(observed[["mutation", "domain", "in_rbd", "count", "frequency",
                      "curated_defining", "indel_ambiguous"]]
            .style.format({"frequency": "{:.3f}"}))
""")

code("""
# Agreement between what we observe at high frequency and the curated catalogue.
rows = []
for variant in ["Alpha", "Delta", "Omicron"]:
    sub = dataset[dataset["variant"] == variant]
    if sub.empty:
        continue
    curated = {m["pos"] for clade in DEFINING.get(variant, {}).values() for m in clade}
    counts = collections.Counter()
    for seq in sub["spike"]:
        for _, pos, _ in align.mutations_vs_reference(seq, REFERENCE_SPIKE):
            counts[pos] += 1
    frequent = {p for p, c in counts.items() if c / len(sub) >= 0.5}
    rows.append({
        "variant": variant, "n_sequences": len(sub),
        "curated_positions": len(curated),
        "positions_at_50pct_or_more": len(frequent),
        "curated_recovered": len(curated & frequent),
        "recovery_rate": round(len(curated & frequent) / len(curated), 3) if curated else np.nan,
        "frequent_not_curated": len(frequent - curated),
    })
display(pd.DataFrame(rows))
print(textwrap.fill(
    "Positions seen at high frequency but absent from the curated set are mostly "
    "sublineage markers: the catalogue describes the parent clade, while the "
    "isolates sampled belong to descendant sublineages that carry additional "
    "mutations. Curated positions that are not recovered are chiefly deletions in "
    "the NTD loop at 141-145, where protein-level alignment cannot resolve "
    "placement; this is explained in Stage 3 and revisited in Stage 8.", 95))
""")

md("""
### Linking the strongest signals to spike function

The domains carrying the most differentiated positions per residue are the
S1/S2 furin cleavage site, the RBD and the NTD, which is what the literature on
spike function would predict.

- **Receptor-binding domain (319–541)** contains the residues that contact
  ACE2. Omicron's substitutions here include K417N, S477N, T478K, E484A,
  Q493R, Q498R, N501Y and Y505H. N501Y increases ACE2 affinity and is shared
  with Alpha; the cluster around 484–505 is associated with escape from
  neutralising antibodies raised against earlier variants.
- **Furin cleavage site (681–685)** sits at the S1/S2 boundary and must be
  cleaved for membrane fusion. Alpha carries P681H and Delta carries P681R,
  both of which increase cleavage efficiency. That two variants independently
  alter the same five-residue window, by different substitutions, is the
  signature of convergent selection rather than shared ancestry.
- **N-terminal domain (14–305)** carries an antigenic supersite. The
  recurrent deletions here, including 69–70, 141–145 and 211–212, remove
  epitopes and recur across unrelated lineages.
- **HR1 and the connector domain** contribute to the fusion machinery.
  Alpha's S982A and D1118H and Omicron's Q954H and N969K fall here, away from
  the receptor-binding surface entirely, which is part of why variant identity
  is not localised to the RBD.
""")

Path("build/cells_part3.json").write_text(json.dumps(CELLS))
print(f"part 3: {len(CELLS)} cells")
