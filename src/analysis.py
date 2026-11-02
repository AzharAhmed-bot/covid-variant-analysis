"""Stage 5: statistical comparison of spike variation between variants.

Two questions are tested.

Per position. For each of the 1273 reference positions, does the proportion of
isolates carrying a non-reference residue differ between Alpha, Delta and
Omicron? Each position gives a 3 x 2 contingency table, tested with a
chi-square test of independence, or Fisher's exact test where expected counts
are too small for the chi-square approximation to hold. Testing 1273 positions
inflates false positives, so p-values are corrected with Benjamini-Hochberg and
reported as q-values.

Mutation burden. Does the number of RBD mutations per isolate differ between
variants? Counts are small integers and far from normal, so the omnibus test is
Kruskal-Wallis, followed by pairwise Mann-Whitney U tests with
Benjamini-Hochberg correction across the three comparisons.

Statistical significance is separated from effect size throughout. With
hundreds of sequences a difference of a few percent reaches significance while
meaning little biologically, so the frequency difference is reported alongside
every q-value.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import align
import spike_reference as ref

AMBIGUOUS_CHAR = "X"


def mutation_matrix(df: pd.DataFrame, reference: str,
                    column: str = "spike") -> pd.DataFrame:
    """Boolean matrix of isolates x positions, True where the residue differs.

    Unresolved (X) positions are recorded as missing rather than as mutations,
    so they are excluded from each position's test instead of being counted as
    changes.
    """
    n, length = len(df), len(reference)
    data = np.zeros((n, length), dtype=object)
    for i, seq in enumerate(df[column].to_numpy()):
        row = np.fromiter((c != r for c, r in zip(seq, reference)),
                          dtype=bool, count=length)
        amb = np.fromiter((c == AMBIGUOUS_CHAR for c in seq),
                          dtype=bool, count=length)
        row_obj = row.astype(object)
        row_obj[amb] = None
        data[i] = row_obj
    return pd.DataFrame(data, index=df.index,
                        columns=pd.RangeIndex(1, length + 1, name="position"))


def per_position_tests(mut: pd.DataFrame, groups: pd.Series,
                       reference: str, min_count: int = 1) -> pd.DataFrame:
    """Test each position for a difference in mutation frequency between groups."""
    labels = sorted(groups.dropna().unique())
    rows = []
    for pos in mut.columns:
        col = mut[pos]
        observed, freqs = [], {}
        for label in labels:
            vals = col[groups == label].dropna()
            mutated = int(vals.sum())
            total = int(len(vals))
            observed.append([mutated, total - mutated])
            freqs[label] = mutated / total if total else np.nan
        table = np.array(observed, dtype=float)
        total_mutated = table[:, 0].sum()
        if total_mutated < min_count or table.sum() == 0:
            continue
        # Chi-square needs adequate expected counts; otherwise use Fisher.
        try:
            chi2, p, _, expected = stats.chi2_contingency(table)
            test = "chi-square"
            if (expected < 5).any():
                raise ValueError
        except Exception:
            try:
                if table.shape[0] == 2:
                    _, p = stats.fisher_exact(table)
                else:
                    p = stats.fisher_exact(table)[1] if table.shape == (2, 2) else \
                        stats.chi2_contingency(table + 0.5)[1]
                test = "Fisher / corrected chi-square"
            except Exception:
                continue
        values = [v for v in freqs.values() if not np.isnan(v)]
        rows.append({
            "position": int(pos),
            "ref_aa": reference[int(pos) - 1],
            "domain": ref.annotate_position(int(pos)),
            "in_rbd": ref.in_rbd(int(pos)),
            "in_rbm": ref.in_rbm(int(pos)),
            **{f"freq_{label}": freqs[label] for label in labels},
            "max_freq_difference": (max(values) - min(values)) if values else np.nan,
            "n_mutated_total": int(total_mutated),
            "test": test,
            "p_value": p,
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_value"] = multipletests(out["p_value"], method="fdr_bh")[1]
    return out.sort_values("q_value").reset_index(drop=True)


def burden_comparison(df: pd.DataFrame, value: str = "n_rbd_mutations",
                      group: str = "variant") -> dict:
    """Kruskal-Wallis omnibus plus BH-corrected pairwise Mann-Whitney tests."""
    labels = sorted(df[group].dropna().unique())
    samples = [df.loc[df[group] == g, value].to_numpy() for g in labels]
    H, p_omnibus = stats.kruskal(*samples)

    pairs, raw_p, effects = [], [], []
    for a, b in itertools.combinations(labels, 2):
        x = df.loc[df[group] == a, value].to_numpy()
        y = df.loc[df[group] == b, value].to_numpy()
        U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        # Rank-biserial correlation: the effect size matching Mann-Whitney.
        rbc = 2 * U / (len(x) * len(y)) - 1
        pairs.append(f"{a} vs {b}")
        raw_p.append(p)
        effects.append({"median_a": float(np.median(x)), "median_b": float(np.median(y)),
                        "rank_biserial": float(rbc), "n_a": len(x), "n_b": len(y)})
    q = multipletests(raw_p, method="fdr_bh")[1] if raw_p else []
    pairwise = pd.DataFrame([
        {"comparison": c, "p_value": p, "q_value": qq, **e}
        for c, p, qq, e in zip(pairs, raw_p, q, effects)
    ])
    return {
        "groups": labels,
        "medians": {g: float(np.median(s)) for g, s in zip(labels, samples)},
        "kruskal_H": float(H),
        "kruskal_p": float(p_omnibus),
        "pairwise": pairwise,
    }


def domain_summary(results: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Significant positions per spike domain, normalised by domain length."""
    sig = results[results["q_value"] < alpha]
    rows = []
    for domain, (start, end) in ref.SPIKE_DOMAINS.items():
        length = end - start + 1
        tested = int((results["domain"] == domain).sum())
        n_sig = int((sig["domain"] == domain).sum())
        rows.append({
            "domain": domain, "start": start, "end": end, "length_aa": length,
            "positions_tested": tested, "significant": n_sig,
            "significant_per_100aa": round(100 * n_sig / length, 2),
        })
    return pd.DataFrame(rows).sort_values("significant_per_100aa", ascending=False)


def top_mutations_by_variant(df: pd.DataFrame, reference: str,
                             variant: str, n: int = 15) -> pd.DataFrame:
    """Most frequent mutations within one variant, with domain annotation."""
    sub = df[df["variant"] == variant]
    counts: dict[tuple, int] = {}
    for seq in sub["spike"]:
        for m in align.mutations_vs_reference(seq, reference):
            counts[m] = counts.get(m, 0) + 1
    rows = [{
        "mutation": align.format_mutation(*m),
        "position": m[1],
        "domain": ref.annotate_position(m[1]),
        "in_rbd": ref.in_rbd(m[1]),
        "count": c,
        "frequency": c / len(sub) if len(sub) else np.nan,
        "indel_ambiguous": align.is_indel_ambiguous(m[1]),
    } for m, c in counts.items()]
    return (pd.DataFrame(rows).sort_values("count", ascending=False)
            .head(n).reset_index(drop=True))
