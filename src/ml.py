"""Stage 7: classification, ablation and clustering.

Three things are modelled, for different reasons.

Geographic origin, within one variant. This is the honest predictive task.
Region is not part of how a lineage is defined, so a model has to find real
signal, and it may fail. The variant is held fixed because region and variant
are confounded in GenBank: Omicron records are overwhelmingly North American,
so a model free to read the variant could infer the region from that instead
of from geography. Holding the variant fixed removes that shortcut.

Variant, as an ablation. Predicting variant from spike is circular, because
lineages are defined by these mutations, and both the model and a trivial
baseline score near 100%. It is still worth running as an ablation: comparing
the full spike against an RBD-masked spike and the RBD alone measures how much
variant signal sits outside the receptor-binding domain, which is a real
quantity rather than a restatement of the labels.

Clustering, unsupervised. Whether structure in spike sequence recovers the
variant labels without being told them, scored with the adjusted Rand index.

Leakage is controlled in three ways: identical mutation profiles are removed
before splitting, so no sequence appears in both training and test; feature
selection and scaling run inside each cross-validation fold through a
pipeline, never on the full dataset; and the test set is evaluated once.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (adjusted_rand_score, balanced_accuracy_score,
                             classification_report, confusion_matrix,
                             roc_auc_score)
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline

import spike_reference as ref

RANDOM_STATE = 42


def mutation_features(df: pd.DataFrame, reference: str,
                      column: str = "spike") -> pd.DataFrame:
    """Binary isolates x positions matrix: 1 where the residue differs.

    Unresolved positions are encoded 0 (treated as "no evidence of a
    mutation") rather than 1, so missing data cannot masquerade as signal.
    """
    length = len(reference)
    rows = np.zeros((len(df), length), dtype=np.int8)
    for i, seq in enumerate(df[column].to_numpy()):
        for j, (c, r) in enumerate(zip(seq, reference)):
            if c != r and c != "X":
                rows[i, j] = 1
    return pd.DataFrame(rows, index=df.index,
                        columns=[f"pos_{p}" for p in range(1, length + 1)])


def mask_region(X: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Zero out a residue range, to remove its contribution to the model."""
    out = X.copy()
    cols = [f"pos_{p}" for p in range(start, end + 1) if f"pos_{p}" in out.columns]
    out[cols] = 0
    return out


def keep_only_region(X: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Keep only a residue range."""
    cols = [f"pos_{p}" for p in range(start, end + 1) if f"pos_{p}" in X.columns]
    return X[cols]


def build_models(k_features: int = 200) -> dict[str, Pipeline]:
    """Pipelines so that feature selection happens inside each fold."""
    return {
        "Majority-class baseline": Pipeline([
            ("clf", DummyClassifier(strategy="most_frequent")),
        ]),
        "Logistic regression": Pipeline([
            ("select", SelectKBest(chi2, k=k_features)),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced",
                                       random_state=RANDOM_STATE)),
        ]),
        "Random forest": Pipeline([
            ("select", SelectKBest(chi2, k=k_features)),
            ("clf", RandomForestClassifier(n_estimators=400, class_weight="balanced",
                                           random_state=RANDOM_STATE, n_jobs=-1)),
        ]),
    }


def evaluate(X: pd.DataFrame, y: pd.Series, label: str,
             k_features: int = 200, test_size: float = 0.25,
             n_splits: int = 5, verbose: bool = True) -> pd.DataFrame:
    """Cross-validate on a training split, then score once on a held-out test set."""
    k_features = min(k_features, X.shape[1])
    counts = y.value_counts()
    if len(counts) < 2 or counts.min() < n_splits:
        if verbose:
            print(f"  [{label}] skipped: smallest class has {counts.min()} members, "
                  f"too few for {n_splits}-fold cross-validation")
        return pd.DataFrame()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=RANDOM_STATE)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_STATE)

    rows = []
    for name, pipe in build_models(k_features).items():
        cv_scores = cross_val_score(pipe, X_train, y_train, cv=cv,
                                    scoring="balanced_accuracy", n_jobs=1)
        pipe.fit(X_train, y_train)
        pred = pipe.predict(X_test)
        bacc = balanced_accuracy_score(y_test, pred)
        auc = np.nan
        if len(counts) == 2 and hasattr(pipe, "predict_proba"):
            try:
                proba = pipe.predict_proba(X_test)[:, 1]
                auc = roc_auc_score((y_test == pipe.classes_[1]).astype(int), proba)
            except Exception:
                pass
        rows.append({
            "task": label, "model": name,
            "cv_balanced_accuracy": round(float(np.mean(cv_scores)), 3),
            "cv_sd": round(float(np.std(cv_scores)), 3),
            "test_balanced_accuracy": round(float(bacc), 3),
            "test_roc_auc": round(float(auc), 3) if auc == auc else np.nan,
            "n_train": len(y_train), "n_test": len(y_test),
            "n_features": X.shape[1],
        })
    out = pd.DataFrame(rows)
    if verbose:
        print(f"\n  [{label}]  n={len(y)}  classes={dict(counts)}")
        print(out[["model", "cv_balanced_accuracy", "cv_sd",
                   "test_balanced_accuracy", "test_roc_auc"]].to_string(index=False))
    return out


def detailed_report(X: pd.DataFrame, y: pd.Series, k_features: int = 200,
                    test_size: float = 0.25) -> tuple[str, pd.DataFrame]:
    """Per-class report and confusion matrix for the stronger model."""
    k_features = min(k_features, X.shape[1])
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=RANDOM_STATE)
    pipe = build_models(k_features)["Random forest"]
    pipe.fit(X_train, y_train)
    pred = pipe.predict(X_test)
    labels = sorted(y.unique())
    cm = pd.DataFrame(confusion_matrix(y_test, pred, labels=labels),
                      index=[f"true {l}" for l in labels],
                      columns=[f"pred {l}" for l in labels])
    return classification_report(y_test, pred, zero_division=0), cm


def rbd_ablation(X: pd.DataFrame, y: pd.Series, k_features: int = 200) -> pd.DataFrame:
    """How much variant signal survives outside the receptor-binding domain."""
    start, end = ref.RBD
    variants = {
        "Full spike": X,
        "RBD masked": mask_region(X, start, end),
        "RBD only": keep_only_region(X, start, end),
        "NTD only": keep_only_region(X, *ref.SPIKE_DOMAINS["NTD"]),
    }
    frames = [evaluate(Xv, y, label, k_features=k_features, verbose=False)
              for label, Xv in variants.items()]
    return pd.concat([f for f in frames if not f.empty], ignore_index=True)


def clustering_report(X: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    """Do unsupervised clusters recover the variant labels?"""
    k = y.nunique()
    rows = []
    for name, model in [
        ("KMeans", KMeans(n_clusters=k, n_init=10, random_state=RANDOM_STATE)),
        ("Agglomerative (Ward)", AgglomerativeClustering(n_clusters=k)),
    ]:
        labels = model.fit_predict(X.to_numpy(dtype=float))
        rows.append({"method": name, "n_clusters": k,
                     "adjusted_rand_index": round(adjusted_rand_score(y, labels), 3)})
    return pd.DataFrame(rows)


def redundancy_curve(X: pd.DataFrame, y: pd.Series, sizes=(1, 2, 3, 5, 10, 20, 50, 100),
                     n_repeats: int = 20, random_state: int = RANDOM_STATE) -> pd.DataFrame:
    """How many spike positions are needed to separate the variants?

    Masking the RBD does not reduce variant accuracy at all, which says the
    signal is not localised: it is repeated across the protein. To measure that
    redundancy, randomly chosen subsets of positions of increasing size are used
    as the only features. The size at which accuracy saturates is the number of
    arbitrary positions sufficient to identify a variant.

    Positions are drawn from those that actually vary, since constant columns
    carry no information and would only dilute the sample.
    """
    rng = np.random.default_rng(random_state)
    varying = [c for c in X.columns if X[c].nunique() > 1]
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)
    rows = []
    for k in sizes:
        if k > len(varying):
            continue
        scores = []
        for _ in range(n_repeats):
            cols = list(rng.choice(varying, size=k, replace=False))
            pipe = Pipeline([("clf", RandomForestClassifier(
                n_estimators=200, class_weight="balanced",
                random_state=random_state, n_jobs=-1))])
            scores.append(cross_val_score(pipe, X[cols], y, cv=cv,
                                          scoring="balanced_accuracy").mean())
        rows.append({
            "n_positions": k,
            "mean_balanced_accuracy": round(float(np.mean(scores)), 3),
            "sd": round(float(np.std(scores)), 3),
            "min": round(float(np.min(scores)), 3),
            "max": round(float(np.max(scores)), 3),
        })
    return pd.DataFrame(rows)


def top_discriminative_positions(X: pd.DataFrame, y: pd.Series,
                                 reference: str, n: int = 15) -> pd.DataFrame:
    """Positions a random forest relies on most, with domain annotation."""
    model = RandomForestClassifier(n_estimators=500, class_weight="balanced",
                                   random_state=RANDOM_STATE, n_jobs=-1)
    model.fit(X, y)
    imp = pd.Series(model.feature_importances_, index=X.columns)
    top = imp.sort_values(ascending=False).head(n)
    rows = []
    for col, importance in top.items():
        pos = int(col.split("_")[1])
        group_freq = {g: float(X.loc[y == g, col].mean()) for g in sorted(y.unique())}
        rows.append({
            "position": pos,
            "ref_aa": reference[pos - 1],
            "domain": ref.annotate_position(pos),
            "in_rbd": ref.in_rbd(pos),
            "importance": round(float(importance), 4),
            **{f"freq_{g}": round(v, 3) for g, v in group_freq.items()},
        })
    return pd.DataFrame(rows)
