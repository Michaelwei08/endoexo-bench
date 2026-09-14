#!/usr/bin/env python3
"""Baselines, the learned discriminator, and the two cross-validation schemes.

The task: among reads whose best competitive hit is the EXOGENOUS reference --
i.e. the reads a detection pipeline would COUNT -- separate those that really
came from the exogenous strain from those that came from an endogenous locus.
That is the false-positive mechanism the whole benchmark is about.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_curve
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def fpr_at_sensitivity(y_true: np.ndarray, score: np.ndarray, target: float = 0.95) -> float:
    """False-positive rate at the loosest threshold still reaching `target` recall.

    This is the operating point a detection assay is actually tuned at; ROC AUC
    averages over thresholds nobody would ever use.
    """
    if y_true.sum() == 0 or (1 - y_true).sum() == 0:
        return float("nan")
    fpr, tpr, _ = roc_curve(y_true, score)
    ok = tpr >= target
    return float(fpr[ok][0]) if ok.any() else float("nan")


def rule_baselines(feats: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """The rules actually in production use, translated to these features.

    unique_best is the AS > XS test used throughout the source pipeline;
    gap_ge_6 is motivated by the observed residual cross-mappings sitting at an
    AS-XS gap of only 3 to 5.
    """
    return {
        "rule_unique_best": (feats["as_minus_xs"] > 0).astype(int),
        "rule_unique_best_gap6": ((feats["as_minus_xs"] >= 6)).astype(int),
        "rule_gap6_and_len100": ((feats["as_minus_xs"] >= 6) &
                                 (feats["aligned_len"] >= 100)).astype(int),
        "rule_no_close_alt": (feats["n_alt_within5"] == 0).astype(int),
    }


def confusion_rates(y_true: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    """(sensitivity, false-positive rate) for a hard 0/1 rule."""
    pos, neg = y_true == 1, y_true == 0
    sens = float(pred[pos].mean()) if pos.any() else float("nan")
    fpr = float(pred[neg].mean()) if neg.any() else float("nan")
    return sens, fpr


def cv_scores(X: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
              scheme: str, model: str, n_splits: int = 5, seed: int = 0) -> np.ndarray:
    """Out-of-fold scores under a grouped or a read-level split."""
    if scheme == "grouped":
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        split = splitter.split(X, y, groups)
    elif scheme == "read_level":
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        split = splitter.split(X, y)
    else:
        raise ValueError(scheme)

    oof = np.full(len(y), np.nan)
    for train, test in split:
        if len(np.unique(y[train])) < 2:
            continue
        clf = _make_model(model, seed)
        clf.fit(X[train], y[train])
        oof[test] = clf.predict_proba(X[test])[:, 1]
    return oof


def _make_model(model: str, seed: int):
    if model == "logreg":
        return make_pipeline(StandardScaler(),
                             LogisticRegression(max_iter=2000, random_state=seed))
    if model == "gbm":
        return HistGradientBoostingClassifier(max_iter=200, random_state=seed)
    raise ValueError(model)


def evaluate_oof(y: np.ndarray, oof: np.ndarray, target: float = 0.95) -> dict[str, float]:
    ok = ~np.isnan(oof)
    if ok.sum() == 0 or len(np.unique(y[ok])) < 2:
        return {"pr_auc": float("nan"), "fpr_at_sens": float("nan"), "n_scored": int(ok.sum())}
    return {
        "pr_auc": float(average_precision_score(y[ok], oof[ok])),
        "fpr_at_sens": fpr_at_sensitivity(y[ok], oof[ok], target),
        "n_scored": int(ok.sum()),
    }
