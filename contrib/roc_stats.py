"""
roc_stats.py  --  ROC, AUC, and bootstrap confidence intervals
==============================================================

Small, dependency-light helpers for turning per-scan (SCR, label) pairs into
ROC curves, AUC point estimates, and bootstrap confidence intervals.

CONVENTION
----------
`label` is 1 for tumour-containing scans (positives) and 0 for healthy scans
(negatives). `scr` is the raw signal-to-clutter ratio in dB (see
`detection_metrics.raw_scr`). The detection threshold is swept over 0..30 dB.
"""
import numpy as np


def auc_from_scores(scr, label, thr=None):
    """Sweep the SCR threshold and integrate the ROC curve to get AUC."""
    if thr is None:
        thr = np.linspace(0, 30, 1000)
    P = label == 1
    N = ~P
    nP, nN = P.sum(), N.sum()
    if nP == 0 or nN == 0:
        return np.nan
    tpr = np.array([((scr >= t) & P).sum() / nP for t in thr])
    fpr = np.array([((scr >= t) & N).sum() / nN for t in thr])
    order = np.argsort(fpr)
    return abs(np.trapezoid(tpr[order], fpr[order]))


def roc_curve_points(scr, label, thr=None):
    """Return (fpr, tpr) arrays for plotting an ROC curve."""
    if thr is None:
        thr = np.linspace(0, 30, 1000)
    P = label == 1
    N = ~P
    tpr = np.array([((scr >= t) & P).sum() / P.sum() for t in thr])
    fpr = np.array([((scr >= t) & N).sum() / N.sum() for t in thr])
    return fpr, tpr


def bootstrap_auc_ci(scr, label, n_boot=2000, seed=0):
    """
    Bootstrap (point estimate, 2.5th percentile, 97.5th percentile) of AUC.

    Resamples scans with replacement n_boot times. With small cohorts the CI is
    wide; report it honestly alongside the point estimate.
    """
    rng = np.random.default_rng(seed)
    n = len(scr)
    point = auc_from_scores(scr, label)
    aucs = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        s, l = scr[i], label[i]
        if (l == 1).any() and (l == 0).any():
            a = auc_from_scores(s, l)
            if not np.isnan(a):
                aucs.append(a)
    aucs = np.array(aucs)
    return point, np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)


def paired_bootstrap_difference(scr_a, lab_a, scr_b, lab_b,
                                n_boot=2000, seed=1):
    """
    Paired bootstrap of AUC(method_a) - AUC(method_b).

    Returns (mean_difference, ci_low, ci_high, prob_a_greater_than_b).
    Inputs must be aligned (same scan order) so the same resample index applies
    to both methods.
    """
    rng = np.random.default_rng(seed)
    n = len(scr_a)
    diffs = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        a = auc_from_scores(scr_a[i], lab_a[i])
        b = auc_from_scores(scr_b[i], lab_b[i])
        if not (np.isnan(a) or np.isnan(b)):
            diffs.append(a - b)
    diffs = np.array(diffs)
    return (diffs.mean(),
            np.percentile(diffs, 2.5),
            np.percentile(diffs, 97.5),
            (diffs > 0).mean())
