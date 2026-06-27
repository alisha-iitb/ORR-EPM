# Contributions — what this fork adds on top of [TysonReimer/ORR-EPM](https://github.com/TysonReimer/ORR-EPM)

This fork uses the upstream ORR / ORR-EPM-T algorithm **as its reconstruction
engine**. It does not modify that algorithm. What it adds is a new **application
(dense breast)** and a new **methodology (reference-free clutter suppression)**,
together with the evaluation code needed to study them. Each contribution below is
mapped to the module that implements it.

---

## 1. Dense-breast (BI-RADS 3) evaluation pipeline

**Module:** `contrib/clutter_suppression.py`, `contrib/detection_metrics.py`

The study is restricted to the dense (BI-RADS 3) cohort of UM-BMID Gen-3 and
includes **healthy (tumour-free) phantoms as negatives**, so that specificity and
threshold-free ROC can be computed (rather than detection on tumour-present scans
alone). Healthy dense scans use the adipose reference so positives and negatives
are processed identically.

---

## 2. The ORR-EPM-T partitioned phase factor

**Module:** `contrib/epmt_phase_factor.py`

The full enhanced-physics model (ORR-EPM-T) requires a **complex, partitioned
phase factor** built from `get_phase_fac_partitioned` before reconstruction:
- the air→breast→air partitioning gives the path-dependent delay (the "-T"),
- the imaginary part of the breast wavenumber gives loss/attenuation (`complex_k`).

`enh_phys` then carries only `['beam','spherical','gain']`. This module assembles
that factor for the dense-cohort experiments. (Empirically, the gradient-descent
step size for this grid/sampling is `3e5`; larger values diverge.)

---

## 3. Reference-free clutter suppression

**Module:** `contrib/clutter_suppression.py`

A single interface, `suppress_clutter(mode, ...)`, provides the strategies
compared in the thesis:
- `ideal` / `adipose` — reference subtraction (baseline; not clinically practical),
- `svd` — single-scan, reference-free,
- `pca` — population-learned clutter subspace,
- plus `avgtrace`, `tdg`.

The reference-free methods are the deployable alternative this work introduces:
none requires a tumour-free copy of the patient's breast.

---

## 4. Leakage-corrected leave-one-out PCA

**Module:** `contrib/leave_one_out_pca.py`

When PCA learns its clutter subspace from healthy scans that are also used as test
negatives, the model has "seen" the data it is scored on (train/test leakage).
`pca_clutter_loo(..., exclude_id=...)` excludes each scored negative from its own
clutter model, giving an honest estimate. In this study LOO **raised** AUC
(e.g. DAS 63.3% → 77.8%), indicating the uncorrected estimate had been suppressing
genuine tumour contrast.

---

## 5. Threshold-free ROC, AUC, and bootstrap statistics

**Module:** `contrib/roc_stats.py`

`auc_from_scores` (threshold sweep), `roc_curve_points`, `bootstrap_auc_ci`
(2000-resample 95% CIs), and `paired_bootstrap_difference` (method-vs-method
significance). All reported AUCs carry bootstrap CIs; significance against the
baseline is assessed by whether the baseline falls outside the interval.

---

## 6. Detection metric (LE / SCR / detection)

**Module:** `contrib/detection_metrics.py`

`evaluate_tumor_detection` returns localization error (LE), signal-to-clutter
ratio (SCR, measured on a 3 cm-diameter region centred on the image maximum, per
the EPM paper), and a detection decision (SCR and LE thresholds). `raw_scr`
returns SCR alone for ROC threshold sweeps. The grid size is passed in explicitly
so scoring always matches the reconstruction grid.

---

## 7. Convergence logging helper

**Module:** `contrib/patch_algos.py`

A small, optional helper used by this pipeline to **log the optimizer's cost
history** so convergence can be verified (the cost should drop substantially
across iterations). It makes `orr_recon` return `(img, cost_history)` and is
written surgically so it touches only `_orr_cpu`/`orr_recon` and leaves the
beamformers untouched. This is a convenience for the evaluation, not a change to
the reconstruction algorithm itself.

---

## Empirical findings (results, not code modules)

Two thesis conclusions are **experimental findings** produced by running the
modules above across the method × clutter-suppression grid, rather than standalone
functions:

- **Enhanced physics helps only when clutter is removed ideally** — ORR-EPM-T is
  strongest under reference/SVD suppression and uniquely detects the smallest
  tumours, but does not benefit from learned PCA.
- **Reconstruction and clutter suppression interact, not compose** — ORR-EPM-T's
  ordering across suppression modes is the reverse of the simpler methods'; in the
  deployable (reference-free) regime, the simpler methods win.

These are documented with their supporting numbers in [`results/`](results/).
