# Results — Dense-Breast Evaluation (DDP Stage 2)

> This folder is **optional and self-contained**. Delete it to remove all
> thesis-specific results and leave a clean code fork. Nothing in `contrib/`
> depends on it.

All results are on **dense (BI-RADS 3) phantoms** from UM-BMID Gen-3 unless
stated otherwise. Dense cohort = **15 tumour + 15 healthy** scans; fatty
(BI-RADS 1) cohort = **35 + 35**.

---

## Table 1 — Detection sensitivity by tumour size (ideal subtraction)

Scans detected / total (3 scans per size).

| Tumour size | DAS | DMAS | ORR | ORR-EPM-T |
|-------------|-----|------|-----|-----------|
| 10 mm | 1/3 | 1/3 | 1/3 | **2/3** |
| 15 mm | 1/3 | 2/3 | 2/3 | 2/3 |
| 20 mm | 3/3 | 3/3 | 3/3 | 3/3 |
| 25 mm | 3/3 | 3/3 | 3/3 | 3/3 |
| 30 mm | 3/3 | 3/3 | 3/3 | 3/3 |

All methods detect ≥20 mm reliably; ORR-EPM-T uniquely catches a 10 mm tumour.

---

## Table 2 — Clutter-suppression AUC matrix (the central result)

AUC (%) on dense phantoms, by reconstruction method × clutter strategy.

| Method | Adipose ref | SVD | PCA | LOO-PCA |
|--------|-------------|-----|-----|---------|
| DAS | 37.8 | 43.3 | 63.3 | **77.8** |
| DMAS | 40.4 | 44.2 | 72.2 | **78.2** |
| ORR | 49.3 | 50.4 | 66.7 | **79.6** |
| ORR-EPM-T | 53.1 | 64.4 | 60.0 | 58.2 |

**Reading it:** performance climbs steadily Adipose → SVD → PCA → LOO-PCA for
DAS/DMAS/ORR. ORR-EPM-T is the exception — flat-to-declining — which is the basis
of the "physics and clutter suppression are partial substitutes" finding.

---

## Table 3 — LOO-PCA AUC with 95% bootstrap CIs (2000 resamples)

| Method | LOO-PCA AUC | 95% CI | Adipose baseline | Baseline below CI? |
|--------|-------------|--------|------------------|--------------------|
| DAS | 77.8% | [58.2, 91.3] | 37.8% | **Yes** |
| DMAS | 78.2% | [59.5, 92.6] | 40.4% | **Yes** |
| ORR | 79.6% | [61.3, 91.6] | 49.3% | **Yes** |

For every method the baseline falls below the 95% CI → the improvement from
reference-free PCA is **statistically significant** on this sample.

---

## Table 4 — Localization error on detected scans

| Method | Mean LE (cm) | n detected |
|--------|--------------|-----------|
| DAS | 0.65 | 11 |
| DMAS | 0.68 | 12 |
| ORR | 0.67 | 13 |
| ORR-EPM-T | 1.14 | 13 |

All methods localize within ~2 cm. ORR-EPM-T's slightly larger LE reflects its
edge-emphasizing path-dependent physics.

---

## Figures

Place the exported PNGs in `results/figures/`:

| File | Description |
|------|-------------|
| `FIG1_clutter_comparison_complete.png` | Clutter-suppression AUC bar chart with 95% CIs (all four methods × four conditions). |
| `FIG2_roc_dense_fatty.png` | ROC curves, dense vs. fatty, under adipose subtraction. |
| `FIG3_representative_recons.png` | Representative ORR reconstructions under adipose / SVD / PCA (cyan = true tumour). |
| `LE_boxplot_dense_detected.png` | Localization-error distribution on detected scans. |

---

## Raw data

CSV outputs (AUC inputs: per-scan SCR + label) can be placed in `results/csv/`.
The master summary is `MASTER_summary_table.csv`.

---

## Honest caveats

- Small cohort (15+15 dense): AUC CIs are wide (~±15 pts); the *comparison to
  baseline* is significant, the point estimates are not precise.
- Phantoms, not patients: rigid phantoms favour clutter methods that assume
  consistent clutter; patient validation is required.
- PCA needs a healthy database; clinical performance with mismatched populations
  is unverified.
- Left/right differential imaging could not be evaluated on dense tissue:
  UM-BMID Gen-3 has only one BI-RADS 3 variant per shell.
