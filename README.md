# Dense-Breast Evaluation & Reference-Free Clutter Suppression for Radar-Based Microwave Imaging

Fork of [TysonReimer/ORR-EPM](https://github.com/TysonReimer/ORR-EPM). The upstream ORR / ORR-EPM-T framework is used **as the reconstruction engine**; this fork adds an evaluation study on top of it. Produced for a Dual Degree Project (DDP) thesis.

## Overview

Microwave breast imaging (MBI) is limited primarily by poor specificity in **dense (BI-RADS 3) breasts** — the population where mammography is weakest. The strongest published MBI results rely on **reference subtraction** (imaging the breast with and without the tumour), which is not obtainable clinically. This work evaluates, on the dense cohort of the open-access **UM-BMID Gen-3** dataset:

1. Whether radar MBI works on dense breast specifically, and
2. Whether it works **without an idealised reference scan**, using reference-free clutter suppression (SVD, PCA, leave-one-out PCA).

All methods are scored with threshold-free ROC/AUC and 2000-resample bootstrap confidence intervals.

## Contributions

- **Dense-breast evaluation pipeline** restricted to BI-RADS 3, with healthy phantoms included as negatives to enable ROC/specificity analysis.
- **Reference-free clutter suppression** (SVD, PCA, and a leakage-corrected leave-one-out PCA) as a clinically deployable alternative to reference subtraction.
- **Leakage correction**: excluding each scored scan from the subspace that cleans it (LOO) *raises* AUC, indicating the uncorrected estimate suppressed genuine tumour contrast.
- **Physics × clutter interaction**: ORR-EPM-T's ordering across suppression modes is the reverse of the simpler methods', showing the two stages interact rather than compose.
- **Evaluation protocol**: threshold-free ROC/AUC with bootstrap CIs; significance assessed against the baseline interval.

## Key findings

Dense (BI-RADS 3) cohort, 15 tumour + 15 healthy scans:

- **Clutter removal dominates the choice of reconstruction algorithm.** Under the realistic adipose reference, all methods sat at or below chance (AUC 37.8–53.1%): it removes skin but leaves tumour-like fibroglandular clutter.
- **Reference-free LOO-PCA recovers usable performance**, lifting DAS/DMAS/ORR to **77.8–79.6% AUC**, with bootstrap intervals excluding the baseline — and requiring no second scan.
- **Leakage correction improves the estimate** (e.g. DAS 63.3% → 77.8%).
- **Enhanced physics helps only under ideal clutter removal.** ORR-EPM-T was strongest under reference/SVD suppression and uniquely detected the smallest tumours, but did not benefit from learned PCA; in the deployable regime the simpler methods win.
- **ORR-EPM-T detects well but localises poorly** (LE 1.14 cm vs ~0.65 cm), consistent with its path-dependent attenuation emphasising the tumour's superficial edge.

**Bottom line:** dense-breast detection no longer requires a clinically-unobtainable reference scan; the gains come from clutter removal rather than from a more detailed reconstruction model.

Full tables, figures, and CSVs are in [`results/`](results/)

## Methods

**Reconstruction.** DAS and DMAS beamformers; ORR (gradient-descent minimisation of the data-fidelity cost); ORR-EPM-T (ORR with beam/gain, spherical spreading, loss via a complex wavenumber, and a path-dependent air→breast→air delay carried in a complex partitioned phase factor).

**Clutter suppression.** Reference subtraction (`ideal`/`adipose`, baseline); SVD (single-scan, reference-free); PCA (clutter subspace learned from same-shell healthy scans); leave-one-out PCA (PCA with each scored negative excluded from its own subspace).

## Repository layout

---

## Repository layout

```
.
├── umbms/                     # upstream ORR/ORR-EPM-T package (unchanged engine)
├── run/                       # upstream run scripts
├── contrib/                   # >>> this work: modular, documented library
│   ├── epmt_phase_factor.py       # build the partitioned ORR-EPM-T phase factor
│   ├── clutter_suppression.py     # selectable strategies: ideal/adipose/svd/pca/...
│   ├── leave_one_out_pca.py       # leakage-free PCA evaluation
│   ├── detection_metrics.py       # LE / SCR / detection scoring
│   ├── roc_stats.py               # ROC, AUC, bootstrap CIs, paired bootstrap
│   └── patch_algos.py             # small helper: log optimizer cost history for convergence checks
├── results/                   # >>> OPTIONAL thesis results (tables, figures, CSVs); deletable
├── CONTRIBUTIONS.md           # what this work adds, mapped to the code
├── PUSH_GUIDE.md              # how to put this on your GitHub
└── README.md
```

`contrib/` and `results/` are self-contained; deleting `results/` leaves a clean
code fork.

---

## Acknowledgements

This work uses the ORR / ORR-EPM-T framework and the UM-BMID dataset of
**T. Reimer, S. Christie, I. Prykhodko, and S. Pistorius** (University of
Manitoba) as its reconstruction engine and data source. The reconstruction
algorithm and dataset are theirs; this fork contributes a dense-breast evaluation
and a reference-free clutter-suppression study built on top of them.

**References**
1. T. Reimer and S. Pistorius, "An Optimization-Based Approach to Radar Image
   Reconstruction in Breast Microwave Sensing," *Sensors*, 21(24):8172, 2021.
2. T. Reimer, S. Christie, I. Prykhodko, S. Pistorius, "Enhanced Physics
   Modelling in Radar-Based Microwave Imaging for Breast Cancer Detection,"
   *IEEE J. Electromagn. RF Microw. Med. Biol.*, 2025.
3. T. Reimer, J. Krenkevich, S. Pistorius, "An open-access experimental dataset
   for breast microwave imaging," *EuCAP*, 2020.
