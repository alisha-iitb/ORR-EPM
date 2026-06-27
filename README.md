# Dense-Breast Evaluation & Reference-Free Clutter Suppression for Radar-Based Microwave Imaging

> **This is a fork of [TysonReimer/ORR-EPM](https://github.com/TysonReimer/ORR-EPM).**
> The upstream repository provides the Optimization-based Radar Reconstruction
> (ORR) framework and its enhanced-physics variant (ORR-EPM-T) from Reimer *et al.*
> That framework is used here **as the reconstruction engine**.
>
> This fork adds a new **application and evaluation study** built on top of it: it
> asks whether radar-based microwave breast imaging can work on **dense (BI-RADS 3)
> breasts** — the case where mammography is weakest — and whether it can do so
> **without the idealised reference scan** that the strongest published results
> rely on but that no real patient can provide. The contribution is therefore one
> of **new problem scope and new methodology** (reference-free clutter suppression),
> not a modification of the upstream algorithm.
>
> Produced as part of a Dual Degree Project (DDP) thesis; published so that anyone
> continuing this line of work has a clean, documented, runnable starting point.

---

## Motivation

Microwave breast imaging (MBI) is a non-ionizing, low-cost modality that exploits
the large dielectric contrast between tumour and healthy tissue. Its persistent
weakness across two decades of literature is **specificity in dense breasts** —
precisely the population for whom mammography is least effective and the clinical
need is greatest.

The strongest published MBI results depend on **reference subtraction**: imaging
the breast with and without the tumour and subtracting the two. That is fine in a
phantom lab but impossible in a clinic, where no tumour-free copy of a patient's
breast exists. This work therefore studies two questions on the dense cohort of
the open-access **UM-BMID Gen-3** dataset:

1. Does radar MBI work on **dense** breast specifically?
2. Can it work **without an idealised reference scan**, using clinically
   deployable, reference-free clutter suppression?

---

## Contributions

This fork contributes, on top of the ORR/ORR-EPM-T framework:

1. **A dense-breast (BI-RADS 3) evaluation pipeline.** The analysis is restricted
   to the dense cohort and includes healthy (tumour-free) phantoms as negatives,
   enabling specificity and threshold-free ROC analysis.

2. **Reference-free clutter suppression as the deployable alternative.** SVD, PCA,
   and a **leakage-corrected leave-one-out PCA** are introduced as substitutes for
   reference subtraction, none of which requires a second patient scan.

3. **A leakage-correction result.** Leaving each scored scan out of the subspace
   that cleans it (LOO) *raises* AUC rather than lowering it, showing the
   uncorrected PCA estimate had been suppressing genuine tumour contrast.

4. **An empirical characterization of the physics × clutter interaction.** The
   reconstruction algorithm and the clutter-suppression stage interact rather than
   compose: ORR-EPM-T is strongest under reference/SVD suppression but does *not*
   benefit from learned PCA, so in the deployable regime the *simpler* methods win.

5. **A rigorous evaluation protocol.** All methods are scored with threshold-free
   ROC/AUC and 2000-resample bootstrap confidence intervals, with significance
   assessed by whether the baseline falls outside the interval.

These are packaged as a modular `contrib/` library (see *Repository layout*) so
each piece can be reused independently.

---

## Key findings (from the thesis)

On dense (BI-RADS 3) phantoms; 15 tumour + 15 healthy scans.

1. **How clutter is removed matters more than which reconstruction algorithm is
   used.** With the realistic adipose reference, every method sat at or below
   chance (AUC 37.8–53.1%): the reference removes skin but leaves the bulky,
   tumour-like fibroglandular clutter.

2. **Reference-free leave-one-out PCA recovers usable performance** — lifting DAS,
   DMAS, and base ORR from 37.8–49.3% to **77.8–79.6% AUC**, with bootstrap
   intervals that exclude the baseline. Needing no second scan, it is clinically
   deployable.

3. **The honest (leakage-corrected) estimator is also the better one.** Correcting
   leakage in plain PCA *raised* AUC (e.g. DAS 63.3% → 77.8%).

4. **Enhanced physics helps only when clutter is removed ideally.** Under
   reference or SVD suppression, ORR-EPM-T was strongest and uniquely detected the
   smallest tumours — reproducing the published small-tumour advantage, now
   isolated to dense breast.

5. **Reconstruction and clutter suppression interact; they do not compose.**
   ORR-EPM-T peaked under SVD (64.4%) and did not benefit from PCA — the reverse
   ordering of the simpler methods. In the deployable regime, the simpler methods
   win.

6. **ORR-EPM-T detects well but localises poorly** (LE 1.14 cm vs ~0.65 cm),
   consistent with its path-dependent attenuation concentrating response at the
   tumour's superficial edge.

> **Bottom line:** accurate dense-breast detection no longer requires a reference
> scan that cannot be obtained in practice; the gains come from better clutter
> removal rather than from an increasingly detailed reconstruction model.

Full tables, figures, and CSVs are in [`results/`](results/), which is optional
and deletable (see below).

---

## Key concepts (quick reference)

**Radar forward model.** Each method relates a reflectivity image `R(r)` to the
measured signal via `S11(ω, rₐ) = ∫ R(r)·exp(−jω·2τ(r,rₐ)) dr`, with
`τ = |r − rₐ|/v`. Reconstruction inverts this.

- **DAS / DMAS** — delay-and-sum and its multiply variant (classical beamformers).
- **ORR** — minimizes the data-fidelity cost by gradient descent.
- **ORR-EPM-T** — ORR with enhanced physics (beam/gain, spherical spreading, loss
  via a complex wavenumber, and a path-dependent air→breast→air delay carried in a
  complex *partitioned phase factor*).

**Clutter suppression** (the air–skin reflection is ~100× the tumour):
- *Reference subtraction* (ideal/adipose) — needs a tumour-free scan; not clinical.
- *SVD* — remove top singular values of one scan; reference-free, single-scan.
- *PCA* — learn the clutter subspace from same-shell healthy scans.
- *Leave-one-out PCA* — PCA with each scored negative excluded from its own model
  (removes train/test leakage).

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
