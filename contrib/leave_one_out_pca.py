"""
leave_one_out_pca.py  --  Leakage-free PCA clutter suppression
==============================================================

WHY
---
The plain PCA clutter-suppression mode (see `clutter_suppression.py`) learns a
clutter subspace from healthy scans of the same shell. When those same healthy
scans are ALSO used as the test negatives in an ROC/AUC evaluation, the model
has effectively "seen" the data it is being scored on -- a form of train/test
leakage that inflates the AUC.

Leave-one-out (LOO) fixes this: when scoring a given healthy (negative) scan,
that scan is excluded from the clutter model. Each negative is therefore scored
by a model that has never seen it, giving an honest estimate of how the method
would perform on a genuinely new scan -- the only number that matters clinically.

Tumour (positive) scans are never used to train the clutter model, so they do
not need to be left out; only the negatives do.

RESULT IN THIS PROJECT
----------------------
On dense (BI-RADS 3) phantoms, LOO did not reduce AUC relative to the leaky
estimate -- it slightly improved it (e.g. DMAS 72.2% -> 78.2%), consistent with
a cleaner, less over-fit clutter subspace. The adipose-reference baseline fell
below the 95% bootstrap CI of the LOO-PCA AUC for every reconstruction method.
"""
import numpy as np
from sklearn.decomposition import PCA


def pca_clutter_loo(idx, e, s11, md, exclude_id=None, target_birads=3):
    """
    PCA clutter suppression with optional leave-one-out.

    Parameters
    ----------
    idx : int            index of the target scan in s11/md
    e : dict             metadata for the target scan
    s11, md :            full dataset
    exclude_id : int     scan id to hold OUT of the clutter model. Pass the
                         target's own id when scoring a healthy NEGATIVE, so it
                         is excluded from its own model (leave-one-out). Pass
                         None for tumour scans (no leakage to remove).
    target_birads : int  density class of healthy scans used for the subspace.

    Returns
    -------
    cleaned (n_freqs, n_antennas) complex array, or None if too few healthy
    same-shell scans remain to fit the subspace.
    """
    tar_fd = s11[idx, :, :]
    F, A = tar_fd.shape
    shell = e["phant_id"].split("F")[0]

    healthy = [s11[i, :, :] for i, h in enumerate(md)
               if h.get("birads") == target_birads
               and np.isnan(h.get("tum_diam", np.nan))
               and h["phant_id"].split("F")[0] == shell
               and h["id"] != exclude_id]          # <-- the leave-one-out

    if len(healthy) < 2:
        return None

    H = np.array(healthy)
    N_h = H.shape[0]
    H_real = np.hstack([np.real(H.reshape(N_h, F * A)),
                        np.imag(H.reshape(N_h, F * A))])
    pca = PCA(n_components=min(3, N_h))
    pca.fit(H_real)

    t_real = np.hstack([np.real(tar_fd.reshape(1, F * A)),
                        np.imag(tar_fd.reshape(1, F * A))])
    clutter = pca.inverse_transform(pca.transform(t_real))
    clutter = clutter[0, :F * A] + 1j * clutter[0, F * A:]
    return tar_fd - clutter.reshape(F, A)
