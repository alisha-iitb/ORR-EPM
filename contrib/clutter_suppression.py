"""
clutter_suppression.py  --  Selectable clutter-suppression strategies
=====================================================================

The air-skin reflection dominates a raw breast-microwave scan by roughly two
orders of magnitude, and must be removed before reconstruction. This module
collects the strategies compared in this project behind a single interface,
`suppress_clutter(mode, ...)`, so the reconstruction engine can switch methods
with one argument.

STRATEGIES
----------
  "ideal"     Subtract the matched fibroglandular (healthy) reference scan.
              Removes skin AND glandular clutter. Best-case, but requires a
              tumour-free copy of the same breast -> NOT clinically practical.

  "adipose"   Subtract the adipose-only reference (skin + fat, no glandular).
              Removes skin but leaves glandular clutter. Also reference-based.

  "avgtrace"  Subtract the per-frequency mean across antennas. Reference-free.

  "svd"       Zero the top singular values of the frequency x antenna matrix
              (the strong, antenna-consistent skin/structure response).
              Reference-free, single-scan -> most clinically deployable.

  "tdg"       Average-trace subtraction followed by time-domain gating with a
              Tukey window. Reference-free.

  "pca"       Learn a clutter subspace from healthy scans of the SAME shell and
              project it out of the target. Removes skin + common glandular
              clutter. Requires a database of healthy scans.

NOTE ON LEAKAGE
---------------
The plain "pca" mode learns the clutter subspace from same-shell healthy scans
that may also be used as test negatives elsewhere -- a form of train/test
leakage. For an honest estimate, use `pca_clutter_loo` (leave-one-out), which
excludes the scan being scored from its own clutter model. See
`contrib/leave_one_out_pca.py`.
"""
import numpy as np
from scipy.signal import windows
from sklearn.decomposition import PCA


def suppress_clutter(clutter_mode, original_index, tar_md, tar_fd,
                     s11, md, expt_ids, target_birads=3):
    """
    Return the clutter-suppressed frequency-domain signal, or None to skip.

    Parameters
    ----------
    clutter_mode : str   one of {"ideal","adipose","avgtrace","svd","tdg","pca"}
    original_index : int index of this scan in s11/md
    tar_md : dict        metadata for this scan
    tar_fd : ndarray     (n_freqs, n_antennas) complex S11 of this scan
    s11, md, expt_ids :  the full dataset arrays / lists
    target_birads : int  density class used to gather healthy scans for PCA
    """
    # ---- reference-based methods ------------------------------------------
    if clutter_mode == "ideal":
        try:
            ref_idx = expt_ids.index(tar_md["fib_ref_id"])
            return tar_fd - s11[ref_idx, :, :]
        except (ValueError, KeyError):
            print(f"  [!] fib_ref_id not found for scan {tar_md.get('id')}; skipping.")
            return None

    if clutter_mode == "adipose":
        try:
            ref_idx = expt_ids.index(tar_md["adi_ref_id"])
            return tar_fd - s11[ref_idx, :, :]
        except (ValueError, KeyError):
            print(f"  [!] adi_ref_id not found for scan {tar_md.get('id')}; skipping.")
            return None

    # ---- reference-free methods -------------------------------------------
    if clutter_mode == "avgtrace":
        # Subtract the mean trace over antennas (the antenna-common response).
        return tar_fd - np.mean(tar_fd, axis=1, keepdims=True)

    if clutter_mode == "svd":
        # Remove the two strongest singular components (skin / bulk structure).
        # NOTE: removing too few leaves clutter; too many erodes the tumour.
        U, S, Vh = np.linalg.svd(tar_fd, full_matrices=False)
        S2 = np.copy(S)
        S2[0:2] = 0.0
        return U @ np.diag(S2) @ Vh

    if clutter_mode == "tdg":
        # Average-trace removal + time-domain Tukey gating.
        x = tar_fd - np.mean(tar_fd, axis=1, keepdims=True)
        td = np.fft.ifft(x, axis=0)
        w = windows.tukey(td.shape[0], alpha=0.3).reshape(td.shape[0], 1)
        return np.fft.fft(td * w, axis=0)

    if clutter_mode == "pca":
        # Population-learned clutter subspace from same-shell healthy scans.
        F, A = tar_fd.shape
        shell = tar_md["phant_id"].split("F")[0]
        healthy = [s11[i, :, :] for i, e in enumerate(md)
                   if e.get("birads") == target_birads
                   and np.isnan(e.get("tum_diam", np.nan))
                   and e["phant_id"].split("F")[0] == shell]
        if len(healthy) < 2:
            print("  [!] Not enough healthy scans of this shell for PCA; skipping.")
            return None
        H = np.array(healthy)
        N_h = H.shape[0]
        # Complex data is split into real+imag so PCA (real-valued) can fit it.
        H_real = np.hstack([np.real(H.reshape(N_h, F * A)),
                            np.imag(H.reshape(N_h, F * A))])
        pca = PCA(n_components=min(3, N_h))
        pca.fit(H_real)
        t_real = np.hstack([np.real(tar_fd.reshape(1, F * A)),
                            np.imag(tar_fd.reshape(1, F * A))])
        clutter = pca.inverse_transform(pca.transform(t_real))
        clutter = clutter[0, :F * A] + 1j * clutter[0, F * A:]
        return tar_fd - clutter.reshape(F, A)

    raise ValueError(f"Unknown clutter_mode: {clutter_mode}")
