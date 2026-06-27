"""
epmt_phase_factor.py  --  Building the partitioned ORR-EPM-T phase factor
========================================================================

ORR-EPM-T's enhanced physics is carried in a COMPLEX, PARTITIONED phase factor
that must be built (via `get_phase_fac_partitioned`) before calling `orr_recon`:
  * the air->breast->air partitioning gives the path-dependent delay (the "-T"),
  * the imaginary part of the breast wavenumber gives loss/attenuation
    (`complex_k`).

`enh_phys` then carries only ['beam','spherical','gain']. This module assembles
that factor for the dense-cohort experiments.

STEP SIZE
---------
For a 100x100 grid with 96 MHz frequency sampling, the gradient-descent step size
must be tuned: an empirical sweep found 3e5 to be the largest stable step (the
larger values implied elsewhere diverge to infinity on this configuration). Use
STEP_EPM_T = 3e5 here.
"""
import gc
import multiprocessing as mp
import numpy as np
from umbms.tdelay.partitioned import get_phase_fac_partitioned


def build_phase_fac_T(ant_rad, adi_rad, recon_fs, m_size, roi_rad,
                      gly_eps, gly_conds, eps_0, c, n_cores,
                      ini_a_ang=-136.0, ant_t_delay=0.19e-9):
    """
    Build the complex, partitioned phase factor for ORR-EPM-T.

    Returns the phase factor array, or None on failure (e.g. unknown breast
    radius). `gly_eps`/`gly_conds` are the (interpolated) glycerin dielectric
    permittivity and conductivity at each reconstruction frequency.
    """
    try:
        workers = mp.Pool(n_cores)
        # Complex breast wavenumber: real part = phase speed, imag part = loss.
        breast_k = (2 * np.pi * recon_fs
                    * np.sqrt(gly_eps - 1j * gly_conds
                              / (2 * np.pi * recon_fs * eps_0)) / c)
        air_k = 2 * np.pi * recon_fs / c

        pf = get_phase_fac_partitioned(
            ant_rho=ant_rad, m_size=m_size, roi_rho=roi_rad,
            air_k=air_k, breast_k=breast_k, ini_a_ang=ini_a_ang,
            adi_rad=adi_rad, phant_x=0, phant_y=0, worker_pool=workers)
        workers.close(); workers.join(); del workers; gc.collect()

        # Apply antenna time delay on top of the partitioned path delay.
        phase_fac_T = (np.exp(-1j * pf)
                       * np.exp(-1j * 2 * np.pi
                                * recon_fs[:, None, None, None] * ant_t_delay))
        del pf; gc.collect()
        return phase_fac_T
    except Exception as exc:
        print(f"  [X] Partitioned phase factor failed: {exc}")
        return None
