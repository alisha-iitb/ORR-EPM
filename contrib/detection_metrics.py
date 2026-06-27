"""
detection_metrics.py  --  Localization error, SCR, and detection scoring
========================================================================

Scores a reconstructed image against the known tumour position. Used to turn
each reconstruction into the (LE, SCR, detection) triple and the raw SCR scores
that feed the ROC analysis.

CONVENTIONS USED HERE
---------------------
1. PIXEL -> COORDINATE MAPPING. Row 0 is the top of the image (+y); the mapping
   is  y = +roi_rad - row*pixel - pixel/2,  x = -roi_rad + col*pixel + pixel/2.
   Getting this convention consistent with the dataset's tumour coordinates is
   essential -- an inconsistent sign mirrors the image and places the scored peak
   far from the true tumour.

2. SCR DEFINITION (per the EPM paper). SCR is the ratio (dB) of the peak inside a
   3 cm-diameter region centred on the IMAGE MAXIMUM to the peak outside it -- not
   a region centred on the true tumour. This is the quantity the detection
   threshold is meant to act on.

3. GRID SIZE. The reconstruction grid size (`m_size`) is passed in explicitly so
   scoring always matches the grid the image was reconstructed on.
"""
import numpy as np


def evaluate_tumor_detection(img_array, tum_x, tum_y, tum_diam,
                             roi_rad, m_size,
                             scr_threshold_db=1.5, le_tolerance_cm=0.5):
    """
    Returns (localization_error_cm, scr_db, detected_str).

    Parameters
    ----------
    img_array : 2D complex/real reconstruction
    tum_x, tum_y : true tumour centre [cm]
    tum_diam : tumour diameter [cm]
    roi_rad : image half-width [cm]
    m_size : pixels per dimension (must match the reconstruction grid!)
    scr_threshold_db : detection threshold on SCR
    le_tolerance_cm : positioning tolerance added to the tumour radius
    """
    img_mag = np.abs(img_array)
    img_mag[~np.isfinite(img_mag)] = 0
    pixel_size = (2 * roi_rad) / m_size

    # Coordinate grids. NOTE the y-axis convention (row 0 = +roi_rad), which is
    # (see convention note 1 in the module docstring)
    Y, X = np.ogrid[:m_size, :m_size]
    y_coords = roi_rad - Y * pixel_size - (pixel_size / 2)
    x_coords = -roi_rad + X * pixel_size + (pixel_size / 2)

    # Location of the image maximum (the candidate tumour).
    max_idx = np.unravel_index(np.argmax(img_mag), img_mag.shape)
    pred_y = roi_rad - max_idx[0] * pixel_size - (pixel_size / 2)
    pred_x = -roi_rad + max_idx[1] * pixel_size + (pixel_size / 2)
    le = np.sqrt((pred_x - tum_x) ** 2 + (pred_y - tum_y) ** 2)

    # SCR: peak inside a 1.5 cm-radius (3 cm-diameter) region centred on the
    # IMAGE MAXIMUM, versus the peak outside it (per the EPM paper).
    dist_from_peak = np.sqrt((x_coords - pred_x) ** 2 + (y_coords - pred_y) ** 2)
    roi_mask = dist_from_peak <= 1.5
    clutter_mask = dist_from_peak > 1.5
    max_signal = np.max(img_mag[roi_mask]) if np.any(roi_mask) else 0
    max_clutter = np.max(img_mag[clutter_mask]) if np.any(clutter_mask) else 0
    scr_linear = max_signal / max_clutter if max_clutter > 0 else 0
    scr_db = 20 * np.log10(scr_linear) if scr_linear > 0 else 0

    # Detection requires high SCR AND the peak being near the true tumour.
    le_threshold = (tum_diam / 2.0) + le_tolerance_cm
    detected = "Yes" if (scr_db >= scr_threshold_db and le <= le_threshold) else "No"
    return round(le, 2), round(scr_db, 2), detected


def raw_scr(img, roi_rad, m_size):
    """SCR in dB only (no detection decision). Used for ROC threshold sweeps."""
    m = np.abs(img)
    m[~np.isfinite(m)] = 0
    ps = (2 * roi_rad) / m_size
    Y, X = np.ogrid[:m_size, :m_size]
    yc = roi_rad - Y * ps - ps / 2
    xc = -roi_rad + X * ps + ps / 2
    pk = np.unravel_index(np.argmax(m), m.shape)
    py = roi_rad - pk[0] * ps - ps / 2
    px = -roi_rad + pk[1] * ps + ps / 2
    d = np.sqrt((xc - px) ** 2 + (yc - py) ** 2)
    sig = m[d <= 1.5].max()
    clut = m[d > 1.5].max() if np.any(d > 1.5) else 0
    return 20 * np.log10(sig / clut) if clut > 0 else 0.0
