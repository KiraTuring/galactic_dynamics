#!/usr/bin/env python3
"""
Fit Multi-Gaussian Expansion (MGE) from an HST WFPC2 image using mgefit.

Pipeline:
    1. Read HST F814W drizzled image (SCI extension)
    2. find_galaxy() — locate center, ellipticity, PA
    3. sectors_photometry() — radial surface brightness in sectors
    4. mge_fit_sectors_regularized() — fit Gaussians with regularization
    5. Convert units: counts → AB mag/arcsec² → L⊙ pc⁻²
    6. Save results + diagnostic plots

Usage:
    python fit_mge_hst.py <galaxy>
    python fit_mge_hst.py NGC4621
    python fit_mge_hst.py NGC4621 --fwhm 0.13 --ngauss 15

Output (data/processed/{galaxy}/):
    mge_f814w.ecsv            — MGE parameters
    mge_f814w_components.png  — Gaussian ellipses + radial profile
    mge_f814w_contours.png    — Contour overlay on image
    mge_f814w_radial.png      — Radial profile vs fit
    mge_f814w_comparison.png  — Comparison with existing MGE (if found)

Dependencies:
    mgefit (Cappellari 2002), astropy, numpy, matplotlib
"""

import argparse, glob, sys
from pathlib import Path

import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from mgefit.find_galaxy import find_galaxy
from mgefit.sectors_photometry import sectors_photometry
from mgefit.mge_fit_sectors_regularized import mge_fit_sectors_regularized


# --- Physical constants ---
C_AA = 2.99792458e18          # speed of light in Å/s
ARCSEC_PER_RAD = 206265.0     # arcsec/rad
AB_ZP = -48.6                 # AB magnitude zeropoint
M_SUN_F814W = 4.56            # Solar absolute AB magnitude in F814W (Willmer 2018)

# --- Paths ---
_JAM_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
DATA_RAW = _JAM_ROOT / "data" / "raw"
DATA_PROC = _JAM_ROOT / "data" / "processed"


def read_hst_image(path):
    """Read SCI extension from HST drizzled FITS, return (data, metadata)."""
    with fits.open(path) as hdu:
        sci = hdu["SCI"]
        data = sci.data.astype(np.float64)

        cd1_1 = sci.header.get("CD1_1")
        cd2_2 = sci.header.get("CD2_2")
        if cd1_1 is not None:
            pixscale = abs(cd1_1) * 3600.0
        elif cd2_2 is not None:
            pixscale = abs(cd2_2) * 3600.0
        else:
            raise KeyError("Cannot determine pixel scale from CD matrix")

        photflam = sci.header.get("PHOTFLAM")
        centrwv = hdu[0].header.get("CENTRWV", 8012.0)

    return data, pixscale, photflam, centrwv


def background_level(image, fraction=0.99):
    """Estimate background threshold from cumulative flux."""
    flat = np.sort(image.flatten())[::-1]
    cum = np.cumsum(flat)
    cum /= cum[-1]
    return flat[cum <= fraction][-1]


def flux_to_abmag(f_nu):
    """Flux density (erg/cm²/s/Hz) → AB mag."""
    return -2.5 * np.log10(f_nu) - 48.6


def abmag_to_solar(mu_ab, m_sun=M_SUN_F814W):
    """AB mag/arcsec² → L⊙ pc⁻²."""
    return ARCSEC_PER_RAD ** 2 * 10 ** (-0.4 * (mu_ab - m_sun))


def plot_components(surf, sigma, q_obs, savepath):
    """MGE Gaussian ellipses + radial profile."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    for s, q in zip(sigma, q_obs):
        ell = Ellipse((0, 0), width=2 * s, height=2 * s * q,
                       angle=0, alpha=0.3, facecolor="C0", edgecolor="C0")
        ax.add_patch(ell)
    ax.set_aspect("equal")
    ax.set_xlim(-max(sigma) * 1.1, max(sigma) * 1.1)
    ax.set_ylim(-max(sigma) * 1.1, max(sigma) * 1.1)
    ax.set_xlabel("x (arcsec)")
    ax.set_ylabel("y (arcsec)")
    ax.set_title("MGE Gaussian components")

    ax = axes[1]
    ax.errorbar(sigma * (q_obs + 1) / 2, surf,
                xerr=[sigma * (1 - q_obs) / 2],
                fmt="o", capsize=4)
    ax.set_xlabel("sigma (arcsec)")
    ax.set_ylabel(r"$I\ (L_\odot\ \mathrm{pc}^{-2})$")
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def fit_mge_hst(galaxy, fwhm=None, ngauss=15, trim_margin=0.05):
    """
    Run full MGE fitting pipeline on an HST F814W image.

    Parameters
    ----------
    galaxy : str
        Galaxy name (e.g. "NGC4621").
    fwhm : float or None
        PSF FWHM in arcsec. Default (None): 0.13 arcsec for WFPC2 F814W.
    ngauss : int
        Number of Gaussian components (default 15).
    trim_margin : float
        Fractional margin to trim from image edges (default 0.05).
    """
    out_dir = DATA_PROC / galaxy
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Locate F814W image ---
    pattern = str(DATA_RAW / galaxy / "hst_*f814w*drz.fits")
    hst_list = sorted(glob.glob(pattern))
    if not hst_list:
        raise FileNotFoundError(
            f"No HST F814W drz FITS found in {DATA_RAW / galaxy}/"
        )
    hst_path = hst_list[0]
    print(f"\n{'=' * 60}")
    print(f"MGE fit: {galaxy}")
    print(f"{'=' * 60}")
    print(f"Image: {hst_path}")

    # --- 2. Read image ---
    image, pixscale, photflam, centrwv = read_hst_image(hst_path)
    ny, nx = image.shape
    print(f"Size: {nx}×{ny}, Scale: {pixscale:.4f} arcsec/pix")
    print(f"PHOTFLAM: {photflam:.4e}, CENTRWV: {centrwv:.1f} Å")

    # Trim
    mx = int(nx * trim_margin)
    my = int(ny * trim_margin)
    image_trim = image[my:ny - my, mx:nx - mx]
    print(f"Trimmed: {image_trim.shape}")

    # Background level
    bg = background_level(image_trim)
    print(f"Background cutoff (top 99% flux): {bg:.3f}")

    # --- 3. Find galaxy ---
    sec = find_galaxy(image_trim, plot=False, fraction=0.1)
    print(f"Center: ({sec.xmed:.1f}, {sec.ymed:.1f}) pix")
    print(f"Ellipticity: {sec.eps:.4f}, PA: {sec.pa:.2f}°")

    # Pixel-to-arcsec offset for plotting
    x_arcsec = (np.arange(image_trim.shape[1]) - sec.ymed) * pixscale
    y_arcsec = (np.arange(image_trim.shape[0]) - sec.xmed) * pixscale

    # --- 4. Sector photometry ---
    pho = sectors_photometry(
        image_trim, sec.eps, sec.theta, sec.xmed, sec.ymed,
        plot=False,
    )
    print(f"Sectors: {len(pho.radius)} points")

    # --- 5. PSF ---
    if fwhm is None:
        fwhm = 0.13  # typical WFPC2 PC F814W
        print(f"PSF FWHM: {fwhm:.2f} arcsec (default for WFPC2 F814W)")
    else:
        print(f"PSF FWHM: {fwhm:.2f} arcsec (user)")
    sigmapsf = fwhm / (2.355 * pixscale)
    print(f"PSF sigma: {sigmapsf:.2f} pix")

    # --- 6. MGE fit ---
    mge = mge_fit_sectors_regularized(
        pho.radius, pho.angle, pho.counts, sec.eps,
        ngauss=ngauss,
        plot=False,
        scale=pixscale,
        sigmapsf=sigmapsf,
        linear=False,
        outer_slope=4,
        bulge_disk=False,
        qbounds=[0.02, 0.999],
    )

    total_counts, sigma_pix, q_obs = mge.sol
    ng = len(total_counts)

    # --- 7. Unit conversion ---
    sigma = sigma_pix * pixscale  # arcsec

    # AB zeropoint of this image: m_AB = -2.5 log10(counts) + ZP
    # f_nu = counts * photflam * centrwv² / c  (erg/cm²/s/Hz)
    # m_AB = -2.5 log10(f_nu) - 48.6
    # So ZP = -2.5 log10(photflam * centrwv² / c) - 48.6
    zp_ab = flux_to_abmag(photflam * centrwv ** 2 / C_AA)

    # Surface brightness per Gaussian (counts/s/arcsec²)
    sb_counts = total_counts / (2 * np.pi * sigma ** 2 * q_obs)
    # → AB mag/arcsec²
    mu_ab = flux_to_abmag(sb_counts * photflam * centrwv ** 2 / C_AA)
    # → L⊙ pc⁻²
    surf = abmag_to_solar(mu_ab, M_SUN_F814W)

    # --- 8. Print results ---
    print(f"\n{'I (L⊙/pc²)':>15s}  {'sigma (arcsec)':>15s}  {'q':>8s}")
    print("-" * 42)
    for i in range(ng):
        print(f"{surf[i]:>15.1f}  {sigma[i]:>15.3f}  {q_obs[i]:>8.3f}")

    # --- 9. Save ECSV ---
    table = Table({
        "I": surf,
        "sigma": sigma,
        "q": q_obs,
        "pa_twist": np.zeros(ng),
    })
    ecsv_path = out_dir / "mge_f814w.ecsv"
    ascii.write(table, ecsv_path, format="ecsv", overwrite=True)
    print(f"\nSaved: {ecsv_path}")

    # --- 10. Diagnostic plots ---

    # 10a. Gaussian components
    plot_components(surf, sigma, q_obs,
                    str(out_dir / "mge_f814w_components.png"))

    # 10b. Contour overlay (mge_print_contours: xc=row, yc=col)
    from mgefit.mge_print_contours import mge_print_contours
    plt.figure(figsize=(8, 8))
    mge_print_contours(
        image_trim, sec.theta, sec.xmed, sec.ymed, mge.sol,
        minlevel=bg * 2,
        sigmapsf=sigmapsf,
        normpsf=[1.0],
        scale=pixscale,
    )
    plt.savefig(str(out_dir / "mge_f814w_contours.png"), dpi=150)
    plt.close()
    print(f"  Saved: {out_dir / 'mge_f814w_contours.png'}")

    # 10c. Radial profile
    pho_radius = pho.radius * pixscale
    pho_sb = pho.counts / pixscale ** 2
    pho_mu = flux_to_abmag(pho_sb * photflam * centrwv ** 2 / C_AA)
    pho_surf = abmag_to_solar(pho_mu, M_SUN_F814W)

    radi = np.logspace(np.log10(pho_radius.min() * 0.5),
                       np.log10(pho_radius.max() * 2), 300)
    mge_sb = np.zeros_like(radi)
    for f, s, q in zip(surf, sigma, q_obs):
        mge_sb += f * np.exp(-radi ** 2 / (2 * s ** 2))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(radi, mge_sb, "k-", label="MGE fit", lw=2)
    ax.scatter(pho_radius, pho_surf, s=5, c="C0", alpha=0.4, label="Data")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(pho_radius.min() * 0.5, pho_radius.max() * 1.5)
    ymin = np.percentile(pho_surf[pho_surf > 0], 5)
    ymax = np.percentile(pho_surf, 99.5)
    ax.set_ylim(ymin * 0.3, ymax * 2)
    ax.set_xlabel("Radius (arcsec)")
    ax.set_ylabel(r"$I\ (L_\odot\ \mathrm{pc}^{-2})$")
    ax.legend()
    plt.tight_layout()
    plt.savefig(str(out_dir / "mge_f814w_radial.png"), dpi=150)
    plt.close()
    print(f"  Saved: {out_dir / 'mge_f814w_radial.png'}")

    # 10d. Compare with existing MGE (different data sources — not directly comparable)
    existing = out_dir / "mge.ecsv"
    if existing.exists():
        old = ascii.read(existing)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.errorbar(sigma * (q_obs + 1) / 2, surf,
                    xerr=sigma * (1 - q_obs) / 2,
                    fmt="o", capsize=4, label=f"F814W (HST PC)")
        ax.errorbar(old["sigma"] * (old["q"] + 1) / 2, old["I"],
                    xerr=old["sigma"] * (1 - old["q"]) / 2,
                    fmt="s", capsize=4, label=f"Existing (SAURON/SDSS r)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        # Set limits to show the overlapping range
        all_i = np.concatenate([surf, old["I"]])
        all_s = np.concatenate([sigma, old["sigma"]])
        ax.set_xlim(all_s.min() * 0.5, all_s.max() * 2)
        ax.set_ylim(all_i.min() * 0.3, all_i.max() * 3)
        ax.set_xlabel("sigma (arcsec)")
        ax.set_ylabel(r"$I\ (L_\odot\ \mathrm{pc}^{-2})$")
        ax.set_title("MGE comparison (different instruments — offset expected)")
        ax.legend()
        plt.tight_layout()
        plt.savefig(str(out_dir / "mge_f814w_comparison.png"), dpi=150)
        plt.close()
        print(f"  Saved: {out_dir / 'mge_f814w_comparison.png'}")

    print(f"\nDone. Results in {out_dir}/")
    return surf, sigma, q_obs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit MGE from HST F814W image"
    )
    parser.add_argument("galaxy", help="Galaxy name (e.g. NGC4621)")
    parser.add_argument("--fwhm", type=float, default=None,
                        help="PSF FWHM in arcsec (default: 0.13 for WFPC2 F814W)")
    parser.add_argument("--ngauss", type=int, default=15,
                        help="Number of Gaussian components (default 15)")
    parser.add_argument("--trim", type=float, default=0.05,
                        help="Edge trim fraction (default 0.05)")
    args = parser.parse_args()

    fit_mge_hst(args.galaxy, fwhm=args.fwhm, ngauss=args.ngauss,
                trim_margin=args.trim)
