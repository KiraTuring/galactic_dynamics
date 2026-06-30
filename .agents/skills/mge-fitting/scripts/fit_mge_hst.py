#!/usr/bin/env python3
"""
Fit Multi-Gaussian Expansion (MGE) from an HST WFPC2 image using mgefit.

Pipeline:
    1. Read HST F814W drizzled image (SCI extension)
    2. find_galaxy() — locate center, ellipticity, PA
    3. sectors_photometry() — radial surface brightness in sectors
    4. mge_fit_sectors_regularized() — fit Gaussians with regularization (Cappellari 2002)
    5. Convert units: counts → AB mag/arcsec² → L⊙ pc⁻²
    6. Save results + diagnostic plots

Usage:
    python fit_mge_hst.py <galaxy>
    python fit_mge_hst.py NGC4621
    python fit_mge_hst.py NGC4621 --fwhm 0.13 --ngauss 15
    python fit_mge_hst.py NGC5845                # uses F555W if F814W unavailable

Output (data/processed/{galaxy}/):
    mge_{filter_name}.ecsv            — MGE parameters
    mge_{filter_name}_components.png  — Gaussian ellipses + radial profile
    mge_{filter_name}_contours.png    — Contour overlay on image
    mge_{filter_name}_radial.png      — Radial profile vs fit
    mge_{filter_name}_comparison.png  — Comparison with existing MGE (if found)

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

# Use JAM's prep.read_hst_image for WHT-masked image
_JAM_PATH = Path(__file__).resolve().parent.parent.parent.parent / "JAM"
sys.path.insert(0, str(_JAM_PATH))
from jam_fit.prep import read_hst_image as read_hst_wht


# --- Physical constants ---
C_AA = 2.99792458e18          # speed of light in Å/s
ARCSEC_PER_RAD = 206265.0     # arcsec/rad
AB_ZP = -48.6                 # AB magnitude zeropoint
M_SUN_AB = {            # Solar absolute AB magnitudes (Willmer 2018)
    "F555W": 4.82,      # V-band
    "F606W": 4.77,      # wide V
    "F702W": 4.71,      # R-band
    "F814W": 4.56,      # I-band
    "F160W": 5.36,      # H-band (NICMOS)
}
M_SUN_DEFAULT = 4.56    # fallback

# --- Paths ---
_JAM_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
DATA_RAW = _JAM_ROOT / "data" / "raw"
DATA_PROC = _JAM_ROOT / "data" / "processed"


def read_phot_cal(path):
    """Read PHOTFLAM / CENTRWV / filter from HST FITS for unit conversion."""
    with fits.open(path) as hdu:
        sci = hdu["SCI"]
        photflam = sci.header.get("PHOTFLAM")
        centrwv = hdu[0].header.get("CENTRWV", 8012.0)
        filt = hdu[0].header.get("FILTNAM1", "").strip().upper()
    return photflam, centrwv, filt


def background_level(image, fraction=0.99):
    """Estimate background threshold from cumulative flux."""
    flat = np.sort(image.flatten())[::-1]
    cum = np.cumsum(flat)
    cum /= cum[-1]
    return flat[cum <= fraction][-1]


def estimate_sky(image, margin=0.1):
    """Estimate sky level from image corners (mode of edge pixels)."""
    ny, nx = image.shape
    my, mx = int(ny * margin), int(nx * margin)
    corners = np.concatenate([
        image[:my, :mx].ravel(),
        image[:my, -mx:].ravel(),
        image[-my:, :mx].ravel(),
        image[-my:, -mx:].ravel(),
    ])
    # Remove outliers
    lo, hi = np.percentile(corners, [5, 95])
    sky = corners[(corners >= lo) & (corners <= hi)].mean()
    return sky


def flux_to_abmag(f_nu):
    """Flux density (erg/cm²/s/Hz) → AB mag."""
    return -2.5 * np.log10(f_nu) - 48.6


def abmag_to_solar(mu_ab, m_sun):
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


def fit_mge_hst(galaxy, fwhm=None, ngauss=20, trim_margin=0.05, sky=None):
    """
    Run full MGE fitting pipeline on an HST F814W image.

    Parameters
    ----------
    galaxy : str
        Galaxy name (e.g. "NGC4621").
    fwhm : float or None
        PSF FWHM in arcsec. Default (None): 0.13 arcsec for WFPC2 F814W.
    ngauss : int
        Number of Gaussian components (default 20).
    trim_margin : float
        Fractional margin to trim from image edges (default 0.05).
    sky : float or None
        Sky background in counts/s/pix. If None, estimate automatically.
        Sky-subtracted image is used for sector photometry.
    """
    out_dir = DATA_PROC / galaxy
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Locate HST image ---
    patterns = ["hst_*f814w*drz.fits", "hst_*f555w*drz.fits",
                "hst_*f702w*drz.fits", "hst_*drz.fits"]
    hst_list = []
    for p in patterns:
        hst_list = sorted(glob.glob(str(DATA_RAW / galaxy / p)))
        if hst_list:
            break
    if not hst_list:
        raise FileNotFoundError(
            f"No HST drz FITS found for {galaxy} in {DATA_RAW / galaxy}/"
        )
    hst_path = hst_list[0]
    print(f"\n{'=' * 60}")
    print(f"MGE fit: {galaxy}")
    print(f"{'=' * 60}")
    print(f"Image: {hst_path}")
    # --- 2. Read image (with WHT mask from prep.py) ---
    image_ma, wht, extent, pixscale = read_hst_wht(hst_path)
    photflam, centrwv, filter_name = read_phot_cal(hst_path)
    ny, nx = image_ma.shape
    print(f"Size: {nx}×{ny}, Scale: {pixscale:.4f} arcsec/pix")
    print(f"PHOTFLAM: {photflam:.4e}, CENTRWV: {centrwv:.1f} Å, Filter: {filter_name}")

    # Select solar magnitude by filter
    m_sun = M_SUN_AB.get(filter_name, M_SUN_DEFAULT)

    # Trim
    mx = int(nx * trim_margin)
    my = int(ny * trim_margin)
    image_trim = image_ma[my:ny - my, mx:nx - mx].data.astype(np.float64)
    wht_trim = wht[my:ny - my, mx:nx - mx]
    print(f"Trimmed: {image_trim.shape}")

    # Bad pixels mask from WHT
    badpixels = wht_trim <= 0
    nbad = badpixels.sum()
    if nbad > 0:
        print(f"Bad pixels masked (WHT <= 0): {nbad} ({100 * nbad / badpixels.size:.2f}%)")

    # --- 3. Sky background ---
    if sky is None:
        sky_est = estimate_sky(image_trim)
        sky_est = max(sky_est, 0.0)
        print(f"Sky background (estimated from corners): {sky_est:.4f} counts/s/pix")
    else:
        sky_est = sky
        print(f"Sky background (user): {sky_est:.4f} counts/s/pix")

    sky_est = float(sky_est)
    image_sky = image_trim - sky_est

    # Background level (for sector photometry cutoff)
    bg = background_level(image_sky)
    print(f"Background cutoff (top 99% flux): {bg:.3f}")

    # --- 4. Find galaxy (use sky-subtracted data; replace bad pixels with 0) ---
    image_clean = np.where(badpixels, 0.0, image_sky)
    sec = find_galaxy(image_clean, plot=False, fraction=0.1)
    print(f"Center: ({sec.xmed:.1f}, {sec.ymed:.1f}) pix")
    print(f"Ellipticity: {sec.eps:.4f}, PA: {sec.pa:.2f}°")

    # Pixel-to-arcsec offset for plotting
    x_arcsec = (np.arange(image_clean.shape[1]) - sec.ymed) * pixscale
    y_arcsec = (np.arange(image_clean.shape[0]) - sec.xmed) * pixscale

    # --- 5. Sector photometry (with badpixels mask) ---
    pho = sectors_photometry(
        image_clean, sec.eps, sec.theta, sec.xmed, sec.ymed,
        badpixels=badpixels,
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

    # --- 6. MGE fit (regularized, follows Cappellari 2002 / 2006) ---
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

    # Solar magnitude for this filter
    m_sun = M_SUN_AB.get(filter_name, M_SUN_DEFAULT)
    print(f"Filter: {filter_name}, CENTRWV={centrwv:.0f} Å, M⊙_AB={m_sun:.2f}")

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
    surf = abmag_to_solar(mu_ab, m_sun)

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
    ecsv_path = out_dir / f"mge_{filter_name}.ecsv"
    ascii.write(table, ecsv_path, format="ecsv", overwrite=True)
    print(f"\nSaved: {ecsv_path}")

    # --- 10. Diagnostic plots ---

    # 10a. Gaussian components
    plot_components(surf, sigma, q_obs,
                    str(out_dir / f"mge_{filter_name}_components.png"))

    # 10b. Contour overlay (full FOV + zoom, manual extent)
    ny_s, nx_s = image_sky.shape
    y_pix = np.arange(ny_s) - sec.xmed
    x_pix = np.arange(nx_s) - sec.ymed
    ext_arcsec = [x_pix[0] * pixscale, x_pix[-1] * pixscale,
                  y_pix[0] * pixscale, y_pix[-1] * pixscale]

    # PSF-convolved MGE model image (rotation matches _gauss2d_mge)
    model = np.zeros_like(image_sky)
    ang_rad = np.radians(sec.theta - 90)
    cos_a, sin_a = np.cos(ang_rad), np.sin(ang_rad)
    yy, xx = np.meshgrid(y_pix, x_pix, indexing='ij')
    # x' = row*cos + col*sin  (major axis),  y' = col*cos - row*sin  (minor axis)
    x_maj = yy * cos_a + xx * sin_a
    x_min = xx * cos_a - yy * sin_a
    for w, s, q in zip(total_counts, sigma_pix, q_obs):
        sx = np.sqrt(s ** 2 + sigmapsf ** 2)
        sy = np.sqrt((s * q) ** 2 + sigmapsf ** 2)
        g = np.exp(-0.5 * ((x_maj / sx) ** 2 + (x_min / sy) ** 2))
        model += w / (2 * np.pi * sx * sy) * g

    peak = image_sky[int(round(sec.xmed)), int(round(sec.ymed))]
    nlevels = 12
    lvls = np.logspace(np.log10(max(bg, 1e-3)),
                       np.log10(peak * 0.9), nlevels)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, zoom in zip(axes, [None, 3]):
        ax.contour(image_sky, levels=lvls, colors='k', linewidths=1,
                    extent=ext_arcsec)
        ax.contour(model, levels=lvls, colors='r', linewidths=1,
                    extent=ext_arcsec, linestyles='--')
        ax.axhline(0, color='gray', lw=0.5, alpha=0.5)
        ax.axvline(0, color='gray', lw=0.5, alpha=0.5)
        ax.set_xlabel('arcsec')
        ax.set_ylabel('arcsec')
        ax.set_aspect('equal')
        if zoom:
            ax.set_xlim(-zoom, zoom)
            ax.set_ylim(-zoom, zoom)
            ax.set_title(f'Zoom {zoom}″ × {zoom}″')
        else:
            ax.set_title(f'Full FOV ({nx_s * pixscale:.0f}″ × {ny_s * pixscale:.0f}″)')
    plt.tight_layout()
    plt.savefig(str(out_dir / f"mge_{filter_name}_contours.png"), dpi=150)
    plt.close()
    print(f"  Saved: {out_dir / f'mge_{filter_name}_contours.png'}")

    # 10c. Radial profile (major + minor axis)
    pho_radius = pho.radius * pixscale
    pho_sb = pho.counts / pixscale ** 2
    pho_mu = flux_to_abmag(pho_sb * photflam * centrwv ** 2 / C_AA)
    pho_surf = abmag_to_solar(pho_mu, m_sun)

    # Split by angle: major axis (|angle| < 15°) vs minor axis (|angle - 90| < 15°)
    angle = np.asarray(pho.angle)
    major_sel = (np.abs(angle) < 15) | (np.abs(angle - 180) < 15)
    minor_sel = np.abs(angle - 90) < 15

    # MGE model along major axis (q=1) and minor axis (q=q_obs)
    radi = np.logspace(np.log10(pho_radius.min() * 0.5),
                       np.log10(pho_radius.max() * 2), 300)
    mge_major = np.zeros_like(radi)
    mge_minor = np.zeros_like(radi)
    for f, s, q in zip(surf, sigma, q_obs):
        mge_major += f * np.exp(-radi ** 2 / (2 * s ** 2))
        mge_minor += f * np.exp(-radi ** 2 / (2 * (s * q) ** 2))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(radi, mge_major, "k-", label="MGE major axis", lw=2)
    ax.plot(radi, mge_minor, "k--", label="MGE minor axis", lw=2)
    ax.scatter(pho_radius[major_sel], pho_surf[major_sel],
               s=8, c="C0", alpha=0.4, label="Data (major)", marker="o")
    ax.scatter(pho_radius[minor_sel], pho_surf[minor_sel],
               s=8, c="C1", alpha=0.4, label="Data (minor)", marker="s")
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
    plt.savefig(str(out_dir / f"mge_{filter_name}_radial.png"), dpi=150)
    plt.close()
    print(f"  Saved: {out_dir / f'mge_{filter_name}_radial.png'}")

    # 10d. Compare with existing MGE (different data sources — not directly comparable)
    existing = out_dir / "mge.ecsv"
    if existing.exists():
        old = ascii.read(existing)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.errorbar(sigma * (q_obs + 1) / 2, surf,
                    xerr=sigma * (1 - q_obs) / 2,
                    fmt="o", capsize=4, label=f"{filter_name.upper()} (HST PC)")
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
        plt.savefig(str(out_dir / f"mge_{filter_name}_comparison.png"), dpi=150)
        plt.close()
        print(f"  Saved: {out_dir / f'mge_{filter_name}_comparison.png'}")

    print(f"\nDone. Results in {out_dir}/")
    return surf, sigma, q_obs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit MGE from HST F814W image"
    )
    parser.add_argument("galaxy", help="Galaxy name (e.g. NGC4621)")
    parser.add_argument("--fwhm", type=float, default=None,
                        help="PSF FWHM in arcsec (default: 0.13 for WFPC2 F814W)")
    parser.add_argument("--ngauss", type=int, default=20,
                        help="Number of Gaussian components (default 20)")
    parser.add_argument("--trim", type=float, default=0.05,
                        help="Edge trim fraction (default 0.05)")
    parser.add_argument("--sky", type=float, default=None,
                        help="Sky background in counts/s/pix (default: auto)")
    args = parser.parse_args()

    fit_mge_hst(args.galaxy, fwhm=args.fwhm, ngauss=args.ngauss,
                trim_margin=args.trim, sky=args.sky)
