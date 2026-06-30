#!/usr/bin/env python3
"""
Fit Multi-Gaussian Expansion (MGE) from an HST WFPC2 image using mgefit.

Pipeline:
    1. find_galaxy() — locate center, ellipticity, PA
    2. sectors_photometry() — radial surface brightness in sectors
    3. mge_fit_sectors_regularized() — fit Gaussians (Cappellari 2002)
    4. Convert units: counts → AB mag/arcsec² → L⊙ pc⁻²
    5. Save results + diagnostic plots

Usage:
    python fit_mge_hst.py <galaxy>
    python fit_mge_hst.py NGC4621
    python fit_mge_hst.py NGC4621 --fwhm 0.13 --ngauss 15
    python fit_mge_hst.py NGC5845

Output (data/processed/{galaxy}/):
    mge_{filter_name}.ecsv            — MGE parameters
    mge_{filter_name}_components.png  — Gaussian ellipses + radial profile
    mge_{filter_name}_contours.png    — Contour overlay on image
    mge_{filter_name}_radial.png      — Radial profile vs fit
    mge_{filter_name}_comparison.png  — Comparison with existing MGE (if found)
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

_JAM_PATH = Path(__file__).resolve().parent.parent.parent.parent / "JAM"
sys.path.insert(0, str(_JAM_PATH))
from jam_fit.prep import read_hst_image as read_hst_wht


# --- Physical constants ---
C_AA = 2.99792458e18
ARCSEC_PER_10PC = 20626.5  # arcsec per pc at 10 pc = 206265/10
AB_ZP = -48.6
M_SUN_AB = {
    "F555W": 4.82, "F606W": 4.77, "F702W": 4.71,
    "F814W": 4.56, "F160W": 5.36,
}
M_SUN_DEFAULT = 4.56

_JAM_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
DATA_RAW = _JAM_ROOT / "data" / "raw"
DATA_PROC = _JAM_ROOT / "data" / "processed"


# ---- helpers ----

def flux_to_abmag(f_nu):
    return -2.5 * np.log10(f_nu) - 48.6


def abmag_to_solar(mu_ab, m_sun):
    return ARCSEC_PER_10PC ** 2 * 10 ** (-0.4 * (mu_ab - m_sun))


def background_level(image, fraction=0.99):
    flat = np.sort(image.flatten())[::-1]
    cum = np.cumsum(flat)
    cum /= cum[-1]
    return flat[cum <= fraction][-1]


def estimate_sky(image, margin=0.1):
    ny, nx = image.shape
    my, mx = int(ny * margin), int(nx * margin)
    corners = np.concatenate([
        image[:my, :mx].ravel(), image[:my, -mx:].ravel(),
        image[-my:, :mx].ravel(), image[-my:, -mx:].ravel(),
    ])
    lo, hi = np.percentile(corners, [5, 95])
    return corners[(corners >= lo) & (corners <= hi)].mean()


# ---- I/O helpers ----

def _read_phot_cal(path):
    with fits.open(path) as hdu:
        sci = hdu["SCI"]
        photflam = sci.header.get("PHOTFLAM")
        centrwv = hdu[0].header.get("CENTRWV", 8012.0)
        filt = hdu[0].header.get("FILTNAM1", "").strip().upper()
    return photflam, centrwv, filt


def _find_hst_image(galaxy):
    patterns = ["hst_*f814w*drz.fits", "hst_*f555w*drz.fits",
                "hst_*f702w*drz.fits", "hst_*drz.fits"]
    for p in patterns:
        lst = sorted(glob.glob(str(DATA_RAW / galaxy / p)))
        if lst:
            return lst[0]
    raise FileNotFoundError(f"No HST drz image for {galaxy} in {DATA_RAW / galaxy}/")


# ---- MGE unit conversion ----

def _mge_counts_to_surf(total_counts, sigma_pix, q_obs, pixscale,
                        photflam, centrwv, m_sun):
    sigma = sigma_pix * pixscale
    sb_counts = total_counts / (2 * np.pi * sigma ** 2 * q_obs)
    f_nu = sb_counts * photflam * centrwv ** 2 / C_AA
    mu_ab = flux_to_abmag(f_nu)
    surf = abmag_to_solar(mu_ab, m_sun)
    return surf, sigma


def _compute_mge_model(total_counts, sigma_pix, q_obs,
                       sec, sigmapsf, ny, nx):
    y_pix = np.arange(ny) - sec.xmed
    x_pix = np.arange(nx) - sec.ymed
    ang_rad = np.radians(sec.theta - 90)
    cos_a, sin_a = np.cos(ang_rad), np.sin(ang_rad)
    yy, xx = np.meshgrid(y_pix, x_pix, indexing='ij')
    x_maj = yy * cos_a + xx * sin_a
    x_min = xx * cos_a - yy * sin_a

    model = np.zeros((ny, nx))
    for w, s, q in zip(total_counts, sigma_pix, q_obs):
        sx = np.sqrt(s ** 2 + sigmapsf ** 2)
        sy = np.sqrt((s * q) ** 2 + sigmapsf ** 2)
        g = np.exp(-0.5 * ((x_maj / sx) ** 2 + (x_min / sy) ** 2))
        model += w / (2 * np.pi * sx * sy) * g
    return model


# ---- plotting ----

def _plot_components(surf, sigma, q_obs, savepath):
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
                xerr=[sigma * (1 - q_obs) / 2], fmt="o", capsize=4)
    ax.set_xlabel("sigma (arcsec)")
    ax.set_ylabel(r"$I\ (L_\odot\ \mathrm{pc}^{-2})$")
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def _plot_contours(image, model, sec, bg, pixscale, savepath):
    ny_s, nx_s = image.shape
    x_pix = (np.arange(nx_s) - sec.ymed) * pixscale
    y_pix = (np.arange(ny_s) - sec.xmed) * pixscale
    ext = [x_pix[0], x_pix[-1], y_pix[0], y_pix[-1]]

    peak = image[int(round(sec.xmed)), int(round(sec.ymed))]
    lvls = np.logspace(np.log10(max(bg, 1e-3)),
                       np.log10(peak * 0.9), 12)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, zoom in zip(axes, [None, 3]):
        ax.contour(image, levels=lvls, colors='k', linewidths=1, extent=ext)
        ax.contour(model, levels=lvls, colors='r', linewidths=1,
                    extent=ext, linestyles='--')
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
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def _plot_radial(pho, surf, sigma, q_obs, pixscale, photflam, centrwv, m_sun, savepath):
    pho_radius = pho.radius * pixscale
    pho_sb = pho.counts / pixscale ** 2
    pho_fnu = pho_sb * photflam * centrwv ** 2 / C_AA
    pho_surf = abmag_to_solar(flux_to_abmag(pho_fnu), m_sun)

    angle = np.asarray(pho.angle)
    major_sel = (np.abs(angle) < 15) | (np.abs(angle - 180) < 15)
    minor_sel = np.abs(angle - 90) < 15

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
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def _plot_comparison(surf, sigma, q_obs, filter_name, out_dir):
    existing = out_dir / "mge.ecsv"
    if not existing.exists():
        return
    old = ascii.read(existing)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(sigma * (q_obs + 1) / 2, surf,
                xerr=sigma * (1 - q_obs) / 2,
                fmt="o", capsize=4, label=f"{filter_name} (HST PC)")
    ax.errorbar(old["sigma"] * (old["q"] + 1) / 2, old["I"],
                xerr=old["sigma"] * (1 - old["q"]) / 2,
                fmt="s", capsize=4, label="Existing (SAURON/SDSS r)")
    ax.set_xscale("log")
    ax.set_yscale("log")
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


# ---- pipeline ----

def fit_mge_hst(galaxy, fwhm=None, ngauss=20, trim_margin=0.05, sky=None):
    out_dir = DATA_PROC / galaxy
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Locate & read image
    hst_path = _find_hst_image(galaxy)
    print(f"\n{'=' * 60}\nMGE fit: {galaxy}\n{'=' * 60}")
    print(f"Image: {hst_path}")
    image_ma, wht, _, pixscale = read_hst_wht(hst_path)
    photflam, centrwv, filter_name = _read_phot_cal(hst_path)
    m_sun = M_SUN_AB.get(filter_name, M_SUN_DEFAULT)
    print(f"Size: {image_ma.shape[1]}×{image_ma.shape[0]}, "
          f"Scale: {pixscale:.4f}\", Filter: {filter_name}")

    # 2. Trim + sky subtract + mask bad pixels
    mx = int(image_ma.shape[1] * trim_margin)
    my = int(image_ma.shape[0] * trim_margin)
    image_trim = image_ma[my:-my, mx:-mx].data.astype(np.float64)
    wht_trim = wht[my:-my, mx:-mx]
    badpixels = wht_trim <= 0
    print(f"Trimmed: {image_trim.shape}", end="")
    if badpixels.any():
        print(f", bad pixels: {badpixels.sum()} ({100 * badpixels.sum() / badpixels.size:.1f}%)")
    else:
        print()

    if sky is None:
        sky = max(estimate_sky(image_trim), 0.0)
        print(f"Sky (auto): {sky:.4f} counts/s/pix")
    else:
        print(f"Sky (user): {sky:.4f} counts/s/pix")
    image_clean = np.where(badpixels, 0.0, image_trim - sky)
    bg = background_level(image_clean)
    print(f"Background cutoff: {bg:.3f}")

    # 3. Find galaxy
    sec = find_galaxy(image_clean, plot=False, fraction=0.1)
    print(f"Center: ({sec.xmed:.1f}, {sec.ymed:.1f}) pix, "
          f"eps={sec.eps:.4f}, PA={sec.pa:.1f}°")

    # 4. Sector photometry
    pho = sectors_photometry(
        image_clean, sec.eps, sec.theta, sec.xmed, sec.ymed,
        badpixels=badpixels, plot=False)
    print(f"Sectors: {len(pho.radius)} points")

    # 5. PSF
    if fwhm is None:
        fwhm = 0.13
        print(f"PSF FWHM: {fwhm:.2f}\" (default)")
    else:
        print(f"PSF FWHM: {fwhm:.2f}\" (user)")
    sigmapsf = fwhm / (2.355 * pixscale)
    print(f"PSF sigma: {sigmapsf:.2f} pix")

    # 6. MGE fit
    mge = mge_fit_sectors_regularized(
        pho.radius, pho.angle, pho.counts, sec.eps,
        ngauss=ngauss, plot=False, scale=pixscale, sigmapsf=sigmapsf,
        linear=False, outer_slope=4, bulge_disk=False,
        qbounds=[0.02, 0.999])
    total_counts, sigma_pix, q_obs = mge.sol
    ng = len(total_counts)

    # 7. Unit conversion
    surf, sigma = _mge_counts_to_surf(
        total_counts, sigma_pix, q_obs, pixscale,
        photflam, centrwv, m_sun)

    print(f"\n{'I (L⊙/pc²)':>15s}  {'sigma (arcsec)':>15s}  {'q':>8s}")
    print("-" * 38)
    for i in range(ng):
        print(f"{surf[i]:>15.1f}  {sigma[i]:>10.3f}  {q_obs[i]:>8.3f}")

    # 8. Save ECSV
    table = Table({"I": surf, "sigma": sigma, "q": q_obs,
                   "pa_twist": np.zeros(ng)})
    ecsv_path = out_dir / f"mge_{filter_name}.ecsv"
    ascii.write(table, ecsv_path, format="ecsv", overwrite=True)
    print(f"\nSaved: {ecsv_path}")

    # 9. Diagnostics
    prefix = out_dir / f"mge_{filter_name}"

    _plot_components(surf, sigma, q_obs, str(prefix) + "_components.png")

    model = _compute_mge_model(
        total_counts, sigma_pix, q_obs, sec, sigmapsf,
        image_clean.shape[0], image_clean.shape[1])
    #_plot_contours(image_clean, model, sec, bg, pixscale,
    #               str(prefix) + "_contours.png")
    # Use mge_print_contours (handles central pixel integration correctly)
    from mgefit.mge_print_contours import mge_print_contours

    # mge_print_contours uses img[xc, yc] convention (xc=row, yc=col)
    r0, c0 = int(round(sec.xmed)), int(round(sec.ymed))
    if image_clean[r0, c0] <= 0:
        roi = image_clean[max(0,r0-3):r0+4, max(0,c0-3):c0+4]
        dr, dc = np.unravel_index(roi.argmax(), roi.shape)
        r0, c0 = max(0,r0-3) + dr, max(0,c0-3) + dc
        print(f"  Center pixel masked, using nearby peak at row={r0} col={c0} = {image_clean[r0, c0]:.2f}")

    def _mge_contour(ax, zoom):
        mge_print_contours(image_clean, sec.theta, r0, c0, mge.sol,
                           minlevel=max(bg, 1e-3), sigmapsf=sigmapsf,
                           normpsf=[1.0], scale=pixscale, magstep=0.5)
        if zoom:
            ax.set_xlim(-zoom, zoom)
            ax.set_ylim(-zoom, zoom)
            ax.set_title(f'Zoom {zoom}″ × {zoom}″')
        else:
            ax.set_title(f'Full FOV ({image_clean.shape[1]*pixscale:.0f}″ × {image_clean.shape[0]*pixscale:.0f}″)')

    fig = plt.figure(figsize=(16, 7))
    ax1 = fig.add_subplot(121)
    plt.sca(ax1)
    _mge_contour(ax1, None)
    ax2 = fig.add_subplot(122)
    plt.sca(ax2)
    _mge_contour(ax2, 3)
    plt.tight_layout()
    plt.savefig(str(prefix) + "_contours.png", dpi=150)
    plt.close()
    print(f"  Saved: {prefix}_contours.png")

    _plot_radial(pho, surf, sigma, q_obs, pixscale, photflam, centrwv, m_sun,
                 str(prefix) + "_radial.png")

    _plot_comparison(surf, sigma, q_obs, filter_name, out_dir)

    print(f"\nDone. Results in {out_dir}/")
    return surf, sigma, q_obs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit MGE from HST image")
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
