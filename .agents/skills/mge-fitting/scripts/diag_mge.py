#!/usr/bin/env python3
"""
Diagnose MGE fit quality along major and minor axes.

Usage:
    python diag_mge.py <galaxy> --filter F814W [--fwhm 0.13] [--tag os2]

Output: tables of r, data, model, residual for major/minor axis.
"""

import argparse, sys
from pathlib import Path

import numpy as np
from astropy.io import fits, ascii

_JAM_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
DATA_RAW = _JAM_ROOT / "data" / "raw"
DATA_PROC = _JAM_ROOT / "data" / "processed"

sys.path.insert(0, str(_JAM_ROOT / "JAM"))
from jam_fit.prep import read_hst_image as read_hst_wht

from mgefit.find_galaxy import find_galaxy
from mgefit.sectors_photometry import sectors_photometry

from mgefit.mge_fit_sectors import mge_fit_sectors


def load_mge_params(path):
    tbl = ascii.read(path)
    # Convert I (L⊙/pc²) back to total_counts (counts/s/pixel integrated)
    # We need pixscale, photflam, etc. to reverse the conversion.
    # Better: just read the ECSV and use it directly for model prediction
    # in physical units, then compare with data in physical units.
    return tbl


def main(galaxy, filt="F814W", fwhm=None, ngauss=20, tag="", outer_slope=4):
    out_dir = DATA_PROC / galaxy

    # 1. Locate image
    pattern = f"hst_*{filt.lower()}*drz.fits"
    import glob
    lst = glob.glob(str(DATA_RAW / galaxy / pattern))
    if not lst:
        print(f"No {filt} image found for {galaxy}")
        return
    hst_path = lst[0]
    print(f"Image: {hst_path}")

    # 2. Read
    image_ma, wht, _, pixscale = read_hst_wht(hst_path)
    if fwhm is None:
        fwhm = 0.13
    sigmapsf = fwhm / (2.355 * pixscale)

    ny, nx = image_ma.shape
    image_raw = image_ma.data.astype(np.float64)
    wht_img = wht
    badpixels = wht_img <= 0
    image_clean = np.where(badpixels, 0.0, image_raw)

    # 3. Find galaxy
    sec = find_galaxy(image_clean, plot=False, fraction=0.1)
    print(f"Center: ({sec.xmed:.1f}, {sec.ymed:.1f}), eps={sec.eps:.4f}, theta={sec.theta:.1f}")

    # 4. Sector photometry
    pho = sectors_photometry(
        image_clean, sec.eps, sec.theta, sec.xmed, sec.ymed,
        badpixels=badpixels, plot=False)
    r_arcsec = pho.radius * pixscale
    angle = np.asarray(pho.angle)
    counts = np.asarray(pho.counts)

    # 5. MGE fit (same params as original; uses mge_fit_sectors for yfit access)
    mge = mge_fit_sectors(
        pho.radius, pho.angle, pho.counts, sec.eps,
        ngauss=ngauss, plot=False, scale=pixscale, sigmapsf=sigmapsf,
        linear=False, outer_slope=outer_slope, bulge_disk=False,
        qbounds=[0.02, 0.999])

    print(f"Fit: {mge.sol.shape[1]} components, chi2={mge.chi2:.1f}")

    # 6. Model prediction at data points
    model_counts = mge.yfit
    residual = 1 - model_counts / counts

    # 7. Split by axis
    major = (np.abs(angle) < 15) | (np.abs(angle - 180) < 15)
    minor = np.abs(angle - 90) < 15

    print(f"\n{'='*70}")
    print(f"  Residuals along MAJOR axis ({major.sum()} pts)")
    print(f"{'='*70}")
    hdr = f"{'r(arcsec)':>12s}  {'data':>10s}  {'model':>10s}  {'data-model':>10s}  {'resid':>7s}"
    print(f"  {hdr}")
    for i in np.where(major)[0]:
        print(f"  {r_arcsec[i]:>7.3f}  {counts[i]:>10.4f}  {model_counts[i]:>10.4f}  "
              f"{counts[i]-model_counts[i]:>10.4f}  {residual[i]:>+7.4f}")

    # Binned summary
    print(f"\n{'='*70}")
    print(f"  Binned residuals — MAJOR AXIS")
    print(f"{'='*70}")
    h2 = f"{'r range':>12s}  {'n':>4s}  {'bias':>8s}  {'stdev':>7s}  {'chi2':>7s}"
    print(f"  {h2}")
    bins = np.logspace(np.log10(r_arcsec[major].min()*0.8),
                       np.log10(r_arcsec[major].max()*1.2), 12)
    for i in range(len(bins)-1):
        lo, hi = bins[i], bins[i+1]
        m = major & (r_arcsec >= lo) & (r_arcsec < hi)
        n = m.sum()
        if n < 2: continue
        e = residual[m]
        print(f"  {lo:5.2f}-{hi:5.2f}  {n:4d}  {e.mean():+8.4f}  {e.std():7.4f}  {e@e:7.2f}")

    print(f"\n{'='*70}")
    print(f"  Residuals along MINOR axis ({minor.sum()} pts)")
    print(f"{'='*70}")
    hdr = f"{'r(arcsec)':>12s}  {'data':>10s}  {'model':>10s}  {'data-model':>10s}  {'resid':>7s}"
    print(f"  {hdr}")
    for i in np.where(minor)[0]:
        print(f"  {r_arcsec[i]:>7.3f}  {counts[i]:>10.4f}  {model_counts[i]:>10.4f}  "
              f"{counts[i]-model_counts[i]:>10.4f}  {residual[i]:>+7.4f}")

    print(f"\n{'='*70}")
    print(f"  Binned residuals — MINOR AXIS")
    print(f"{'='*70}")
    h2 = f"{'r range':>12s}  {'n':>4s}  {'bias':>8s}  {'stdev':>7s}  {'chi2':>7s}"
    print(f"  {h2}")
    for i in range(len(bins)-1):
        lo, hi = bins[i], bins[i+1]
        m = minor & (r_arcsec >= lo) & (r_arcsec < hi)
        n = m.sum()
        if n < 2: continue
        e = residual[m]
        print(f"  {lo:5.2f}-{hi:5.2f}  {n:4d}  {e.mean():+8.4f}  {e.std():7.4f}  {e@e:7.2f}")

    print(f"\n  Overall chi2 = {mge.chi2:.1f}, ncomp = {mge.sol.shape[1]}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagnose MGE fit residuals")
    parser.add_argument("galaxy", help="Galaxy name")
    parser.add_argument("--filter", default="F814W", help="Filter (default F814W)")
    parser.add_argument("--fwhm", type=float, default=None, help="PSF FWHM")
    parser.add_argument("--ngauss", type=int, default=20, help="N Gaussians")
    parser.add_argument("--tag", default="", help="Output tag")
    parser.add_argument("--outer-slope", type=float, default=4, help="outer slope")
    args = parser.parse_args()
    main(args.galaxy, filt=args.filter, fwhm=args.fwhm,
         ngauss=args.ngauss, tag=args.tag, outer_slope=args.outer_slope)
