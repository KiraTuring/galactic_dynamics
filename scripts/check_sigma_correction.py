#!/usr/bin/env python3
"""Check SAURON sigma x1.05 correction by convolving OASIS sigma to SAURON resolution.

Convolves OASIS sigma with the differential PSF (SAURON minus OASIS in quadrature)
then compares 1D major-axis slit profiles:
  - OASIS sigma (original)
  - OASIS sigma (PSF-convolved)
  - SAURON sigma
  - SAURON sigma x 1.05

Output: results/JAM/{galaxy}/sigma_offset_check.pdf
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml
from pathlib import Path

# Resolve paths
_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = _SCRIPT_DIR.parent
_JAM = ROOT / "JAM"

sys.path.insert(0, str(_JAM))

from jam_fit.prep import read_oasis_kin, read_sauron_kin, vel_pa_corr, get_rms
from jam_fit.psf import bins_convolve


def diff_psf_sigma(oasis_psf, sauron_psf):
    """Differential PSF sigma using core (sigma1) components only.

    The broad Gaussian wings have negligible effect on sigma profiles.
    diff = sqrt(sauron_sigma1² - oasis_sigma1²)
    """
    s2_o = oasis_psf["sigma1"]**2
    s2_s = sauron_psf["sigma1"]**2
    return np.sqrt(max(0, s2_s - s2_o))


def slit_binnum(xbin, ybin, angle=0, slit=0.3, rmax=None):
    """Extract bin indices along a slit of given angle and width.

    Parameters
    ----------
    xbin, ybin : ndarray
        Bin centroid coordinates (PA-aligned, arcsec).
    angle : float
        Position angle of the slit in degrees.
    slit : float
        Half-width of the slit in arcsec.
    rmax : float or None
        Maximum radius to include.

    Returns
    -------
    slit_bin : ndarray
        Indices of bins within the slit, sorted by radius along slit.
    rbin : ndarray
        Radius along the slit for each selected bin.
    """
    theta = np.deg2rad(angle)
    mask = np.abs(ybin * np.cos(theta) - xbin * np.sin(theta)) < slit / 2
    slit_idx = np.where(mask)[0]
    rbin = xbin[slit_idx] * np.cos(theta) + ybin[slit_idx] * np.sin(theta)
    order = np.argsort(rbin)
    slit_idx = slit_idx[order]
    rbin = rbin[order]

    if rmax is not None:
        keep = np.abs(rbin) < rmax
        slit_idx = slit_idx[keep]
        rbin = rbin[keep]

    return slit_idx, rbin


def check_galaxy(galaxy, pa_oasis, release):
    """Run the sigma offset check for one galaxy."""
    print(f"\n{'='*50}\n  {galaxy}\n{'='*50}")

    # ---- 1. Load data ----
    raw = ROOT / "data" / "raw" / galaxy

    oasis_file = str(raw / "oasis" / f"kinematics_oasis_{galaxy}.fits")
    if not os.path.exists(oasis_file):
        oasis_file += ".gz"
    cube_file = str(raw / "sauron" / f"MS_{galaxy}_{release}_C2D.fits")
    kin_file  = str(raw / "sauron" / f"{galaxy}_{release}_idl.fits.gz")

    print(f"  OASIS: {oasis_file}")
    print(f"  SAURON cube: {cube_file}")
    print(f"  SAURON kin:  {kin_file}")

    # OASIS
    xb0_o, yb0_o, kins0_o, dkins_o, xpix_o, ypix_o, binnum_o, flux_o = read_oasis_kin(oasis_file)
    xbin_o, ybin_o, kins_o, pa_fit_o = vel_pa_corr(xb0_o, yb0_o, kins0_o, dkins_o, binnum_o, flux_o, pa=pa_oasis)
    rms_o, erms_o = get_rms(kins_o, dkins_o)

    # SAURON
    xb0_s, yb0_s, kins0_s, dkins0_s, xpix_s, ypix_s, binnum_s, flux_s = read_sauron_kin(cube_file, kin_file)
    xbin_s, ybin_s, kins_s, pa_fit_s = vel_pa_corr(xb0_s, yb0_s, kins0_s, dkins0_s, binnum_s, flux_s, pa=pa_oasis)
    rms_s, erms_s = get_rms(kins_s, dkins0_s)

    # ---- 2. Load PSF ----
    oasis_psf = yaml.safe_load(open(
        ROOT / "results" / "JAM" / galaxy / "psf_fit" / "psf_result.yaml"))
    sauron_psf = yaml.safe_load(open(
        ROOT / "results" / "JAM" / galaxy / "psf_fit_sauron" / "psf_result.yaml"))

    print(f"  OASIS PSF:  sigma=[{oasis_psf['sigma1']:.3f}, {oasis_psf['sigma2']:.3f}], w={oasis_psf['weight1']:.3f}")
    print(f"  SAURON PSF: sigma=[{sauron_psf['sigma1']:.3f}, {sauron_psf['sigma2']:.3f}], w={sauron_psf['weight1']:.3f}")

    diff_sigma = diff_psf_sigma(oasis_psf, sauron_psf)
    print(f"  Differential PSF sigma (single-Gaussian): {diff_sigma:.4f} arcsec")

    # ---- 3. Convolve OASIS V and Vrms², then reconstruct sigma ----
    v_conv_bin, _ = bins_convolve(
        xpix_o, ypix_o, binnum_o,
        kins_o[0][binnum_o],  # V (first moment)
        flux_o,
        [diff_sigma], [1.0]
    )
    rms2_conv_bin, _ = bins_convolve(
        xpix_o, ypix_o, binnum_o,
        rms_o[binnum_o] ** 2,  # Vrms²
        flux_o,
        [diff_sigma], [1.0]
    )
    sig_conv_bin = np.sqrt(np.maximum(0, rms2_conv_bin - v_conv_bin ** 2))

    # ---- 4. 1D slits ----
    slit_o_bins, r_o = slit_binnum(xbin_o, ybin_o, angle=0, slit=0.3, rmax=None)
    slit_s_bins, r_s = slit_binnum(xbin_s, ybin_s, angle=0, slit=1.0, rmax=10)

    # ---- 5. Plot ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1.5]})

    offset = 1.05

    # Top: sigma vs radius
    ax1.errorbar(r_o, kins_o[1][slit_o_bins], dkins_o[1][slit_o_bins],
                 fmt="o", ms=4, alpha=0.3, color="C0", label="OASIS sigma")
    ax1.errorbar(r_o, sig_conv_bin[slit_o_bins],
                 fmt="-", lw=1.5, color="k", label="OASIS sigma (PSF-convolved)")
    ax1.errorbar(r_s, kins_s[1][slit_s_bins], dkins0_s[1][slit_s_bins],
                 fmt="s", ms=4, alpha=0.4, color="gray", label="SAURON sigma")
    ax1.errorbar(r_s, kins_s[1][slit_s_bins] * offset, dkins0_s[1][slit_s_bins],
                 fmt="s", ms=4, alpha=0.7, color="C3", label=f"SAURON sigma x {offset:.2f}")

    ax1.set_ylabel("sigma (km/s)")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.set_title(galaxy)
    ax1.minorticks_on()

    # Bottom: ratio
    # Interpolate convolved OASIS sigma to SAURON radii
    from scipy.interpolate import interp1d
    interp_conv = interp1d(r_o, sig_conv_bin[slit_o_bins],
                           bounds_error=False, fill_value=np.nan)

    conv_at_s = interp_conv(r_s)

    ratio_nocorr = kins_s[1][slit_s_bins] / conv_at_s
    ratio_corr = kins_s[1][slit_s_bins] * offset / conv_at_s

    # Summary statistics
    med_nocorr = np.nanmedian(ratio_nocorr)
    med_corr = np.nanmedian(ratio_corr)
    scatter_nocorr = 1.5 * np.nanmedian(np.abs(ratio_nocorr - med_nocorr))
    scatter_corr = 1.5 * np.nanmedian(np.abs(ratio_corr - med_corr))
    print(f"  Ratio SAURON / OASIS_conv:  median={med_nocorr:.4f}, scatter={scatter_nocorr:.4f}")
    print(f"  Ratio SAURON x {offset} / OASIS_conv: median={med_corr:.4f}, scatter={scatter_corr:.4f}")

    # ---- 6. All-bin comparison ----
    out_dir = ROOT / "results" / "JAM" / galaxy
    from scipy.interpolate import LinearNDInterpolator
    interp = LinearNDInterpolator(np.column_stack([xbin_o, ybin_o]), sig_conv_bin)
    sig_conv_at_s = interp(xbin_s, ybin_s)
    valid = np.isfinite(sig_conv_at_s) & (sig_conv_at_s > 0) & (kins_s[1] > 0)

    all_ratio_nocorr = kins_s[1][valid] / sig_conv_at_s[valid]
    all_ratio_corr = kins_s[1][valid] * offset / sig_conv_at_s[valid]

    med_all_nocorr = np.nanmedian(all_ratio_nocorr)
    med_all_corr = np.nanmedian(all_ratio_corr)
    scatter_all_nocorr = 1.5 * np.nanmedian(np.abs(all_ratio_nocorr - med_all_nocorr))
    scatter_all_corr = 1.5 * np.nanmedian(np.abs(all_ratio_corr - med_all_corr))

    n_valid = np.sum(valid)
    print(f"\n  All-bin ({n_valid} SAURON bins):")
    print(f"    SAURON / OASIS_conv:  median={med_all_nocorr:.4f}, scatter={scatter_all_nocorr:.4f}")
    print(f"    SAURON x {offset} / OASIS_conv: median={med_all_corr:.4f}, scatter={scatter_all_corr:.4f}")

    # Figure 2: all-bin scatter + ratio histogram
    fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))
    ax3.scatter(sig_conv_at_s[valid], kins_s[1][valid], c='gray', alpha=0.3, s=4, label='SAURON')
    ax3.scatter(sig_conv_at_s[valid], kins_s[1][valid] * offset, c='C3', alpha=0.4, s=4,
                label=f'SAURON x {offset}')
    lims = [min(sig_conv_at_s[valid].min(), kins_s[1][valid].min()),
            max(sig_conv_at_s[valid].max(), (kins_s[1][valid]*offset).max())]
    ax3.plot(lims, lims, 'k--', alpha=0.3)
    ax3.set_xlabel("OASIS sigma (PSF-convolved) [km/s]")
    ax3.set_ylabel("SAURON sigma [km/s]")
    ax3.legend(fontsize=8)
    ax3.set_title(f"{galaxy}: all-bin comparison")

    bins = np.linspace(0.6, 1.4, 41)
    ax4.hist(all_ratio_nocorr, bins=bins, alpha=0.4, color='gray',
             label=f'no corr (med={med_all_nocorr:.3f})')
    ax4.hist(all_ratio_corr, bins=bins, alpha=0.6, color='C3',
             label=f'x {offset} (med={med_all_corr:.3f})')
    ax4.axvline(1.0, color='k', ls='--', alpha=0.3)
    ax4.set_xlabel("SAURON sigma / OASIS sigma (PSF-convolved)")
    ax4.set_ylabel("N bins")
    ax4.legend(fontsize=8)
    ax4.set_title(f"{galaxy}: ratio distribution")
    fig2.tight_layout()

    out_all = out_dir / "sigma_offset_check_allbins.pdf"
    fig2.savefig(out_all, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"  -> {out_all}")

    ax2.axhline(1.0, color="k", ls="--", alpha=0.3)
    ax2.errorbar(r_s, ratio_nocorr, fmt="s", ms=3, alpha=0.4, color="gray",
                 label=f"SAURON / OASIS_conv (med={med_nocorr:.3f})")
    ax2.errorbar(r_s, ratio_corr, fmt="s", ms=3, alpha=0.7, color="C3",
                 label=f"SAURON x {offset} / OASIS_conv (med={med_corr:.3f})")

    ax2.set_xlabel("Radius along major axis (arcsec)")
    ax2.set_ylabel("Ratio")
    ax2.legend(loc="upper right", fontsize=7)
    ax2.minorticks_on()
    ax2.set_ylim(0.7, 1.3)

    fig.tight_layout()

    out_file = out_dir / "sigma_offset_check.pdf"
    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_file}")


if __name__ == "__main__":
    GALAXIES = {
        "NGC4552": {"pa": 122.5, "release": "r1"},
        "NGC5846": {"pa": 80.0,  "release": "r5"},
        "NGC5813": {"pa": 146.0, "release": "r3"},
    }
    for g, info in GALAXIES.items():
        check_galaxy(g, info["pa"], info["release"])
    print("\nDone.")
