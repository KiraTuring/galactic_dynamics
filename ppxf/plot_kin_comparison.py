import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from plotbin.display_pixels import display_pixels
from ppxfprep import read_oasis_cube, vorbin_spectrum, binNum_gen


def load_ref_kinematics(ref_file):
    hdu = fits.open(ref_file)[1].data
    xpix = hdu['XPIX'][0]
    ypix = hdu['YPIX'][0]
    binnum = hdu['BINNUM'][0]
    xbin = hdu['XBIN'][0]
    ybin = hdu['YBIN'][0]
    kin = {
        'v': hdu['VEL'][0],
        'sigma': hdu['SIG'][0],
        'h3': hdu['H3'][0],
        'h4': hdu['H4'][0],
        'dv': hdu['DVEL'][0],
        'dsigma': hdu['DSIG'][0],
        'dh3': hdu['DH3'][0],
        'dh4': hdu['DH4'][0],
    }
    if 'H5' in hdu.columns.names:
        kin['h5'] = hdu['H5'][0]
        kin['h6'] = hdu['H6'][0]
        kin['dh5'] = hdu['DH5'][0]
        kin['dh6'] = hdu['DH6'][0]
    return xbin, ybin, kin, xpix, ypix, binnum


def load_our_kinematics(ecsv_file, cube_file, target_sn):
    t = ascii.read(ecsv_file)
    xbin = np.array(t['xbin'])
    ybin = np.array(t['ybin'])
    kin = {
        'v': np.array(t['v']),
        'sigma': np.array(t['sigma']),
        'h3': np.array(t['h3']),
        'h4': np.array(t['h4']),
        'dv': np.array(t['dv']),
        'dsigma': np.array(t['dsigma']),
        'dh3': np.array(t['dh3']),
        'dh4': np.array(t['dh4']),
    }
    if 'h5' in t.colnames:
        kin['h5'] = np.array(t['h5'])
        kin['h6'] = np.array(t['h6'])
        kin['dh5'] = np.array(t['dh5'])
        kin['dh6'] = np.array(t['dh6'])

    _, _, x_pix, y_pix, _, _, _ = read_oasis_cube(cube_file)
    binnum = binNum_gen(x_pix, y_pix, xbin, ybin)

    return xbin, ybin, kin, x_pix, y_pix, binnum


def bin_to_pixel(binnum, bin_vals):
    return bin_vals[binnum]


def plot_kin_comparison(xbin1, ybin1, kin1, xpix1, ypix1, binnum1,
                        xbin2, ybin2, kin2, xpix2, ypix2, binnum2,
                        galaxy, systemic_vel=None, output_file=None):
    has_h6 = 'h5' in kin1 and 'h5' in kin2
    if has_h6:
        quantities = ['v', 'sigma', 'h3', 'h4', 'h5', 'h6']
        labels = [r'$V - V_{\rm sys}$', r'$\sigma$', r'$h_3$', r'$h_4$', r'$h_5$', r'$h_6$']
        units = ['km/s', 'km/s', '', '', '', '']
        ncols = 6
        figsize = (30, 10)
    else:
        quantities = ['v', 'sigma', 'h3', 'h4']
        labels = [r'$V - V_{\rm sys}$', r'$\sigma$', r'$h_3$', r'$h_4$']
        units = ['km/s', 'km/s', '', '']
        ncols = 4
        figsize = (22, 10)

    if systemic_vel is None:
        systemic_vel = np.median(kin1['v'])

    vlims = {}
    for q in quantities:
        if q == 'v':
            vals1 = kin1[q] - systemic_vel
            vals2 = kin2[q] - systemic_vel
        else:
            vals1 = kin1[q]
            vals2 = kin2[q]
        lo = min(np.percentile(vals1, 2), np.percentile(vals2, 2))
        hi = max(np.percentile(vals1, 98), np.percentile(vals2, 98))
        margin = 0.05 * (hi - lo) if hi > lo else 1.0
        if q != 'sigma':
            absmax = max(abs(lo), abs(hi))
            vlims[q] = (-absmax, absmax)
        else:
            vlims[q] = (lo - margin, hi + margin)

    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import Normalize

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, ncols, figure=fig, wspace=0.3, hspace=0.25,
                  left=0.03, right=0.97, top=0.92, bottom=0.06)
    suffix = ' (h6)' if has_h6 else ''
    fig.suptitle(f'{galaxy} — Kinematic Comparison{suffix}', fontsize=16, y=0.97)

    row_labels = ['This work (pPXF)', 'Reference (OASIS)']

    for row, (xpix, ypix, binnum, bin_kin, row_label) in enumerate([
        (xpix1, ypix1, binnum1, kin1, row_labels[0]),
        (xpix2, ypix2, binnum2, kin2, row_labels[1]),
    ]):
        for col, (q, lbl, unit) in enumerate(zip(quantities, labels, units)):
            bin_vals = bin_kin[q] - systemic_vel if q == 'v' else bin_kin[q]
            pixel_vals = bin_to_pixel(binnum, bin_vals)
            ax = fig.add_subplot(gs[row, col])
            plt.sca(ax)
            display_pixels(xpix, ypix, pixel_vals,
                           vmin=vlims[q][0], vmax=vlims[q][1],
                           check_grid=False, colorbar=False)
            ax.set_title(lbl if row == 0 else '', fontsize=14)
            if col == 0:
                ax.set_ylabel(row_label, fontsize=12)

            norm = Normalize(vmin=vlims[q][0], vmax=vlims[q][1])
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='sauron'),
                              ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
            if unit:
                cb.set_label(unit, fontsize=9)

    if output_file:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    plt.close(fig)

    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Plot kinematic comparison between pPXF results and reference OASIS kinematics")
    parser.add_argument("galaxy", help="Galaxy name (e.g. NGC5813)")
    parser.add_argument("--ecsv", default=None,
                        help="Our pPXF ECSV kinematics file (default: output/{GALAXY}/gauss_hermite_kins_h4.ecsv)")
    parser.add_argument("--ref", default=None,
                        help="Reference OASIS kinematics FITS (default: ../data/raw/{GALAXY}/oasis/kinematics_oasis_{GALAXY}.fits)")
    parser.add_argument("--cube", default=None,
                        help="OASIS FITS cube for pixel coordinates (default: auto-detect)")
    parser.add_argument("--target-sn", type=float, default=40,
                        help="Voronoi target S/N used (for binning reconstruction)")
    parser.add_argument("--systemic-vel", type=float, default=None,
                        help="Systemic velocity for centering V map (default: median of our V)")
    parser.add_argument("--output", default=None,
                        help="Output PNG file (default: output/{GALAXY}/kin_comparison.png)")
    args = parser.parse_args()

    import glob
    galaxy = args.galaxy
    ecsv_file = args.ecsv or f"output/{galaxy}/gauss_hermite_kins_h4.ecsv"
    ref_file = args.ref or f"../data/raw/{galaxy}/oasis/kinematics_oasis_{galaxy}.fits"
    output_file = args.output or f"output/{galaxy}/kin_comparison.png"
    cube_file = args.cube
    if cube_file is None:
        candidates = sorted(glob.glob(f"../data/raw/{galaxy}/oasis/MS_{galaxy}_oas_r*_E3D.fits*"))
        cube_file = candidates[0] if candidates else None

    print(f"=== Kinematic comparison for {galaxy} ===")
    print(f"  Our data:  {ecsv_file}")
    print(f"  Reference: {ref_file}")
    print(f"  Cube:      {cube_file}")

    xbin1, ybin1, kin1, xpix1, ypix1, binnum1 = load_our_kinematics(
        ecsv_file, cube_file, args.target_sn)
    xbin2, ybin2, kin2, xpix2, ypix2, binnum2 = load_ref_kinematics(ref_file)

    systemic_vel = args.systemic_vel or np.median(kin1['v'])
    print(f"  Systemic V: {systemic_vel:.1f} km/s")
    print(f"  Our bins: {len(xbin1)}, Ref bins: {len(xbin2)}")

    plot_kin_comparison(xbin1, ybin1, kin1, xpix1, ypix1, binnum1,
                        xbin2, ybin2, kin2, xpix2, ypix2, binnum2,
                        galaxy, systemic_vel=systemic_vel, output_file=output_file)


if __name__ == "__main__":
    main()
