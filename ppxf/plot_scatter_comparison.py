import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii


def main():
    parser = argparse.ArgumentParser(description="Bin-by-bin scatter comparison of pPXF vs reference kinematics")
    parser.add_argument("galaxy", help="Galaxy name (e.g. NGC5813)")
    parser.add_argument("--ecsv", default=None,
                        help="Our pPXF ECSV file (default: output/{GALAXY}_refbins/gauss_hermite_kins_h4.ecsv)")
    parser.add_argument("--ref", default=None,
                        help="Reference OASIS kinematics FITS (default: ../data/raw/{GALAXY}/oasis/kinematics_oasis_{GALAXY}.fits)")
    parser.add_argument("--output", default=None,
                        help="Output PNG file (default: output/{GALAXY}_refbins/scatter_comparison.png)")
    args = parser.parse_args()

    galaxy = args.galaxy
    ecsv_file = args.ecsv or f"output/{galaxy}_refbins/gauss_hermite_kins_h4.ecsv"
    ref_file = args.ref or f"../data/raw/{galaxy}/oasis/kinematics_oasis_{galaxy}.fits"
    output_file = args.output or f"output/{galaxy}_refbins/scatter_comparison.png"

    t = ascii.read(ecsv_file)
    d = fits.open(ref_file)[1].data

    has_h6 = 'h5' in t.colnames

    ours = {
        'v': np.array(t['v']),
        'sigma': np.array(t['sigma']),
        'h3': np.array(t['h3']),
        'h4': np.array(t['h4']),
    }
    ref = {
        'v': d['VEL'][0],
        'sigma': d['SIG'][0],
        'h3': d['H3'][0],
        'h4': d['H4'][0],
    }
    if has_h6:
        ours['h5'] = np.array(t['h5'])
        ours['h6'] = np.array(t['h6'])
        ref['h5'] = d['H5'][0]
        ref['h6'] = d['H6'][0]

    if has_h6:
        quantities = ['v', 'sigma', 'h3', 'h4', 'h5', 'h6']
        labels = [r'$V - V_{\rm sys}$ [km/s]', r'$\sigma$ [km/s]',
                  r'$h_3$', r'$h_4$', r'$h_5$', r'$h_6$']
        nrows, ncols = 2, 3
        figsize = (18, 12)
    else:
        quantities = ['v', 'sigma', 'h3', 'h4']
        labels = [r'$V - V_{\rm sys}$ [km/s]', r'$\sigma$ [km/s]', r'$h_3$', r'$h_4$']
        nrows, ncols = 2, 2
        figsize = (12, 12)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)

    for ax, q, lbl in zip(axes.flat, quantities, labels):
        x = ref[q]
        y = ours[q]
        if q == 'v':
            x = x - np.median(x)
            y = y - np.median(y)
        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        margin = 0.05 * (hi - lo) if hi > lo else 1.0
        lo -= margin
        hi += margin

        ax.scatter(x, y, s=8, alpha=0.5, edgecolors='none')
        ax.plot([lo, hi], [lo, hi], 'k--', lw=0.8)

        if q != 'sigma' and lo < 0 < hi:
            ax.axhline(0, color='grey', ls='--', lw=0.6)
            ax.axvline(0, color='grey', ls='--', lw=0.6)

        diff = y - x
        r = np.corrcoef(x, y)[0, 1]
        med_diff = np.median(diff)
        med_abs_diff = np.median(np.abs(diff))

        if q in ('v', 'sigma'):
            text = f'$r$ = {r:.3f}\nmedian $\\Delta$ = {med_diff:+.1f} km/s\nmedian $|\\Delta|$ = {med_abs_diff:.1f} km/s'
        else:
            text = f'$r$ = {r:.3f}\nmedian $\\Delta$ = {med_diff:+.4f}\nmedian $|\\Delta|$ = {med_abs_diff:.4f}'
        ax.text(0.05, 0.95, text, transform=ax.transAxes, va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect('equal')
        ax.set_xlabel(f'Reference {lbl}', fontsize=12)
        ax.set_ylabel(f'pPXF {lbl}', fontsize=12)
        ax.set_title(lbl, fontsize=14)

    suffix = ' (h6)' if has_h6 else ''
    fig.suptitle(f'{galaxy} — Bin-by-bin Kinematic Comparison{suffix}', fontsize=16, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close(fig)


if __name__ == "__main__":
    main()
