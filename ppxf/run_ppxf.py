import argparse
import glob
import os
import numpy as np
from datetime import datetime
from time import perf_counter as clock
import yaml

from ppxfprep import (
    read_oasis_cube,
    vorbin_spectrum,
    bin_spectrum,
    read_miles_lib,
    log_rebin_and_normalize,
    prepare_miles_templates,
    run_ppxf,
    run_ppxf_bootstrap,
    save_fits_oasis,
    save_fits_sauron,
    save_kinematics_ecsv,
    mask_to_intervals,
)
import ppxf.ppxf_util as util
from astropy.io import fits


def main():
    parser = argparse.ArgumentParser(
        description="pPXF kinematic extraction pipeline for OASIS IFU datacubes")
    parser.add_argument("galaxy", help="Galaxy name (e.g. NGC4621)")
    parser.add_argument("--cube", default=None,
                        help="OASIS FITS cube path (default: auto-detect from ../data/raw/{GALAXY}/oasis/)")
    parser.add_argument("--ref-bins", default=None, const="auto", nargs="?",
                        help="Use reference OASIS kinematics FITS for binning (default: auto-detect; skip to use Voronoi)")
    parser.add_argument("--miles", default="./miles_lib/MILES_library_v9.1_FITS",
                        help="MILES library directory")
    parser.add_argument("--target-sn", type=float, default=60,
                        help="Voronoi binning target S/N (default: 60, ignored if --ref-bins)")
    parser.add_argument("--redshift", type=float, required=True,
                        help="Galaxy redshift")
    parser.add_argument("--fwhm-gal", type=float, default=5.4,
                        help="Instrumental FWHM in Angstrom (default: 5.4)")
    parser.add_argument("--moments", type=int, default=6,
                        help="GH moments for pPXF fit: 4 or 6 (default: 6)")
    parser.add_argument("--bias", type=float, default=0.5,
                        help="pPXF bias parameter (default: 0.5)")
    parser.add_argument("--degree", type=int, default=4,
                        help="Additive polynomial degree (default: 4)")
    parser.add_argument("--bootstrap", type=int, default=100,
                        help="Bootstrap iterations for error estimation (0=skip, default: 100)")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory (default: output/{GALAXY}[_refbins])")
    parser.add_argument("--trial-moments", type=int, nargs='+', default=[4, 6],
                        help="Moments to test on total spectrum (default: 4 6)")
    args = parser.parse_args()

    galaxy = args.galaxy
    use_ref_bins = args.ref_bins is not None
    if use_ref_bins and args.ref_bins == "auto":
        candidates = sorted(glob.glob(f"../data/raw/{galaxy}/oasis/kinematics_oasis_{galaxy}.fits*"))
        if candidates:
            args.ref_bins = candidates[0]
        else:
            print(f"  WARNING: No reference kinematics found, falling back to Voronoi binning")
            use_ref_bins = False
            args.ref_bins = None
    if use_ref_bins:
        args.output_dir = args.output_dir or f"output/{galaxy}_refbins"
    else:
        args.output_dir = args.output_dir or f"output/{galaxy}"
    os.makedirs(args.output_dir, exist_ok=True)

    if args.cube is None:
        candidates = sorted(glob.glob(f"../data/raw/{galaxy}/oasis/MS_{galaxy}_oas_r*_E3D.fits*"))
        if candidates:
            args.cube = candidates[0]
        else:
            args.cube = f"../data/raw/{galaxy}/oasis/MS_{galaxy}_oas_r3_E3D.fits"

    print(f"=== pPXF pipeline for {galaxy} ===")
    print(f"  Cube:      {args.cube}")
    print(f"  Redshift:  {args.redshift}")
    print(f"  FWHM_gal:  {args.fwhm_gal} A")
    print(f"  Moments:   {args.moments}")
    print(f"  Bias:      {args.bias}")
    print(f"  Bootstrap: {args.bootstrap}")
    if use_ref_bins:
        print(f"  Ref bins:  {args.ref_bins}")
    print(f"  Output:    {args.output_dir}/")

    print("\n--- Step 1: Read OASIS cube ---")
    lamRange, spectrum, x, y, flux, sn, variance = read_oasis_cube(args.cube)
    print(f"  Spectrum shape: {spectrum.shape}, N_spaxels: {len(x)}")

    if use_ref_bins:
        print(f"\n--- Step 2: Use reference binning from {args.ref_bins} ---")
        ref = fits.open(args.ref_bins)[1].data
        x_gen = ref['XBIN'][0]
        y_gen = ref['YBIN'][0]
        nbins_ref = len(x_gen)
        binNum = ref['BINNUM'][0]
        weighted_spectrum, binFlux = bin_spectrum(spectrum, binNum)
        weighted_variance, _ = bin_spectrum(variance, binNum)

        if weighted_spectrum.shape[0] < nbins_ref:
            n_missing = nbins_ref - weighted_spectrum.shape[0]
            weighted_spectrum = np.vstack([weighted_spectrum, np.zeros((n_missing, weighted_spectrum.shape[1]))])
            weighted_variance = np.vstack([weighted_variance, np.zeros((n_missing, weighted_variance.shape[1]))])
            binFlux = np.append(binFlux, np.zeros(n_missing))

        good = binFlux > 0
        weighted_spectrum = weighted_spectrum[good]
        weighted_variance = weighted_variance[good]
        x_gen = x_gen[good]
        y_gen = y_gen[good]
        binFlux = binFlux[good]
        print(f"  N_bins (reference): {len(x_gen)}")
    else:
        print(f"\n--- Step 2: Voronoi binning (target_sn={args.target_sn}) ---")
        binNum, x_gen, y_gen, weighted_spectrum, binFlux = vorbin_spectrum(
            spectrum, x, y, flux, sn, target_sn=args.target_sn)
        weighted_variance, _ = bin_spectrum(variance, binNum)
        print(f"  N_bins (Voronoi): {len(x_gen)}")

    print("\n--- Step 3: Log-rebin + normalize ---")
    rebin_spectrum, rebin_variance, ln_lam, velscale = log_rebin_and_normalize(
        lamRange, weighted_spectrum, weighted_variance, binFlux)
    print(f"  Vel_scale: {velscale:.2f} km/s, Rebin shape: {rebin_spectrum.shape}")

    print("\n--- Step 4: Prepare MILES templates ---")
    lamRange_miles, spectrum_miles = read_miles_lib(dirname=args.miles)
    templates, ln_lam_temp, lam_temp = prepare_miles_templates(
        lamRange_miles, spectrum_miles, args.fwhm_gal, velscale)
    print(f"  Templates shape: {templates.shape}")

    print("\n--- Step 5: Trial fit on total spectrum ---")
    total_spectrum = np.atleast_2d(np.sum(rebin_spectrum, axis=1)).T
    total_noise = np.atleast_2d(np.sqrt(np.sum(rebin_variance, axis=1))).T

    pp_trial = None
    for m in args.trial_moments:
        _, _, pps = run_ppxf(templates, total_spectrum, velscale, ln_lam, lam_temp,
                             args.redshift, noise=total_noise, plot=False,
                             quiet=True, moments=m, bias=args.bias)
        print(f"  moments={m}: chi2={pps[0].chi2:.2f}, sol={pps[0].sol}")
        if m == args.moments:
            pp_trial = pps[0]

    if pp_trial is None:
        _, _, pps = run_ppxf(templates, total_spectrum, velscale, ln_lam, lam_temp,
                             args.redshift, noise=total_noise, plot=False,
                             quiet=True, moments=args.moments, bias=args.bias)
        pp_trial = pps[0]

    besttemp = templates @ pp_trial.weights.copy()
    print(f"  Best template computed from weights")

    print(f"\n--- Step 6: pPXF fit on all bins (moments={args.moments}) ---")
    t = clock()
    kin_list, dkin_list, pp_list = run_ppxf(
        besttemp, rebin_spectrum, velscale, ln_lam, lam_temp, args.redshift,
        noise=np.sqrt(rebin_variance), quiet=True, moments=args.moments,
        bias=args.bias, degree=args.degree)
    print(f'  Elapsed: {clock() - t:.1f} s, N_bins fitted: {len(kin_list)}')

    if args.bootstrap > 0:
        print(f"\n--- Step 7: Bootstrap error estimation (nrand={args.bootstrap}) ---")
        t = clock()
        kin_list, dkin_list, pp_list = run_ppxf_bootstrap(
            besttemp, rebin_spectrum, velscale, ln_lam, lam_temp, args.redshift,
            noise=np.sqrt(rebin_variance), quiet=True, moments=args.moments,
            bias=args.bias, degree=args.degree, nrand=args.bootstrap)
        print(f'  Elapsed: {clock() - t:.1f} s')

    kins = np.array(kin_list)
    dkins = np.array(dkin_list)

    print("\n--- Step 8: Save results ---")
    npz_file = os.path.join(args.output_dir, f"ppxf_kins_mileslib_h{args.moments}_b{args.bias}.npz")
    np.savez(npz_file, kin_list=kin_list, dkin_list=dkin_list,
             x_gen=x_gen, y_gen=y_gen, velscale=velscale, ln_lam=ln_lam)
    print(f"  Saved: {npz_file}")

    fits_file = os.path.join(args.output_dir, f"ppxf_kins_{galaxy}_h{args.moments}.fits")
    save_fits_oasis(x_gen, y_gen, kins, dkins, x, y, binNum,
                    flux, fits_file, moments=args.moments)
    print(f"  Saved: {fits_file}")

    ecsv_file = os.path.join(args.output_dir, f"gauss_hermite_kins_h{args.moments}.ecsv")
    save_kinematics_ecsv(kin_list, dkin_list, x_gen, y_gen, ecsv_file,
                         args.moments, galaxy, binNum=binNum, x_pix=x, y_pix=y)
    print(f"  Saved: {ecsv_file}")

    lam_range_temp = [np.min(lam_temp), np.max(lam_temp)]
    mask = util.determine_mask(ln_lam, lam_range_temp, args.redshift)
    intervals = mask_to_intervals(np.exp(ln_lam), ~mask)
    print(f"\n  Masked wavelength intervals: {intervals}")

    desc_file = os.path.join(args.output_dir, "description.yaml")
    desc = {}
    if os.path.exists(desc_file):
        with open(desc_file) as f:
            desc = yaml.safe_load(f) or {}
    if 'runs' not in desc:
        desc['runs'] = []
    run_record = {
        'timestamp': datetime.now().isoformat(timespec='seconds'),
        'moments': args.moments,
        'bias': args.bias,
        'bootstrap': args.bootstrap,
        'degree': args.degree,
        'fwhm_gal': args.fwhm_gal,
        'redshift': args.redshift,
        'n_bins': len(x_gen),
        'binning': 'ref_bins' if use_ref_bins else f'voronoi_sn{args.target_sn}',
        'output_files': [
            os.path.basename(npz_file),
            os.path.basename(fits_file),
            os.path.basename(ecsv_file),
        ],
    }
    desc['runs'].append(run_record)
    with open(desc_file, 'w') as f:
        yaml.dump(desc, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    print(f"  Updated: {desc_file}")

    print(f"\n=== Done! Kinematics for {galaxy} saved to {args.output_dir}/ ===")


if __name__ == "__main__":
    main()
