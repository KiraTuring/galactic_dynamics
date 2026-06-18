import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from plotbin.display_pixels import display_pixels
import ppxf.ppxf_util as util
from ppxf.ppxf import ppxf


def read_oasis_cube(fits_file, verbose=False, plot=False):
    hdu = fits.open(fits_file)
    if verbose:
        print("FITS file structure:")
        hdu.info()
        print("Header of the first extension:")
        print(repr(hdu[0].header))
        print("Header of the second extension:")
        print(repr(hdu[1].header))
        print("Header of the fourth extension:")
        print(repr(hdu[3].header))

    header = hdu[1].header
    data = hdu[1].data

    spectra = data["data_spe"]
    variance = data["stat_spe"]
    x = data["xpos"]
    y = data["ypos"]

    npix = spectra.shape[1]
    lam = header['CRVALS'] + header['CDELTS']*np.arange(npix)
    lam_range = np.array([lam[0], lam[-1]])

    flux = np.mean(spectra, axis=1)
    noise = np.mean(np.sqrt(variance), axis=1)
    sn = flux / noise

    if plot:
        plt.figure(1)
        display_pixels(x, y, np.log(np.mean(spectra, 1)))
        plt.xlabel("arcsec")
        plt.ylabel("arcsec")

        plt.figure(2, figsize=(20, 4))
        plt.errorbar(lam, spectra[1, :], yerr=np.sqrt(variance[1, :]))
        plt.xlabel(r"wavelength $\AA$")
        plt.ylabel("Flux")

    return lam_range, spectra, x, y, flux, sn, variance


def bin_spectrum(spectrum, binNum):
    weighted_sums = np.zeros((np.max(binNum)+1, spectrum.shape[1]))
    for i in range(spectrum.shape[1]):
        weighted_sums[:, i] = np.bincount(binNum, weights=spectrum[:, i])

    binFlux = np.bincount(binNum)
    weighted_spectrum = weighted_sums / binFlux.reshape(-1, 1)

    return weighted_spectrum, binFlux


def vorbin_spectrum(spectrum, x, y, flux, sn, target_sn=60):
    binNum, x_gen, y_gen, x_bar, y_bar, sn_new, nPixels, scale = voronoi_2d_binning(
        x, y, flux, flux/sn, target_sn, plot=1, quiet=1)

    weighted_spectrum, binFlux = bin_spectrum(spectrum, binNum)

    return binNum, x_gen, y_gen, weighted_spectrum, binFlux


def read_miles_lib(dirname='./miles_lib/MILES_library_v9.1_FITS'):
    snum = 985
    lamRange0, spectrum0 = _read_single_miles_spectrum(1, dirname=dirname)
    miles_spectrum = np.zeros((snum, spectrum0.shape[1]))
    for sid in range(1, snum+1):
        lamRange1, spectrum1 = _read_single_miles_spectrum(sid, dirname=dirname)
        if (lamRange1 != lamRange0).any():
            raise ValueError(
                f"Different wavelength ranges for MILES spectrum {sid}")
        if spectrum1.shape != spectrum0.shape:
            raise ValueError(f"Different shapes for MILES spectrum {sid}")
        miles_spectrum[sid-1] = spectrum1

    return lamRange0, miles_spectrum


def _read_single_miles_spectrum(sid, dirname):
    path = f'{dirname}/s{sid:04d}.fits'
    hdu = fits.open(path)
    spectrum = hdu[0].data
    header = hdu[0].header
    lam_range = header['CRVAL1'] + \
        np.array([0., header['CDELT1'] * (header['NAXIS1'] - 1)])
    return lam_range, spectrum


def log_rebin_and_normalize(lamRange, spectrum, variance, binFlux, velscale=None):

    rebin_spectrum, ln_lam, velscale = util.log_rebin(
        lamRange, binFlux * spectrum.T, velscale=velscale)
    rebin_variance, _, _ = util.log_rebin(
        lamRange, binFlux * variance.T, velscale=velscale)

    norm = np.median(rebin_spectrum)
    rebin_spectrum /= norm
    rebin_variance /= norm**2

    return rebin_spectrum, rebin_variance, ln_lam, velscale


def prepare_miles_templates(lamRange_miles, spectrum_miles, fwhm_gal, velscale):

    fwhm_miles = 2.5

    fwhm_diff2 = np.clip(fwhm_gal**2 - fwhm_miles**2, 0, None)
    sigma = np.sqrt(fwhm_diff2) / np.sqrt(4 * np.log(4))

    lam_miles = np.linspace(
        lamRange_miles[0], lamRange_miles[-1], spectrum_miles.shape[1])
    spectrum_miles_con = util.varsmooth(lam_miles, spectrum_miles.T, sigma)

    templates, ln_lam_temp = util.log_rebin(
        lamRange_miles, spectrum_miles_con, velscale=velscale)[:2]
    lam_temp = np.exp(ln_lam_temp)

    return templates, ln_lam_temp, lam_temp


def run_ppxf(templates, galaxy, velscale, ln_lam=None, lam_temp=None, redshift=0, goodPixels=None,
             noise=None, noise_level=0.0047, plot=False, quiet=False, bias=None, moments=6, degree=4):

    if noise is None:
        noise = np.full_like(galaxy, noise_level)

    if goodPixels is None:
        lam_range_temp = [np.min(lam_temp), np.max(lam_temp)]
        goodPixels = util.determine_goodpixels(
            ln_lam, lam_range_temp, redshift)

    c = 299792.458
    vel = c*np.log(1 + redshift)
    start = [vel, 200.]

    lam = None
    if ln_lam is not None:
        lam = np.exp(ln_lam)
    kin_list = []
    dkin_list = []
    pp_list = []
    nbins = galaxy.shape[1]
    report_every = max(1, min(20, nbins // 5))
    show_progress = nbins > 3
    for i, (gal, err) in enumerate(zip(galaxy.T, noise.T)):
        pp = ppxf(templates, gal, err, velscale, start,
                  goodpixels=goodPixels, moments=moments, lam=lam,
                  lam_temp=lam_temp, degree=degree, quiet=quiet, bias=bias)
        if bias is None and quiet is False:
            print(f'PPXF bias={pp.bias}')
        if plot:
            pp.plot()
            plt.show()
        kin_list.append(pp.sol)
        dkin_list.append(pp.error*np.sqrt(pp.chi2))
        pp_list.append(pp)
        start = pp.sol.copy()
        if show_progress and ((i + 1) % report_every == 0 or i + 1 == nbins):
            print(f"  Fitted {i + 1}/{nbins} bins")

    return kin_list, dkin_list, pp_list


def bootstrap_residuals(model, resid, wild=True):
    if wild:
        eps = resid * (2 * np.random.randint(2, size=resid.size) - 1)
    else:
        eps = np.random.choice(resid, size=resid.size)

    return model + eps


def run_ppxf_bootstrap(templates, galaxy, velscale, ln_lam=None, lam_temp=None, redshift=0,
                       goodPixels=None, noise=None, noise_level=0.0047, plot=False, quiet=False,
                       bias=None, moments=6, degree=4, nrand=9):
    if noise is None:
        noise = np.full_like(galaxy, noise_level)

    if goodPixels is None:
        lam_range_temp = [np.min(lam_temp), np.max(lam_temp)]
        goodPixels = util.determine_goodpixels(ln_lam, lam_range_temp, redshift)

    c = 299792.458
    vel = c * np.log(1 + redshift)
    start = [vel, 200.]

    lam = None
    if ln_lam is not None:
        lam = np.exp(ln_lam)
    kin_list = []
    dkin_list = []
    pp_list = []
    nbins = galaxy.shape[1]
    report_every = max(1, min(20, nbins // 5))
    show_progress = nbins > 3

    for i, (gal, err) in enumerate(zip(galaxy.T, noise.T)):
        pp = ppxf(templates, gal, err, velscale, start,
                  goodpixels=goodPixels, moments=moments, lam=lam,
                  lam_temp=lam_temp, degree=degree, quiet=quiet, bias=bias)
        if bias is None and quiet is False:
            print(f'PPXF bias={pp.bias}')
        if plot:
            pp.plot()
            plt.show()

        pp_list.append(pp)
        bestfit0 = pp.bestfit.copy()
        resid = gal - bestfit0
        start = pp.sol.copy()
        sol_list = np.zeros((nrand, moments))
        for j in range(nrand):
            galaxy1 = bootstrap_residuals(bestfit0, resid)
            pp = ppxf(templates, galaxy1, err, velscale, start,
                      goodpixels=goodPixels, moments=moments, lam=lam,
                      lam_temp=lam_temp, degree=degree, quiet=quiet, bias=0)
            sol_list[j] = pp.sol
        kin_list.append(np.mean(sol_list, axis=0))
        dkin_list.append(np.std(sol_list, axis=0))
        if show_progress and ((i + 1) % report_every == 0 or i + 1 == nbins):
            print(f"  Bootstrap {i + 1}/{nbins} bins")

    return kin_list, dkin_list, pp_list


def save_fits_sauron(kins, dkins, x, y, flux, filename, moments=4):
    cols = [
        fits.Column(name='VPXF', array=kins[:, 0], format='E'),
        fits.Column(name='SPXF', array=kins[:, 1], format='E'),
        fits.Column(name='H3PXF', array=kins[:, 2], format='E'),
        fits.Column(name='H4PXF', array=kins[:, 3], format='E'),
        fits.Column(name='EVPXF', array=dkins[:, 0], format='E'),
        fits.Column(name='ESPXF', array=dkins[:, 1], format='E'),
        fits.Column(name='EH3PXF', array=dkins[:, 2], format='E'),
        fits.Column(name='EH4PXF', array=dkins[:, 3], format='E'),
    ]
    if moments >= 6:
        cols.extend([
            fits.Column(name='H5PXF', array=kins[:, 4], format='E'),
            fits.Column(name='H6PXF', array=kins[:, 5], format='E'),
            fits.Column(name='EH5PXF', array=dkins[:, 4], format='E'),
            fits.Column(name='EH6PXF', array=dkins[:, 5], format='E'),
        ])
    cols.extend([
        fits.Column(name='XS', array=x, format='E'),
        fits.Column(name='YS', array=y, format='E'),
        fits.Column(name='FLUX', array=flux, format='E'),
        fits.Column(name='NO', array=range(len(kins)), format='J'),
    ])
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto(filename, overwrite=True)


def save_fits_oasis(xbin, ybin, kins, dkins, x, y, binnum, flux, filename, moments=4):
    nbins = kins.shape[0]
    npixs = x.size

    kin_names = ['VEL', 'SIG', 'H3', 'H4', 'H5', 'H6'][:moments]
    dkin_names = ['DVEL', 'DSIG', 'DH3', 'DH4', 'DH5', 'DH6'][:moments]

    cols = [
        fits.Column(name='XBIN', array=[xbin], format=f'{nbins}E'),
        fits.Column(name='YBIN', array=[ybin], format=f'{nbins}E'),
    ]
    for i, (kn, dn) in enumerate(zip(kin_names, dkin_names)):
        cols.append(fits.Column(name=kn, array=[kins[:, i]], format=f'{nbins}E'))
        cols.append(fits.Column(name=dn, array=[dkins[:, i]], format=f'{nbins}E'))
    cols.extend([
        fits.Column(name='XPIX', array=[x], format=f'{npixs}E'),
        fits.Column(name='YPIX', array=[y], format=f'{npixs}E'),
        fits.Column(name='BINNUM', array=[binnum], format=f'{npixs}J'),
        fits.Column(name='SURF', array=[flux], format=f'{npixs}E'),
    ])
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto(filename, overwrite=True)


def save_kinematics_ecsv(kin_list, dkin_list, x_gen, y_gen, filepath, moments, galaxy,
                         binNum=None, x_pix=None, y_pix=None):
    from astropy import table as atable

    kins = np.array(kin_list)
    dkins = np.array(dkin_list)
    nbins = len(kins)

    col_names = ['v', 'sigma', 'h3', 'h4', 'h5', 'h6'][:moments]
    dcol_names = ['dv', 'dsigma', 'dh3', 'dh4', 'dh5', 'dh6'][:moments]

    t = atable.Table()
    t['binID_dynamite'] = np.arange(1, nbins + 1)
    for i, name in enumerate(col_names):
        t[name] = kins[:, i]
    for i, name in enumerate(dcol_names):
        t[name] = dkins[:, i]
    t['xbin'] = x_gen
    t['ybin'] = y_gen
    t['n_gh'] = np.full(nbins, moments, dtype=int)
    t['is_good'] = np.ones(nbins, dtype=bool)

    if binNum is not None:
        t.meta['galaxy'] = galaxy
        t.meta['moments'] = moments
    if x_pix is not None and y_pix is not None:
        t.meta['n_pix'] = len(x_pix)

    t.write(filepath, overwrite=True, format='ascii.ecsv')
    return t


def mask_to_intervals(x, mask):
    mask = np.asarray(mask)
    x = np.asarray(x)

    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []

    breaks = np.where(np.diff(idx) > 1)[0]

    starts = np.insert(idx[breaks + 1], 0, idx[0])
    ends = np.append(idx[breaks], idx[-1])

    return [(x[s], x[e]) for s, e in zip(starts, ends)]


def binNum_gen(x, y, xBin, yBin):
    if x.size < 1e4:
        binNum = np.argmin((x[:, None] - xBin)**2 + (y[:, None] - yBin)**2, axis=1)
    else:
        binNum = np.zeros(x.size, dtype=int)
        for j, (xj, yj) in enumerate(zip(x, y)):
            binNum[j] = np.argmin((xj - xBin)**2 + (yj - yBin)**2)
    return binNum
