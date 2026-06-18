import corner
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import ndimage, signal
from ppxf.ppxf import ppxf, rebin
from tqdm import tqdm
from time import perf_counter as clock


def read_sauron_spectcube(fits_file, verbose=False, plot=True):
    with fits.open(fits_file) as hdul:
        if verbose:
            print("FITS file structure:")
            hdul.info()
            print("Header of the first extension:")
            print(repr(hdul[0].header))

        header = hdul[0].header
        spectra = hdul[0].data
        variance = hdul[1].data
        lam_range = header['CRVAL1'] + \
            np.array([0., header['CDELT1'] * (header['NAXIS1'] - 1)])

        print(header['NAXIS1'], spectra.T.shape[0])

        table = hdul[2].data
        x = table["A"]
        y = table["D"]
        flux = table["FLUX"]
        sn = table["SN"]

    if plot:
        plt.figure(1)
        from plotbin.display_pixels import display_pixels
        display_pixels(x, y, np.log(np.mean(spectra, 1)))
        plt.xlabel("arcsec")
        plt.ylabel("arcsec")

        lam = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        plt.figure(2)
        plt.plot(lam, spectra[1, :])
        plt.xlabel(r"wavelength $\AA$")
        plt.ylabel("Flux")

    return lam_range, spectra, x, y, flux, sn, variance


def mock_spectrum(templates, vel, sigma, h3, h4, sn, h5=0, h6=0, factor=10, rng=None, ln_lam=None, ln_lam_temp=None):

    dx = int(abs(vel) + 5*sigma)
    x = np.linspace(-dx, dx, 2*dx*factor + 1)
    w = (x - vel)/sigma
    w2 = w**2
    gauss = np.exp(-0.5*w2)
    gauss /= np.sum(gauss)
    h3poly = w*(2.*w2 - 3.)/np.sqrt(3.)
    h4poly = (w2*(4.*w2 - 12.) + 3.)/np.sqrt(24.)
    h5poly = w*(w2*(4*w2 - 20.) + 15.)/np.sqrt(60.)
    h6poly = (w2*(w2*(w2 - 15.) + 45.) - 15.)/np.sqrt(720.)
    losvd = gauss * (1. + h3*h3poly + h4*h4poly + h5*h5poly + h6*h6poly)

    galaxy = signal.fftconvolve(templates, losvd, mode="same")
    if factor is not None:
        galaxy = rebin(galaxy, factor)
    if ln_lam is not None and ln_lam_temp is not None:
        from scipy.interpolate import interp1d
        interp_func = interp1d(
            ln_lam_temp, galaxy, kind='linear', bounds_error=False, fill_value="extrapolate")
        galaxy = interp_func(ln_lam)

    noise = np.clip(galaxy, 1, None)/sn
    if rng is None:
        try:
            galaxy = np.random.normal(galaxy, noise)
        except:
            print("Error adding noise to the galaxy spectrum. Check the dimensions.")
            print(f"Galaxy shape: {galaxy.shape}, Noise shape: {noise.shape}")
            print(np.min(galaxy), np.max(galaxy))
            print(np.min(noise), np.max(noise))
            raise
    else:
        galaxy = rng.normal(galaxy, noise)

    return galaxy, noise, dx


def ppxf_mc_sim(templates, velscale, nruns=1000, h3=0.1, h4=0.1, sn=60, moments=4, bias=0.5, factor=10, seed=123, sigma_range=(0.5, 4)):

    rng = np.random.default_rng(seed)
    starNew = ndimage.zoom(templates, factor, order=3)
    star = rebin(starNew, factor)
    m = nruns

    velV = rng.uniform(size=m)
    sigmaV = np.linspace(sigma_range[0], sigma_range[1], m)
    result = np.zeros((m, moments))

    t = clock()
    for j, (vel, sigma) in tqdm(enumerate(zip(velV, sigmaV))):

        quiet = False if j < 3 else True

        galaxy, noise, dx = mock_spectrum(
            starNew, vel, sigma, h3, h4, sn, factor=factor, rng=rng)
        start = np.array([vel + rng.uniform(-1, 1), sigma *
                          rng.uniform(0.8, 1.2)])*velscale

        pp = ppxf(star, galaxy, noise, velscale, start,
                  goodpixels=np.arange(dx, galaxy.size - dx),
                  plot=False, moments=moments, bias=bias, quiet=quiet)
        result[j, :] = pp.sol

    print('Calculation time: %.2f s' % (clock()-t))
    return velV, sigmaV, result


def plot_percentiles(x0, y0, nbins=10, sigma=1):

    n = x0.size//nbins
    x = x0[:nbins*n]
    y = y0[:nbins*n]

    j = np.argsort(x)
    x = x[j]
    y = y[j]

    x = x.reshape(-1, n)
    y = y.reshape(-1, n)
    xMean = np.mean(x, axis=1)

    from scipy import stats
    p = 200*stats.norm.cdf(sigma) - 100

    perc = (100 + np.array([-p, 0, p]))/2
    yMin, yMedian, yMax = np.percentile(y, perc, axis=1)

    plt.plot(x0, y0, '.', zorder=0)
    plt.plot(xMean, yMedian, 'lime')
    plt.fill_between(xMean, yMin, yMax, facecolor='gold')


def plot_bias(velV, sigmaV, velscale, h3, h4, result, bias, h5=0, h6=0):

    ngh = result.shape[1]
    gh_list = np.zeros(ngh-2)
    gh_list[0] = h3
    gh_list[1] = h4

    nrow = ngh//2
    ncolumn = 2
    plt.clf()
    fig, axes = plt.subplots(nrow, ncolumn, figsize=(10, nrow*4))
    plt.subplot(nrow, ncolumn, 1)
    plot_percentiles(sigmaV*velscale, result[:, 0] - velV*velscale)
    plt.axhline(0, color='r')
    plt.axvline(velscale, linestyle='dashed')
    plt.axvline(2*velscale, linestyle='dashed')
    plt.ylim(-60, 60)
    plt.xlabel(r'$\sigma_{\rm in}\ (km\ s^{-1})$')
    plt.ylabel(r'$V - V_{\rm in}\ (km\ s^{-1})$')
    plt.text(2.05*velscale, -15, r'2$\times$velscale')

    plt.subplot(nrow, ncolumn, 2)
    plot_percentiles(sigmaV*velscale, result[:, 1] - sigmaV*velscale)
    plt.axhline(0, color='r')
    plt.axvline(velscale, linestyle='dashed')
    plt.axvline(2*velscale, linestyle='dashed')
    plt.ylim(-60, 60)
    plt.xlabel(r'$\sigma_{in}\ (km\ s^{-1})$')
    plt.ylabel(r'$\sigma - \sigma_{\rm in}\ (km\ s^{-1})$')
    plt.text(2.05*velscale, -15, r'2$\times$velscale')

    for i, h in enumerate(gh_list, start=3):
        plt.subplot(nrow, ncolumn, i)
        plot_percentiles(sigmaV*velscale, result[:, i-1])
        plt.axhline(h, color='r')
        plt.axhline(0, linestyle='dotted', color='limegreen')
        plt.axvline(velscale, linestyle='dashed')
        plt.axvline(2*velscale, linestyle='dashed')
        plt.ylim(-0.2+h, 0.2+h)
        plt.xlabel(r'$\sigma_{\rm in}\ (km\ s^{-1})$')
        plt.ylabel(f'$h_{i}$')
        plt.text(2.05*velscale, h - 0.15, r'2$\times$velscale')

    plt.suptitle(f'bias={bias}')
    plt.tight_layout()


def plot_corner(result, labels=None, fontsize=10, threshold=0.1, **kwargs):

    ndim = result.shape[1]
    if labels is None:
        labels = [f'h{i+1}' for i in range(ndim)]
        labels[0] = "vel (km/s)"
        labels[1] = "sig (km/s)"
    fig = corner.corner(result, labels=labels, show_titles=True, **kwargs)
    axes = fig.get_axes()
    for i in range(ndim):
        for j in range(i):
            ax = axes[i * ndim + j]

            x = result[:, j]
            y = result[:, i]
            r = np.corrcoef(x, y)[0, 1]
            if abs(r) >= threshold:
                ax.text(
                    0.05, 0.95,
                    rf"$r$ = {r:.2f}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=fontsize,
                )
    return fig

