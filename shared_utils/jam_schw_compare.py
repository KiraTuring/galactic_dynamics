import numpy as np
import matplotlib.pyplot as plt

from jam_fit.utils import logistic, sigma_ratio_logistic, stellar_mass, rotate_moments
from schwarzpy.model import read_nn_mer, read_ap, plot_vel_ellipsoid
from jam_fit import JamInput, jam_model


def plot_anisotropy(nn_mer_file, km_in_arcsec=0.05, direction='theta',
                    fig=None, ax=None):
    """Plot Schwarzschild anisotropy profile (sigma_r/sigma_t) in the
    meridional plane.

    Parameters
    ----------
    nn_mer_file : str
        Path to nn_mer.out from Schwarzschild.
    km_in_arcsec : float
        Length conversion: km per arcsec.
    direction : str
        'theta', 'phi', or 'both' (tangential average).
    fig, ax : optional
        Existing Figure/Axes.

    Returns
    -------
    rplot : ndarray
        Radius array in arcsec.
    avgplot : ndarray
        Mass-weighted average sigma_ratio vs radius.
    """
    from schwarzpy.model import read_nn_mer

    rvals, thvals, nnmer_out = read_nn_mer(nn_mer_file)
    nEner, nTheta = len(rvals), len(thvals)

    rplot = 10 ** rvals / km_in_arcsec

    sigmar2 = nnmer_out.T[4].reshape(nEner, nTheta)
    if direction == 'theta':
        sigmat2 = nnmer_out.T[5].reshape(nEner, nTheta)
    elif direction == 'phi':
        sigmat2 = (nnmer_out.T[8] - nnmer_out.T[7] ** 2).reshape(nEner, nTheta)
    elif direction == 'both':
        sigmat2_theta = nnmer_out.T[5].reshape(nEner, nTheta)
        sigmat2_phi = (nnmer_out.T[8] - nnmer_out.T[7] ** 2).reshape(nEner, nTheta)
        sigmat2 = (sigmat2_theta + sigmat2_phi) / 2
    else:
        raise ValueError("Invalid direction, should be 'theta', 'phi' or 'both'")

    weights = nnmer_out.T[3].reshape(nEner, nTheta)
    sigma_ratio = np.sqrt(sigmar2 / sigmat2)

    sigmar2_radial_avg = np.average(sigmar2, weights=weights, axis=1)
    sigmat2_radial_avg = np.average(sigmat2, weights=weights, axis=1)
    avgplot = np.sqrt(sigmar2_radial_avg / sigmat2_radial_avg)

    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=200)
    for i in range(nTheta):
        ax.plot(rplot, sigma_ratio[:, i],
                label=r'$\theta = %.f^\circ$' % np.degrees(thvals[i]),
                color=f'C{i}', alpha=0.5)
    ax.plot(rplot, avgplot, '-', color='black', label='Average')
    ax.axhline(1, color='k', alpha=0.5)
    ax.set_xscale('log')
    ax.set_xlabel('Radius (arcsec)')
    if direction == 'theta':
        ax.set_ylabel(r'$\sigma_r/\sigma_\theta$')
    elif direction == 'phi':
        ax.set_ylabel(r'$\sigma_r/\sigma_\phi$')
    elif direction == 'both':
        ax.set_ylabel(r'$\sigma_r/\sigma_t$')
    ax.legend(ncol=2)

    return rplot, avgplot


def plot_2D_anisotropy(nn_mer_path, aperture_path, distMpc=17.78,
                       aperture=None, fig=None, ax=None, color='r', alpha=0.7):
    """Plot 2D velocity ellipsoids from Schwarzschild orbit data in the
    meridional plane.

    Parameters
    ----------
    nn_mer_path : str
        Path to nn_mer.out.
    aperture_path : str
        Path to aperture.dat.
    distMpc : float
        Distance in Mpc.
    aperture : float or None
        Aperture size in arcsec (half-width).  None reads from file.
    fig, ax : optional
        Existing Figure/Axes.

    Returns
    -------
    fig, ax, rplot, thvals
    """

    lcorner, rcorner, rotation, xlen, ylen, pixelscale = read_ap(aperture_path)
    if aperture is None:
        aperture = 0.5 * xlen * pixelscale

    rvals, thvals, orbit_data = read_nn_mer(nn_mer_path)

    pc_in_arcsec = distMpc / 206265 * 1e6
    km_in_arcsec = pc_in_arcsec * 3.086e13
    rplot = 10 ** rvals / km_in_arcsec

    nEner, nTheta = len(rvals), len(thvals)
    orbit_data = orbit_data.reshape(nEner, nTheta, -1)
    scale = np.pi * rplot / 2 / nTheta

    R = rplot[orbit_data[:, :, 0].astype(int) - 1]
    TH = thvals[orbit_data[:, :, 1].astype(int) - 1]
    R3 = R ** 3

    DeltaR3 = np.concatenate((R3[[0], :], np.diff(R3, axis=0)), axis=0)
    cosTH = -np.cos(np.linspace(0, 0.5 * np.pi, nTheta + 1))
    DeltacosTH = np.broadcast_to(np.diff(cosTH), DeltaR3.shape)

    X = R * np.sin(TH)
    Y = R * np.cos(TH)

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    for i in range(nEner):
        for j in range(nTheta):
            index_r = int(orbit_data[i, j, 0] - 1)
            index_th = int(orbit_data[i, j, 1] - 1)
            matrix = np.array([[orbit_data[i, j, 4], orbit_data[i, j, 6]],
                               [orbit_data[i, j, 6], orbit_data[i, j, 5]]])
            phi_ratio = np.sqrt(
                (orbit_data[i, j, 8] - orbit_data[i, j, 7] ** 2)
                / orbit_data[i, j, 4])
            plot_vel_ellipsoid(
                matrix, scale[i],
                center=(rplot[index_r] * np.sin(thvals[index_th]),
                        rplot[index_r] * np.cos(thvals[index_th])),
                aperture=(aperture, aperture),
                theta=thvals[index_th], show_axes=False, color=color,
                ax=ax, alpha=alpha)

    obs_density = orbit_data[:, :, 2] / DeltaR3 / DeltacosTH
    mod_density = orbit_data[:, :, 3] / DeltaR3 / DeltacosTH
    vmin, vmax = np.percentile(np.log10(obs_density), [1, 99])
    levels = np.linspace(vmin, vmax, 30)
    ax.contour(X, Y, np.log10(obs_density), colors='k', levels=levels,
               linestyles='solid')
    ax.contour(X, Y, np.log10(mod_density), colors=color, levels=levels,
               linestyles='solid', alpha=alpha)
    ax.set_xlabel('R (arcsec)')
    ax.set_ylabel('Z (arcsec)')

    return fig, ax, rplot, thvals


def plot_2D_anisotropy_logi(r, th, ratio0, ratio1, r0, alpha,
                            fig=None, ax=None, color='b', alpha_plot=0.7):
    """Plot 2D velocity ellipsoids from a JAM logistic-anisotropy profile.

    Parameters
    ----------
    r : ndarray
        Radial grid in arcsec.
    th : ndarray
        Angular grid in radians.
    ratio0, ratio1, r0, alpha : float
        Logistic anisotropy parameters (see sigma_ratio_logistic).
    fig, ax : optional
        Existing Figure/Axes.

    Returns
    -------
    fig, ax
    """

    scale = np.pi * r / 2 / len(th)
    logi_anisotropy = sigma_ratio_logistic(r, ratio0, ratio1, r0, alpha)

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    for i, r_val in enumerate(r):
        for j, th_val in enumerate(th):
            matrix = np.array([[logi_anisotropy[i], 0], [0, 1]])
            plot_vel_ellipsoid(
                matrix, scale[i],
                center=(r_val * np.sin(th_val),
                        r_val * np.cos(th_val)),
                theta=th_val,
                show_axes=False, color=color, ax=ax, alpha=alpha_plot)

    ax.set_xlabel('R (arcsec)')
    ax.set_ylabel('Z (arcsec)')

    return fig, ax


def sigma_ratio_const(r, ratio, thvals):
    r2, th2, rth = rotate_moments(ratio**2, 1, 0, thvals)
    ratio_proj = np.sqrt(r2/th2)
    return ratio_proj


def plot_2D_anisotropy_const(r, th, ratio, align='cyl',
                             fig=None, ax=None, color='b', alpha=0.7):
    """Plot 2D velocity ellipsoids from a JAM logistic-anisotropy profile.

    Parameters
    ----------
    r : ndarray
        Radial grid in arcsec.
    th : ndarray
        Angular grid in radians.
    ratio : float
        Constant anisotropy ratio sigma_r/sigma_theta.
    fig, ax : optional
        Existing Figure/Axes.

    Returns
    -------
    fig, ax
    """

    scale = np.pi * r / 2 / len(th)
    r_mesh, th_mesh = np.meshgrid(r, th, indexing='ij')

    if align == 'cyl':
        TH = th_mesh.flatten()-np.pi/2
    elif align == 'sph':
        TH = np.zeros_like(r_mesh.flatten())
    else:
        raise ValueError("Invalid align, should be 'cyl' or 'sph'")
    r2, th2, rth = rotate_moments(ratio**2, 1, 0, TH)
    r2 = r2.reshape(r_mesh.shape)
    th2 = th2.reshape(r_mesh.shape)
    rth = rth.reshape(r_mesh.shape)

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    for i, r_val in enumerate(r):
        for j, th_val in enumerate(th):
            matrix = np.array([[r2[i, j], rth[i, j]], [rth[i, j], th2[i, j]]])
            plot_vel_ellipsoid(
                matrix, scale[i],
                center=(r_val * np.sin(th_val),
                        r_val * np.cos(th_val)),
                theta=th_val,
                show_axes=False, color=color, ax=ax, alpha=alpha)

    ax.set_xlabel('R (arcsec)')
    ax.set_ylabel('Z (arcsec)')

    return fig, ax


def plot_2D_density(nn_mer_path, aperture_path, distMpc=17.78,
                    aperture=None):
    """Plot 2D intrinsic density from Schwarzschild in the meridional plane.

    Parameters
    ----------
    nn_mer_path : str
        Path to nn_mer.out.
    aperture_path : str
        Path to aperture.dat.
    distMpc : float
        Distance in Mpc.
    aperture : float or None
        Half-width in arcsec.  None reads from file.
    """
    from schwarzpy.model import read_nn_mer, read_ap

    lcorner, rcorner, rotation, xlen, ylen, pixelscale = read_ap(aperture_path)
    aperture = 0.5 * xlen * pixelscale if aperture is None else aperture
    rvals, thvals, orbit_data = read_nn_mer(nn_mer_path)

    pc_in_arcsec = distMpc / 206265 * 1e6
    km_in_arcsec = pc_in_arcsec * 3.086e13
    rplot = 10 ** rvals / km_in_arcsec

    nEner, nTheta = len(rvals), len(thvals)
    orbit_data = orbit_data.reshape(nEner, nTheta, -1)
    scale = np.pi * rplot / 2 / nTheta

    R = rplot[orbit_data[:, :, 0].astype(int) - 1]
    TH = thvals[orbit_data[:, :, 1].astype(int) - 1]
    R3 = R ** 3

    DeltaR3 = np.concatenate((R3[[0], :], np.diff(R3, axis=0)), axis=0)
    cosTH = -np.cos(np.linspace(0, 0.5 * np.pi, nTheta + 1))
    DeltacosTH = np.broadcast_to(np.diff(cosTH), DeltaR3.shape)

    X = R * np.sin(TH)
    Y = R * np.cos(TH)

    obs_density = orbit_data[:, :, 4] / DeltaR3 / DeltacosTH
    mod_density = orbit_data[:, :, 4] / DeltaR3 / DeltacosTH
    levels = np.linspace(np.nanmin(np.log10(obs_density)),
                         np.nanmax(np.log10(obs_density)), 30)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.contour(X, Y, np.log10(obs_density), colors='k', levels=levels,
               linestyles='solid')
    ax.contour(X, Y, np.log10(mod_density), colors='r', levels=levels,
               linestyles='solid')
    im = ax.contourf(X, Y, np.log10(obs_density))
    plt.colorbar(im, label='log10(mod/obs)')
    ax.set_xlabel('R (arcsec)')
    ax.set_ylabel('Z (arcsec)')
    ax.set_title('Intrinsic Density Comparison')
    ax.set_aspect('equal')
    ax.set_xlim(0, aperture)
    ax.set_ylim(0, aperture)


def plot_2D_moment(nn_mer_path, aperture_path, aperture=None,
                   jam_input=None, pars_dict=None,
                   moment_names=None, show_ratio=False):
    """Plot 2D comparison of Schwarzschild vs JAM intrinsic moments in the
    meridional plane.

    Side-by-side panels: [AxiSchw | JAM | JAM - Schw] for each moment.

    Parameters
    ----------
    nn_mer_path : str
        Path to nn_mer.out.
    aperture_path : str
        Path to aperture.dat.
    aperture : float or None
        Half-width in arcsec.  None reads from file.
    jam_input : JamInput
        JAM configuration object.
    pars_dict : dict
        JAM parameter values keyed by name.
    moment_names : list of str
        Subset of: 'r2', 'th2', 'phi2', 't2', 'total2', 'density', etc.

    Returns
    -------
    R, TH, schw_moments, jam_moments, R_jam, TH_jam, jam_moments_intp
    """
    from schwarzpy.model import read_nn_mer, read_ap
    from jam_fit import jam_model

    if moment_names is None:
        moment_names = ['r2', 'th2', 'phi2']

    lcorner, rcorner, rotation, xlen, ylen, pixelscale = read_ap(aperture_path)
    aperture = 0.5 * xlen * pixelscale if aperture is None else aperture

    rvals, thvals, orbit_data = read_nn_mer(nn_mer_path)

    if jam_input is None or pars_dict is None:
        raise ValueError("jam_input and pars_dict must be provided")
    pc_in_arcsec = jam_input.dist / 206265 * 1e6
    km_in_arcsec = pc_in_arcsec * 3.086e13
    rplot = 10 ** rvals / km_in_arcsec

    orbit_data = orbit_data.reshape(len(rvals), len(thvals), -1)
    R = rplot[orbit_data[:, :, 0].astype(int) - 1]
    TH = thvals[orbit_data[:, :, 1].astype(int) - 1]

    X = R * np.sin(TH)
    Y = R * np.cos(TH)

    R3 = R ** 3
    DeltaR3 = np.concatenate((R3[[0], :], np.diff(R3, axis=0)), axis=0)
    cosTH = -np.cos(np.linspace(0, 0.5 * np.pi, len(thvals) + 1))
    DeltacosTH = np.broadcast_to(np.diff(cosTH), DeltaR3.shape)

    schw_moments = {}
    schw_moments['mass'] = orbit_data[:, :, 3]
    schw_moments['density'] = orbit_data[:, :, 3] / DeltaR3 / DeltacosTH
    schw_moments['r2'] = orbit_data[:, :, 4]
    schw_moments['th2'] = orbit_data[:, :, 5]
    schw_moments['rth'] = orbit_data[:, :, 6]
    schw_moments['phi1'] = orbit_data[:, :, 7]
    schw_moments['phi2'] = orbit_data[:, :, 8]
    schw_moments['phi2sig'] = schw_moments['phi2'] - schw_moments['phi1'] ** 2
    schw_moments['t2'] = (schw_moments['th2'] + schw_moments['phi2sig']) / 2
    schw_moments['total2'] = (schw_moments['r2'] + schw_moments['th2']
                              + schw_moments['phi2sig']) / 3

    jam_kwargs = jam_input.get_kwargs()
    out_model_intr = jam_model(pars_dict, Rbin=X.flatten(), zbin=Y.flatten(),
                               type='intr', **jam_kwargs)

    jam_moments = {}
    flux = out_model_intr.flux.reshape(X.shape)
    jam_moments['mass'] = flux * 10 ** pars_dict['lg_ml'] * DeltaR3 * DeltacosTH
    jam_moments['density'] = flux * 10 ** pars_dict['lg_ml']
    jam_r2, jam_th2, jam_phi2sig, jam_phi2 = out_model_intr.model

    if jam_kwargs['align'] == 'cyl':
        jam_r2, jam_th2, _ = rotate_moments(jam_th2, jam_r2, np.zeros_like(jam_r2), TH.flatten())

    jam_moments['r2'] = jam_r2.reshape(X.shape)
    jam_moments['th2'] = jam_th2.reshape(X.shape)
    jam_moments['phi2sig'] = jam_phi2sig.reshape(X.shape)
    jam_moments['phi2'] = jam_phi2.reshape(X.shape)
    jam_moments['phi1'] = np.sqrt(jam_moments['phi2'] - jam_moments['phi2sig'])
    jam_moments['t2'] = (jam_moments['th2'] + jam_moments['phi2sig']) / 2
    jam_moments['total2'] = (jam_moments['r2'] + jam_moments['th2']
                             + jam_moments['phi2sig']) / 3

    r_jam = np.geomspace(np.min(rplot), np.max(rplot), 50)
    theta_jam = np.linspace(0, 0.5 * np.pi, 50)
    R_jam, TH_jam = np.meshgrid(r_jam, theta_jam, indexing='ij')
    X_jam = R_jam * np.sin(TH_jam)
    Y_jam = R_jam * np.cos(TH_jam)
    out_model_intr_intp = jam_model(pars_dict, Rbin=X_jam.flatten(),
                                    zbin=Y_jam.flatten(),
                                    type='intr', **jam_kwargs)

    jam_moments_intp = {}
    flux_intp = out_model_intr_intp.flux.reshape(X_jam.shape)
    jam_moments_intp['density'] = flux_intp * 10 ** pars_dict['lg_ml']
    jam_r2_i, jam_th2_i, jam_phi2sig_i, jam_phi2_i = out_model_intr_intp.model

    if jam_kwargs['align'] == 'cyl':
        jam_r2_i, jam_th2_i, _ = rotate_moments(jam_th2_i, jam_r2_i, np.zeros_like(jam_r2_i), TH_jam.flatten())

    jam_moments_intp['r2'] = jam_r2_i.reshape(X_jam.shape)
    jam_moments_intp['th2'] = jam_th2_i.reshape(X_jam.shape)
    jam_moments_intp['phi2sig'] = jam_phi2sig_i.reshape(X_jam.shape)
    jam_moments_intp['phi2'] = jam_phi2_i.reshape(X_jam.shape)
    jam_moments_intp['t2'] = (jam_moments_intp['th2']
                              + jam_moments_intp['phi2sig']) / 2
    jam_moments_intp['total2'] = (jam_moments_intp['r2']
                                  + jam_moments_intp['th2']
                                  + jam_moments_intp['phi2sig']) / 3

    from plotbin.sauron_colormap import register_sauron_colormap
    register_sauron_colormap()

    fig, axes = plt.subplots(len(moment_names), 3,
                             figsize=(9, 3 * len(moment_names)),
                             sharex=True, sharey=True, squeeze=False)

    vmin, vmax = np.percentile(np.sqrt(schw_moments['r2']), [0, 100])
    if show_ratio:
        vmin, vmax = 0.0, 2.0

    dmin, dmax = np.percentile(np.log10(schw_moments['density']), [1, 99])
    levels = 30
    cmap = 'sauron'

    for i, moment_name in enumerate(moment_names):
        schw_value = np.sqrt(schw_moments[moment_name])
        jam_value = np.sqrt(jam_moments[moment_name])
        jam_value_intp = np.sqrt(jam_moments_intp[moment_name])

        if show_ratio:
            # schw_rms = np.log10(np.sqrt(schw_moments['r2'])/schw_rms)
            # jam_rms = np.log10(np.sqrt(jam_moments['r2'])/jam_rms)
            # jam_rms_intp = np.log10(np.sqrt(jam_moments_intp['r2'])/jam_rms_intp)
            schw_value = np.sqrt(schw_moments['r2'])/schw_value
            jam_value = np.sqrt(jam_moments['r2'])/jam_value
            jam_value_intp = np.sqrt(jam_moments_intp['r2'])/jam_value_intp

        ax = axes[i, 0]
        im = ax.contourf(X, Y, schw_value, levels=levels,
                         cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax, label='km/s')
        ax.scatter(X.flatten(), Y.flatten(), c='k', s=5, zorder=10)
        ax.contour(X, Y, np.log10(schw_moments['density']),
                   colors='k', alpha=0.5, linestyles='solid',
                   levels=levels, vmin=dmin, vmax=dmax)
        ax.set_aspect('equal')
        ax.set_xlim(0, aperture)
        ax.set_ylim(0, aperture)
        if show_ratio:
            ax.set_ylabel(r'$\sigma_{\rm r2}/\sigma_{\rm %s}$' % moment_name)
        else:
            ax.set_ylabel(r'Intrinsic <V$_{\rm %s}$>' % moment_name)
        if i == 0:
            ax.set_title('AxiSchw')

        ax = axes[i, 1]
        im = ax.contourf(X_jam, Y_jam,
                         jam_value_intp,
                         levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax, label='km/s')
        ax.scatter(X.flatten(), Y.flatten(), c='k', s=5, zorder=10)
        ax.contour(X_jam, Y_jam,
                   np.log10(jam_moments_intp['density']),
                   colors='k', alpha=0.5, linestyles='solid',
                   levels=levels, vmin=dmin, vmax=dmax)
        ax.set_aspect('equal')
        ax.set_xlim(0, aperture)
        ax.set_ylim(0, aperture)
        if i == 0:
            ax.set_title('JAM')

        ax = axes[i, 2]
        max_val = np.max(np.abs(jam_value - schw_value))
        im = ax.contourf(X, Y, jam_value - schw_value,
                         levels=levels, cmap=cmap,
                         vmin=-max_val, vmax=max_val)
        plt.colorbar(im, ax=ax)
        ax.scatter(X.flatten(), Y.flatten(), c='k', s=5, zorder=10)
        ax.contour(X, Y, np.log10(jam_moments['density']),
                   colors='k', alpha=0.5, linestyles='solid',
                   levels=levels, vmin=dmin, vmax=dmax)
        ax.set_aspect('equal')
        ax.set_xlim(0, aperture)
        ax.set_ylim(0, aperture)
        if i == 0:
            ax.set_title('JAM - Schw')

    return R, TH, schw_moments, jam_moments, R_jam, TH_jam, jam_moments_intp
