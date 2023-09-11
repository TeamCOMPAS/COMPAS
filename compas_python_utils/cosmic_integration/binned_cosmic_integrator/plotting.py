import numpy as np
import matplotlib.pyplot as plt
from corner import corner
from typing import List, Tuple
import warnings
from .conversions import m1_m2_to_eta_chirp_mass

CMAP = 'inferno'


def plot_detection_rate_matrix(
        detection_rate: np.ndarray,
        chirp_masses: np.array,
        redshifts: np.array,
        normalise: bool = True,
        annotation: str = None,
) -> plt.Figure:
    """
    Plot the detection rate matrix as a 2D heatmap with marginal histograms.
    """
    warnings.filterwarnings("ignore")
    assert detection_rate.shape == (len(chirp_masses), len(redshifts))
    z, mc, rate2d = redshifts, chirp_masses, detection_rate
    low_mc, high_mc = np.min(mc), np.max(mc)
    low_z, high_z = np.min(z), np.max(z)

    chirp_mass_rate = np.sum(detection_rate, axis=1)
    redshift_rate = np.sum(detection_rate, axis=0)
    norm_kwgs = dict()
    if normalise:
        rate2d = np.exp(np.log(rate2d) - np.max(np.log(rate2d)))
        rate2d = rate2d / np.max(rate2d)
        quantiles = {
            '1σ': np.quantile(rate2d, [0.997])[0],
            '2σ': np.quantile(rate2d, [0.95])[0],
            '3σ': np.quantile(rate2d, [0.68])[0]
        }
        min_q, max_q = '3σ', '2σ'
        rate2d[rate2d < quantiles[min_q]] = quantiles[min_q]
        rate2d[rate2d > quantiles[max_q]] = quantiles[max_q]
        norm_kwgs = dict(vmax=quantiles[max_q], vmin=quantiles[min_q])

    mc_range = [low_mc, high_mc]
    z_range = [low_z, high_z]

    fig = plt.figure(figsize=(5, 5))
    gs = plt.GridSpec(4, 4)

    ax_2d = fig.add_subplot(gs[1:4, 0:3])
    ax_top = fig.add_subplot(gs[0, 0:3])
    ax_right = fig.add_subplot(gs[1:4, 3])

    zz, mcc = np.meshgrid(z, mc)
    cbar = ax_2d.pcolormesh(
        zz,
        mcc,
        rate2d,
        cmap=CMAP,
        norm="linear",
        **norm_kwgs
    )

    ax_2d.set_xlabel("Redshift")
    ax_2d.set_ylabel("Chirp mass ($M_{\odot}$)")
    ax_2d.set_facecolor("black")

    if annotation is None:
        n_events_per_year = np.nansum(detection_rate)
        annotation = f"Grid: {rate2d.T.shape}\nN det: {n_events_per_year:.2f}/yr"
    ax_2d.annotate(
        annotation,
        xy=(1, 0),
        xycoords="axes fraction",
        xytext=(-5, 5),
        textcoords="offset points",
        ha="right",
        va="bottom",
        color="white",
    )
    kwgs = dict(color="red", lw=1)
    ax_right.step(
        chirp_mass_rate,
        chirp_masses,
        **kwgs,
    )
    ax_top.step(
        redshifts,
        redshift_rate,
        **kwgs,
    )
    ax_right.axis("off")
    ax_top.axis("off")

    ax_right.set_ylim(*mc_range)
    ax_2d.set_ylim(*mc_range)
    ax_top.set_xlim(*z_range)
    ax_2d.set_xlim(*z_range)

    # remove space between subplots
    fig.subplots_adjust(hspace=0, wspace=0)

    cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.6])
    fig.colorbar(
        cbar,
        cax=cbar_ax,
        orientation="vertical",
        label="Rate (yr$^{-1}$)",
        aspect=5,
    )
    cbar_ax.tick_params(labelsize=8, length=0)
    cbar_ax.yaxis.set_ticks_position("left")
    if normalise:
        cbar_ax.set_yticklabels([min_q, max_q])
        cbar_ax.set_yticks([quantiles[min_q], quantiles[max_q]])

    # set all cbar_ax spline color to white
    for spine in cbar_ax.spines.values():
        spine.set_edgecolor("white")

    fig.subplots_adjust(hspace=0, wspace=0)
    return fig


def plot_sfr_and_metallicity(
        redshift: np.array, sfr: np.array,
        metallicities: np.array, dPdlogZ: np.ndarray,
        p_draw_metallicity: np.array,
        metallicity_label: str,
        sf_label: str,
        redshift_range: List, logZ_range: List,
) -> plt.Figure:
    fig, axes = plt.subplots(3, 1, figsize=(5, 8))
    ax = axes[0]
    ax.plot(redshift, sfr, label="SFR")
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"SFR [$M_{\odot}/\rm{yr}/\rm{Mpc}^3$]")
    ax.set_yscale("log")
    ax.set_xlim(*redshift_range)
    txt = sf_label.replace(" ", "\n")
    ax.text(0.05, 0.5, txt, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', horizontalalignment='left')

    ax = axes[1]
    dPdlogZ = dPdlogZ / np.max(dPdlogZ)
    im = ax.imshow(dPdlogZ, extent=[*redshift_range, *logZ_range], aspect="auto", cmap=CMAP)
    ax.set_xlabel("Redshift")
    ax.set_ylabel("logZ")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\frac{dP}{d\log Z}$", rotation=0, labelpad=-10)
    # set the color bar ticks to be 0 and 1 and label them min and max
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels([r"min", r"max"])
    txt = metallicity_label.split(";")
    txt = [t.replace("N", "\mathcal{N}") for t in txt]
    txt = "\n".join([f"${t}$" for t in txt])
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', horizontalalignment='left', bbox=props)

    ax = axes[2]
    # plot the probability of drawing a metallicity
    ps = [p_draw_metallicity] * len(metallicities)
    ax.fill_between(metallicities, ps, 0, alpha=0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel("logZ")
    ax.set_ylabel("p(Z)")
    fig.tight_layout()
    return fig


def plot_snr_grid(
        snr_grid_at_1Mpc: np.ndarray,
        m1: np.array,
        m2: np.array,
        snr: np.array,
        pdetection: np.array,
        snr_threshold: float,
) -> plt.Figure:
    fig, axes = plt.subplots(3, 1, figsize=(5, 9))
    ax = axes[0]
    im = ax.imshow(
        snr_grid_at_1Mpc,
        extent=[m1.min(), m1.max(), m2.min(), m2.max()],
        aspect="auto",
        origin="lower",
        cmap=CMAP,
        vmin=snr_threshold,
    )
    ax.set_xlabel(r"$m_1$ [$M_{\odot}$]")
    ax.set_ylabel(r"$m_2$ [$M_{\odot}$]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"SNR @ 1 Mpc")

    ax = axes[1]
    eta, mc = m1_m2_to_eta_chirp_mass(m1, m2)
    im = ax.imshow(
        snr_grid_at_1Mpc,
        extent=[eta.min(), eta.max(), mc.min(), mc.max()],
        aspect="auto",
        origin="lower",
        cmap=CMAP,
        vmin=snr_threshold,
    )
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel(r"$\mathcal{M}_{\rm chirp}$ [$M_{\odot}$]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"SNR @ 1 Mpc")

    ax = axes[2]
    ax.plot(snr, pdetection)
    ax.set_xlabel(r"SNR")
    ax.set_ylabel(r"$P({\rm detection})$")
    # determine pdet at snr threshold
    pdet_at_thresh = pdetection[snr == snr_threshold]
    snr_at_pdet_max = snr[pdetection >= 0.985][0]
    ax.set_xlim(left=0, right=snr_at_pdet_max)
    ax.set_ylim(bottom=pdet_at_thresh, top=1)
    fig.tight_layout()
    return plt.gcf()


def plot_bbh_population(
        data: np.ndarray, params: List[str]
) -> plt.Figure:
    n_sys = len(data)
    # mask out the outliers
    mask = np.array([True] * len(data))
    for i in range(data.shape[1]):
        mask_up = data[:, i] < np.quantile(data[:, i], 0.975)
        mask_down = data[:, i] > np.quantile(data[:, i], 0.025)
        mask = mask & mask_up & mask_down
    fig = corner(
        data[mask],
        bins=25,
        labels=params,
        color="tab:orange",
        max_n_ticks=3,
        plot_datapoints=False,
        plot_density=True,
        plot_contours=False,
        fill_contours=False,
        pcolor_kwargs=dict(edgecolors="tab:gray", alpha=1, linewidths=0.05),
        quiet=True,
    )
    axes = fig.get_axes()
    # for the diagonal axes
    for i in range(len(params)):
        ax = axes[i * (len(params) + 1)]
        ax.set_title(params[i])
        for splines in ax.spines.values():
            splines.set_visible(False)
        ax.spines["bottom"].set_visible(True)

    fig.suptitle(f"BBH Population ({n_sys:,} BBHs)")
    return fig
