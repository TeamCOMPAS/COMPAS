from scipy.interpolate import interp1d
from .bbh_population import BBHPopulation
from .cosmological_model import CosmologicalModel
from .snr_grid import SNRGrid
import matplotlib.pyplot as plt
from tqdm.auto import trange
from .gpu_utils import xp




def binned_detection_rates_from_compas_output(
        compas_path, cosmological_parameters={}, max_detectable_redshift=1.0,
        chirp_mass_bins=100, redshift_bins=100,
):
    bbh_population = BBHPopulation(compas_path)
    cosmological_model = CosmologicalModel(**cosmological_parameters)
    snr_grid = SNRGrid()

    sorted_idx = xp.argsort(bbh_population.chirp_mass)
    redshift = cosmological_model.redshift


    if chirp_mass_bins is None:
        chirp_mass_bins = bbh_population.chirp_mass[sorted_idx]
    elif isinstance(chirp_mass_bins, int):
        chirp_mass_bins = xp.linspace(3, 40, chirp_mass_bins)

    if redshift_bins is None:
        redshift_bins = redshift[redshift < max_detectable_redshift]
    elif isinstance(redshift_bins, int):
        redshift_bins = xp.linspace(0, max_detectable_redshift, redshift_bins)

    grid = _compute_binned_detection_rates(
        bbh_population, cosmological_model, snr_grid,
        max_detectable_redshift=max_detectable_redshift,
        chirp_mass_bins=chirp_mass_bins,
        redshift_bins=redshift_bins,
    )
    return grid, chirp_mass_bins, redshift_bins


def _compute_binned_detection_rates(
        bbh_population: BBHPopulation,
        cosmological_model: CosmologicalModel,
        snr_grid: SNRGrid,
        chirp_mass_bins: xp.ndarray,
        redshift_bins: xp.ndarray,
        max_detectable_redshift: float = 1.0,
):
    n_formed = cosmological_model.sfr / bbh_population.avg_sf_mass_needed  # Divide the star formation rate density by the representative SF mass

    # calculate the formation and merger rates using what we computed above
    detection_rate = _find_merger_rates(
        redshifts=cosmological_model.redshift,
        times=cosmological_model.time,
        time_first_SF=cosmological_model.time_first_star_formation,
        n_formed=n_formed,
        dPdlogZ=cosmological_model.dPdlogZ,
        metallicities=cosmological_model.metallicities,
        p_draw_metallicity=cosmological_model.p_draw_metallicity,
        COMPAS_metallicites=bbh_population.z_zams,
        COMPAS_delay_times=bbh_population.t_delay,
        COMPAS_Mc=bbh_population.chirp_mass,
        COMPAS_eta=bbh_population.eta,
        distances=cosmological_model.distance,
        snr_grid_at_1Mpc=snr_grid.snr_grid_at_1Mpc,
        detection_probability_from_snr=snr_grid.pdetection,
        snr_Mc_bins=snr_grid.chirp_mass,
        snr_eta_bins=snr_grid.eta,
        snr_bins=snr_grid.snr,
        comoving_vol=cosmological_model.comoving_volume,
        max_detectable_redshift=max_detectable_redshift,
        redshift_bins=redshift_bins,
        chirp_mass_bins=chirp_mass_bins,
    )

    return detection_rate


def _find_merger_rates(
        # COMPAS parameters
        COMPAS_metallicites: xp.array,
        COMPAS_delay_times: xp.array,
        COMPAS_Mc: xp.array,
        COMPAS_eta: xp.array,
        # cosmological parameters
        redshifts:xp.array,
        times:xp.array,
        time_first_SF:float,
        n_formed:int,
        dPdlogZ:xp.ndarray,
        metallicities:xp.array,
        p_draw_metallicity:float,
        distances:xp.array,
        comoving_vol: xp.array,
        # detection parameters
        snr_grid_at_1Mpc:xp.ndarray,
        detection_probability_from_snr:xp.array,
        snr_Mc_bins:xp.array,
        snr_eta_bins:xp.array,
        snr_bins:xp.array,
        max_detectable_redshift:float,
        # binning parameters
        chirp_mass_bins:xp.array,
        redshift_bins:xp.array,
)->xp.ndarray:
    n_binaries = len(COMPAS_Mc)

    # Create empty arrays to store the results
    formation_rate = xp.zeros(len(redshifts))
    idx_max_z = xp.digitize(max_detectable_redshift, redshifts)
    merger_rate = xp.zeros(idx_max_z)
    detection_probability = xp.zeros(idx_max_z)
    det_matrix = xp.zeros(shape=(len(COMPAS_Mc), len(redshifts[:idx_max_z])))
    # det_matrix = xp.zeros(shape=(len(chirp_mass_bins), len(redshift_bins)))

    # create time->redshift interpolator
    times_to_redshifts = interp1d(times, redshifts)

    # sort the COMPAS data by chirp mass
    sorted_idx = xp.argsort(COMPAS_Mc)
    COMPAS_Mc = COMPAS_Mc[sorted_idx]
    COMPAS_eta = COMPAS_eta[sorted_idx]
    COMPAS_metallicites = COMPAS_metallicites[sorted_idx]
    COMPAS_delay_times = COMPAS_delay_times[sorted_idx]


    # go through each binary in the COMPAS data
    for i in trange(n_binaries, desc='Computing detection rates'):
        # calculate formation rate (see Neijssel+19 Section 4) - note this uses dPdlogZ for *closest* metallicity
        observed_metallicity_mask = xp.digitize(COMPAS_metallicites[i], metallicities)
        formation_rate[:] = n_formed * dPdlogZ[:, observed_metallicity_mask] / p_draw_metallicity

        # calculate the time at which the binary formed if it merges at this redshift
        time_of_formation = times - COMPAS_delay_times[i]

        # we have only calculated formation rate up to z=max(redshifts),
        # so we need to only find merger rates for formation times at z<max(redshifts)
        # first locate the index above which the binary would have formed before z=max(redshifts)
        idx_earliest_tform = xp.digitize(time_first_SF, time_of_formation)
        if idx_earliest_tform > idx_max_z:
            idx_earliest_tform = idx_max_z + 1

        will_merge = idx_earliest_tform > 0

        if will_merge:
            # work out the redshift at the time of formation
            z_of_formation = times_to_redshifts(time_of_formation[:idx_earliest_tform - 1])
            z_of_formation_index = xp.digitize(z_of_formation, redshifts)
            # set the merger rate at z (with z<10) to the formation rate at z_form
            if idx_earliest_tform > idx_max_z:
                idx_earliest_tform = idx_max_z + 1
            merger_rate[:idx_earliest_tform - 1] = formation_rate[z_of_formation_index][:idx_max_z]
        else:
            merger_rate[:] = 0

        # get bin indices for Mc and eta
        eta_index = xp.digitize(COMPAS_eta[i], snr_eta_bins)
        Mc_index = xp.digitize(COMPAS_Mc[i] * (1 + redshifts[:idx_max_z]), snr_Mc_bins)

        # lookup SNRs using the eta and Mc indices
        snrs = snr_grid_at_1Mpc[eta_index, Mc_index] / distances[:idx_max_z]
        idx = xp.clip(xp.digitize(snrs, snr_bins), 0, len(snr_bins) - 1)
        detection_probability[:] = detection_probability_from_snr[idx]

        det_matrix[i, :] = merger_rate * detection_probability * comoving_vol[:idx_max_z]
    return det_matrix


def plot_detection_rate_matrix(
        detection_rate,
        chirp_masses,
        redshifts,
):
    z, mc, rate2d = redshifts, chirp_masses, detection_rate
    low_mc, high_mc = xp.min(mc), xp.max(mc)
    low_z, high_z = xp.min(z), xp.max(z)
    n_events_per_year = xp.nansum(detection_rate)
    chirp_mass_rate = xp.sum(detection_rate, axis=1)
    redshift_rate = xp.sum(detection_rate, axis=0)
    rate2d = xp.exp(xp.log(rate2d) - xp.max(xp.log(rate2d)))
    rate2d = rate2d / xp.max(rate2d)
    quantiles = {
        '1sig': xp.quantile(rate2d, [0.997])[0],
        '2sig': xp.quantile(rate2d, [0.95])[0],
        '3sig': xp.quantile(rate2d, [0.68])[0]
    }
    min_q, max_q = '3sig', '2sig'
    rate2d[rate2d < quantiles[min_q]] = quantiles[min_q]
    rate2d[rate2d > quantiles[max_q]] = quantiles[max_q]

    mc_range = [low_mc, high_mc]
    z_range = [low_z, high_z]

    fig = plt.figure(figsize=(5, 5))
    gs = plt.GridSpec(4, 4)

    ax_2d = fig.add_subplot(gs[1:4, 0:3])
    ax_top = fig.add_subplot(gs[0, 0:3])
    ax_right = fig.add_subplot(gs[1:4, 3])

    zz, mcc = xp.meshgrid(z, mc)
    cbar = ax_2d.pcolormesh(
        zz,
        mcc,
        rate2d,
        cmap='viridis',
        norm="linear",
        vmax=quantiles[max_q],
        vmin=quantiles[min_q],
    )

    ax_2d.set_xlabel("Redshift")
    ax_2d.set_ylabel("Chirp mass ($M_{\odot}$)")
    ax_2d.set_facecolor("black")
    annote = f"Grid: {rate2d.T.shape}\nN det: {n_events_per_year:.2f}/yr"
    ax_2d.annotate(
        annote,
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
    cbar_ax.set_yticklabels([min_q, max_q])
    cbar_ax.set_yticks([quantiles[min_q], quantiles[max_q]])

    # set all cbar_ax spline color to white
    for spine in cbar_ax.spines.values():
        spine.set_edgecolor("white")

    fig.tight_layout()
    return fig
