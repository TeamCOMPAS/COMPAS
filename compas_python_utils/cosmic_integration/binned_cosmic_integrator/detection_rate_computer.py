from tqdm.auto import trange
from scipy.interpolate import interp1d
import numpy as np

from .bbh_population import BBHPopulation
from .cosmological_model import CosmologicalModel
from .snr_grid import SNRGrid
from .gpu_utils import xp
from .bin_2d_data import bin_2d_data


def compute_binned_detection_rates(
        bbh_population: BBHPopulation,
        cosmological_model: CosmologicalModel,
        snr_grid: SNRGrid,
        chirp_mass_bins: np.ndarray,
        redshift_bins: np.ndarray,
        max_detectable_redshift: float = 1.0,
        verbose=True
) -> np.ndarray:
    """
    Compute the detection rate matrix for a given BBH population and cosmological model.

    If the GPU is available, this function moves the np.ndarray objects to the GPU and perform the computation there.
    If the GPU is not available, this function will perform the computation on the CPU.

    """
    n_formed = cosmological_model.sfr / bbh_population.avg_sf_mass_needed
    # Divide the star formation rate density by the representative SF mass

    # calculate the formation and merger rates using what we computed above
    # return as np array
    return _find_merger_rates(
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
        verbose=verbose
    )


def _find_merger_rates(
        # COMPAS parameters
        COMPAS_metallicites: xp.array,
        COMPAS_delay_times: xp.array,
        COMPAS_Mc: xp.array,
        COMPAS_eta: xp.array,
        # cosmological parameters
        redshifts: xp.array,
        times: xp.array,
        time_first_SF: float,
        n_formed: int,
        dPdlogZ: xp.ndarray,
        metallicities: xp.array,
        p_draw_metallicity: float,
        distances: xp.array,
        comoving_vol: xp.array,
        # detection parameters
        snr_grid_at_1Mpc: xp.ndarray,
        detection_probability_from_snr: xp.array,
        snr_Mc_bins: xp.array,
        snr_eta_bins: xp.array,
        snr_bins: xp.array,
        max_detectable_redshift: float,
        # binning parameters
        chirp_mass_bins: xp.array,
        redshift_bins: xp.array,
        #
        verbose: bool = True
) -> xp.ndarray:
    """
    Compute the detection rate matrix for a given BBH population and cosmological model.
    :return:
        np.ndarray of shape (len(chirp_mass_bins), len(redshift_bins))
    """
    n_binaries = len(COMPAS_Mc)
    # sort the COMPAS data by chirp mass
    sorted_idx = xp.argsort(COMPAS_Mc)
    COMPAS_Mc = COMPAS_Mc[sorted_idx]
    COMPAS_eta = COMPAS_eta[sorted_idx]
    COMPAS_metallicites = COMPAS_metallicites[sorted_idx]
    COMPAS_delay_times = COMPAS_delay_times[sorted_idx]

    # Create empty arrays to store the results
    formation_rate = xp.zeros(len(redshifts))
    idx_max_z = xp.digitize(max_detectable_redshift, redshifts, right=True)
    merger_rate = xp.zeros(idx_max_z)
    detection_probability = xp.zeros(idx_max_z)
    det_matrix = xp.zeros(shape=(len(chirp_mass_bins), idx_max_z))

    # create time->redshift interpolator
    times_to_redshifts = interp1d(times, redshifts)

    # go through each binary in the COMPAS data
    for i in trange(
            n_binaries,
            desc='Computing detection rates',
            disable=not verbose
    ):
        mc_bin_idx = xp.digitize(COMPAS_Mc[i], chirp_mass_bins, right=True)
        if mc_bin_idx == len(chirp_mass_bins):
            # skip binaries with chirp mass above the highest bin
            continue

        # calculate formation rate (see Neijssel+19 Section 4) - note this uses dPdlogZ for *closest* metallicity
        observed_metallicity_mask = xp.digitize(COMPAS_metallicites[i], metallicities, right=True)
        formation_rate[:] = n_formed * dPdlogZ[:, observed_metallicity_mask] / p_draw_metallicity

        # calculate the time at which the binary formed if it merges at this redshift
        time_of_formation = times - COMPAS_delay_times[i]

        # we have only calculated formation rate up to z=max(redshifts),
        # so we need to only find merger rates for formation times at z<max(redshifts)
        # first locate the index above which the binary would have formed before z=max(redshifts)
        idx_earliest_tform = xp.digitize(time_first_SF, time_of_formation, right=True)
        if idx_earliest_tform > idx_max_z:
            idx_earliest_tform = idx_max_z + 1

        will_merge = idx_earliest_tform > 0

        if will_merge:
            # work out the redshift at the time of formation
            z_of_formation = times_to_redshifts(time_of_formation[:idx_earliest_tform - 1])
            z_of_formation_index = xp.digitize(z_of_formation, redshifts, right=True)
            # set the merger rate at z (with z<10) to the formation rate at z_form
            if idx_earliest_tform > idx_max_z:
                idx_earliest_tform = idx_max_z + 1
            merger_rate[:idx_earliest_tform - 1] = formation_rate[z_of_formation_index][:idx_max_z]
        else:
            merger_rate[:] = 0

        # get bin indices for Mc and eta
        eta_index = xp.digitize(COMPAS_eta[i], snr_eta_bins, right=True)
        Mc_index = xp.digitize(COMPAS_Mc[i] * (1 + redshifts[:idx_max_z]), snr_Mc_bins, right=True)

        # lookup SNRs using the eta and Mc indices
        snrs = snr_grid_at_1Mpc[eta_index, Mc_index] / distances[:idx_max_z]
        idx = xp.clip(xp.digitize(snrs, snr_bins, right=True), 0, len(snr_bins) - 1)
        detection_probability[:] = detection_probability_from_snr[idx]

        # add the detection rate to the detection matrix
        det_matrix[mc_bin_idx, :] += merger_rate * detection_probability * comoving_vol[:idx_max_z]

    # bin the detection matrix in redshift
    det_matrix = bin_2d_data(det_matrix, redshifts[:idx_max_z], redshift_bins, axis=1)

    return det_matrix
