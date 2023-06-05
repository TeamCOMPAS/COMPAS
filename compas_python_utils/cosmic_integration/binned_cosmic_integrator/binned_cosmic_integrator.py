from typing import List, Optional, Union
from .bbh_population import BBHPopulation
from astropy.cosmology import WMAP9 as cosmology
from .. import selection_effects
from scipy.interpolate import interp1d
import numpy as np

def get_binned_detection_rates(
        compas_output_path: str,
        m1_min: Optional[float] = 5, m1_max: Optional[float] = 150,
        m2_min: Optional[float] = 0.1,
        binary_fraction: Optional[float] = 1.0,
        chirp_mass_bins: Union[int, List[float]] = 100,
        redshift_bins: Union[int, List[float]] = 100,
):
    bbh_population = BBHPopulation(compas_output_path, m1_min, m1_max, m2_min, binary_fraction)

    # calculate the redshifts array and its equivalents
    redshifts, n_redshifts_detection, times, time_first_SF, distances, shell_volumes = calculate_redshift_related_params(
        max_redshift, max_redshift_detection, redshift_step, z_first_SF)

    # find the star forming mass per year per Gpc^3 and convert to total number formed per year per Gpc^3
    sfr = find_sfr(redshifts, a=aSF, b=bSF, c=cSF, d=dSF)  # functional form from Madau & Dickinson 2014

    # Calculate the representative SF mass
    Average_SF_mass_needed = (bbh_population.mass_evolved_per_binary * bbh_population.n_systems)
    print('Average_SF_mass_needed = ',
          Average_SF_mass_needed)  # print this, because it might come in handy to know when writing up results :)
    n_formed = sfr / Average_SF_mass_needed  # Divide the star formation rate density by the representative SF mass

    # work out the metallicity distribution at each redshift and probability of drawing each metallicity in COMPAS
    dPdlogZ, metallicities, p_draw_metallicity = find_metallicity_distribution(
        redshifts, min_logZ_COMPAS=np.log(
        np.min(COMPAS.initialZ)),
           max_logZ_COMPAS=np.log(
               np.max(COMPAS.initialZ)),
           mu0=mu0, muz=muz, sigma_0=sigma0,
           sigma_z=sigmaz, alpha=alpha,
           min_logZ=min_logZ, max_logZ=max_logZ,
           step_logZ=step_logZ)

    # calculate the formation and merger rates using what we computed above
    formation_rate, merger_rate = find_formation_and_merger_rates(n_binaries, redshifts, times, time_first_SF, n_formed,
                                                                  dPdlogZ,
                                                                  metallicities, p_draw_metallicity,
                                                                  COMPAS.metallicitySystems,
                                                                  COMPAS.delayTimes, COMPAS.sw_weights)

    # create lookup tables for the SNR at 1Mpc as a function of the masses and the probability of detection as a function of SNR
    snr_grid_at_1Mpc, detection_probability_from_snr = compute_snr_and_detection_grids(sensitivity, snr_threshold,
                                                                                       Mc_max, Mc_step,
                                                                                       eta_max, eta_step, snr_max,
                                                                                       snr_step)

    # use lookup tables to find the probability of detecting each binary at each redshift
    detection_probability = find_detection_probability(chirp_masses, etas, redshifts, distances, n_redshifts_detection,
                                                       n_binaries,
                                                       snr_grid_at_1Mpc, detection_probability_from_snr, Mc_step,
                                                       eta_step, snr_step)

    # finally, compute the detection rate using Neijssel+19 Eq. 2
    detection_rate = np.zeros(shape=(n_binaries, n_redshifts_detection))
    detection_rate = merger_rate[:, :n_redshifts_detection] * detection_probability \
                     * shell_volumes[:n_redshifts_detection] / (1 + redshifts[:n_redshifts_detection])

    return detection_rate




def find_formation_and_merger_rates(n_binaries, redshifts, times, time_first_SF, n_formed, dPdlogZ, metallicities,
                                    p_draw_metallicity,
                                    COMPAS_metallicites, COMPAS_delay_times, COMPAS_weights=None):
    """
        Find both the formation and merger rates for each binary at each redshift

        Args:
            n_binaries          --> [int]            Number of DCO binaries in the arrays
            redshifts           --> [list of floats] Redshifts at which to evaluate the rates
            times               --> [list of floats] Equivalent of the redshifts in terms of age of the Universe
            n_formed            --> [float]          Binary formation rate (number of binaries formed per year per cubic Gpc) represented by each simulated COMPAS binary
            dPdlogZ             --> [2D float array] Probability of getting a particular logZ at a certain redshift
            metallicities       --> [list of floats] Metallicities at which dPdlogZ is evaluated; if this is None, assume that metallicity weighting is 1 (corresponds to all SFR happening at one fixed metallicity)
            p_draw_metallicity  --> [float]          Probability of drawing a certain metallicity in COMPAS (float because assuming uniform)
            COMPAS_metallicites --> [list of floats] Metallicity of each binary in COMPAS data
            COMPAS_delay_times  --> [list of floats] Delay time of each binary in COMPAS data
            COMPAS_weights      --> [list of floats] Adaptive sampling weights for each binary in COMPAS data (defaults to all 1s for unweighted samples)

        Returns:
            formation_rate      --> [2D float array] Formation rate for each binary at each redshift
            merger_rate         --> [2D float array] Merger rate for each binary at each redshift
    """
    # check if weights were provided, if not use uniform weights
    if COMPAS_weights is None:
        COMPAS_weights = np.ones(n_binaries)

    # initalise rates to zero
    n_redshifts = len(redshifts)
    redshift_step = redshifts[1] - redshifts[0]
    formation_rate = np.zeros(shape=(n_binaries, n_redshifts))
    merger_rate = np.zeros(shape=(n_binaries, n_redshifts))

    # interpolate times and redshifts for conversion
    times_to_redshifts = interp1d(times, redshifts)

    # make note of the first time at which star formation occured
    age_first_sfr = time_first_SF

    # go through each binary in the COMPAS data
    for i in range(n_binaries):
        # if metallicities array is None, assume all SFR happened at one fixed metallicity
        if metallicities is None:
            formation_rate[i, :] = n_formed * COMPAS_weights[i]
        # calculate formation rate (see Neijssel+19 Section 4) - note this uses dPdlogZ for *closest* metallicity
        else:
            formation_rate[i, :] = n_formed * dPdlogZ[:,
                                              np.digitize(COMPAS_metallicites[i], metallicities)] / p_draw_metallicity * \
                                   COMPAS_weights[i]

        # calculate the time at which the binary formed if it merges at this redshift
        time_of_formation = times - COMPAS_delay_times[i]

        # we have only calculated formation rate up to z=max(redshifts), so we need to only find merger rates for formation times at z<max(redshifts)
        # first locate the index above which the binary would have formed before z=max(redshifts)
        first_too_early_index = np.digitize(age_first_sfr, time_of_formation)

        # include the whole array if digitize returns end of array and subtract one so we don't include the time past the limit
        first_too_early_index = first_too_early_index + 1 if first_too_early_index == n_redshifts else first_too_early_index

        # as long as that doesn't preclude the whole range
        if first_too_early_index > 0:
            # work out the redshift at the time of formation
            z_of_formation = times_to_redshifts(time_of_formation[:first_too_early_index - 1])

            # calculate which index in the redshift array these redshifts correspond to
            z_of_formation_index = np.ceil(z_of_formation / redshift_step).astype(int)

            # set the merger rate at z (with z<10) to the formation rate at z_form
            merger_rate[i, :first_too_early_index - 1] = formation_rate[i, z_of_formation_index]
    return formation_rate, merger_rate


def compute_snr_and_detection_grids(sensitivity="O1", snr_threshold=8.0, Mc_max=300.0, Mc_step=0.1,
                                    eta_max=0.25, eta_step=0.01, snr_max=1000.0, snr_step=0.1):
    """
        Compute a grid of SNRs and detection probabilities for a range of masses and SNRs

        These grids are computed to allow for interpolating the values of the snr and detection probability. This function
        combined with find_detection_probability() could be replaced by something like
            for i in range(n_binaries):
                detection_probability = selection_effects.detection_probability(COMPAS.mass1[i],COMPAS.mass2[i],
                                            redshifts, distances, GWdetector_snr_threshold, GWdetector_sensitivity)
        if runtime was not important.

        Args:
            sensitivity                    --> [string]         Which detector sensitivity to use: one of ["design", "O1", "O3"]
            snr_threshold                  --> [float]          What SNR threshold required for a detection
            Mc_max                         --> [float]          Maximum chirp mass in grid
            Mc_step                        --> [float]          Step in chirp mass to use in grid
            eta_max                        --> [float]          Maximum symmetric mass ratio in grid
            eta_step                       --> [float]          Step in symmetric mass ratio to use in grid
            snr_max                        --> [float]          Maximum snr in grid
            snr_step                       --> [float]          Step in snr to use in grid

        Returns:
            snr_grid_at_1Mpc               --> [2D float array] The snr of a binary with masses (Mc, eta) at a distance of 1 Mpc
            detection_probability_from_snr --> [list of floats] A list of detection probabilities for different SNRs
    """
    # get interpolator given sensitivity
    interpolator = selection_effects.SNRinterpolator(sensitivity)

    # create chirp mass and eta arrays
    Mc_array = np.arange(Mc_step, Mc_max + Mc_step, Mc_step)
    eta_array = np.arange(eta_step, eta_max + eta_step, eta_step)

    # convert to total, primary and secondary mass arrays
    Mt_array = Mc_array / eta_array[:, np.newaxis] ** 0.6
    M1_array = Mt_array * 0.5 * (1. + np.sqrt(1. - 4 * eta_array[:, np.newaxis]))
    M2_array = Mt_array - M1_array

    # interpolate to get snr values if binary was at 1Mpc
    snr_grid_at_1Mpc = interpolator(M1_array, M2_array)

    # precompute a grid of detection probabilities as a function of snr
    snr_array = np.arange(snr_step, snr_max + snr_step, snr_step)
    detection_probability_from_snr = selection_effects.detection_probability_from_snr(snr_array, snr_threshold)

    return snr_grid_at_1Mpc, detection_probability_from_snr


def find_detection_probability(Mc, eta, redshifts, distances, n_redshifts_detection, n_binaries, snr_grid_at_1Mpc,
                               detection_probability_from_snr,
                               Mc_step=0.1, eta_step=0.01, snr_step=0.1):
    """
        Compute the detection probability given a grid of SNRs and detection probabilities with masses

        Args:
            Mc                             --> [list of floats] Chirp mass of binaries in COMPAS
            eta                            --> [list of floats] Symmetric mass ratios of binaries in COMPAS
            redshifts                      --> [list of floats] List of redshifts
            distances                      --> [list of floats] List of distances corresponding to redshifts
            n_redshifts_detection          --> [int]            Index (in redshifts) to which we evaluate detection probability
            n_binaries                     --> [int]            Number of merging binaries in the COMPAS file
            snr_grid_at_1Mpc               --> [2D float array] The snr of a binary with masses (Mc, eta) at a distance of 1 Mpc
            detection_probability_from_snr --> [list of floats] A list of detection probabilities for different SNRs
            Mc_step                        --> [float]          Step in chirp mass to use in grid
            eta_step                       --> [float]          Step in symmetric mass ratio to use in grid
            snr_step                       --> [float]          Step in snr to use in grid
    """
    # by default, set detection probability to one
    detection_probability = np.ones(shape=(n_binaries, n_redshifts_detection))

    # for each binary in the COMPAS file
    for i in range(n_binaries):
        # shift frames for the chirp mass
        Mc_shifted = Mc[i] * (1 + redshifts[:n_redshifts_detection])

        # work out the closest index to the given values of eta and Mc
        eta_index = np.round(eta[i] / eta_step).astype(int) - 1
        Mc_index = np.round(Mc_shifted / Mc_step).astype(int) - 1

        # lookup values for the snr (but make sure you don't go over the top of the array)
        snrs = np.ones(n_redshifts_detection) * 0.00001
        Mc_below_max = Mc_index < snr_grid_at_1Mpc.shape[1]
        snrs[Mc_below_max] = snr_grid_at_1Mpc[eta_index, Mc_index[Mc_below_max]]

        # convert these snr values to the correct distances
        snrs = snrs / distances[:n_redshifts_detection]

        # lookup values for the detection probability (but make sure you don't go over the top of the array)
        detection_list_index = np.round(snrs / snr_step).astype(int) - 1
        snr_below_max = detection_list_index < len(detection_probability_from_snr)
        snr_below_min = detection_list_index < 0

        # remember we set probability = 1 by default? Because if we don't set it here, we have snr > max snr
        # which is 1000 by default, meaning very detectable
        detection_probability[i, snr_below_max] = detection_probability_from_snr[detection_list_index[snr_below_max]]
        # on the other hand, if SNR is too low, the detection probability is effectively zero
        detection_probability[i, snr_below_min] = 0

    return detection_probability




