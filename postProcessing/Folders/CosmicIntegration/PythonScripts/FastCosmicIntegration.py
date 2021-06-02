import numpy as np
import h5py  as h5
import os
import time
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmology
from scipy.interpolate import interp1d
import scipy
import ClassCOMPAS
import selection_effects
import warnings
import astropy.units as u
import argparse

import astropy.units as u

def calculate_redshift_related_params(max_redshift=10.0, max_redshift_detection=1.0, redshift_step=0.001, z_first_SF = 10.0):
    """ 
        Given limits on the redshift, create an array of redshifts, times, distances and volumes

        Args:
            max_redshift           --> [float]          Maximum redshift to use for calculations
            max_redshift_detection --> [float]          Maximum redshift to calculate detection rates (must be <= max_redshift)
            redshift_step          --> [float]          size of step to take in redshift
            z_first_SF             --> [float]          redshift of first star formation

        Returns:
            redshifts              --> [list of floats] List of redshifts between limits supplied
            n_redshifts_detection  --> [int]            Number of redshifts in list that should be used to calculate detection rates
            times                  --> [list of floats] Equivalent of redshifts but converted to age of Universe
            distances              --> [list of floats] Equivalent of redshifts but converted to luminosity distances
            shell_volumes          --> [list of floats] Equivalent of redshifts but converted to shell volumes
    """
    # create a list of redshifts and record lengths
    redshifts = np.arange(0, max_redshift + redshift_step, redshift_step)
    n_redshifts_detection = int(max_redshift_detection / redshift_step)

    # convert redshifts to times and ensure all times are in Myr
    times = cosmology.age(redshifts).to(u.Myr).value

    # and time of first Sf
    time_first_SF = cosmology.age(z_first_SF).to(u.Myr).value

    # convert redshifts to distances and ensure all distances are in Mpc (also avoid D=0 because division by 0)
    distances = cosmology.luminosity_distance(redshifts).to(u.Mpc).value
    distances[0] = 0.001

    # convert redshifts to volumnes and ensure all volumes are in Gpc^3
    volumes = cosmology.comoving_volume(redshifts).to(u.Gpc**3).value

    # split volumes into shells and duplicate last shell to keep same length
    shell_volumes = np.diff(volumes)
    shell_volumes = np.append(shell_volumes, shell_volumes[-1])

    return redshifts, n_redshifts_detection, times, time_first_SF, distances, shell_volumes

def find_sfr(redshifts):
    """
        Calculate the star forming mass per unit volume per year using
        Neijssel+19 Eq. 6

        Args:
            redshifts --> [list of floats] List of redshifts at which to evaluate the sfr

        Returns:
            sfr       --> [list of floats] Star forming mass per unit volume per year for each redshift
    """
    # use Neijssel+19 to get value in mass per year per cubic Mpc and convert to per cubic Gpc then return
    sfr = 0.01 * ((1+redshifts)**2.77) / (1 + ((1+redshifts)/2.9)**4.7) * u.Msun / u.yr / u.Mpc**3
    return sfr.to(u.Msun / u.yr / u.Gpc**3).value


def find_metallicity_distribution(redshifts, min_logZ_COMPAS, max_logZ_COMPAS,
                                  mu0=0.035, muz=-0.23, sigma_0=0.39, sigma_z=0.0, alpha =0.0,
                                  min_logZ  =-12.0, max_logZ  =0.0, step_logZ = 0.01):
                                 
    """
    Calculate the distribution of metallicities at different redshifts using a log skew normal distribution
    the log-normal distribution is a special case of this log skew normal distribution distribution, and is retrieved by setting 
    the skewness to zero (alpha = 0). 
    Based on the method in Neijssel+19. Default values of mu0=0.035, muz=-0.23, sigma_0=0.39, sigma_z=0.0, alpha =0.0, 
    retrieve the dP/dZ distribution used in Neijssel+19

    NOTE: This assumes that metallicities in COMPAS are drawn from a flat in log distribution!

    Args:
        max_redshift       --> [float]          max redshift for calculation
        redshift_step      --> [float]          step used in redshift calculation
        min_logZ_COMPAS    --> [float]          Minimum logZ value that COMPAS samples
        max_logZ_COMPAS    --> [float]          Maximum logZ value that COMPAS samples
        
        mu0    = 0.025    --> [float]           location (mean in normal) at redshift 0
        muz    = -0.05    --> [float]           redshift scaling/evolution of the location
        sigma_0 = 1.25     --> [float]          Scale (variance in normal) at redshift 0
        sigma_z = 0.05     --> [float]          redshift scaling of the scale (variance in normal)
        alpha   = -1.77    --> [float]          shape (skewness, alpha = 0 retrieves normal dist)

        min_logZ           --> [float]          Minimum logZ at which to calculate dPdlogZ (influences normalization)
        max_logZ           --> [float]          Maximum logZ at which to calculate dPdlogZ (influences normalization)
        step_logZ          --> [float]          Size of logZ steps to take in finding a Z range

    Returns:
        dPdlogZ            --> [2D float array] Probability of getting a particular logZ at a certain redshift
        metallicities      --> [list of floats] Metallicities at which dPdlogZ is evaluated
        p_draw_metallicity --> float            Probability of drawing a certain metallicity in COMPAS (float because assuming uniform)
    """
    import scipy

    ##################################
    # the PDF of a standard normal distrtibution
    def normal_PDF(x):
        return 1./(np.sqrt(2* np.pi)) * np.exp(-(1./2) * (x)**2 )

    ##################################
    # the CDF of a standard normal distrtibution
    def normal_CDF(x):
        return 1./2. * (1 + scipy.special.erf(x/np.sqrt(2)) )
        
    ##################################
    # Log-Linear redshift dependence of sigma
    sigma = sigma_0* 10**(sigma_z*redshifts)
    
    ##################################
    # Follow Langer & Norman 2007? in assuming that mean metallicities evolve in z as:
    mean_metallicities = mu0 * 10**(muz * redshifts) 
        
    # Now we re-write the expected value of ou log-skew-normal to retrieve mu
    beta = alpha/(np.sqrt(1 + (alpha)**2))
    PHI  = normal_CDF(beta * sigma) 
    mu_metallicities = np.log(mean_metallicities/2. * 1./(np.exp(0.5*sigma**2) * PHI )  ) 

    ##################################
    # create a range of metallicities (thex-values, or random variables)
    log_metallicities = np.arange(min_logZ, max_logZ + step_logZ, step_logZ)
    metallicities = np.exp(log_metallicities)


    ##################################
    # probabilities of log-skew-normal (without the factor of 1/Z since this is dp/dlogZ not dp/dZ)
    dPdlogZ = 2./(sigma[:,np.newaxis]) * normal_PDF((log_metallicities -  mu_metallicities[:,np.newaxis])/sigma[:,np.newaxis]) * normal_CDF(alpha * (log_metallicities -  mu_metallicities[:,np.newaxis])/sigma[:,np.newaxis] )

    ##################################
    # normalise the distribution over al metallicities
    norm = dPdlogZ.sum(axis=-1) * step_logZ
    dPdlogZ = dPdlogZ /norm[:,np.newaxis]

    ##################################
    # assume a flat in log distribution in metallicity to find probability of drawing Z in COMPAS
    p_draw_metallicity = 1 / (max_logZ_COMPAS - min_logZ_COMPAS)
    
    return dPdlogZ, metallicities, p_draw_metallicity




def find_formation_and_merger_rates(n_binaries, redshifts, times, time_first_SF, n_formed, dPdlogZ, metallicities, p_draw_metallicity,
                                    COMPAS_metallicites, COMPAS_delay_times, COMPAS_weights=None):
    """
        Find both the formation and merger rates for each binary at each redshift

        Args:
            n_binaries          --> [int]            Number of DCO binaries in the arrays
            redshifts           --> [list of floats] Redshifts at which to evaluate the rates
            times               --> [list of floats] Equivalent of the redshifts in terms of age of the Universe
            n_formed            --> [float]          Binary formation rate (number of binaries formed per year per cubic Gpc) represented by each simulated COMPAS binary
            dPdlogZ             --> [2D float array] Probability of getting a particular logZ at a certain redshift
            metallicities       --> [list of floats] Metallicities at which dPdlogZ is evaluated
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
        # calculate formation rate (see Neijssel+19 Section 4) - note this uses dPdlogZ for *closest* metallicity
        formation_rate[i, :] = n_formed * dPdlogZ[:, np.digitize(COMPAS_metallicites[i], metallicities)] / p_draw_metallicity * COMPAS_weights[i]

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
    Mt_array = Mc_array / eta_array[:,np.newaxis]**0.6
    M1_array = Mt_array * 0.5 * (1. + np.sqrt(1. - 4 * eta_array[:,np.newaxis]))
    M2_array = Mt_array - M1_array

    # interpolate to get snr values if binary was at 1Mpc
    snr_grid_at_1Mpc = interpolator(M1_array, M2_array)

    # precompute a grid of detection probabilities as a function of snr
    snr_array = np.arange(snr_step, snr_max + snr_step, snr_step)
    detection_probability_from_snr = selection_effects.detection_probability_from_snr(snr_array, snr_threshold)

    return snr_grid_at_1Mpc, detection_probability_from_snr

def find_detection_probability(Mc, eta, redshifts, distances, n_redshifts_detection, n_binaries, snr_grid_at_1Mpc, detection_probability_from_snr,
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
        #on the other hand, if SNR is too low, the detection probability is effectively zero
        detection_probability[i, snr_below_min] = 0

    return detection_probability

def find_detection_rate(path, filename="COMPAS_Output.h5", dco_type="BBH", weight_column=None,
                        merges_hubble_time=True, pessimistic_CEE=True, no_RLOF_after_CEE=True,
                        max_redshift=10.0, max_redshift_detection=1.0, redshift_step=0.001,
                        z_first_SF = 10, m1_min=5 * u.Msun, m1_max=150 * u.Msun, m2_min=0.1 * u.Msun, fbin=0.7,
                        Z0=0.035, alpha=-0.23, sigma=0.39, step_logZ=0.01,
                        sensitivity="O1", snr_threshold=8, 
                        Mc_max=300.0, Mc_step=0.1, eta_max=0.25, eta_step=0.01,
                        snr_max=1000.0, snr_step=0.1):
    """
        The main function of this file. Finds the detection rate, formation rate and merger rate for each
        binary in a COMPAS file at a series of redshifts defined by intput. Also returns relevant COMPAS
        data.

        NOTE: This code assumes that assumes that metallicities in COMPAS are drawn from a flat in log distribution

        Args:
            ===================================================
            == Arguments for finding and masking COMPAS file ==
            ===================================================
            path                   --> [string] Path to the COMPAS file that contains the output
            filename               --> [string] Name of the COMPAS file
            dco_type               --> [string] Which DCO type to calculate rates for: one of ["all", "BBH", "BHNS", "BNS"]
            weight_column          --> [string] Name of column in "DoubleCompactObjects" file that contains adaptive sampling weights
                                                    (Leave this as None if you have unweighted samples)
            merges_in_hubble_time  --> [bool]   whether to mask binaries that don't merge in a Hubble time
            no_RLOF_after_CEE      --> [bool]   whether to mask binaries that have immediate RLOF after a CCE
            pessimistic_CEE        --> [bool]   whether to mask binaries that go through Optimistic CE scenario

            ===========================================
            == Arguments for creating redshift array ==
            ===========================================
            max_redshift           --> [float]  Maximum redshift to use in array
            max_redshift_detection --> [float]  Maximum redshift to calculate detection rates (must be <= max_redshift)
            redshift_step          --> [float]  Size of step to take in redshift

            ====================================================================
            == Arguments for determining star forming mass per sampled binary ==
            ====================================================================
            m1_min                 --> [float]  Minimum primary mass sampled by COMPAS
            m1_max                 --> [float]  Maximum primary mass sampled by COMPAS
            m2_min                 --> [float]  Minimum secondary mass sampled by COMPAS
            fbin                   --> [float]  Binary fraction used by COMPAS

            =======================================================================
            == Arguments for creating metallicity distribution and probabilities ==
            =======================================================================
            Z0                     --> [float]  Parameter used for calculating mean metallicity (see Neijssel+19 Eq.7-9)
            alpha                  --> [float]  Parameter used for calculating mean metallicity (see Neijssel+19 Eq.7-9)
            sigma                  --> [float]  Parameter used for calculating mean metallicity (see Neijssel+19 Eq.7-9)
            step_logZ              --> [float]  Size of logZ steps to take in finding a Z range

            =======================================================
            == Arguments for determining detection probabilities ==
            =======================================================
            sensitivity            --> [string] Which detector sensitivity to use: one of ["design", "O1", "O3"]
            snr_threshold          --> [float]  What SNR threshold required for a detection
            Mc_max                 --> [float]  Maximum chirp mass in grid
            Mc_step                --> [float]  Step in chirp mass to use in grid
            eta_max                --> [float]  Maximum symmetric mass ratio in grid
            eta_step               --> [float]  Step in symmetric mass ratio to use in grid
            snr_max                --> [float]  Maximum snr in grid
            snr_step               --> [float]  Step in snr to use in grid

        Returns:
            detection_rate         --> [2D float array] Detection rate for each binary at each redshift in 1/yr
            formation_rate         --> [2D float array] Formation rate for each binary at each redshift in 1/yr/Gpc^3
            merger_rate            --> [2D float array] Merger rate for each binary at each redshift in 1/yr/Gpc^3
            redshifts              --> [list of floats] List of redshifts
            COMPAS                 --> [Object]         Relevant COMPAS data in COMPASData Class
    """

    # assert that input will not produce errors
    assert max_redshift_detection <= max_redshift, "Maximum detection redshift cannot be below maximum redshift"
    assert m1_min <= m1_max, "Minimum sampled primary mass cannot be above maximum sampled primary mass"
    assert np.logical_and(fbin >= 0.0, fbin <= 1.0), "Binary fraction must be between 0 and 1"
    assert Mc_step < Mc_max, "Chirp mass step size must be less than maximum chirp mass"
    assert eta_step < eta_max, "Symmetric mass ratio step size must be less than maximum symmetric mass ratio"
    assert snr_step < snr_max, "SNR step size must be less than maximum SNR"

    nonnegative_args = [(max_redshift, "max_redshift"), (max_redshift_detection, "max_redshift_detection"), (m1_min.value, "m1_min"), (m1_max.value, "m1_max"),
                        (m2_min.value, "m2_min"), (Z0, "Z0"), (sigma, "sigma"), (step_logZ, "step_logZ"), (snr_threshold, "snr_threshold"), (Mc_max, "Mc_max"),
                        (Mc_step, "Mc_step"), (eta_max, "eta_max"), (eta_step, "eta_step"), (snr_max, "snr_max"), (snr_step, "snr_step")]
    for arg, arg_str in nonnegative_args:
        assert arg >= 0.0, "{} must be nonnegative".format(arg_str)

    # warn if input is not advisable
    if redshift_step > max_redshift_detection:
        warnings.warn("Redshift step is greater than maximum detection redshift", stacklevel=2)
    if Mc_step > 1.0:
        warnings.warn("Chirp mass step is greater than 1.0, large step sizes can produce unpredictable results", stacklevel=2)
    if eta_step > 0.1:
        warnings.warn("Symmetric mass ratio step is greater than 0.1, large step sizes can produce unpredictable results", stacklevel=2)
    if snr_step > 1.0:
        warnings.warn("SNR step is greater than 1.0, large step sizes can produce unpredictable results", stacklevel=2)

    # start by getting the necessary data from the COMPAS file
    COMPAS = ClassCOMPAS.COMPASData(path, fileName=filename, Mlower=m1_min, Mupper=m1_max, m2_min=m2_min, binaryFraction=fbin, suppress_reminder=True)
    COMPAS.setCOMPASDCOmask(types=dco_type, withinHubbleTime=merges_hubble_time, pessimistic=pessimistic_CEE, noRLOFafterCEE=no_RLOF_after_CEE)
    COMPAS.setCOMPASData()
    COMPAS.set_sw_weights(weight_column)
    COMPAS.find_star_forming_mass_per_binary_sampling()

    # compute the chirp masses and symmetric mass ratios only for systems of interest
    chirp_masses = (COMPAS.mass1*COMPAS.mass2)**(3/5) / (COMPAS.mass1 + COMPAS.mass2)**(1/5)
    etas = COMPAS.mass1 * COMPAS.mass2 / (COMPAS.mass1 + COMPAS.mass2)**2
    n_binaries = len(chirp_masses)
    # another warning on poor input
    if max(chirp_masses)*(1+max_redshift_detection) < Mc_max:
        warnings.warn("Maximum chirp mass used for detectability calculation is below maximum binary chirp mass * (1+maximum redshift for detectability calculation)", stacklevel=2)

    # calculate the redshifts array and its equivalents
    redshifts, n_redshifts_detection, times, time_first_SF, distances, shell_volumes = calculate_redshift_related_params(max_redshift, max_redshift_detection, redshift_step, z_first_SF)

    # find the star forming mass per year per Gpc^3 and convert to total number formed per year per Gpc^3
    sfr = find_sfr(redshifts)
    n_formed = sfr / (COMPAS.mass_evolved_per_binary.value * COMPAS.n_systems)

    # work out the metallicity distribution at each redshift and probability of drawing each metallicity in COMPAS
    dPdlogZ, metallicities, p_draw_metallicity = find_metallicity_distribution(redshifts, np.log(np.min(COMPAS.initialZ)),
                                                                                np.log(np.max(COMPAS.initialZ)), Z0, alpha, sigma, step_logZ)

    # calculate the formation and merger rates using what we computed above
    formation_rate, merger_rate = find_formation_and_merger_rates(n_binaries, redshifts, times, time_first_SF, n_formed, dPdlogZ,
                                                                    metallicities, p_draw_metallicity, COMPAS.metallicitySystems,
                                                                    COMPAS.delayTimes, COMPAS.sw_weights)

    # create lookup tables for the SNR at 1Mpc as a function of the masses and the probability of detection as a function of SNR
    snr_grid_at_1Mpc, detection_probability_from_snr = compute_snr_and_detection_grids(sensitivity, snr_threshold, Mc_max, Mc_step,
                                                                                    eta_max, eta_step, snr_max, snr_step)

    # use lookup tables to find the probability of detecting each binary at each redshift
    detection_probability = find_detection_probability(chirp_masses, etas, redshifts, distances, n_redshifts_detection, n_binaries,
                                                        snr_grid_at_1Mpc, detection_probability_from_snr, Mc_step, eta_step, snr_step)

    # finally, compute the detection rate using Neijssel+19 Eq. 2
    detection_rate = np.zeros(shape=(n_binaries, n_redshifts_detection))
    detection_rate = merger_rate[:, :n_redshifts_detection] * detection_probability \
                    * shell_volumes[:n_redshifts_detection] / (1 + redshifts[:n_redshifts_detection])

    return detection_rate, formation_rate, merger_rate, redshifts, COMPAS

def plot_rates(formation_rate, merger_rate, detection_rate, redshifts, chirp_masses):
    """
        Show a summary plot of the results, it also returns the summaries that it computes

        Args:
            detection_rate            --> [2D float array] Detection rate for each binary at each redshift in 1/yr
            formation_rate            --> [2D float array] Formation rate for each binary at each redshift in 1/yr/Gpc^3
            merger_rate               --> [2D float array] Merger rate for each binary at each redshift in 1/yr/Gpc^3
            redshifts                 --> [list of floats] List of redshifts
            COMPAS                    --> [Object]         Relevant COMPAS data in COMPASData Class

        Returns:
            total_formation_rate      --> [list of floats] Total formation rate over each redshift
            total_merger_rate         --> [list of floats] Total merger rate over each redshift
            total_detection_rate      --> [list of floats] Total detection rate over each redshift
            cumulative_detection_rate --> [list of floats] Cumulative detection rate over each redshift
            detection_rate_by_binary  --> [list of floats] Detection rate for each binary over every redshift
    """
    # sum things up across binaries
    total_formation_rate = np.sum(formation_rate, axis=0)
    total_merger_rate = np.sum(merger_rate, axis=0)
    total_detection_rate = np.sum(detection_rate, axis=0)
    
    # and across redshifts
    cumulative_detection_rate = np.cumsum(total_detection_rate)
    detection_rate_by_binary = np.sum(detection_rate, axis=1)

    # set some constants for the plots
    plt.rc('font', family='serif')
    fs = 20
    lw = 2

    fig, axes = plt.subplots(2, 2, figsize=(20, 20))

    axes[0,0].plot(redshifts, total_formation_rate, lw=lw)
    axes[0,0].set_xlabel('Redshift', fontsize=fs)
    axes[0,0].set_ylabel(r'Formation rate $[\rm \frac{\mathrm{d}N}{\mathrm{d}Gpc^3 \mathrm{d}yr}]$', fontsize=fs)

    axes[0,1].plot(redshifts, total_merger_rate, lw=lw)
    axes[0,1].set_xlabel('Redshift', fontsize=fs)
    axes[0,1].set_ylabel(r'Merger rate $[\rm \frac{\mathrm{d}N}{\mathrm{d}Gpc^3 \mathrm{d}yr}]$', fontsize=fs)

    axes[1,0].plot(redshifts[:len(cumulative_detection_rate)], cumulative_detection_rate, lw=lw)
    axes[1,0].set_xlabel('Redshift', fontsize=fs)
    axes[1,0].set_ylabel(r'Cumulative detection rate $[\rm \frac{\mathrm{d}N}{\mathrm{d}yr}]$', fontsize=fs)

    axes[1,1].hist(chirp_masses, weights=detection_rate_by_binary, bins=25, range=(0, 50))
    axes[1,1].set_xlabel(r'Chirp mass, $\mathcal{M}_c$', fontsize=fs)
    axes[1,1].set_ylabel(r'Mass distrbution of detections $[\rm \frac{\mathrm{d}N}{\mathrm{d}\mathcal{M}_c \mathrm{d}yr}]$', fontsize=fs)

    for ax in axes.flatten():
        ax.tick_params(labelsize=0.7*fs)

    plt.show()


##################################################################
### 
### Run it!
###
##################################################################
if __name__ == "__main__":
    # Define command line options for the most commonly varied options
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", dest= 'path',  help="Path to the COMPAS file that contains the output",type=str, default = './')
    parser.add_argument("--filename", dest= 'fname',  help="Name of the COMPAS file",type=str, default = "COMPAS_Output.h5")

    # Parameters determining dP/dZ and SFR(z), default options from Neijssel 2019
    parser.add_argument("--mu0", dest= 'mu0',  help="mean metallicity at redshhift 0",type=float, default=0.035)
    parser.add_argument("--muz", dest= 'muz',  help="redshift evolution of mean metallicity, dPdlogZ",type=float, default=-0.23)
    parser.add_argument("--sigma0", dest= 'sigma0',  help="variance in metallicity density distribution, dPdlogZ",type=float, default=0.39)
    parser.add_argument("--sigmaz", dest= 'sigmaz',  help="redshift evolution of variance, dPdlogZ",type=float, default=0.0)
    parser.add_argument("--alpha", dest= 'alpha',  help="skewness of mtallicity density distribution, dPdlogZ",type=float, default=0.0)
    parser.add_argument("--aSF", dest= 'aSF',  help="Parameter for shape of SFR(z)",type=float, default=0.01) 
    parser.add_argument("--bSF", dest= 'bSF',  help="Parameter for shape of SFR(z)",type=float, default=2.77)
    parser.add_argument("--cSF", dest= 'cSF',  help="Parameter for shape of SFR(z)",type=float, default=2.90)
    parser.add_argument("--dSF", dest= 'dSF',  help="Parameter for shape of SFR(z)",type=float, default=4.70)
 
    args = parser.parse_args()


    start_CI = time.time()
    # Do the cosmic integration
    detection_rate, formation_rate, merger_rate, redshifts, COMPAS = find_detection_rate(args.path, filename=args.fname, dco_type=Dco_type, weight_column="mixture_weight",
                            max_redshift=Max_redshift, max_redshift_detection=Max_redshift_detection, redshift_step=Redshift_step, z_first_SF= z_first_SF,
                            m1_min=M1_min, m1_max=M1_max, m2_min=M2_min, fbin=Fbin,
                            aSF = args.aSF, bSF = args.bSF, cSF = args.cSF, dSF = args.dSF, 
                            mu0=args.mu0, muz=args.muz, sigma0=args.sigma0, sigmaz=args.sigmaz, alpha=args.alpha, 
                            min_logZ=-12.0, max_logZ=0.0, step_logZ=0.01,
                            sensitivity=Sensitivity, snr_threshold=Snr_threshold, 
                            Mc_min=0.1, Mc_max=300.0, Mc_step=0.1, eta_min=0.01, eta_max=0.25, eta_step=0.01,
                            snr_min=0.1, snr_max=1000.0, snr_step=0.1)
    end_CI = time.time()

    # If you want, you can also plot your result
    start_plot = time.time()
    chirp_masses = COMPAS.chirp_masses()[COMPAS.DCO_masks[Dco_type]].value
    plot_rates(args.path, formation_rate, merger_rate, detection_rate, redshifts, chirp_masses, show_plot = False, mu0=args.mu0, muz=args.muz, sigma0=args.sigma0, sigmaz=args.sigmaz, alpha=args.alpha)
    end_plot = time.time()

    print('CI took ', end_CI - start_CI, 's')
    print('plot took ', end_plot - start_plot, 's')

