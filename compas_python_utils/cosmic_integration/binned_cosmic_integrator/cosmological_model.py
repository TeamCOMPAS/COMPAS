from astropy import units as u
from astropy.cosmology import Planck15 as cosmology
import numpy as np
from scipy.stats import norm as NormDist

from ..FastCosmicIntegration import find_sfr, find_metallicity_distribution



MAX_REDSHIFT = 10.0
REDSHIFT_STEP = 0.001

class CosmologicalModel:
    def __init__(
            self, sf_params, mu0, sigma0,
            min_observed_log_metallicity=-4, max_n_observed_log_metallicity=0,
    ):
        self.aSF = sf_params[0]
        self.bSF = sf_params[1]
        self.cSF = sf_params[2]
        self.dSF = sf_params[3]
        self.mu0 = mu0
        self.sigma0 = sigma0


        self.redshifts = np.arange(0, MAX_REDSHIFT + REDSHIFT_STEP, REDSHIFT_STEP)
        self.sfr = find_sfr(self.redshifts, a=self.aSF, b=self.bSF, c=self.cSF, d=self.dSF)
        dPdlogZ, metallicities, p_draw_metallicity =  find_metallicity_distribution(

        )



    def plot(self):
        pass






def find_sfr(redshifts, a=0.01, b=2.77, c=2.90, d=4.70):
    """
        Calculate the star forming mass per unit volume per year following
        Neijssel+19 Eq. 6, using functional form of Madau & Dickinson 2014

        Args:
            redshifts --> [list of floats] List of redshifts at which to evaluate the sfr

        Returns:
            sfr       --> [list of floats] Star forming mass per unit volume per year for each redshift
    """
    # get value in mass per year per cubic Mpc and convert to per cubic Gpc then return
    sfr = a * ((1 + redshifts) ** b) / (1 + ((1 + redshifts) / c) ** d) * u.Msun / u.yr / u.Mpc ** 3
    return sfr.to(u.Msun / u.yr / u.Gpc ** 3).value


def find_metallicity_distribution(redshifts, min_logZ_COMPAS, max_logZ_COMPAS,
                                  mu0=0.035, muz=-0.23, sigma_0=0.39, sigma_z=0.0, alpha=0.0,
                                  min_logZ=-12.0, max_logZ=0.0, step_logZ=0.01):
    """
    Calculate the distribution of metallicities at different redshifts using a log skew normal distribution
    the log-normal distribution is a special case of this log skew normal  distribution, and is retrieved by setting
    the skewness to zero (alpha = 0).
    Based on the method in Neijssel+19. Default values of mu0=0.035, muz=-0.23, sigma_0=0.39, sigma_z=0.0, alpha =0.0,
    retrieve the dP/dZ distribution used in Neijssel+19.  See van Son+2022 for skewed log normal distribution.

    NOTE: This assumes that metallicities in COMPAS are drawn from a flat in log distribution!

    Args:
        max_redshift       --> [float]          max redshift for calculation
        redshift_step      --> [float]          step used in redshift calculation
        min_logZ_COMPAS    --> [float]          Minimum logZ value that COMPAS samples
        max_logZ_COMPAS    --> [float]          Maximum logZ value that COMPAS samples

        mu0    =  0.035    --> [float]          location (mean in normal) at redshift 0
        muz    = -0.23    --> [float]           redshift scaling/evolution of the location
        sigma_0 = 0.39     --> [float]          Scale (variance in normal) at redshift 0
        sigma_z = 0.00     --> [float]          redshift scaling of the scale (variance in normal)
        alpha   = 0.00    --> [float]           shape (skewness, alpha = 0 retrieves normal dist)

        min_logZ           --> [float]          Minimum logZ at which to calculate dPdlogZ (influences normalization)
        max_logZ           --> [float]          Maximum logZ at which to calculate dPdlogZ (influences normalization)
        step_logZ          --> [float]          Size of logZ steps to take in finding a Z range

    Returns:
        dPdlogZ            --> [2D float array] Probability of getting a particular logZ at a certain redshift
        metallicities      --> [list of floats] Metallicities at which dPdlogZ is evaluated
        p_draw_metallicity --> float            Probability of drawing a certain metallicity in COMPAS (float because assuming uniform)
    """
    ##################################
    # Log-Linear redshift dependence of sigma
    sigma = sigma_0 * 10 ** (sigma_z * redshifts)

    ##################################
    # Follow Langer & Norman 2006 in assuming that mean metallicities evolve in z as:
    mean_metallicities = mu0 * 10 ** (muz * redshifts)

    # Now we re-write the expected value of the log-skew-normal to retrieve mu
    beta = alpha / (np.sqrt(1 + (alpha) ** 2))
    PHI = NormDist.cdf(beta * sigma)
    mu_metallicities = np.log(mean_metallicities / 2. * 1. / (np.exp(0.5 * sigma ** 2) * PHI))

    ##################################
    # create a range of metallicities (the x-values, or random variables)
    log_metallicities = np.arange(min_logZ, max_logZ + step_logZ, step_logZ)
    metallicities = np.exp(log_metallicities)

    ##################################
    # probabilities of log-skew-normal (without the factor of 1/Z since this is dp/dlogZ not dp/dZ)
    dPdlogZ = 2. / (sigma[:, np.newaxis]) * NormDist.pdf(
        (log_metallicities - mu_metallicities[:, np.newaxis]) / sigma[:, np.newaxis]) * NormDist.cdf(
        alpha * (log_metallicities - mu_metallicities[:, np.newaxis]) / sigma[:, np.newaxis])

    ##################################
    # normalise the distribution over all metallicities; this choice of normalisation assumes that metallicities outside the COMPAS range have yields of zero
    norm = dPdlogZ.sum(axis=-1) * step_logZ
    dPdlogZ = dPdlogZ / norm[:, np.newaxis]

    ##################################
    # assume a flat in log distribution in sampled metallicity to find probability of drawing Z in COMPAS
    p_draw_metallicity = 1 / (max_logZ_COMPAS - min_logZ_COMPAS)

    return dPdlogZ, metallicities, p_draw_metallicity