import numpy as np
from scipy.integrate import quad
import h5py as h5       # for reading in data
import functools

@functools.lru_cache()
def __get_imf_normalisation_values(m1=0.01, m2=0.08, m3=0.5, m4=200.0, a12=0.3, a23=1.3, a34=2.3):
    b1 = 1 / (
            (m2 ** (1 - a12) - m1 ** (1 - a12)) / (1 - a12)
            + m2 ** (-(a12 - a23)) * (m3 ** (1 - a23) - m2 ** (1 - a23)) / (1 - a23)
            + m2 ** (-(a12 - a23)) * m3 ** (-(a23 - a34)) * (m4 ** (1 - a34) - m3 ** (1 - a34)) / (1 - a34)
    )
    b2 = b1 * m2 ** (-(a12 - a23))
    b3 = b2 * m3 ** (-(a23 - a34))
    return b1, b2, b3

@np.vectorize
def IMF(m, m1=0.01, m2=0.08, m3=0.5, m4=200.0, a12=0.3, a23=1.3, a34=2.3):
    """Calculate the fraction of stellar mass between m and m + dm for a three part broken power law.

    Default values follow Kroupa (2001)
            zeta(m) ~ m^(-a_ij)
    Parameters
    ----------
    m : `float` or `np.ndarray`
        Mass at which to evaluate
    mi : float, optional
        masses at which to transition the slope
    aij : float, optional
        slope of the IMF between mi and mj
    Returns
    -------
    imf_vals
        IMF evaluated at the given masses
    """
    # calculate normalisation constants that ensure the IMF is continuous
    b1, b2, b3 = __get_imf_normalisation_values(m1, m2, m3, m4, a12, a23, a34)

    # evaluate IMF either at a point or for a list of points
    if m1 <= m < m2:
        return b1 * m ** (-a12)
    elif m2 <= m < m3:
        return b2 * m ** (-a23)
    elif m3 <= m < m4:
        return b3 * m ** (-a34)
    else:
        return 0.0


def get_COMPAS_fraction(m1_low, m1_upp, m2_low, f_bin, mass_ratio_pdf_function=lambda q: 1,
                        m1=0.01, m2=0.08, m3=0.5, m4=200.0, a12=0.3, a23=1.3, a34=2.3):
    """Calculate the fraction of mass in a COMPAS population relative to the total Universal population. This
    can be used to normalise the rates of objects from COMPAS simulations.

    Parameters
    ----------
    m1_low : `float`
        Lower limit on the sampled primary mass
    m1_upp : `float`
        Upper limit on the sampled primary mass
    m2_low : `float`
        Lower limit on the sampled secondary mass
    f_bin : `float`
        Binary fraction
    mass_ratio_pdf_function : `function`, optional
        Function to calculate the mass ratio PDF, by default a uniform mass ratio distribution
    mi, aij : `float`
        Settings for the IMF choice, see `IMF` for details, by default follows Kroupa (2001)

    Returns
    -------
    fraction
        The fraction of mass in a COMPAS population relative to the total Universal population
    """ 
    # first, for normalisation purposes, we can find the integral with no COMPAS cuts
    def full_integral(mass, m1, m2, m3, m4, a12, a23, a34):
        primary_mass = IMF(mass, m1, m2, m3, m4, a12, a23, a34) * mass
        
        # find the expected companion mass given the mass ratio pdf function
        expected_secondary_mass = quad(lambda q: q * mass_ratio_pdf_function(q), 0, 1)[0] * primary_mass
        
        single_stars = (1 - f_bin) * primary_mass
        binary_stars = f_bin * (primary_mass + expected_secondary_mass)
        return single_stars + binary_stars
    full_mass = quad(full_integral, m1, m4, args=(m1, m2, m3, m4, a12, a23, a34))[0]
    
    # now we do a similar integral but for the COMPAS regime
    def compas_integral(mass, m2_low, f_bin, m1, m2, m3, m4, a12, a23, a34):
        # define the primary mass in the same way
        primary_mass = IMF(mass, m1, m2, m3, m4, a12, a23, a34) * mass
        
        # find the fraction that are below the m2 mass cut
        f_below_m2low = quad(mass_ratio_pdf_function, 0, m2_low / mass)[0]
        
        # expectation value of the secondary mass given the m2 cut and mass ratio pdf function
        expected_secondary_mass = quad(lambda q: q * mass_ratio_pdf_function(q), m2_low / mass, 1)[0] * primary_mass
        
        # return total mass of binary stars that have m2 above the cut
        return f_bin * (1 - f_below_m2low) * (primary_mass + expected_secondary_mass)
    
    compas_mass = quad(compas_integral, m1_low, m1_upp, args=(m2_low, f_bin, m1, m2, m3, m4, a12, a23, a34))[0]
    return compas_mass / full_mass


def retrieveMassEvolvedPerZ(path):
    with h5.File(path, 'r') as f:
        allSystems = f['BSE_System_Parameters']
        metals = (allSystems['Metallicity@ZAMS(1)'])[()]
        m1s = (allSystems['Mass@ZAMS(1)'])[()]
        m2s = (allSystems['Mass@ZAMS(2)'])[()]
        unique_metals = np.unique(metals)
        total = np.zeros(len(unique_metals))
        for i, Z in enumerate(unique_metals):
            mask = metals == Z
            total[i] = np.sum(m1s[mask]) + np.sum(m2s[mask])
    return total


def totalMassEvolvedPerZ(path, Mlower, Mupper, m2_low, binaryFraction, mass_ratio_pdf_function=lambda q: 1,
                         m1=0.01, m2=0.08, m3=0.5, m4=200., a12=-0.3, a23=-1.3, a34=-2.3):

    # calculate the fraction of mass in the COMPAS simulation vs. the real population without sample cuts
    fraction = get_COMPAS_fraction(m1_low=Mlower, m1_upp=Mupper, m2_low=m2_low, f_bin=binaryFraction,
                                   mass_ratio_pdf_function=mass_ratio_pdf_function,
                                   m1=m1, m2=m2, m3=m3, m4=m4, a12=a12, a23=a23, a34=a34)
    multiplicationFactor = 1 / fraction

    # get the mass evolved for each metallicity bin and convert to a total mass using the fraction
    MassEvolvedPerZ = retrieveMassEvolvedPerZ(path)
    totalMassEvolvedPerMetallicity = MassEvolvedPerZ / fraction

    return multiplicationFactor, totalMassEvolvedPerMetallicity


def analytical_star_forming_mass_per_binary_using_kroupa_imf(
        m1_max, m1_min, m2_min, fbin=1., imf_mass_bounds=[0.01,0.08,0.5,200]
):
    """
    Analytical computation of the mass of stars formed per binary star formed within the
    [m1 min, m1 max] and [m2 min, ..] rage,
    using the Kroupa IMF:

        p(M) \propto M^-0.3 for M between m1 and m2
        p(M) \propto M^-1.3 for M between m2 and m3;
        p(M) = alpha * M^-2.3 for M between m3 and m4;

    @Ilya Mandel's derivation
    """
    m1, m2, m3, m4 = imf_mass_bounds
    alpha = (-(m4**(-1.3)-m3**(-1.3))/1.3 - (m3**(-0.3)-m2**(-0.3))/(m3*0.3) + (m2**0.7-m1**0.7)/(m2*m3*0.7))**(-1)
    # average mass of stars (average mass of all binaries is a factor of 1.5 larger)
    m_avg = alpha * (-(m4**(-0.3)-m3**(-0.3))/0.3 + (m3**0.7-m2**0.7)/(m3*0.7) + (m2**1.7-m1**1.7)/(m2*m3*1.7))
    # fraction of binaries that COMPAS simulates
    fint = -alpha / 1.3 * (m1_max ** (-1.3) - m1_min ** (-1.3)) + alpha * m2_min / 2.3 * (m1_max ** (-2.3) - m1_min ** (-2.3))
    # mass represented by each binary simulated by COMPAS
    m_rep = (1/fint) * m_avg * (1.5 + (1-fbin)/fbin)
    return m_rep