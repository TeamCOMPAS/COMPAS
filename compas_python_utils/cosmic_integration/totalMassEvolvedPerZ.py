import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import h5py as h5
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
    https://arxiv.org/abs/astro-ph/0009005
    Equation 1-2

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



def get_COMPAS_fraction(m1_low, m1_upp, m2_low, f_bin=None,
                        mass_ratio_pdf_function=lambda q: 1,
                        m1=0.01, m2=0.08, m3=0.5, m4=200.0,
                        a12=0.3, a23=1.3, a34=2.3):
    """
    Calculate the fraction of mass in a COMPAS population relative to the total Universal population.
    Can be used to normalise the rates of objects from COMPAS simulations.

    Parameters
    ----------
    m1_low, m1_upp : float
        Primary mass cuts in COMPAS simulation
    m2_low : float
        Secondary mass cutoff
    f_bin : float or None
        Binary fraction. If None, use a stepwise mass-dependent binary fraction.
    mass_ratio_pdf_function : function
        PDF of mass ratio q
    mi, aij : float
        IMF breakpoints and slopes
    """
    binary_bin_edges = [m1, 0.08, 0.5, 1, 10, m4]
    
    def get_binary_fraction(mass):
        binaryFractions = [0.1, 0.225, 0.5, 0.8, 1.0]
        for i in range(len(binary_bin_edges) - 1):
            if binary_bin_edges[i] <= mass < binary_bin_edges[i + 1]:
                return binaryFractions[i]
        return 0  # catch-all

    def IMF_mass(mass):
        return IMF(mass, m1, m2, m3, m4, a12, a23, a34) * mass

    def integrand_full(mass, f_bin):
        local_f_bin = get_binary_fraction(mass) if f_bin is None else f_bin
        expected_q = quad(lambda q: q * mass_ratio_pdf_function(q), 0, 1)[0]
        return (1 - local_f_bin + local_f_bin * (1 + expected_q)) * IMF_mass(mass)

    def integrand_compas(mass, f_bin):
        local_f_bin = get_binary_fraction(mass) if f_bin is None else f_bin
        q_min = m2_low / mass
        if q_min >= 1:
            return 0  # No valid secondaries
        f_q = quad(mass_ratio_pdf_function, q_min, 1)[0]
        expected_q = quad(lambda q: q * mass_ratio_pdf_function(q), q_min, 1)[0]
        return local_f_bin * f_q * (1 + expected_q) * IMF_mass(mass)

    # split integral at binary fraction steps if f_bin is None (i.e. variable and like a step function)
    def split_integral(func, a, b, f_bin):
        total = 0
        for edge_start, edge_end in zip(binary_bin_edges[:-1], binary_bin_edges[1:]):
            left = max(a, edge_start)
            right = min(b, edge_end)
            if left < right:
                result, _ = quad(func, left, right, args=(f_bin,))
                total += result
        return total

    if f_bin is None:
        full_mass = split_integral(integrand_full, m1, m4, f_bin)
        compas_mass = split_integral(integrand_compas, m1_low, m1_upp, f_bin)
    else:
        full_mass = quad(integrand_full, m1, m4, args=(f_bin,))[0]
        compas_mass = quad(integrand_compas, m1_low, m1_upp, args=(f_bin,))[0]

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
                         m1=0.01, m2=0.08, m3=0.5, m4=200., a12=0.3, a23=1.3, a34=2.3):
    """
    Calculate the total mass evolved per metallicity as a function of redshift in a COMPAS simulation.
    """
    # calculate the fraction of mass in the COMPAS simulation vs. the real population without sample cuts
    fraction = get_COMPAS_fraction(m1_low=Mlower, m1_upp=Mupper, m2_low=m2_low, f_bin=binaryFraction,
                                   mass_ratio_pdf_function=mass_ratio_pdf_function,
                                   m1=m1, m2=m2, m3=m3, m4=m4, a12=a12, a23=a23, a34=a34)
    multiplicationFactor = 1 / fraction

    # Warning: This is slow and error prone! esp if you sample metallicities smoothly
    # get the mass evolved for each metallicity bin and convert to a total mass using the fraction
    MassEvolvedPerZ = retrieveMassEvolvedPerZ(path)

    totalMassEvolvedPerMetallicity = MassEvolvedPerZ / fraction 

    return multiplicationFactor, totalMassEvolvedPerMetallicity


def star_forming_mass_per_binary(
        path,
        Mlower, Mupper, m2_low, binaryFraction, mass_ratio_pdf_function=lambda q: 1,
        m1=0.01, m2=0.08, m3=0.5, m4=200., a12=0.3, a23=1.3, a34=2.3):
    """
    Calculate the total mass of stars formed per binary star formed within the COMPAS simulation.
    """
    fraction = get_COMPAS_fraction(m1_low=Mlower,m1_upp=Mupper,m2_low=m2_low,
                                   f_bin=binaryFraction,mass_ratio_pdf_function=mass_ratio_pdf_function,
                                   m1=m1, m2=m2, m3=m3, m4=m4,
                                   a12=a12, a23=a23, a34=a34)

    # get the total mass in COMPAS and number of binaries
    with h5.File(path, 'r') as f:
        allSystems = f['BSE_System_Parameters']
        m1s = (allSystems['Mass@ZAMS(1)'])[()]
        m2s = (allSystems['Mass@ZAMS(2)'])[()]
        n_binaries = len(m1s)
        total_star_forming_mass_in_COMPAS = sum(m1s) + sum(m2s)

    total_star_forming_mass = total_star_forming_mass_in_COMPAS / fraction
    return total_star_forming_mass / n_binaries


def inverse_sample_IMF(
        n_samples = int(1e5),
        m_min=0.01, m_max=200,
        m1=0.01, m2=0.08, m3=0.5, m4=200., a12=0.3, a23=1.3, a34=2.3,
        cdf_pts=int(1e4)
        ):
    m = np.linspace(m_min, m_max, cdf_pts)
    imf_values = IMF(m, m1, m2, m3, m4, a12, a23, a34)
    cumulative = np.cumsum(imf_values)
    cumulative -= cumulative.min()
    f = interp1d(cumulative/cumulative.max(), m)
    return f(np.random.random(n_samples))

def draw_samples_from_kroupa_imf(
        Mlower, Mupper, m2_low,
        m1=0.01, m2=0.08, m3=0.5, m4=200., a12=0.3, a23=1.3, a34=2.3,
        n_samples = int(1e5)
):
    """
    Draw samples from the Kroupa IMF
    """
    m1_samples = inverse_sample_IMF(n_samples=n_samples,
        m_min=Mlower, m_max=Mupper,
        m1=m1, m2=m2, m3=m3, m4=m4, a12=a12, a23=a23, a34=a34
    )
    m2_samples = m1_samples * np.random.random(n_samples)
    mask = (Mlower < m1_samples) & (m1_samples <= Mupper) & (m2_low < m2_samples)
    return m1_samples[mask] , m2_samples[mask]





##################################################
# Old version of analytical calculation
##################################################
def old_analytical_star_forming_mass_per_binary_using_kroupa_imf(
        m1_min, m1_max, m2_min, fbin=1., imf_mass_bounds=[0.01,0.08,0.5,200]
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
    if m1_min < m3:
        raise ValueError(f"This analytical derivation requires IMF break m3  < m1_min ({m3} !< {m1_min})")
    alpha = (-(m4**(-1.3)-m3**(-1.3))/1.3 - (m3**(-0.3)-m2**(-0.3))/(m3*0.3) + (m2**0.7-m1**0.7)/(m2*m3*0.7))**(-1)
    # average mass of stars (average mass of all binaries is a factor of 1.5 larger)
    m_avg = alpha * (-(m4**(-0.3)-m3**(-0.3))/0.3 + (m3**0.7-m2**0.7)/(m3*0.7) + (m2**1.7-m1**1.7)/(m2*m3*1.7))
    # fraction of binaries that COMPAS simulates
    fint = -alpha / 1.3 * (m1_max ** (-1.3) - m1_min ** (-1.3)) + alpha * m2_min / 2.3 * (m1_max ** (-2.3) - m1_min ** (-2.3))
    # mass represented by each binary simulated by COMPAS
    m_rep = (1/fint) * m_avg * (1.5 + (1-fbin)/fbin)
    return m_rep


###################################################
def analytical_star_forming_mass_per_binary_using_kroupa_imf(
        m1_min, m1_max, m2_min, fbin=1., imf_mass_bounds=[0.01,0.08,0.5,200]
):
    """
    Analytical computation of the mass of stars formed per binary star formed within the
    [m1 min, m1 max] and [m2 min, ..] rage, using the Kroupa IMF:

        p(M) \propto M^-0.3 for M between m1 and m2
        p(M) \propto M^-1.3 for M between m2 and m3;
        p(M) = alpha * M^-2.3 for M between m3 and m4;

    m1_min, m1_max are the min and max sampled primary masses
    m2_min is the min sampled secondary mass

    This function further assumes a flat mass ratio distribution with qmin = m2_min/m1, and  m2_max = m1_max
    Lieke base on Ilya Mandel's derivation
    """
    # Kroupa IMF 
    m1, m2, m3, m4 = imf_mass_bounds
    continuity_constants = [1./(m2*m3), 1./(m3), 1.0]  
    IMF_powers = [-0.3, -1.3, -2.3]  

    if m1_min < m3:
        raise ValueError(f"This analytical derivation requires IMF break m3  < m1_min ({m3} !< {m1_min})")
    if m1_min > m1_max:
        raise ValueError(f"Minimum sampled primary mass cannot be above maximum sampled primary mass: m1_min ({m1_min} !<  m1_max {m1_max})")
    if m1_max > m4:
        raise ValueError(f"Maximum sampled primary mass cannot be above maximum mass of Kroupa IMF:  m1_max ({m1_max} !<  m4 {m4})")
    
    # normalize IMF over the complete mass range:
    alpha = (-(m4**(-1.3)-m3**(-1.3))/1.3 - (m3**(-0.3)-m2**(-0.3))/(m3*0.3) + (m2**0.7-m1**0.7)/(m2*m3*0.7))**(-1)

    # we want to compute M_stellar_sys_in_universe / N_binaries_in_COMPAS
    #  = N_binaries_in_universe/N_binaries_in_COMPAS * N_stellar_sys_in_universe/N_binaries_in_universe * M_stellar_sys_in_universe/N_stellar_sys_in_universe
    #  = 1/fint * 1/fbin * average mass of a stellar system in the Universe

    # fint =  N_binaries_in_COMPAS/N_binaries_in_universe: fraction of binaries that COMPAS simulates
    fint = -alpha / 1.3 * (m1_max ** (-1.3) - m1_min ** (-1.3)) + alpha * m2_min / 2.3 * (m1_max ** (-2.3) - m1_min ** (-2.3))

    # Next for N_stellar_sys_in_universe/N_binaries_in_universe * M_stellar_sys_in_universe/N_stellar_sys_in_universe
    # N_stellar_sys_in_universe/N_binaries_in_universe = the binary fraction 
    # fbin edges and values are chosen to approximately follow Figure 1 from Offner et al. (2023)
    binary_bin_edges = [m1, 0.08, 0.5, 1, 10, m4]    
    if fbin == None:
        # use a binary fraction that varies with mass
        binaryFractions = [0.1, 0.225, 0.5, 0.8, 1.0] 
    else:
        # otherwise use a constant binary fraction
        binaryFractions = [fbin] * 5

    # M_stellar_sys_in_universe/N_stellar_sys_in_universe = average mass of a stellar system in the Universe,
    # we are computing 1/fbin * M_stellar_sys_in_universe/N_stellar_sys_in_universe, skipping steps this leads to:
    # int_A^B (1/fb(m1) + 0.5) m1 P(m1) dm1. 
    # This is a double piecewise integral, i.e. pieces over the binary fraction bins and IMF mass bins.
    piece_wise_integral = 0

    # For every binary fraction bin
    for i in range(len(binary_bin_edges) - 1):
        fbin = binaryFractions[i] # Binary fraction for this range

        # And every piece of the Kroupa IMF
        for j in range(len(imf_mass_bounds) - 1):
            exponent = IMF_powers[j] # IMF exponent for these masses

            # Check if the binary fraction bin overlaps with the IMF mass bin
            if binary_bin_edges[i + 1] <= imf_mass_bounds[j] or binary_bin_edges[i] >= imf_mass_bounds[j + 1]:
                continue  # No overlap

            # Integrate from the most narrow range
            m_start = max(binary_bin_edges[i], imf_mass_bounds[j])
            m_end = min(binary_bin_edges[i + 1], imf_mass_bounds[j + 1])

            # Compute the definite integral:
            integral = ( m_end**(exponent + 2) - m_start**(exponent + 2) ) / (exponent + 2) * continuity_constants[j]

            # Compute the sum term
            sum_term = (1 /fbin + 0.5) * integral
            piece_wise_integral += sum_term

    # combining them:
    Average_mass_stellar_sys_per_fbin = alpha * piece_wise_integral

    # Now compute the average mass per binary in COMPAS M_stellar_sys_in_universe / N_binaries_in_COMPAS
    M_sf_Univ_per_N_binary_COMPAS = (1/fint) * Average_mass_stellar_sys_per_fbin

    return M_sf_Univ_per_N_binary_COMPAS