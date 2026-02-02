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
    fbinary_bin_edges = [m1, 0.08, 0.5, 1, 10, m4]
    
    def get_binary_fraction(mass):
        binaryFractions = [0.1, 0.225, 0.5, 0.8, 1.0]
        for i in range(len(fbinary_bin_edges) - 1):
            if mass < fbinary_bin_edges[0]:  
                # Mass below lowest binary fraction bin edge (shouldn't happen)
                return binaryFractions[0]
            if fbinary_bin_edges[i] <= mass < fbinary_bin_edges[i + 1]:
                return binaryFractions[i]
            if mass >= fbinary_bin_edges[-1]:
                # Mass above highest binary fraction bin edge (shouldn't happen)
                return binaryFractions[-1]

    def integrand_full(mass, f_bin):
        local_f_bin = get_binary_fraction(mass) if f_bin is None else f_bin
        expected_q = quad(lambda q: q * mass_ratio_pdf_function(q), 0, 1)[0]
        # mass of single stars = (1 - f_bin) * m1 
        # mass of binaries = f_bin * (1 + <q>) * m1 
        expected_mass_all_stellar_sys =(1 + local_f_bin * expected_q) * mass * IMF(mass, m1, m2, m3, m4, a12, a23, a34) 
        return expected_mass_all_stellar_sys

    def integrand_compas(mass, f_bin):
        local_f_bin = get_binary_fraction(mass) if f_bin is None else f_bin
        # Only binaries contribute in COMPAS population
        # Integrand is (1 + q) * f_bin * m1 * P(m1) * P(q),
        q_min = m2_low / mass
        if q_min >= 1:
            return 0  # No valid secondaries
        # Integrate (1 + q)P(q) dq over q from q_min to 1, 
        # we get p(q)dq:
        p_qdq = quad(mass_ratio_pdf_function, q_min, 1)[0]
        # and q P(q) dq (= expected_q)
        expected_q = quad(lambda q: q * mass_ratio_pdf_function(q), q_min, 1)[0]

        expected_mass_compas_binaries = (p_qdq + expected_q) * local_f_bin * mass * IMF(mass, m1, m2, m3, m4, a12, a23, a34)
        return expected_mass_compas_binaries


    # split integral at binary fraction steps if f_bin is None (i.e. variable and like a step function)
    def split_integral(func, a, b, f_bin):
        total = 0
        for edge_start, edge_end in zip(fbinary_bin_edges[:-1], fbinary_bin_edges[1:]):
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

    fraction = compas_mass / full_mass
    return fraction


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



###################################################
# Analytical calculation of star forming mass per binary
###################################################
def analytical_star_forming_mass_per_binary_using_kroupa_imf(
    m1_min, m1_max, m2_min, fbin=1.0, imf_mass_bounds=(0.01, 0.08, 0.5, 200.0)):
    """
    Takes: 
        m1_min, m1_max, m2_min: COMPAS mass ranges [Msun]
        fbin: binary fraction (if None, use piecewise constant fbin(m1))
        imf_mass_bounds: Kroupa IMF mass bounds [Msun]
    Computes:
      N_bin_in_COMPAS 
      average_stellar_mass_sys = M_sys,Univ / N_sys,Univ (blue)
      M_sf_Univ_per_N_binary_COMPAS = average_stellar_mass_sys / N_bin_in_COMPAS

    Assumes:
      Kroupa IMF: P(m1) = alpha * C_cont,i * m1^{gamma_i}  on IMF segment i
      Piecewise constant on binary-fraction: f_bin(m1) with bins j
      Flat mass ratio: P(q)=U(0,1)  =>  P(m2>m2_min | m1) = 1 - m2_min/m1  (for m1 >= m2_min)
    """
    # -------------------------
    # Kroupa IMF 
    # -------------------------
    m1, m2, m3, m4 = imf_mass_bounds
    continuity_constants = [1.0 / (m2 * m3), 1.0 / m3, 1.0]
    IMF_powers = [-0.3, -1.3, -2.3]

    if m1_min > m1_max:
        raise ValueError(f"Require m1_min <= m1_max, got {m1_min} > {m1_max}.")
    if m1_min < m1:
        raise ValueError(f"Require m1_min >= {m1} (Universe lower IMF bound).")
    if m1_max > m4:
        raise ValueError(f"Require m1_max <= {m4} (Universe upper IMF bound).")

    # Normalization alpha over full IMF range [m1, m4]
    alpha = (
        - (m4 ** (-1.3) - m3 ** (-1.3)) / 1.3
        - (m3 ** (-0.3) - m2 ** (-0.3)) / (m3 * 0.3)
        + (m2 ** (0.7) - m1 ** (0.7)) / (m2 * m3 * 0.7)
    ) ** (-1)

    # -------------------------
    # Binary-fraction bins
    # -------------------------
    fbinary_bin_edges = [m1, 0.08, 0.5, 1.0, 10.0, m4]
    if fbin is None:
        binaryFractions = [0.1, 0.225, 0.5, 0.8, 1.0]
    else:
        binaryFractions = [float(fbin)] * (len(fbinary_bin_edges) - 1)

    # -------------------------
    # Helpers: overlaps and power integrals
    # -------------------------
    def overlap(lo1, hi1, lo2, hi2):
        lo = max(lo1, lo2)
        hi = min(hi1, hi2)
        return lo, hi

    def int_power(lo, hi, power):
        """∫_lo^hi m^{power} dm, for power != -1."""
        if hi <= lo:
            return 0.0
        return (hi ** (power + 1) - lo ** (power + 1)) / (power + 1)

    # -------------------------
    # Average_stellar_mass_sys = ∫_Univ (1+0.5 fbin) m1 P(m1) dm1
    # -------------------------
    av_Mstar_integral_sum = 0.0
    for i in range(len(IMF_powers)):  # IMF segment index i
        gamma_i = IMF_powers[i]
        C_cont_i = continuity_constants[i]
        imf_lo, imf_hi = imf_mass_bounds[i], imf_mass_bounds[i + 1]

        for j in range(len(binaryFractions)):  # fbin bin index j
            fbin_j = binaryFractions[j]
            fb_lo, fb_hi = fbinary_bin_edges[j], fbinary_bin_edges[j + 1]

            A, B = overlap(imf_lo, imf_hi, fb_lo, fb_hi)
            if B <= A:
                continue

            # integrand: (1 + 0.5 fbin_j) * m1 * (alpha*C_cont_i*m1^{gamma_i})
            # => alpha * (1+0.5 fbin_j) * C_cont_i * ∫ m1^{gamma_i+1} dm1
            av_Mstar_integral_sum += (1.0 + 0.5 * fbin_j) * C_cont_i * int_power(A, B, gamma_i + 1)

    average_stellar_mass_sys = alpha * av_Mstar_integral_sum

    # -------------------------
    # N_bin_in_COMPAS = ∫_{m1_min}^{m1_max} P(m1) fbin(m1) (1 - m2_min/m1) dm1
    # -------------------------
    N_bin_compas_sum = 0.0
    for i in range(len(IMF_powers)):  # IMF segment index i
        gamma_i = IMF_powers[i]
        C_cont_i = continuity_constants[i]
        imf_lo, imf_hi = imf_mass_bounds[i], imf_mass_bounds[i + 1]

        for j in range(len(binaryFractions)):  # fbin bin index j
            fbin_j = binaryFractions[j]
            fb_lo, fb_hi = fbinary_bin_edges[j], fbinary_bin_edges[j + 1]

            # overlap additionally with COMPAS m1-range
            A, B = overlap(imf_lo, imf_hi, fb_lo, fb_hi)
            A, B = overlap(A, B, m1_min, m1_max)
            if B <= A:
                continue

            # integrand: alpha*C_cont_i*m^{gamma_i} * fbin_j * (1 - m2_min/m)
            # => alpha*fbin_j*C_cont_i * ∫ (m^{gamma_i} - m2_min*m^{gamma_i-1}) dm
            term_main = int_power(A, B, gamma_i)              # ∫ m^{gamma_i} dm
            term_m2   = int_power(A, B, gamma_i - 1)          # ∫ m^{gamma_i-1} dm
            N_bin_compas_sum += fbin_j * C_cont_i * (term_main - m2_min * term_m2)

    N_bin_in_COMPAS = alpha * N_bin_compas_sum

    if N_bin_in_COMPAS <= 0.0:
        raise ValueError("Computed N_bin_in_COMPAS <= 0; check bounds and m2_min.")

    M_sf_Univ_per_N_binary_COMPAS = average_stellar_mass_sys / N_bin_in_COMPAS

    return M_sf_Univ_per_N_binary_COMPAS

