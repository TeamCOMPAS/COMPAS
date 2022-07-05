#include "WhiteDwarfs.h"

/* Pass LogMassRate, MassRate in units of Msun/yr. Logic of eta values follows the one of Claeys2014.
 * Biggest differences are using Piersanti2014 instead of KatoHachisu 2004, and using eta = 1 for accumulation in helium regime.
 * Also, I assumed that mild and strong regimes lead to change in accretion efficiency.
 */

/* Calculate eta_hydrogen from Claeys+ 2014, appendix B. We have changed the mass accretion limits for
 * Nomoto+ 2007 ones, after applying a quadratic fit to cover the low-mass end.
 *
 * double CalculateetaH(const double p_LogMassRate)
 *
 * @param   [IN]    p_LogMassRate        Logarithm of the mass transfer rate (Msun/yr)
 * @return                               eta_hydrogen, "hydrogen accretion efficiency"
 */

double WhiteDwarfs::CalculateetaH(const double p_LogMassRate) {
    double etaH = 0.0;
    double MdotCritH = -0.98023471 * m_Mass * m_Mass + 2.88247131 * m_Mass - 8.33017155;
    double MdotLowH = -1.2137735 * m_Mass * m_Mass + 3.57319872 * m_Mass - 9.21757267;
    if (utils::Compare(p_LogMassRate, MdotCritH) > -1) {
        etaH = PPOW(10, MdotCritH - p_LogMassRate);
    } else if ((utils::Compare(p_LogMassRate, MdotCritH) == -1) && (utils::Compare(p_LogMassRate, MdotLowH) > -1)) {
        etaH = 1.0;
    }
    return etaH;
}

/* Calculate eta_helium from Claeys+ 2014, appendix B. We have changed the mass accretion limits for
 * Piersanti+ 2014 ones. The different flashes regimes from Piersanti+ 2014 have been merged into one,
 * and the accumulation regime has been change so we can get double detonations. Finally, eta_KH04 has
 * also been updated with the accretion efficiency values from Piersanti+ 2014.
 *
 * double CalculateetaHe(const double p_LogMassRate)
 *
 * @param   [IN]    p_LogMassRate        Logarithm of the mass transfer rate (Msun/yr)
 * @return                               eta_hydrogen, "helium accretion efficiency"
 */

double WhiteDwarfs::CalculateetaHe(const double p_LogMassRate) {
    double etaHe = 0.0;
    double MdotCritHe = -6.84 + 1.349 * m_Mass;
    double MdotLowHe = -8.115 + 2.29 * m_Mass;
    double MdotAccumulation =  -8.313 + 1.018 * m_Mass;

    if (utils::Compare(p_LogMassRate, MdotCritHe) > -1) {
        etaHe = PPOW(10, MdotCritHe - p_LogMassRate);
    } else if ((utils::Compare(p_LogMassRate, MdotCritHe) == -1) && (utils::Compare(p_LogMassRate, MdotLowHe) > -1)) {
        etaHe = 1.0;
    } else if ((utils::Compare(p_LogMassRate, MdotLowHe) == -1) && (utils::Compare(p_LogMassRate, MdotAccumulation) > -1)) {
        etaHe = CalculateetaPTY(PPOW(10, p_LogMassRate));
    } else {
        etaHe = 1.0; // Modified so we can have double detonations
    }
    return etaHe;
}



/* Calculate accretion efficiency as indicated in Piersanti+ 2014. Their recipe works
 * for specific mass values, so a better implementation requires interpolation and
 * extrapolation (specially towards the low-mass end). Right now, we just adopt a
 * piece-wise approach.
 *
 * double CalculateetaPTY(const double p_MassRate)
 *
 * @param   [IN]    p_MassRate           Mass transfer rate (Msun/yr)
 * @return                               etaPTY, accretion efficency during flashes as per Piersanti+ 2014
 */
double WhiteDwarfs::CalculateetaPTY(const double p_MassRate) {
    double etaPTY;
    if (utils::Compare(m_Mass, 0.6) < 1) {
        etaPTY = 6e-3 + 5.1e-2*p_MassRate + 8.3e-3*PPOW(p_MassRate, 2) - 3.317e-4*PPOW(p_MassRate,3);
    } else if ((m_Mass <= 0.7) && (m_Mass > 0.6)) {
        etaPTY = -3.5e-2 + 7.5e-2*p_MassRate - 1.8e-3*PPOW(p_MassRate, 2) + 3.266e-5*PPOW(p_MassRate,3);
    } else if ((utils::Compare(m_Mass, 0.81) < 1) && (utils::Compare(m_Mass, 0.7) > 0)) {
        etaPTY = 9.3e-2 + 1.8e-2*p_MassRate + 1.6e-3*PPOW(p_MassRate, 2) - 4.111e-5*PPOW(p_MassRate,3);
    } else if ((utils::Compare(m_Mass, 0.92) < 1) && (utils::Compare(m_Mass, 0.81) > 0)) {
        etaPTY = -7.59e-2 + 1.54e-2*p_MassRate + 4e-4*PPOW(p_MassRate, 2) - 5.905e-6*PPOW(p_MassRate,3);
    } else {
        etaPTY = -0.323 + 4.1e-2*p_MassRate - 7e-4*PPOW(p_MassRate, 2) + 4.733e-6*PPOW(p_MassRate,3);
    }
    return etaPTY;
}


/*
 * Calculate the luminosity of a White Dwarf as it cools
 *
 * Hurley et al. 2000, eq 90
 *
 *
 * double CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity, const double p_BaryonNumber)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time since White Dwarf formation in Myr
 * @param   [IN]    P_Metallicity               Metallicity of White Dwarf
 * @param   [IN]    p_BaryonNumber              Baryon number - differs per White Dwarf type (HeWD, COWD, ONeWD)
 * @return                                      Luminosity of a White Dwarf in Lsol
 */
double WhiteDwarfs::CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity, const double p_BaryonNumber) {
    return (635.0 * p_Mass * PPOW(p_Metallicity, 0.4)) / PPOW(p_BaryonNumber * (p_Time + 0.1), 1.4);
}


/* Calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the retention efficiency parameter
 *
 *
 * For a given mass transfer rate, this function computes the amount of mass a WD would retain after
 * flashes, as given by appendix B of Claeys+ 2014.
 *
 *
 * DBL_DBL CalculateWDMassAcceptanceRate(const double p_LogDonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_LogDonorMassRate          Logarithm of the mass transfer rate of the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and Retention Efficiency Parameter
 */
DBL_DBL WhiteDwarfs::CalculateWDMassAcceptanceRate(const double p_LogDonorMassRate, const bool p_IsHeRich) {

    double acceptanceRate   = 0.0;                                                          // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                          // accretion fraction - default=0.0
    double thisMassRate;
    double DonorMassRate = PPOW(10,p_LogDonorMassRate);

    if (p_IsHeRich) {
        thisMassRate = DonorMassRate * CalculateetaHe(p_LogDonorMassRate);
    } else {
        thisMassRate = DonorMassRate * CalculateetaHe(p_LogDonorMassRate) * CalculateetaH(p_LogDonorMassRate);
    }
    acceptanceRate   = thisMassRate;
    fractionAccreted = acceptanceRate / DonorMassRate;
    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/*
 * Calculate the radius of a white dwarf - good for all types of WD
 *
 * Hurley et al. 2000, eq 91 (from Tout et al. 1997)
 *
 *
 * double CalculateRadiusOnPhase_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius of a White Dwarf in Rsol (since WD is ~ Earth sized, expect answer around 0.009)
 */
double WhiteDwarfs::CalculateRadiusOnPhase_Static(const double p_Mass) {
    double MCH_Mass_one_third  = std::cbrt(MCH / p_Mass); 
    double MCH_Mass_two_thirds = MCH_Mass_one_third* MCH_Mass_one_third;
    return std::max(NEUTRON_STAR_RADIUS, 0.0115 * std::sqrt((MCH_Mass_two_thirds - 1.0/MCH_Mass_two_thirds )));
}



/* Increase shell size after mass transfer episode. Hydrogen and helium shells are kept separately.
 *
 * void IncrementShell(const double p_AccretedMass, bool p_HeRich) {
 *
 * @param   [IN]    p_AccretedMass              Mass accreted
 * @param   [IN]    p_HeRich                    Material is He-rich or not
 */
void WhiteDwarfs::IncrementShell(const double p_AccretedMass, const bool p_HeRich) {
	if (p_HeRich) {
		m_HeShell += p_AccretedMass;
	} else {
		m_HShell += p_AccretedMass;
	}
}
