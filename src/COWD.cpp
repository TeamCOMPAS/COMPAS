#include "COWD.h"


/*
 * Calculate the luminosity of a CARBON-OXYGEN White Dwarf as it cools
 *
 * Hurley et al. 2000, eq 90
 *
 *
 * double CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double P_Metallicity)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time since White Dwarf formation in Myr
 * @param   [IN]    P_Metallicity               Metallicity of White Dwarf
 * @return                                      Luminosity of a White Dwarf in Lsol
 */
double COWD::CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity) {
    return (635.0 * p_Mass * PPOW(p_Metallicity, 0.4)) / PPOW(WD_Baryon_Number.at(STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF) * (p_Time + 0.1), 1.4);
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE COWD::EvolveToNextPhase() {

    m_Mass       = 0.0;
    m_Radius     = 0.0;
    m_Luminosity = 0.0;
    m_Age        = 0.0;

    return STELLAR_TYPE::MASSLESS_REMNANT;
}
