#include "ONeWD.h"


/*
 * Calculate the luminosity of an OXYGEN-NEON White Dwarf as it cools
 *
 * Hurley et al. 2000, eq 90
 *
 *
 * double CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time since White Dwarf formation in Myr
 * @param   [IN]    P_Metallicity               Metallicity of White Dwarf
 * @return                                      Luminosity of a White Dwarf in Lsol
 */
double ONeWD::CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity) {
    return (635.0 * p_Mass * PPOW(p_Metallicity, 0.4)) / PPOW(WD_Baryon_Number.at(STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF) * (p_Time + 0.1), 1.4);
}
