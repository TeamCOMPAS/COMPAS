#include "WhiteDwarfs.h"

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
    return std::max(NEUTRON_STAR_RADIUS, 0.0115 * sqrt(PPOW(MCH / p_Mass, 2.0 / 3.0) - PPOW(p_Mass / MCH, 2.0 / 3.0)));
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * Hurley et al. 2000, just after eq 106
 *
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double WhiteDwarfs::CalculateMassLossRateHurley() {
    return (utils::Compare(m_Luminosity, 4.0E3) > 0) ? CalculateMassLossRateNieuwenhuijzenDeJager() : 0.0;
}


/*
 * Calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the accretion efficiency parameter
 *
 *
 * The maximum acceptance rate of the accretor star during mass transfer is based on stellar type: this function
 * is for white dwarfs.
 *
 * Mass transfer is assumed Eddington limited for BHs and NSs.  The formalism of Nomoto/Claeys is used for WDs.
 *
 * For non compact objects:
 *
 *    1) Kelvin-Helmholtz (thermal) timescale if THERMAL (thermally limited) mass transfer efficiency
 *    2) Choose a fraction of the mass rate that will be effectively accreted for FIXED fraction mass transfer (as in StarTrack)
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_AccretorMassRate          Thermal mass loss rate of the accretor (this star) - ignored here
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
DBL_DBL WhiteDwarfs::CalculateMassAcceptanceRate(const double p_DonorMassRate,  const double p_AccretorMassRate) {

    double acceptanceRate   = 0.0;                                                          // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                          // accretion fraction - default=0.0

    double mDotEddington    = CalculateEddingtonCriticalRate();

    acceptanceRate   = std::min(mDotEddington, p_DonorMassRate);
    fractionAccreted = acceptanceRate / p_DonorMassRate;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/*
 * Choose timestep for evolution
 *
 * Can obviously do this your own way
 * Given in the discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double WhiteDwarfs::ChooseTimestep(const double p_Time) const {

    double dtk = std::min(std::max(0.1, 10.0 * std::max(0.1, 10.0 * p_Time)), 5.0E2);
    double dte = dtk;

    return std::max(dte, NUCLEAR_MINIMUM_TIMESTEP);
}