#include "HeWD.h"

// RTW: What papers should be referenced here? Piersanti+ 2014 seems to only cover He accretion
/* For HeWD, calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the retention efficiency parameter
 *
 *
 * For a given mass transfer rate, this function computes the amount of mass a WD would retain after
 * flashes, as given by appendix B of Claeys+ 2014. 
 * https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract 
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_LogDonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_LogDonorMassRate          Logarithm of the mass transfer rate of the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
DBL_DBL HeWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {

    m_AccretionRegime = DetermineAccretionRegime(p_IsHeRich, p_DonorMassRate); // Check if accretion leads to stage switch for WDs and returns retention efficiency as well.
                                                                               
    double acceptanceRate   = 0.0;                                                          // acceptance mass rate - default = 0.0
    double fractionAccreted = 1.0;                                                          // accretion fraction - default=1.0
    // RTW: I've moved the discussion on the correct fractionAccreted here.
    // RTW: What should the fraction be here? Is this missing a CalculateMassAcceptance? A fraction of 1 seems suspiciously high                                                                                           
    /* NRS: the literature has not explored what happenes to HeWDs accreting as much as the COWD case. 
     * Using the COWD prescription with ONe WDs is "okay" as they cover a similar range in mass, 
     * but all HeWDs are low mass and the only prescriptions I could find were from StarTrack, which basically motivate all of the following.
     * One of my side-projects is to create a better prescription, but for now this should be used in the same way as the COWD function (I see that you are not returning the fraction value)
     */
    // RTW: I've moved around the functions for clarity. There is now a CalculateMassAcceptanceRate and a DetermineAccretionRegime in both HeWD and COWD (where ONeWD takes from COWD). 
    // The accretion fraction should now be correctly outputed here, though I am still skeptical about setting the accretion fraction to 1. 
    // I realize that this is currently possible for the other WD types as well, but I can't imagine that you would get fully conservative MT if,
    // e.g, you had a giant donor losing its entire envelope on the giant's thermal timescale. 
                                                                               
    if (m_AccretionRegime == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES) { // Flashes restrict accumulation
        fractionAccreted = 0.0;
    }

    acceptanceRate = p_DonorMassRate * fractionAccreted;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/* 
 * Calculate the WD accretion regime based on the MT rate and whether the donor is He rich 
 *
 * The accretion regime is one of the following:
 *
 * Helium accretion that could lead to Sub-Chandrasekhar SN
 * Helium accretion that could lead to helium ignition and rejuvenation
 * Hydrogen Flashes
 * Hydrogen Accumulation
 * Note that we have merged the different flashes regimes from Piersanti+ 2014 into a single regime.
 *
 * ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich, const double p_DonorThermalMassLossRate) 
 *
 * @param   [IN]    p_HeRich                        Whether the accreted material is helium-rich or not
 * @param   [IN]    p_DonorThermalMassLossRate      Donor thermal mass loss rate, in units of Msol / Myr
 * @return                                          Current WD accretion regime
 */

ACCRETION_REGIME HeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorThermalMassLossRate) {
    double Mdot = p_DonorThermalMassLossRate / MYR_TO_YEAR; // Accreted mass rate (M_sun/yr)
    double logMdot = log10(Mdot);                           // Logarithm of the accreted mass rate (M_sun/yr)
    ACCRETION_REGIME regime;
    if (p_HeRich) {
        if (utils::Compare(logMdot, HELIUM_WHITE_DWARF_MCRIT) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR; // Could lead to Sub CH SN Ia
            double massSubCh = -4e8 * Mdot + 1.34; // Minimum mass for Sub-Ch Mass detonation. Eq 62, Belczynski+ 2008.
            if (utils::Compare(m_Mass, massSubCh) >= 0 ) {
                m_IsSubChandrasekharTypeIa = true;
            }
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION; // Could lift degeneracy and evolve into He MS. Requires minimum mass ! on top of the shell size
            if (utils::Compare(m_Mass, MASS_HELIUM_BURN) >= 0) {
                if (utils::Compare(Mdot, 1.64e-6) < 0) { // Accretion limit from eq 61, Belczynski+ 2008.
                    double shellCrit = -7.8e-4 * Mdot + 1.34; // Minimum shell mass of He for detonation. Eq 61, Belczynski+ 2008. This helium should not be burnt, but not implemented this yet. Ruiter+ 2014.
                    if (utils::Compare(m_HeShell, shellCrit) >= 0) {
                        m_ShouldRejuvenate = true;
                    }
                } 
                else {
                    m_ShouldRejuvenate = true;
                }
            }
        }
    } 
    else {
        double Mcrit = log10(m_l0Ritter * PPOW(m_Mass, m_lambdaRitter) / (m_XRitter * 6e18)); // Eq. 60 in Belczynski+ 2008. 6e18 is the energy yield of H burning in ergs/g.
        if (utils::Compare(logMdot, Mcrit) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES; // Flashes restrict accumulation
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION; //Material piles up on the WD. Leads to merger or CEE.
        }
    }

    return regime;
}


/*
 * Allow the evolution towards an HeMS or SN.
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                               Whether the WD should evolve on phase or towards an HeMS/SN.
 */
bool HeWD::ShouldEvolveOnPhase() {
    if (m_ShouldRejuvenate) {
        return false;
    }
    else {
        return !IsSupernova();
    }
}

/*
 * List all conditions for SN (AIC or SN Ia) for HeWD WD. 
 * Each condition should also be a separate clause in EvolveToNextPhase.
 *
 * bool IsSupernova()
 *
 * @return                               Whether WD should undergo AIC or SN Ia
 */
bool HeWD::IsSupernova() const {
    return m_IsSubChandrasekharTypeIa;                                           // Go supernova if mass and He shell are large enough
}

/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */
STELLAR_TYPE HeWD::EvolveToNextPhase() {
    if (m_ShouldRejuvenate) {
        return STELLAR_TYPE::NAKED_HELIUM_STAR_MS; 
    }
    else if (m_IsSubChandrasekharTypeIa) {         // Currently, assume a Type Ia from a HeWD is the same as other WDs. May want to vary in the future
        return ResolveSNIa();
    }
    else {                                         // Should not occur
        SHOW_WARN(ERROR::WARNING, "HeWD told to evolve, but not how.");                                          // show warning
        return ResolveAIC();
    }
}


