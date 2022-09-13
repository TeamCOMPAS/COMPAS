#include "COWD.h"

// RTW: Does this only apply for COWDs and ONeWDs? If so, should move into COWDs.
// NRS: So far, yes. The only exeption in HeWDs (where the default is to accrete everything) is to accrete nothing when ACCRETION_REGIME is HELIUM_WHITE_DWARF_HYDROGEN_FLASHES. I want to explore this as a side project, but for now that is all we have.
/* For COWDs and ONeWDs, calculate:
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
DBL_DBL COWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {

    m_AccretionRegime = DetermineAccretionRegime(p_IsHeRich, p_DonorMassRate); // Check if accretion leads to stage switch for WDs and returns retention efficiency as well.
                                                                               
    double acceptanceRate   = 0.0;                                                          // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                          // accretion fraction - default=0.0
    double logDonorMassRate = log10(p_DonorMassRate);

    if (p_IsHeRich) {
        acceptanceRate = p_DonorMassRate * CalculateetaHe(logDonorMassRate);
    } else {
        acceptanceRate = p_DonorMassRate * CalculateetaHe(logDonorMassRate) * CalculateetaH(logDonorMassRate); 
    }
    fractionAccreted = acceptanceRate / p_DonorMassRate;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/* 
 * Calculate the WD accretion regime based on the MT rate and whether the donor is He rich 
 *
 * The accretion regime is one of the following:
 *
 * Helium Accumulation
 * Helium Flashes
 * Helium Stable Burning
 * Helium Optically-Thick Winds
 * Hydrogen Flashes
 * Hydrogen Stable Burning
 * Hydrogen Optically-Thick Winds
 *
 * Note that we have merged the different flashes regimes from Piersanti+ 2014 into a single regime.
 *
 * ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich, const double p_DonorThermalMassLossRate) 
 *
 * @param   [IN]    p_HeRich                        Whether the accreted material is helium-rich or not
 * @param   [IN]    p_DonorThermalMassLossRate      Donor thermal mass loss rate, in units of Msol / Myr
 * @return                                          Current WD accretion regime
 */

ACCRETION_REGIME COWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorThermalMassLossRate) {
    double logMdot = log10(p_DonorThermalMassLossRate / MYR_TO_YEAR); // Logarithm of the accreted mass (M_sun/yr)
    ACCRETION_REGIME regime; // RTW: Do we want a default?

    if (p_HeRich) {
        // The following coefficients in massTransfer limits come from table A1 in Piersanti+ 2014.
        // RTW: Add "log" and "threshold" to the variable names
        double massTransferCrit = MT_LIMIT_CRIT_PIERSANTI_0 + MT_LIMIT_CRIT_PIERSANTI_1 * m_Mass;
        double massTransferStable = MT_LIMIT_STABLE_PIERSANTI_0 + MT_LIMIT_STABLE_PIERSANTI_1 * m_Mass; // Piersanti+2014 has several Flashes regimes. Here we group them into one.
        double massTransferDetonation = MT_LIMIT_DET_PIERSANTI_0 + MT_LIMIT_DET_PIERSANTI_1 * m_Mass; // Critical value for double detonation regime in Piersanti+ 2014
        if (utils::Compare(logMdot, massTransferStable) < 0) {
            if (utils::Compare(logMdot, massTransferDetonation) > 0) {
                regime = ACCRETION_REGIME::HELIUM_FLASHES;
            } 
            else {
                regime = ACCRETION_REGIME::HELIUM_ACCUMULATION;
                m_DoubleDetonation = true;
                if ((utils::Compare(m_Mass, MASS_DOUBLE_DETONATION_CO) >= 0) && (utils::Compare(m_HeShell, SHELL_CRIT) >= 0)) {
                    // RTW: TODO?
                }
            }
        } 
        else if (utils::Compare(logMdot, massTransferCrit) > 0) {
            regime = ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS;
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_STABLE_BURNING;
            if ((utils::Compare(logMdot, MDOT_OFF_C) > 0) && (utils::Compare(m_Mass, 1.33) > 0)) { // The 1.33 Msol value in the comparison is taken from Wang, Podsiadlowski & Han (2017), sect 3.2.
                m_OffCenterIgnition = true;
            }
        }
    } 
    else {
        // The following coefficients in massTransfer limits come from quadratic fits to Nomoto+ 2007 results (table 5) in Mass vs log10 Mdot space, to cover the low-mass end.
        double massTransferCrit = MT_LIMIT_CRIT_NOMOTO_0 +  MT_LIMIT_CRIT_NOMOTO_1 * m_Mass +  MT_LIMIT_CRIT_NOMOTO_2 * m_Mass * m_Mass;
        double massTransferStable =  MT_LIMIT_STABLE_NOMOTO_0 +  MT_LIMIT_STABLE_NOMOTO_1 * m_Mass +  MT_LIMIT_STABLE_NOMOTO_2 * m_Mass * m_Mass;
        if (utils::Compare(logMdot, massTransferStable) < 0) {
            regime = ACCRETION_REGIME::HYDROGEN_FLASHES;
        } 
        else if (utils::Compare(logMdot, massTransferCrit) > 0) {
            regime = ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS;
        } 
        else {
            regime = ACCRETION_REGIME::HYDROGEN_STABLE_BURNING;
        }
    }

    return regime;
}

/*
 * Allow the evolution towards an ONe WD . From https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.1593W/abstract around
 * the end of section 3.2. Also, allows SN.
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                               Whether the WD should evolve on phase or towards an ONeWD/SN.
 */

bool COWD::ShouldEvolveOnPhase() {
    if (m_OffCenterIgnition) {
        return false;
    }
    else {
        return (m_Mass <= MCH);
    }
}

/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */

STELLAR_TYPE COWD::EvolveToNextPhase() {
    if (m_OffCenterIgnition) {
        return STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    }
    else {
        m_Mass       = 0.0;
        m_Radius     = 0.0;
        m_Luminosity = 0.0;
        m_Age        = 0.0;
        return STELLAR_TYPE::MASSLESS_REMNANT;
    }
}
