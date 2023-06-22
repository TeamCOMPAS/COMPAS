#include "COWD.h"

/* For COWDs, calculate:
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
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate from the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
DBL_DBL COWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {

    m_AccretionRegime = DetermineAccretionRegime(p_IsHeRich, p_DonorMassRate); 
                                                                               
    double acceptanceRate   = 0.0;                                                       // Acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                       // Accretion fraction - default = 0.0

    acceptanceRate = p_DonorMassRate * CalculateEtaHe(p_DonorMassRate);
    if (!p_IsHeRich) acceptanceRate *= CalculateEtaH(p_DonorMassRate);
    fractionAccreted = acceptanceRate / p_DonorMassRate;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/* 
 * Determine the WD accretion regime based on the MT rate and whether the donor is He rich. Also,
 * initialize He-Shell detonation or Off-center ignition when necessary, by changing the value
 * of m_HeShellDetonation or m_OffCenterIgnition (respectively).
 *
 * The accretion regime is one of the options listed in enum ACCRETION_REGIME (constants.h)
 *
 * Note that we have merged the different flashes regimes from Piersanti+ 2014 into a single regime.
 *
 * ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich, const double p_DonorMassLossRate) 
 *
 * @param   [IN]    p_HeRich                 Whether the accreted material is helium-rich or not
 * @param   [IN]    p_DonorMassLossRate      Donor mass loss rate, in units of Msol / Myr
 * @return                                   Current WD accretion regime
 */
ACCRETION_REGIME COWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorMassLossRate) {

    double logMdot          = log10(p_DonorMassLossRate / MYR_TO_YEAR);                                                     // Logarithm of the accreted mass (M_sun/yr)
    ACCRETION_REGIME regime = ACCRETION_REGIME::NONE;

    if (p_HeRich) {
        // The following coefficients in logMassTransfer limits come from table A1 in Piersanti+ 2014.
        double logMassTransferCrit       = WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_0 + WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_1 * m_Mass;
        double logMassTransferStable     = WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_0 + WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_1 * m_Mass;  // Piersanti+2014 has several Flashes regimes. Here we group them into one.
        double logMassTransferDetonation = WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_0 + WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_1 * m_Mass;  // Critical value for double detonation regime in Piersanti+ 2014
        if (utils::Compare(logMdot, logMassTransferStable) < 0) {
            if (utils::Compare(logMdot, logMassTransferDetonation) > 0) {
                regime = ACCRETION_REGIME::HELIUM_FLASHES;
            } 
            else {
                regime = ACCRETION_REGIME::HELIUM_ACCUMULATION;
                if ((utils::Compare(m_Mass, MASS_DOUBLE_DETONATION_CO) >= 0) && (utils::Compare(m_HeShell, WD_HE_SHELL_MCRIT_DETONATION) >= 0)) {
                    m_HeShellDetonation = true;                                                                             // JR: Question: should this be set false if the condition is not satisfied?
                }
            }
        } 
        else if (utils::Compare(logMdot, logMassTransferCrit) > 0) {
            regime = ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS;
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_STABLE_BURNING;
            if ((utils::Compare(logMdot, COWD_LOG_MDOT_MIN_OFF_CENTER_IGNITION) > 0) && (utils::Compare(m_Mass, COWD_MASS_MIN_OFF_CENTER_IGNITION) > 0)) {
                m_OffCenterIgnition = true;                                                                                 // JR: Question: should this be set false if the condition is not satisfied?
            }
        }
    } 
    else {
        // The following coefficients in logMassTransfer limits come from quadratic fits to Nomoto+ 2007 results (table 5) in Mass vs log10 Mdot space, to cover the low-mass end.
        double m_Mass_2 = m_Mass * m_Mass;
        double logMassTransferCrit   = WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_0 + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_1 * m_Mass + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_2 * m_Mass_2;
        double logMassTransferStable = WD_LOG_MT_LIMIT_NOMOTO_STABLE_0   + WD_LOG_MT_LIMIT_NOMOTO_STABLE_1   * m_Mass + WD_LOG_MT_LIMIT_NOMOTO_STABLE_2   * m_Mass_2;

        if (utils::Compare(logMdot, logMassTransferStable) < 0) {
            regime = ACCRETION_REGIME::HYDROGEN_FLASHES;
        } 
        else if (utils::Compare(logMdot, logMassTransferCrit) > 0) {
            regime = ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS;
        } 
        else {
            regime = ACCRETION_REGIME::HYDROGEN_STABLE_BURNING;
        }
    }

    return regime;
}


/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */

STELLAR_TYPE COWD::EvolveToNextPhase() {

    STELLAR_TYPE stellarType;

    if (m_OffCenterIgnition) {
        stellarType = STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    }
    else {                                         
        stellarType = ResolveSNIa(); 
    }
    return stellarType;
}
