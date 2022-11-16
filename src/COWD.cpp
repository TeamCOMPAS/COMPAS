#include "COWD.h"

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
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate from the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
DBL_DBL COWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {

    m_AccretionRegime = DetermineAccretionRegime(p_IsHeRich, p_DonorMassRate); // Check if accretion leads to stage switch for WDs and returns retention efficiency as well.
                                                                               
    double acceptanceRate   = 0.0;                                                       // Acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                       // Accretion fraction - default = 0.0
    double massIntakeRate = std::max(p_DonorMassRate, CalculateEddingtonCriticalRate()); // The incoming matter stream from the donor is Eddington limited

    if (p_IsHeRich) {
        acceptanceRate = massIntakeRate * CalculateEtaHe(massIntakeRate);
    } else {
        acceptanceRate = massIntakeRate * CalculateEtaHe(massIntakeRate) * CalculateEtaH(massIntakeRate); 
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
 * ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich, const double p_DonorMassLossRate) 
 *
 * @param   [IN]    p_HeRich                 Whether the accreted material is helium-rich or not
 * @param   [IN]    p_DonorMassLossRate      Donor mass loss rate, in units of Msol / Myr
 * @return                                   Current WD accretion regime
 */

ACCRETION_REGIME COWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorMassLossRate) {
    double logMdot = log10(p_DonorMassLossRate / MYR_TO_YEAR); // Logarithm of the accreted mass (M_sun/yr)
    ACCRETION_REGIME regime = ACCRETION_REGIME::NONE;

    if (p_HeRich) {
        // The following coefficients in logMassTransfer limits come from table A1 in Piersanti+ 2014.
        double logMassTransferCrit       = MT_LIMIT_PIERSANTI_RG_SS_0 + MT_LIMIT_PIERSANTI_RG_SS_1 *m_Mass;
        double logMassTransferStable     = MT_LIMIT_PIERSANTI_SS_MF_0 + MT_LIMIT_PIERSANTI_SS_MF_1 *m_Mass; // Piersanti+2014 has several Flashes regimes. Here we group them into one.
        double logMassTransferDetonation = MT_LIMIT_PIERSANTI_SF_Dt_0 + MT_LIMIT_PIERSANTI_SF_Dt_1 *m_Mass; // Critical value for double detonation regime in Piersanti+ 2014
        if (utils::Compare(logMdot, logMassTransferStable) < 0) {
            if (utils::Compare(logMdot, logMassTransferDetonation) > 0) {
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
        else if (utils::Compare(logMdot, logMassTransferCrit) > 0) {
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
        // The following coefficients in logMassTransfer limits come from quadratic fits to Nomoto+ 2007 results (table 5) in Mass vs log10 Mdot space, to cover the low-mass end.
        double logMassTransferCrit   = MT_LIMIT_NOMOTO_REDGIANT_0 + MT_LIMIT_NOMOTO_REDGIANT_1 *m_Mass + MT_LIMIT_NOMOTO_REDGIANT_2 *m_Mass*m_Mass;
        double logMassTransferStable =  MT_LIMIT_NOMOTO_STABLE_0  + MT_LIMIT_NOMOTO_STABLE_1   *m_Mass + MT_LIMIT_NOMOTO_STABLE_2   *m_Mass*m_Mass;
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
        return !IsSupernova();
    }
}

/*
 * List all conditions for SN (AIC or SN Ia) for COWD WD. 
 * Each condition should also be a separate clause in EvolveToNextPhase.
 *
 * bool IsSupernova()
 *
 * @return                               Whether WD should undergo AIC or SN Ia
 */

bool COWD::IsSupernova() const {
    return m_DoubleDetonation || IsMassAboveChandrasekhar();      
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
    else if (m_DoubleDetonation) {
        return ResolveSNIa(); // RTW: Is this correct? 
    }
    else if (IsMassAboveChandrasekhar()) {
        return ResolveSNIa(); // RTW: Is this correct? 
    }
    else {                                         // Should not occur
        SHOW_WARN(ERROR::WARNING, "COWD told to evolve, but not how.");                                          // show warning
        return ResolveAIC();
    }
}
