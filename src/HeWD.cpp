#include "HeWD.h"
#include "HeMS.h"


/* For HeWD, calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the retention efficiency parameter
 *
 *
 * For He WDs, we calculate the mass accretion rate following the 
 * StarTrack prescription (Belczynski+ 2008, sect 5.7.1).
 * https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_LogDonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_LogDonorMassRate          Logarithm of the mass transfer rate of the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
DBL_DBL HeWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {

    m_AccretionRegime       = DetermineAccretionRegime(p_IsHeRich, p_DonorMassRate);                                    // Check if accretion leads to stage switch for WDs and returns retention efficiency as well.
                                                                               
    double acceptanceRate   = 0.0;                                                                                      // acceptance mass rate - default = 0.0
    double fractionAccreted = m_AccretionRegime == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES ? 0.0 : 1.0;   // accretion fraction - default = 1.0, but flashes restrict accumulation

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/* 
 * Determine the WD accretion regime based on the MT rate and whether the donor is He rich. Also,
 * initialize Sub-Chandrasekhar SN Ia or rejuvenation (evolution into HeMS) when necessary, by 
 * changing the value of m_IsSubChandrasekharTypeIa or m_ShouldRejuvenate (respectively). 
 *
 * The accretion regime is one of the options listed in enum ACCRETION_REGIME (constants.h)
 *
 * ACCRETION_REGIME DetermineAccretionRegime(const bool p_HeRich, const double p_DonorMassLossRate) 
 *
 * @param   [IN]    p_HeRich             Whether the accreted material is helium-rich or not
 * @param   [IN]    p_DonorMassRate      Donor mass loss rate, in units of Msol / Myr
 * @return                               Current WD accretion regime
 */
ACCRETION_REGIME HeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorMassRate) {
    double Mdot = p_DonorMassRate / MYR_TO_YEAR;                                                        // Accreted mass rate (M_sun/yr)
    ACCRETION_REGIME regime;
    if (p_HeRich) {
        if (utils::Compare(Mdot, HEWD_HE_MDOT_CRIT) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR;                     // Could lead to Sub CH SN Ia
            double massSubCh = -4e8 * Mdot + 1.34;                                                      // Minimum mass for Sub-Ch Mass detonation. Eq 62, Belczynski+ 2008.
            if (utils::Compare(m_Mass, massSubCh) >= 0 ) {
                m_IsSubChandrasekharTypeIa = true;                                                      // JR: Question: should this be set false if the condition is not satisfied?
            }
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION;                              // Could lift degeneracy and evolve into He MS. Requires minimum mass ! on top of the shell size
            if (utils::Compare(m_Mass, HEWD_MINIMUM_MASS_IGNITION) >= 0) {
                if (utils::Compare(Mdot, 1.64e-6) < 0) {                                                // Accretion limit from eq 61, Belczynski+ 2008.
                    double mCritHeShell = -7.8e-4 * Mdot + 1.34;                                        // Minimum shell mass of He for detonation. Eq 61, Belczynski+ 2008. This helium should not be burnt, but not implemented this yet. Ruiter+ 2014.
                    if (utils::Compare(m_HeShell, mCritHeShell) >= 0) {
                        m_ShouldRejuvenate = true;                                                      // JR: Question: should this be set false if the condition is not satisfied?
                    }
                } 
                else {
                    m_ShouldRejuvenate = true;
                }
            }
        }
    } 
    else {
        double Mcrit = m_l0Ritter * PPOW(m_Mass, m_lambdaRitter) / (m_XRitter * Q_HYDROGEN_BURNING);    // Eq. 60 in Belczynski+ 2008. 6e18 is the energy yield of H burning in ergs/g.
        if (utils::Compare(Mdot, Mcrit) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES;                             // Flashes restrict accumulation
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION;                        // Material piles up on the WD. Leads to merger or CEE.
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
STELLAR_TYPE HeWD::EvolveToNextPhase() {

    STELLAR_TYPE stellarType;

    if (m_ShouldRejuvenate) {
        m_CoreMass   = m_Mass;
        m_Radius     = HeMS::CalculateRadiusAtZAMS_Static(m_CoreMass);
        m_Luminosity = HeMS::CalculateLuminosityAtZAMS_Static(m_CoreMass);
        m_Tau        = 0;
        stellarType  = STELLAR_TYPE::NAKED_HELIUM_STAR_MS; 
    }
    else {                                         
        stellarType  = ResolveSNIa();       // Currently, assume a Type Ia from a HeWD is the same as other WDs. May want to vary in the future
    }
    return stellarType;
}
