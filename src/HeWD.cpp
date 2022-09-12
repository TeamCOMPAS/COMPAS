#include "HeWD.h"

// RTW: TODO Header
/* Calculate:
 *
 *     (a) Mass fraction retained after accretion episodes
 *     (b) Accretion regime
 *
 *
 * For the helium WD case, the implementation is based on Belczynski+ 2008.
 *
 * The accretion regime is one of the following:
 *
 * Helium accretion that could lead to Sub-Chandrasekhar SN
 * Helium accretion that could lead to helium ignition and rejuvenation
 * Hydrogen Flashes
 * Hydrogen Accumulation
 *
 *
 * std::tuple<double,int> DetermineAccretionRegime(bool p_HeRich, const double p_DonorThermalMassLossRate)
 *
 * @param   [IN]    p_HeRich                        Whether the accreted material is helium-rich or not
 * @param   [IN]    p_DonorThermalMassLossRate      Donor thermal mass loss rate, in units of Msol / Myr
 * @return                                          Tuple containing fraction of mass that should be retained and accretion regime
 */

ACCRETION_REGIME HeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorThermalMassLossRate) {
    double Mdot = p_DonorThermalMassLossRate / MYR_TO_YEAR; // Accreted mass rate (M_sun/yr)
    double logMdot = log10(Mdot);                           // Logarithm of the accreted mass rate (M_sun/yr)
    double fraction = 1.0;      // RTW: What should the fraction be here? Is this missing a CalculateMassAcceptance? A fraction of 1 seems suspiciously high
    /* NRS: the literature has not explored what happenes to HeWDs accreting as much as the COWD case. 
     * Using the COWD prescription with ONe WDs is "okay" as they cover a similar range in mass, 
     * but all HeWDs are low mass and the only prescriptions I could find were from StarTrack, which basically motivate all of the following.
     * One of my side-projects is to create a better prescription, but for now this should be used in the same way as the COWD function (I see that you are not returning the fraction value)
     */
    ACCRETION_REGIME regime;
    if (p_HeRich) {
        if (utils::Compare(logMdot, HELIUM_WHITE_DWARF_MCRIT) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR; // Could lead to Sub CH SN Ia
            double massSubCh = -4e8 * Mdot + 1.34; // Minimum mass for Sub-Ch Mass detonation. Eq 62, Belczynski+ 2008.
            if (utils::Compare(m_Mass, massSubCh) >= 0 ) {
                m_SubChandrasekhar = true;
            }
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION; // Could lift degeneracy and evolve into He MS. Requires minimum mass ! on top of the shell size
            if (utils::Compare(m_Mass, MASS_HELIUM_BURN) >= 0) {
                if (utils::Compare(Mdot, 1.64e-6) < 0) { // Accretion limit from eq 61, Belczynski+ 2008.
                    double shellCrit = -7.8e-4 * Mdot + 1.34; // Minimum shell mass of He for detonation. Eq 61, Belczynski+ 2008. This helium should not be burnt, but not implemented this yet. Ruiter+ 2014.
                    if (utils::Compare(m_HeShell, shellCrit) >= 0) {
                        m_Rejuvenate = true;
                    }
                } 
                else {
                    m_Rejuvenate = true;
                }
            }
        }
    } 
    else {
        double Mcrit = log10(m_l0Ritter * PPOW(m_Mass, m_lambdaRitter) / (m_XRitter * 6e18)); // Eq. 60 in Belczynski+ 2008. 6e18 is the energy yield of H burning in ergs/g.
        if (utils::Compare(logMdot, Mcrit) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES; // Flashes restrict accumulation
            fraction = 0.0;
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
    if (m_SubChandrasekhar || m_Rejuvenate) {
        return false;
    } else {
        return true;
    }
}

/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar type of the upcoming stage.
 */

STELLAR_TYPE HeWD::EvolveToNextPhase() {
    if (m_Rejuvenate) {
        return STELLAR_TYPE::NAKED_HELIUM_STAR_MS;
    }
    else {
        m_Mass       = 0.0;
        m_Radius     = 0.0;
        m_Luminosity = 0.0;
        m_Age        = 0.0;
        return STELLAR_TYPE::MASSLESS_REMNANT;
    }
}


