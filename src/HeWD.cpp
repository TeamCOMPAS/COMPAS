#include "HeWD.h"

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

std::tuple<double,ACCRETION_REGIME> HeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorThermalMassLossRate) {
    double logMdot = log10(p_DonorThermalMassLossRate / MYR_TO_YEAR); // Logarithm of the accreted mass (M_sun/yr)
    double fraction = 1.0;
    ACCRETION_REGIME regime;
    if (p_HeRich) {
        if (utils::Compare(logMdot, HELIUM_WHITE_DWARF_MCRIT) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR; // Could lead to Sub CH SN Ia
        } else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION; // Could lift degeneracy and evolve into He MS. Requires minimum mass ! on top of the shell size
        }
    } else {
        double Mcrit = log10(m_l0Ritter * PPOW(m_Mass, m_lambdaRitter) / (m_XRitter * 6e18)); // Eq. 60 in Belczynski+ 2008. 6e18 is the energy yield of H burning in ergs/g.
        if (utils::Compare(logMdot, Mcrit) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES; // Flashes restrict accumulation
            fraction = 0.0;
        } else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION; //Material piles up on the WD. Leads to merger or CEE.
        }
    }
    return std::make_tuple(fraction, regime);
}

/* Resolve what happens when the star goes through accretion regimes that might enable a change of phase.
 * Flags are activated and the actual change happens elsewhere.
 *
 * void ResolveAccretionRegime(const int p_Regime, const double p_AccretedMass, const double p_Dt)
 *
 * @param   [IN]    p_Regime                        ACCRETION_REGIME value
 * @param   [IN]    p_DonorThermalMassLossRate      Donor thermal mass loss rate, in units of Msol / Myr
 */

void HeWD::ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_DonorThermalMassLossRate) {
    double Mdot = p_DonorThermalMassLossRate / MYR_TO_YEAR;
    double massSubCh = -4e8 * Mdot + 1.34; // Minimum mass for Sub-Ch Mass detonation. Eq 62, Belczynski+ 2008.
    double shellCrit = -7.8e-4 * Mdot + 1.34; // Minimum shell mass of He for detonation. Eq 61, Belczynski+ 2008. This helium should not be burnt, but not implemented this yet. Ruiter+ 2014.
    if (p_Regime == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR) {
        if (utils::Compare(m_Mass, massSubCh) >= 0 ) {
            m_SubChandrasekhar = true;
        }
    } else if (p_Regime == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION) {
        if (utils::Compare(m_Mass, MASS_HELIUM_BURN) >= 0) {
            if (utils::Compare(Mdot, 1.64e-6) < 0) { // Accretion limit from eq 61, Belczynski+ 2008.
                if (utils::Compare(m_HeShell, shellCrit) >= 0) {
                    m_Rejuvenate = true;
                }
            } else {
                m_Rejuvenate = true;
            }
        }
    }
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

// RTW: is this different to COWD.cpp?


