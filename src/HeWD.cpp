#include "HeWD.h"

/* Calculate:
 *
 *     (a) Mass fraction retained after accretion episodes
 *     (b) Accretion regime
 *
 *
 * For the helium WD case, the implementation is based on Belczynski+ 2008.
 *
 * The accretion regime is just an integer which follows:
 *
 * 0.- Helium Accumulation
 * 1.- Helium Flashes
 * 2.- Helium Stable Burning
 * 3.- Helium Optically-Thick Winds
 * 4.- Hydrogen Flashes
 * 5.- Hydrogen Stable Burning
 * 6.- Hydrogen Optically-Thick Winds
 *
 *
 * std::tuple<double,int> DetermineAccretionRegime(bool p_HeRich, const double p_AccretedMass, const double p_Dt)
 *
 * @param   [IN]    p_HeRich             Whether the accreted material is helium-rich or not
 * @param   [IN]    p_AccretedMass       Total mass accreted
 * @param   [IN]    p_Dt                 Size of the timestep, assumed to be the duration of this particular mass transfer episode
 * @return                               Tuple containing fraction of mass that should be retained and accretion regime
 */

std::tuple<double,ACCRETION_REGIME> HeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_AccretedMass, const double p_Dt) {
    double logMdot = log10(p_AccretedMass / p_Dt) - 6; // Logarithm of the accreted mass (M_sun/yr)
    double fraction = 1.0;
    ACCRETION_REGIME regime;
    if (p_HeRich) {
        if (utils::Compare(logMdot, HELIUM_WHITE_DWARF_MCRIT) < 1) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR; // Could lead to Sub CH SN Ia
        } else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION; // Could lift degeneracy and evolve into He MS. Requires minimum mass ! on top of the shell size
        }
    } else {
        double Mcrit = log10(m_l0Ritter * PPOW(m_Mass, m_lambdaRitter) / (m_XRitter * 6e18));
        if (utils::Compare(logMdot, Mcrit) < 1) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES; // Flashes restrict accumulation
            fraction = 0.0;
        } else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION; //Material piles up on the WD leading to mass loss from the system. Leads to merger or CEE.
        }
    }
    return std::make_tuple(fraction, regime);
}

/* Resolve what happens when the star goes through accretion regimes that might enable a change pf phase.
 * Flags are activated and the actual change happens elsewhere.
 *
 * void ResolveAccretionRegime(const int p_Regime, const double p_AccretedMass, const double p_Dt)
 *
 * @param   [IN]    p_Regime             Integer related to the current accretion regime
 * @param   [IN]    p_AccretedMass       Total mass accreted
 * @param   [IN]    p_Dt                 Size of the timestep, assumed to be the duration of this particular mass transfer episode
 */

void HeWD::ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_AccretedMass, const double p_Dt) {
    double massTransfer = p_AccretedMass / (p_Dt * 1e6);
    double massSubCh = -4e8 * massTransfer + 1.34; // Minimum mass for Sub-Ch Mass detonation.
    double shellCrit = -7.8e-4 * massTransfer + 1.34; // Minimum shell mass of He for detonation. Should not be burnt, but not implemented this yet. Ruiter+ 2014.
    if (p_Regime == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR) {
        if (utils::Compare(m_Mass, massSubCh) > -1 ) {
            m_SubChandrasekhar = true;
        }
    } else if (p_Regime == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION) {
        if ((utils::Compare(m_Mass, MASS_HELIUM_BURN) > -1) && (utils::Compare(m_HeShell, shellCrit) > -1) && (utils::Compare(massTransfer, 1.64e-6) < 0)) {
            m_Rejuvenate = true;
        } else if ((utils::Compare(m_Mass, MASS_HELIUM_BURN) > -1) && (utils::Compare(massTransfer, 1.64e-6) > -1)) {
            m_Rejuvenate = true;
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



