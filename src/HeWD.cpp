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

std::tuple<double,int> HeWD::DetermineAccretionRegime(bool p_HeRich, const double p_AccretedMass, const double p_Dt) {
    double LogMdot = log10(p_AccretedMass / p_Dt) - 6; // Logarithm of the accreted mass (M_sun/yr)
    if (p_HeRich) {
        if (LogMdot <= HELIUM_WHITE_DWARF_MCRIT) {
            return std::make_tuple(1.0, 7); // Could lead to Sub CH SN Ia
        } else {
            return std::make_tuple(1.0, 8); // Could lift degeneracy and evolve into He MS. Requires minimum mass ! on top of the shell size
        }
    } else {
        double Mcrit = log10(m_l0 * PPOW(m_Mass, m_lambda) / (m_X * 6e18));
        if (LogMdot <= Mcrit) {
            return std::make_tuple(0.0, 9); // Flashes restrict accumulation
        } else {
            return std::make_tuple(1.0, 10); //Material piles up on the WD leading to mass loss from the system. Leads to merger or CEE.
        }
    }
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

void HeWD::ResolveAccretionRegime(const int p_Regime, const double p_AccretedMass, const double p_Dt) {
    double Mdot = p_AccretedMass / (p_Dt * 1e6);
    double MassSubCh = -4e8 * Mdot + 1.34; // Minimum mass for Sub-Ch Mass detonation.
    double ShellCrit = -7.8e-4 * Mdot + 1.34; // Minimum shell mass of He for detonation. Should not be burnt, but not implemented this yet. Ruiter+ 2014.
    if (p_Regime == 7) {
        if (m_Mass >= MassSubCh) {
            m_SubChandrasekhar = true;
        }
    } else if (p_Regime == 8) {
        if ((m_Mass >= M_HE_BUR) && (m_HeShell >= ShellCrit) && (Mdot < 1.64e-6)) {
            m_Rejuvenate = true;
        } else if ((m_Mass >= M_HE_BUR) && (Mdot >= 1.64e-6)) {
            m_Rejuvenate = true;
        }
    }
}

/*
 * Allow the evolution towards an HeMS or SN.
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
 */

STELLAR_TYPE HeWD::EvolveToNextPhase() {
    if (m_Rejuvenate) {
        return STELLAR_TYPE::NAKED_HELIUM_STAR_MS;
    }
    else {
        m_Mass = m_Radius = m_Luminosity = m_Age = 0.0;
        return STELLAR_TYPE::MASSLESS_REMNANT;
    }
}



