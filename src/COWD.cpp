#include "COWD.h"

/* Calculate:
 *
 *     (a) Mass fraction retained after accretion episodes
 *     (b) Accretion regime
 *
 *
 * Mass retention is based on appendix B of Claeys+ 2014, but critical accretion rates have been updated using
 * Nomoto+ 2007 for hydrogen and Piersanti+ 2014 for helium. Details in CalculateWDMassAcceptanceRate().
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
 * std::tuple<double,int> DetermineAccretionRegime(bool p_HeRich, const double p_AccretedMass, const double p_Dt)
 *
 * @param   [IN]    p_HeRich             Whether the accreted material is helium-rich or not
 * @param   [IN]    p_AccretedMass       Total mass accreted
 * @param   [IN]    p_Dt                 Size of the timestep, assumed to be the duration of this particular mass transfer episode
 * @return                               Tuple containing fraction of mass that should be retained and accretion regime
 */

std::tuple<double,ACCRETION_REGIME> COWD::DetermineAccretionRegime(const bool p_HeRich, const double p_AccretedMass, const double p_Dt) {
        double logMdot = log10(p_AccretedMass / p_Dt) - 6; // Logarithm of the accreted mass (M_sun/yr)
        double fraction;
        ACCRETION_REGIME regime;
        std::tie(std::ignore, fraction) = CalculateWDMassAcceptanceRate(logMdot, p_HeRich);

        if (p_HeRich) {
                double massTransferCrit = -6.84 + 1.349 * m_Mass;
                double massTransferStable = -8.115 + 2.29 * m_Mass; // Piersanti+2014 has several Flashes regimes. Here we group them into one.
                double massTransferDetonation = -8.313 + 1.018 * m_Mass; // Critical value for double detonation regime in Piersanti+ 2014
                if (utils::Compare(logMdot, massTransferStable) < 0) {
                    if (utils::Compare(logMdot, massTransferDetonation) > 0) {
                        regime = ACCRETION_REGIME::HELIUM_FLASHES;
                    } else {
                        regime = ACCRETION_REGIME::HELIUM_ACCUMULATION;
                        }
                } else if (utils::Compare(logMdot, massTransferCrit) > 0) {
                    regime = ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS;
                } else {
                    regime = ACCRETION_REGIME::HELIUM_STABLE_BURNING;
                }
        } else {
                double massTransferCrit = -0.98023471 * m_Mass * m_Mass + 2.88247131 * m_Mass - 8.33017155; // Quadratic fits to Nomoto+ 2007 (Mass vs log10 Mdot space), to cover the low-mass end.
                double massTransferStable = -1.2137735 * m_Mass * m_Mass + 3.57319872 * m_Mass - 9.21757267;
                if (utils::Compare(logMdot, massTransferStable) < 0) {
                    regime = ACCRETION_REGIME::HYDROGEN_FLASHES;
                } else if (utils::Compare(logMdot, massTransferCrit) > 0) {
                    regime = ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS;
                } else {
                    regime = ACCRETION_REGIME::HYDROGEN_STABLE_BURNING;
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

void COWD::ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_AccretedMass, const double p_Dt) {
    double logMdot = log10(p_AccretedMass / p_Dt) - 6; // Logarithm of the accreted mass (M_sun/yr)
    if (p_Regime == ACCRETION_REGIME::HELIUM_ACCUMULATION) {
        if ((utils::Compare(m_Mass, MASS_DOUBLE_DETONATION_CO) > -1) && (utils::Compare(m_HeShell, SHELL_CRIT) > -1)) {
            m_DoubleDetonation = true;
        }
    } else if (p_Regime == ACCRETION_REGIME::HELIUM_STABLE_BURNING) {
        if ((utils::Compare(logMdot, MDOT_OFF_C) > 0) && (utils::Compare(m_Mass, 1.3) > 0)) {
            m_OffCenterIgnition = true;
        }
    }
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
