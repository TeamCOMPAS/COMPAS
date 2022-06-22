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
 * Note that we have merged the different flashes regimes from Piersanti+ 2014 into a single regime.
 *
 * std::tuple<double,int> DetermineAccretionRegime(bool p_HeRich, const double p_AccretedMass, const double p_Dt)
 *
 * @param   [IN]    p_HeRich             Whether the accreted material is helium-rich or not
 * @param   [IN]    p_AccretedMass       Total mass accreted
 * @param   [IN]    p_Dt                 Size of the timestep, assumed to be the duration of this particular mass transfer episode
 * @return                               Tuple containing fraction of mass that should be retained and accretion regime
 */

std::tuple<double,int> COWD::DetermineAccretionRegime(bool p_HeRich, const double p_AccretedMass, const double p_Dt) {
        double LogMdot = log10(p_AccretedMass / p_Dt) - 6; // Logarithm of the accreted mass (M_sun/yr)
        double Fraction;
        std::tie(std::ignore, Fraction) = CalculateWDMassAcceptanceRate(LogMdot, p_HeRich);

        if (p_HeRich) {
                double MdotCrit = -6.84 + 1.349 * m_Mass;
                double MdotStable = -8.115 + 2.29 * m_Mass; // Piersanti+2014 has several Flashes regimes. Here we group them into one.
                double MdotDet = -8.313 + 1.018 * m_Mass; // Critical value for double detonation regime in Piersanti+ 2014
                if (LogMdot < MdotStable) {
                    if (LogMdot > MdotDet) {
                        return std::make_tuple(Fraction, 1);
                    } else {
                        return std::make_tuple(Fraction, 0);
                        }
                } else if (LogMdot > MdotCrit) {
                    return std::make_tuple(Fraction, 3);
                } else {
                    return std::make_tuple(Fraction, 2);
                }
        } else {
                double MdotCrit = -0.98023471 * m_Mass * m_Mass + 2.88247131 * m_Mass - 8.33017155; // Quadratic fits to Nomoto+ 2007 (Mass vs log10 Mdot space), to cover the low-mass end.
                double MdotStable = -1.2137735 * m_Mass * m_Mass + 3.57319872 * m_Mass - 9.21757267;
                if (LogMdot < MdotStable) {
                    return std::make_tuple(Fraction, 4);
                } else if (LogMdot > MdotCrit) {
                    return std::make_tuple(Fraction, 6);
                } else {
                    return std::make_tuple(Fraction, 5);
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

void COWD::ResolveAccretionRegime(const int p_Regime, const double p_AccretedMass, const double p_Dt) {
    double LogMdot = log10(p_AccretedMass / p_Dt) - 6; // Logarithm of the accreted mass (M_sun/yr)
    if (p_Regime == 0) {
        if ((m_Mass >= M_DETONATION_CO) && (m_HeShell >= SHELL_CRIT)) {
            m_DoubleDetonation = true;
        }
    } else if (p_Regime == 2) {
        if ((LogMdot > MDOT_OFF_C) && (m_Mass > 1.3)) {
            m_OffCenterIgnition = true;
        }
    }
}



/*
 * Allow the evolution towards an ONe WD . From https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.1593W/abstract around
 * the end of section 3.2. Also, allows SN.
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
 */

STELLAR_TYPE COWD::EvolveToNextPhase() {
    if (m_OffCenterIgnition) {
        return STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    }
    else {
        m_Mass = m_Radius = m_Luminosity = m_Age = 0.0;
        return STELLAR_TYPE::MASSLESS_REMNANT;
    }
}
