#include "ONeWD.h"

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
 * @param   [IN]    p_AccretedMass       Total mass accreted in M_Sun
 * @param   [IN]    p_Dt                 Size of the timestep in Myr, assumed to be the duration of this particular mass transfer episode
 * @return                               Tuple containing fraction of mass that should be retained and accretion regime
 */

std::tuple<double,ACCRETION_REGIME> ONeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_AccretedMass, const double p_Dt) {
    double logMdot = log10(p_AccretedMass / (p_Dt * MYR_TO_YEAR)); // Logarithm of the accreted mass (M_sun/yr)
    double fraction;
    ACCRETION_REGIME regime;
    std::tie(std::ignore, fraction) = CalculateWDMassAcceptanceRate(logMdot, p_HeRich); // Assume everything works the same way as a CO WD
    if (p_HeRich) {
        // The following coefficients in massTransfer limits come from table A1 in Piersanti+ 2014.
            double massTransferCrit = MT_LIMIT_CRIT_PIERSANTI_0 + MT_LIMIT_CRIT_PIERSANTI_1 * m_Mass;
            double massTransferStable = MT_LIMIT_STABLE_PIERSANTI_0 + MT_LIMIT_STABLE_PIERSANTI_1 * m_Mass; // Piersanti+2014 has several Flashes regimes. Here we group them into one.
            double massTransferDetonation = MT_LIMIT_DET_PIERSANTI_0 + MT_LIMIT_DET_PIERSANTI_1 * m_Mass; // Critical value for double detonation regime in Piersanti+ 2014
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
        // The following coefficients in massTransfer limits come from quadratic fits to Nomoto+ 2007 results (table 5) in Mass vs log10 Mdot space, to cover the low-mass end.
            double massTransferCrit = MT_LIMIT_CRIT_NOMOTO_0 +  MT_LIMIT_CRIT_NOMOTO_1 * m_Mass +  MT_LIMIT_CRIT_NOMOTO_2 * m_Mass * m_Mass;
            double massTransferStable =  MT_LIMIT_STABLE_NOMOTO_0 +  MT_LIMIT_STABLE_NOMOTO_1 * m_Mass +  MT_LIMIT_STABLE_NOMOTO_2 * m_Mass * m_Mass;
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
