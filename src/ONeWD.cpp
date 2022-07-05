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

std::tuple<double,ACCRETION_REGIME> ONeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_AccretedMass, const double p_Dt) {
        double logMdot = log10(p_AccretedMass / p_Dt) - 6; // Logarithm of the accreted mass (M_sun/yr)
        double fraction;
        ACCRETION_REGIME regime;
        std::tie(std::ignore, fraction) = CalculateWDMassAcceptanceRate(logMdot, p_HeRich); // Assume everything works the same way as a CO WD
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
