#include "ONeWD.h"
#include "NS.h"

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

std::tuple<double,ACCRETION_REGIME> ONeWD::DetermineAccretionRegime(const bool p_HeRich, const double p_DonorThermalMassLossRate) {
    double logMdot = log10(p_DonorThermalMassLossRate / MYR_TO_YEAR); // Logarithm of the accreted mass (M_sun/yr)
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

/*
 * Resolve Accretion-Induced Collapse of a WD
 *
 * Following Hurley et al. 2000, Section 6.2.1
 *
 * An AIC of a ONeWD results in a NS, which we are 
 * here assuming to have a low mass equal to the ECSN
 * remnant NS mass, and no natal kick. 
 *
 * STELLAR_TYPE ONeWD::ResolveAIC() 
 *
 * @return                                      Stellar type of remnant (STELLAR_TYPE::NEUTRON_STAR if SN, otherwise current type)
 */

STELLAR_TYPE ONeWD::ResolveAIC() { 

    if (!IsSupernova()) return m_StellarType;                                           // shouldn't be here if no SN

    m_Mass       = MECS_REM;                                                            // defined in constants.h
    m_Radius     = NS::CalculateRadiusOnPhase_Static(m_Mass);                           // neutronStarRadius in Rsol
    m_Luminosity = NS::CalculateLuminosityOnPhase_Static(m_Mass, m_Age);
    
    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;

    SetSNCurrentEvent(SN_EVENT::AIC);                                                  // AIC happening now
    SetSNPastEvent(SN_EVENT::AIC);                                                     // ... and will be a past event

    return STELLAR_TYPE::NEUTRON_STAR;
}
