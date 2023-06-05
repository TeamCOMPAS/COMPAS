#include "WhiteDwarfs.h"
#include "NS.h"

/* Calculate eta_hydrogen from Claeys+ 2014, appendix B. This parameter depends 
 * on three regimes for the mass transfer rate, which here are distinguished by the 
 * thresholds logMdotUppH and logMdotLowH. In Claeys+ 2014, the mass transfer rate is
 * \dot{M}_{tr} and the thresholds are \dot{M}_{cr,H} and \dot{M}_{cr,H}/8, respectively. 
 *
 * However, we have used improved thresholds from Nomoto+ 2007, in which the 
 * lower boundary is \dot{M}_{stable} and the upper boundary is \dot{M}_{RG}. 
 * More precisely, we implemented quadratic fits to the values in Nomoto+ 2007,
 * table 5, as described in Rodriguez+ (in prep). 
 *
 * double CalculateEtaH(const double p_MassTransferRate)
 *
 * @param   [IN]    p_MassTransferRate     Mass transfer rate onto the WD surface (Msun/yr)
 * @return                                 Hydrogen accretion efficiency
 */
double WhiteDwarfs::CalculateEtaH(const double p_MassTransferRate) {

    double etaH = 0.0;                                      // default return value

    double logMassTransferRate = log10(p_MassTransferRate);
    double m_Mass_2            = m_Mass * m_Mass;

    // The following coefficients come from quadratic fits to Nomoto+ 2007 results (table 5) in Mass vs log10 Mdot space, to cover the low-mass end.
    double logMdotUppH = WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_0 + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_1 * m_Mass + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_2 * m_Mass_2; 
    double logMdotLowH = WD_LOG_MT_LIMIT_NOMOTO_STABLE_0   + WD_LOG_MT_LIMIT_NOMOTO_STABLE_1   * m_Mass + WD_LOG_MT_LIMIT_NOMOTO_STABLE_2   * m_Mass_2;
    
    if (utils::Compare(logMassTransferRate, logMdotUppH) >= 0) {
        etaH = PPOW(10, logMdotUppH - logMassTransferRate);
    } 
    else if (utils::Compare(logMassTransferRate, logMdotLowH) >= 0) {
        etaH = 1.0;
    } 

    return etaH;
}


/* Calculate eta_helium from Claeys+ 2014, appendix B. Similarly to CalculateEtaH
 * above, this parameter depends on four regimes for the mass transfer rate, distinguished
 * here by logMdotUppHe, logMdotMidHe, and logMdotLowHe. In Claeys+ 2014, these thresholds
 * are \dot{M}_{up}, \dot{M}_{cr,He}, and \dot{M}_{low}, respectively. 
 *
 * However, we have again updated the thresholds to those described in Piersanti+ 2014,
 * table A1. The thresholds here are named by the boundaries RG/SS, SS/MF, and SF/Dt, 
 * respectively (see text for details). Note that the different flashes regimes from 
 * Piersanti+ 2014 have been merged into one, i.e we omit the MF/SF boundary, and 
 * the accumulation regime has been change so we can get double detonations. Finally, 
 * eta_KH04 has also been updated with the accretion efficiency values from Piersanti+ 2014.
 *
 * double CalculateEtaHe(const double p_MassTransferRate)
 *
 * @param   [IN]    p_MassTransferRate     Mass transfer rate onto the WD surface (Msun/yr)
 * @return                                 Helium accretion efficiency
 */
double WhiteDwarfs::CalculateEtaHe(const double p_MassTransferRate) {

    double etaHe = 1.0;                                     // default return value - so we can have double detonations
    
    double logMassTransferRate = log10(p_MassTransferRate);

    // The following coefficients in massTransfer limits come from table A1 in Piersanti+ 2014.
    double logMdotUppHe = WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_0 + WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_1 * m_Mass;
    double logMdotMidHe = WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_0 + WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_1 * m_Mass;
    double logMdotLowHe = WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_0 + WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_1 * m_Mass;

    if (utils::Compare(logMassTransferRate, logMdotUppHe) >= 0) {
        etaHe = PPOW(10, logMdotUppHe - logMassTransferRate);
    } 
    else if (utils::Compare(logMassTransferRate, logMdotMidHe) >= 0) {  // JR: do we need this since it's the default?  Or may it change here?  Or just here for clarity?
        etaHe = 1.0;
    } 
    else if (utils::Compare(logMassTransferRate, logMdotLowHe) >= 0) {
        etaHe = CalculateEtaPTY(p_MassTransferRate);
    } 

    return etaHe;
}


/* Calculate accretion efficiency as indicated in Piersanti+ 2014, section A3. Their recipe works
 * for specific mass and Mdot values, so a better implementation requires interpolation and
 * extrapolation (specially towards the low-mass end). Right now, we just adopt a
 * piece-wise approach. Note that the authors also specify that this is based on the first
 * strong flash only, but we use it for all episodes.
 *
 * double CalculateEtaPTY(const double p_MassTransferRate)
 *
 * @param   [IN]    p_MassTransferRate     Mass transfer rate onto the WD surface (Msun/yr)
 * @return                                 Accretion efficiency during the first stron helium flash, Piersanti+ 2014
 */
double WhiteDwarfs::CalculateEtaPTY(const double p_MassTransferRate) {

    double etaPTY = 0.0;                        // default return value

    double massRate   = p_MassTransferRate;
    double massRate_2 = massRate * massRate;
    double massRate_3 = massRate_2 * massRate;

    // Limits on each conditional statement come from masses from each model in Piersanti+ 2014. The final etaPTY value is based on table A3.
    if (utils::Compare(m_Mass, 0.6) <= 0) {
        etaPTY = 6.0e-3   + 5.1e-2  * massRate + 8.3e-3 * massRate_2 - 3.317e-4 * massRate_3;
    } 
    else if  (utils::Compare(m_Mass, 0.7) <= 0) {
        etaPTY = -3.5e-2  + 7.5e-2  * massRate - 1.8e-3 * massRate_2 + 3.266e-5 * massRate_3;
    } 
    else if (utils::Compare(m_Mass, 0.81) <= 0) {
        etaPTY = 9.3e-2   + 1.8e-2  * massRate + 1.6e-3 * massRate_2 - 4.111e-5 * massRate_3;
    } 
    else if (utils::Compare(m_Mass, 0.92) <= 0) { 
        etaPTY = -7.59e-2 + 1.54e-2 * massRate + 4.0e-4 * massRate_2 - 5.905e-6 * massRate_3;
    } 
    else {
        etaPTY = -0.323   + 4.1e-2  * massRate - 7.0e-4 * massRate_2 + 4.733e-6 * massRate_3;
    }

    return etaPTY;
}


/*
 * Calculate the luminosity of a White Dwarf as it cools
 *
 * Hurley et al. 2000, eq 90
 *
 *
 * double CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity, const double p_BaryonNumber)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time since White Dwarf formation in Myr
 * @param   [IN]    P_Metallicity               Metallicity of White Dwarf
 * @param   [IN]    p_BaryonNumber              Baryon number - differs per White Dwarf type (HeWD, COWD, ONeWD)
 * @return                                      Luminosity of a White Dwarf in Lsol
 */
double WhiteDwarfs::CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time, const double p_Metallicity, const double p_BaryonNumber) {
    return (635.0 * p_Mass * PPOW(p_Metallicity, 0.4)) / PPOW(p_BaryonNumber * (p_Time + 0.1), 1.4);
}


/*
 * Calculate the radius of a white dwarf - good for all types of WD
 *
 * Hurley et al. 2000, eq 91 (from Tout et al. 1997)
 *
 *
 * double CalculateRadiusOnPhase_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius of a White Dwarf in Rsol (since WD is ~ Earth sized, expect answer around 0.009)
 */
double WhiteDwarfs::CalculateRadiusOnPhase_Static(const double p_Mass) {
    double MCH_Mass_one_third  = std::cbrt(MCH / p_Mass); 
    double MCH_Mass_two_thirds = MCH_Mass_one_third * MCH_Mass_one_third;
    return std::max(NEUTRON_STAR_RADIUS, 0.0115 * std::sqrt((MCH_Mass_two_thirds - 1.0 / MCH_Mass_two_thirds)));
}


/* 
 * Increase shell size after mass transfer episode. Hydrogen and helium shells are kept separately.
 * Only applies the full mass increase from accretion to one of the shells. Does not account for, e.g,
 * the H layer burning and building up the He layer, which may be desired in the future. - RTW 9/14/22
 *
 * void ResolveShellChange(const double p_AccretedMass)
 *
 * @param   [IN]    p_AccretedMass              Mass accreted
 */
void WhiteDwarfs::ResolveShellChange(const double p_AccretedMass) {
    
    switch (m_AccretionRegime) {

        case ACCRETION_REGIME::HELIUM_ACCUMULATION:
        case ACCRETION_REGIME::HELIUM_FLASHES:
        case ACCRETION_REGIME::HELIUM_STABLE_BURNING:
        case ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION:
	        m_HeShell += p_AccretedMass;
            break;

        case ACCRETION_REGIME::HYDROGEN_FLASHES:
        case ACCRETION_REGIME::HYDROGEN_STABLE_BURNING:
        case ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION:
	        m_HShell += p_AccretedMass;
            break;

        default:
            SHOW_WARN(ERROR::WARNING, "Accretion Regime not set for WD, no mass added to shell.");  // show warning 
    }
}


/*
 * Resolve Accretion-Induced Collapse of a WD
 *
 * Following Hurley et al. 2000, Section 6.2.1
 *
 * An AIC of a WD results in a NS, which we are 
 * here assuming to have a low mass equal to the ECSN
 * remnant NS mass, and no natal kick. 
 *
 * STELLAR_TYPE ResolveAIC() 
 *
 * @return                                      Stellar type of remnant (STELLAR_TYPE::NEUTRON_STAR if SN, otherwise current type)
 */
STELLAR_TYPE WhiteDwarfs::ResolveAIC() { 

    if (!IsSupernova()) return m_StellarType;                                           // shouldn't be here if no SN

    m_Mass                                = MECS_REM;                                   // defined in constants.h
    
    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;

    SetSNCurrentEvent(SN_EVENT::AIC);                                                   // AIC happening now
    SetSNPastEvent(SN_EVENT::AIC);                                                      // ... and will be a past event

    return STELLAR_TYPE::NEUTRON_STAR;
}


/*
 * Resolve Type 1a Supernova 
 *
 * A Type 1a SN results in a massless remnant.
 *
 * STELLAR_TYPE ResolveSNIa() 
 *
 * @return                                      Stellar type of remnant (STELLAR_TYPE::MASSLESS_REMNANT if SN, otherwise current type)
 */
STELLAR_TYPE WhiteDwarfs::ResolveSNIa() { 

    if (!IsSupernova()) return m_StellarType;                                           // shouldn't be here if no SN

    m_Mass       = 0.0;
    m_Radius     = 0.0;
    m_Luminosity = 0.0;
    m_Age        = 0.0;
    
    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;

    SetSNCurrentEvent(SN_EVENT::SNIA);                                                  // SN Type Ia happening now
    SetSNPastEvent(SN_EVENT::SNIA);                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}


/*
 * Resolve Double Detonation     
 *
 * A double detonation results in a massless remnant.
 *
 * STELLAR_TYPE ResolveHeSD() 
 *
 * @return                                      Stellar type of remnant (STELLAR_TYPE::MASSLESS_REMNANT if SN, otherwise current type)
 */
STELLAR_TYPE WhiteDwarfs::ResolveHeSD() { 

    if (!IsSupernova()) return m_StellarType;                                           // shouldn't be here if no SN

    m_Mass       = 0.0;
    m_Radius     = 0.0;
    m_Luminosity = 0.0;
    m_Age        = 0.0;
    
    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;

    SetSNCurrentEvent(SN_EVENT::HeSD);                                                  // SN Type Ia (HeSD) happening now
    SetSNPastEvent(SN_EVENT::HeSD);                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}
