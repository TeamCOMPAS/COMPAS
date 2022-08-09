#include "ONeWD.h"
#include "NS.h"

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
