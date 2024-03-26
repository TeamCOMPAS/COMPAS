#ifndef __NS_h__
#define __NS_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"
#include "BH.h"


class BaseStar;
class Remnants;

class NS: virtual public BaseStar, public Remnants {

public:

    NS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    NS& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically
    static  DBL_DBL_DBL     CalculateCoreCollapseSNParams_Static(const double p_Mass);

    static  double          CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time);

    static  double          CalculatePulsarBirthSpinPeriod();

    static  double          CalculateRadiusOnPhaseInKM_Static(const double p_Mass);                                                                                 // Radius on phase in km
    static  double          CalculateRadiusOnPhase_Static(const double p_Mass)      { return CalculateRadiusOnPhaseInKM_Static(p_Mass) * KM_TO_RSOL; }              // Radius on phase in Rsol

    static  double          CalculateRemnantMass_Static(const double p_COCoreMass)  { return 1.17 + (0.09 * p_COCoreMass); }                                        // Hurley et al., eq 92


protected:
    
    void Initialise() {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                                 // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        
        //Set internal properties to zero to avoid meaningless values
        m_Age        = 0.0;
        m_COCoreMass = 0.0;
        m_HeCoreMass = 0.0;
        m_CoreMass   = 0.0;
        m_Mass0      = 0.0;
        
        m_Radius     = CalculateRadiusOnPhase_Static(m_Mass);                                                                                                       // Set the NS radius, in Rsol
        m_Luminosity = CalculateLuminosityOnPhase_Static(m_Mass, m_Age);                                                                                            // Set the NS luminosity

        CalculateAndSetPulsarParameters();
    }

    void FastForward() {                                                                                                                                                        // Set stellar attributes for stars initialized to this stellar type

        m_Radius                                   = CalculateRadiusOnPhase();
        m_Luminosity                               = CalculateLuminosityOnPhase();
    
        m_InitialLuminosity                        = m_Luminosity;
        m_InitialRadius                            = m_Radius;
        m_InitialStellarType                       = m_StellarType;
        m_StellarTypePrev                          = m_StellarType;
    }

    double m_AngularMomentum_CGS;                                                                                                                                   // Current angular momentum in CGS - only required in NS class
    double m_MomentOfInertia_CGS;                                                                                                                                   // MoI in CGS - only required in NS class


    // member functions - alphabetically
    void            CalculateAndSetPulsarParameters();

    double          CalculateLuminosityOnPhase() const                                      { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age); }            // Use class member variables

    double          CalculateMassLossRate()                                                 { return 0.0; }                                                         // Ensure that NSs don't lose mass in winds
    
    double          CalculateMomentOfInertiaCGS() const;                                                                                                            // MoI in CGS
    double          CalculateMomentOfInertia() const                                        { return CalculateMomentOfInertiaCGS() / MSOL_TO_G / RSOL_TO_CM / RSOL_TO_CM; } // MoI (default is solar units)

    double          CalculatePulsarBirthMagneticField();

    double          CalculateRadiusOnPhase() const                                          { return CalculateRadiusOnPhase_Static(m_Mass); }                       // Use class member variables - returns radius in Rsol

    double          CalculateSpinDownRate(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius) const;

    double          ChooseTimestep(const double p_Time) const;

    STELLAR_TYPE    EvolveToNextPhase()                                                     { return STELLAR_TYPE::BLACK_HOLE; }
    
    bool            ShouldEvolveOnPhase() const                                             { return (m_Mass <= OPTIONS->MaximumNeutronStarMass()); }               // Evolve as a neutron star unless mass > maximum neutron star mass (e.g. through accretion)
    void            SpinDownIsolatedPulsar(const double p_Stepsize);
    void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                               const bool   p_RecycledNS,
                                               const double p_Stepsize,
                                               const double p_MassGainPerTimeStep,
                                               const double p_Epsilon);

};

#endif // __NS_h__