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

    static  double          CalculatePulsarBirthSpinPeriod_Static();

    static  double          CalculateRadiusOnPhaseInKM_Static(const double p_Mass);                                                                                 // Radius on phase in km
    static  double          CalculateRadiusOnPhase_Static(const double p_Mass)      { return CalculateRadiusOnPhaseInKM_Static(p_Mass) * KM_TO_RSOL; }              // Radius on phase in Rsol

    static  double          CalculateRemnantMass_Static(const double p_COCoreMass)  { return 1.17 + (0.09 * p_COCoreMass); }                                        // Hurley et al., eq 92


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                                 // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        
        //Set internal properties to zero to avoid meaningless values
        m_Age = 0.0;
        m_COCoreMass  = 0.0;
        m_HeCoreMass  = 0.0;
        m_CoreMass    = 0.0;
        m_Mass0       = 0.0;
        
        m_Radius      = NS::CalculateRadiusOnPhase_Static(m_Mass);                                                                                                  // Set the NS radius, in Rsol
        m_Luminosity  = NS::CalculateLuminosityOnPhase_Static(m_Mass, m_Age);                                                                                       // Set the NS luminosity

        CalculateAndSetPulsarParameters();
    }

    double m_MomentOfInertia;                                                                                                                                       // in CGS g cm^2
    double m_AngularMomentum;                                                                                                                                       // Current angular momentum in (Msol AU^2 yr-1)
    
    // member functions - alphabetically
            void            CalculateAndSetPulsarParameters();

            double          CalculateLuminosityOnPhase() const                      { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age); }                    // Use class member variables

            double          CalculateMassLossRate()                                 { return 0.0; }                                                                 // Ensure that NSs don't lose mass in winds
    
    static  double          CalculateMomentOfInertia_Static(const double p_Mass, const double p_Radius);

    static  double          CalculatePulsarBirthMagneticField_Static();

            double          CalculateRadiusOnPhase() const                          { return CalculateRadiusOnPhase_Static(m_Mass); }                               // Use class member variables - returns radius in Rsol

    static  double          CalculateSpinDownRate_Static(const double p_Omega,
                                                         const double p_MomentOfInteria,
                                                         const double p_MagField,
                                                         const double p_Radius,
                                                         double const p_Alpha);
            STELLAR_TYPE    EvolveToNextPhase()                                     { return STELLAR_TYPE::BLACK_HOLE; }
    
            bool            ShouldEvolveOnPhase() const                             { return (m_Mass <= OPTIONS->MaximumNeutronStarMass()); }                       // Evolve as a neutron star unless mass > maximum neutron star mass (e.g. through accretion)

            void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                                       const bool   p_RecycledNS,
                                                       const double p_Stepsize,
                                                       const double p_MassGainPerTimeStep,
                                                       const double p_Epsilon);
};

#endif // __NS_h__
