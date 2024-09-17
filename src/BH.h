#ifndef __BH_h__
#define __BH_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"

class BaseStar;
class NS;

class BH: virtual public BaseStar, public Remnants {
    
public:
    
    BH() { m_StellarType = STELLAR_TYPE::BLACK_HOLE; };
    
    BH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::BLACK_HOLE;                                                                                                       // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                 // Initialise if required
    }
    
    BH* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        BH* clone = new BH(*this, p_Initialise);
        clone->SetPersistence(p_Persistence);
        return clone;
    }
    
    static BH* Clone(BH& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        BH* clone = new BH(p_Star, p_Initialise);
        clone->SetPersistence(p_Persistence);
        return clone;
    }
    
    // member functions - alphabetically
    static  DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass);
    static  double      CalculateLuminosityOnPhase_Static()                         { return 1.0E-10; }                                                 // Hurley et al. 2000, eq 96
    static  double      CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass);
    static  double      CalculateRadiusOnPhase_Static(const double p_Mass)          { return 4.24E-6 * p_Mass; }                                        // Schwarzschild radius of Black Hole - Hurley et al. 2000, eq 94
     
    static  double      ReweightSupernovaKickByMass_Static(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass);
   
    
protected:
    
    void Initialise() {
        CalculateTimescales();                                                                                                                          // Initialise timescales
        m_Age = 0.0;                                                                                                                                    // Set age appropriately
    }
    
    
    // member functions - alphabetically
    double  CalculateConvergedMassStepZetaNuclear() const                           { return 0.0; }
    double  CalculateEddingtonCriticalRate() const                                  { return 2.6E-8 * m_Mass * MYR_TO_YEAR; }                           // E.g., Marchant+, 2017, Eq. 3, assuming accretion efficiency of 10%
    double  CalculateLuminosityOnPhase() const                                      { return CalculateLuminosityOnPhase_Static(); }
    double  CalculateMassLossRate()                                                 { return 0.0; }                                                     // Ensure that BHs don't lose mass in winds
    double  CalculateMomentOfInertia() const                                        { return (2.0 / 5.0) * m_Mass * m_Radius * m_Radius; }
    double  CalculateRadiusOnPhase() const                                          { return CalculateRadiusOnPhase_Static(m_Mass); }                   // Use class member variables - returns radius in Rsol
};

#endif // __BH_h__
