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

    BH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    BH& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically
    static  DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass);

    static  double      CalculateLuminosityOnPhase_Static()                         { return 1.0E-10; }                                                 // Hurley et al. 2000, eq 96

    static  double      CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass);

    static  double      CalculateRadiusOnPhase_Static(const double p_Mass)          { return 4.24E-6 * p_Mass; }                                        // Schwarzschild radius of Black Hole - Hurley et al. 2000, eq 94

                                                                                                                                          
protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::BLACK_HOLE;                                                                                                       // Set stellar type
        CalculateTimescales();                                                                                                                          // Initialise timescales
        m_Age = 0.0;                                                                                                                                    // Set age appropriately
   }


    // member functions - alphabetically
    double  CalculateConvergedMassStepZetaNuclear() const                           { return 0.0; }

    double  CalculateEddingtonCriticalRate() const                                  { return 2.6E-8 * m_Mass * MYR_TO_YEAR; }   // E.g., Marchant+, 2017, Eq. 3, assuming accretion efficiency of 10%

    double  CalculateGyrationRadius() const                                         { return 0.0; }                                                     // No tidal coupling to a BH

    double  CalculateLuminosityOnPhase() const                                      { return CalculateLuminosityOnPhase_Static(); }

    double  CalculateMomentOfInertia(const double p_RemnantRadius = 0.0) const      { return (2.0 / 5.0) * m_Mass * m_Radius * m_Radius; }
    double  CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0) const    { return CalculateMomentOfInertia(p_RemnantRadius * RSOL_TO_AU) * RSOL_TO_AU * RSOL_TO_AU; }

    double  CalculateRadiusOnPhase() const                                          { return CalculateRadiusOnPhase_Static(m_Mass); }                   // Use class member variables - returns radius in Rsol
};

#endif // __BH_h__
