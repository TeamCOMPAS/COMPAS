#ifndef __BH_h__
#define __BH_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "NS.h"

class BaseStar;
class NS;

class BH: virtual public BaseStar, public NS {

public:

    BH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), NS(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    BH& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically
    static  DBL_DBL_DBL     CalculateCoreCollapseSNParams_Static(const double p_Mass);

    static  double          CalculateLuminosityOnPhase_Static()                             { return 1.0E-10; }                                                     // Hurley et al. 2000, eq 96

    static  double          CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass);

    static  double          CalculateRadiusOnPhase_Static(const double p_Mass)              { return 4.24E-6 * p_Mass; }                                            // Schwarzschild radius of Black Hole - Hurley et al. 2000, eq 94

    virtual void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                                       const bool   p_RecycledNS,
                                                       const double p_Stepsize,
                                                       const double p_MassGainPerTimeStep,
                                                       const double p_Epsilon)              { return BaseStar::UpdateMagneticFieldAndSpin(p_CommonEnvelope, 
                                                                                                                                          p_RecycledNS,
                                                                                                                                          p_Stepsize, 
                                                                                                                                          p_MassGainPerTimeStep, 
                                                                                                                                          p_Epsilon); }             // Use BaseStar   

                                                                                                                                          
protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::BLACK_HOLE;                                                                                                                   // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                                                                                // Set age appropriately
   }


    // member functions - alphabetically
            double          CalculateConvergedMassStepZetaNuclear()                         { return 0.0; }

            double          CalculateEddingtonCriticalRate()                                { return 1.5E-8 * (m_Radius * RSOL_TO_KM / 10.0) * MYR_TO_YEAR; }       // Sluys 2013 ("Binary Evolution in a Nutshell"), eq 70

            double          CalculateGyrationRadius() const                                 { return 0.0; }                                                         // No tidal coupling to a BH

            double          CalculateLuminosityOnPhase()                                    { return CalculateLuminosityOnPhase_Static(); }

            double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0)    { return (2.0 / 5.0) * m_Mass * m_Radius * m_Radius; }
            double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0)  { return CalculateMomentOfInertia(p_RemnantRadius * RSOL_TO_AU) * RSOL_TO_AU * RSOL_TO_AU; }

            double          CalculateRadiusOnPhase()                                        { return CalculateRadiusOnPhase_Static(m_Mass); }                       // Use class member variables - returns radius in Rsol

            bool            ShouldEvolveOnPhase()                                           { return true; }                                                        // Always
            bool            ShouldSkipPhase()                                               { return false; }                                                       // Don't skip
};

#endif // __BH_h__
