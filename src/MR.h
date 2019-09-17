#ifndef __MR_h__
#define __MR_h__

#include "constants.h"
#include "typedefs.h"

#include "BH.h"

class BaseStar;
class BH;

class MR: virtual public BaseStar, public BH {

public:

    MR(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), BH(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    MR& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::MASSLESS_REMNANT;                                             // Set stellar type

        // ensure it's a massless remnant...
        m_Age         = 0.0;
        m_Mass        = 0.0;
        m_COCoreMass  = 0.0;
        m_HeCoreMass  = 0.0;
        m_CoreMass    = 0.0;
        m_Mass0       = 0.0;
        m_EnvMass     = 0.0;
        m_CoreRadius  = 0.0;
        m_Luminosity  = 0.0;
        m_Radius      = 0.0;
        m_Temperature = 0.0;
    }


    // member functions
   	 double CalculateMomentOfInertia(const double p_RemnantRadius = 0.0)        { return 0.0; }     // no moment of inertia for massless remnants - use 0.0
   	 double CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0)      { return 0.0; }     // no moment of inertia for massless remnants - use 0.0

     void   SetPulsarParameters() { }                                                               // NO-OP

     bool   ShouldEvolveOnPhase()                                               { return true; }    // Always
     bool   ShouldSkipPhase()                                                   { return false; }   // Don't skip
};

#endif // __MR_h__
