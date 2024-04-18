#ifndef __MR_h__
#define __MR_h__

#include "constants.h"
#include "typedefs.h"

#include "Remnants.h"

class BaseStar;
class Remnants;

class MR: virtual public BaseStar, public Remnants {

public:

    MR(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    MR& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::MASSLESS_REMNANT;                                                     // Set stellar type

        // ensure it's a massless remnant...
        m_Age         = 0.0;
        m_Mass        = 0.0;
        m_COCoreMass  = 0.0;
        m_HeCoreMass  = 0.0;
        m_CoreMass    = 0.0;
        m_Mass0       = 0.0;
        m_Luminosity  = 0.0;
        m_Radius      = 0.0;
        m_Temperature = 0.0;
    }


    // member functions
   	 double     CalculateMomentOfInertia() const        { return 0.0; }                                     // No moment of inertia for massless remnants - use 0.0
   	 double     CalculateMomentOfInertiaAU() const      { return 0.0; }                                     // No moment of inertia for massless remnants - use 0.0

     void       SetPulsarParameters() const { }                                                             // NO-OP

     bool       ShouldEvolveOnPhase() const             { return true; }                                    // Always
     bool       ShouldSkipPhase() const                 { return false; }                                   // Don't skip
};

#endif // __MR_h__
