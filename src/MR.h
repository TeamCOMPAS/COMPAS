#ifndef __MR_h__
#define __MR_h__

#include "constants.h"
#include "typedefs.h"
#include "limits.h"

#include "Remnants.h"

class BaseStar;
class Remnants;

class MR: virtual public BaseStar, public Remnants {

public:

    MR() { m_StellarType = STELLAR_TYPE::MASSLESS_REMNANT; };
    
    MR(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::MASSLESS_REMNANT;                                                     // Set stellar type
        if (p_Initialise) Initialise();                                                                     // Initialise if required
    }

    MR* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MR* clone = new MR(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static MR* Clone(MR& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MR* clone = new MR(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
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
    
     double     ChooseTimestep(const double p_Time) const { return std::numeric_limits<double>::max(); }                                                              // Can take arbitrarily long time steps for massless remnants -- nothing is happening

     void       SetPulsarParameters() const { }                                                             // NO-OP

     bool       ShouldEvolveOnPhase() const             { return true; }                                    // Always
     bool       ShouldSkipPhase() const                 { return false; }                                   // Don't skip
};

#endif // __MR_h__
