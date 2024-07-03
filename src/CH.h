#ifndef __CH_h__
#define __CH_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MS_gt_07.h"

class BaseStar;
class MS_gt_07;

class CH: virtual public BaseStar, public MS_gt_07 {

public:

    CH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MS_gt_07(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS;                                                                                                                       // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                             // Initialise if required
    }

    CH* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CH* clone = new CH(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static CH* Clone(CH& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CH* clone = new CH(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                                                                                                // Set age appropriately
        m_CHE = true;                                                                                                                                                               // initially for CH stars                                                                                                                                                            // Set age appropriately
    }

    // member functions
    double          CalculateRadiusOnPhase() const      { return m_InitialRadius; }                                                                                                         // Constant from birth
    double          CalculateRadiusAtPhaseEnd() const   { return CalculateRadiusOnPhase(); }                                                                                        // Same as on phase

    STELLAR_TYPE    EvolveToNextPhase();

    bool            ShouldEvolveOnPhase() const         { return m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)] && (OPTIONS->OptimisticCHE() || m_Omega >= m_OmegaCHE); }   // Evolve on CHE phase if age in MS timescale and spinning at least as fast as CHE threshold

};

#endif // __CH_h__
