#ifndef __CHE_h__
#define __CHE_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "BaseStar.h"
#include "MS_gt_07.h"

class BaseStar;
class MS_gt_07;

class CHE: virtual public BaseStar, public MS_gt_07 {

public:

    CHE(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MS_gt_07(p_BaseStar) {
        if (p_Initialise) Initialise();
    }

    CHE& operator = (const BaseStar &baseStar) {
        static_cast<BaseStar&>(*this) = baseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS;                                                                                                                   // Set stellar type
        CalculateTimescales();                                                                                                                                                  // Initialise timescales
        m_Age = 0.0;                                                                                                                                                            // Set age appropriately
    }

    // member functions
    double          CalculateRadiusOnPhase()        { return m_RZAMS; }                                                                                                         // Constant from birth
    double          CalculateRadiusAtPhaseEnd()     { return CalculateRadiusOnPhase(); }                                                                                        // Same as on phase

    STELLAR_TYPE    EvolveToNextPhase();

    bool            ShouldEvolveOnPhase()           { return m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)] && (OPTIONS->OptimisticCHE() || m_Omega >= m_OmegaCHE); }   // Evolve on CHE phase if age in MS timescale and spinning at least as fast as CHE threshold

};

#endif // __CHE_h__
