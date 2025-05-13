#ifndef __MS_lte_07_h__
#define __MS_lte_07_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_lte_07: virtual public BaseStar, public MainSequence {

public:

    MS_lte_07() { m_StellarType = STELLAR_TYPE::MS_LTE_07; };
    
    MS_lte_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::MS_LTE_07;                                                                                // Set stellar type
        if (p_Initialise) Initialise();                                                                                         // Initialise if required
    }

    MS_lte_07* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_lte_07* clone = new MS_lte_07(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static MS_lte_07* Clone(MS_lte_07& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_lte_07* clone = new MS_lte_07(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                  // Initialise timescales
        m_Age = 0.0;                                                                                                            // Set age appropriately
    }


    // member functions - alphabetically
    double      CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const ;
    double      CalculateCriticalMassRatioHurleyHjellmingWebbink() const                    { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_LTE_07; }
    double      CalculateMassTransferRejuvenationFactor()                                   { return 1.0; }

    ENVELOPE    DetermineEnvelopeType() const                                               { return ENVELOPE::CONVECTIVE; }    // Always CONVECTIVE
};

#endif // __MS_lte_07_h__
