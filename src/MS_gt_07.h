#ifndef __MS_gt_07_h__
#define __MS_gt_07_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_gt_07: virtual public BaseStar, public MainSequence {

public:

    MS_gt_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::MS_GT_07;                                                                                                         // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                 // Initialise if required
    }

    MS_gt_07* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_gt_07* clone = new MS_gt_07(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static MS_gt_07* Clone(MS_gt_07& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_gt_07* clone = new MS_gt_07(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                          // Initialise timescales
        // JR: Age for MS_GT_07 is carried over from CH stars switching to MS after spinning down, so not set to 0.0 here
    }


    // member functions - alphabetically

    double      CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const ;
    double      CalculateCriticalMassRatioHurleyHjellmingWebbink() const                      { return 0.33; }                                          // As coded in BSE. Using the inverse owing to how qCrit is defined in COMPAS. See Hurley et al. 2002 sect. 2.6.1 for additional details.
    double      CalculateMassLossRateHurley();
    double      CalculateMassTransferRejuvenationFactor() const;

    ENVELOPE    DetermineEnvelopeType() const;
};

#endif // __MS_gt_07_h__
