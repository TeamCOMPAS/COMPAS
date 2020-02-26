#ifndef __MS_gt_07_h__
#define __MS_gt_07_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_gt_07: virtual public BaseStar, public MainSequence {

public:

    MS_gt_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {
        if (p_Initialise) Initialise();
    }

    MS_gt_07& operator = (const BaseStar &baseStar) {
        static_cast<BaseStar&>(*this) = baseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::MS_GT_07;                                                                                                         // Set stellar type
        CalculateTimescales();                                                                                                                          // Initialise timescales
        // JR: Age for MS_GT_07 is carried over from CH stars switching to MS after spinning down, so not set to 0.0 here
    }


    // member functions - alphabetically
    double   CalculateMassLossRateHurley();
    double   CalculateMassTransferRejuvenationFactor();

    ENVELOPE DetermineEnvelopeType()                        { return ENVELOPE::RADIATIVE; }                                                             // Always RADIATIVE
    ENVELOPE DetermineEnvelopeTypeHurley2002()              { return utils::Compare(m_Mass, 1.25) < 0 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE; }   // Sometimes CONVECTIVE... JR: todo: why is this different?

    bool     IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate);
};

#endif // __MS_gt_07_h__
