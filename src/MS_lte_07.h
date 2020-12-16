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

    MS_lte_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {
        if (p_Initialise) Initialise();
    }

    MS_lte_07& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::MS_LTE_07;                                                    // Set stellar type
        CalculateTimescales();                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                // Set age appropriately
    }


    // member functions - alphabetically
    double      CalculateMassTransferRejuvenationFactor() const;

    ENVELOPE    DetermineEnvelopeType() const { return ENVELOPE::CONVECTIVE; }                      // Always CONVECTIVE

    bool        IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) const;
};

#endif // __MS_lte_07_h__
