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
        if (p_Initialise) Initialise();
    }

    CH& operator = (const BaseStar &baseStar) {
        static_cast<BaseStar&>(*this) = baseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS;                                                                                                                       // Set stellar type
        CalculateTimescales();                                                                                                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                                                                                                // Set age appropriately
        m_CHE = true;                                                                                                                                                               // initially for CH stars                                                                                                                                                            // Set age appropriately
    }

    // member functions
    double          CalculateRadiusOnPhase() const      { return m_RZAMS; }                                                                                                         // Constant from birth
    double          CalculateRadiusAtPhaseEnd() const   { return CalculateRadiusOnPhase(); }                                                                                        // Same as on phase

    // Luminosity
    double          CalculateLogLuminosityRatio(const double p_Mass, const double p_Tau) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const                                   { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                      // Use class member variables

    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const;
    double          CalculateLuminosityOnPhase() const                                      { return CalculateLuminosityOnPhase(m_Age, m_Mass0, m_LZAMS0); }        // Use class member variables

    // Lifetime
    double          CalculateLogLifetimeRatio(const double p_Mass) const;
    double          CalculateLifetimeRatio(const double p_Mass) const;

    void            UpdateAgeAfterMassLoss();                                                                                                                       // Per Hurley et al. 2000, section 7.1

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    // Mass loss rate
    double          CalculateMassLossRateVink();
    double          CalculateMassLossRateWeightOB(const double p_tau);

    STELLAR_TYPE    EvolveToNextPhase();

    bool            ShouldEvolveOnPhase() const         { return m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)] && (OPTIONS->OptimisticCHE() || m_Omega >= m_OmegaCHE); }   // Evolve on CHE phase if age in MS timescale and spinning at least as fast as CHE threshold

};

#endif // __CH_h__
