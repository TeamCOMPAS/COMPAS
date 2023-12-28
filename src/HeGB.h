#ifndef __HeGB_h__
#define __HeGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeHG.h"


class BaseStar;
class HeHG;

class HeGB: virtual public BaseStar, public HeHG {

public:

    HeGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeHG(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    HeGB& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically
    static  double      CalculateAgeOnPhase_Static(const double p_Mass, const double p_CoreMass, const double p_tHeMS, const DBL_VECTOR &p_GBParams);

    static  double      CalculateCoreMassOnPhase_Static(const double p_Mass, const double p_Time, const double p_tHeMS, const DBL_VECTOR &p_GBParams);

    static  double      CalculateLuminosityOnPhase_Static(const double p_CoreMass, const double p_GBPB, const double p_GBPD);

    static  DBL_DBL     CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity);


protected:

    void Initialise() {
        STELLAR_TYPE previousStellarType = m_StellarType;;
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH;                                                                                                   // Set stellar type
        CalculateTimescales();                                                                                                                                          // Initialise timescales
        if (previousStellarType != STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP)                                                                                     // If not evolving from HeHG...
            m_Age = CalculateAgeOnPhase_Static(m_Mass, m_COCoreMass, m_Timescales[static_cast<int>(TIMESCALE::tHeMS)], m_GBParams);                                     // ... Set age appropriately
    }


    // member functions - alphabetically
    double      CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const ;
    double      CalculateCriticalMassRatioHurleyHjellmingWebbink() const                                            { return 1.28; }                                        // From BSE. Using the inverse owing to how qCrit is defined in COMPAS. See Hurley et al. 2002 sect. 2.6.1 for additional details.

    double      CalculateLuminosityOnPhase(const double p_CoreMass, const double p_GBPB, const double p_GBPD) const { return CalculateLuminosityOnPhase_Static(p_CoreMass, p_GBPB, p_GBPD); }
    double      CalculateLuminosityOnPhase() const                                                                  { return CalculateLuminosityOnPhase(m_CoreMass, m_GBParams[static_cast<int>(GBP::B)], m_GBParams[static_cast<int>(GBP::D)]); }

    double      CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const;
    double      CalculateRadiusOnPhase() const                                                                      { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }

    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const;
    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                                  { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }
            
    ENVELOPE    DetermineEnvelopeType() const                                                                       { return ENVELOPE::CONVECTIVE; }                        // Always CONVECTIVE
};

#endif // __HeGB_h__
