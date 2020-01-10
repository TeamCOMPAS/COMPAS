#ifndef __HeGB_h__
#define __HeGB_h__

#include "constants.h"
#include "typedefs.h"
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
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH;                                                                                                           // Set stellar type
        CalculateTimescales();                                                                                                                                                  // Initialise timescales
        if (previousStellarType != STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP)                                                                                             // If not evolving from HeHG...
            m_Age = CalculateAgeOnPhase_Static(m_Mass, m_COCoreMass, m_Timescales[static_cast<int>(TIMESCALE::tHeMS)], m_GBParams);                                             // ... Set age appropriately
    }


    // member functions - alphabetically
            double      CalculateGyrationRadius()                                                                       { return 0.1; }                                         // Hurley et al., 2000, after eq 109 for giants. Single number approximation.

            double      CalculateLuminosityOnPhase(const double p_CoreMass, const double p_GBPB, const double p_GBPD)   { return CalculateLuminosityOnPhase_Static(p_CoreMass, p_GBPB, p_GBPD); }
            double      CalculateLuminosityOnPhase()                                                                    { return CalculateLuminosityOnPhase(m_CoreMass, m_GBParams[static_cast<int>(GBP::B)], m_GBParams[static_cast<int>(GBP::D)]); }

            double      CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity);
            double      CalculateRadiusOnPhase()                                                                        { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }

            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity);
            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase()                                    { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }
            
            ENVELOPE    DetermineEnvelopeType()                                                                         { return ENVELOPE::CONVECTIVE; }                        // Always CONVECTIVE
            ENVELOPE    DetermineEnvelopeTypeHurley2002()                                                               { return ENVELOPE::CONVECTIVE; }                        // Always CONVECTIVE

            bool        IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate);
};

#endif // __HeGB_h__
