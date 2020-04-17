#ifndef __TPAGB_h__
#define __TPAGB_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "EAGB.h"


class BaseStar;
class EAGB;

class TPAGB: virtual public BaseStar, public EAGB {

public:

    TPAGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), EAGB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    TPAGB& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH;                                                                                                                            // Set stellar type
        CalculateTimescales();                                                                                                                                                                              // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tP)];                                                                                                                                              // Set age appropriately
   }


   // member functions - alphabetically
            double          CalculateCOCoreMassAtPhaseEnd()                                                  { return (utils::Compare(m_COCoreMass, m_GBParams[static_cast<int>(GBP::McSN)]) >= 0 && utils::Compare(m_COCoreMass, m_Mass) < 0) ? m_COCoreMass : m_Mass; }
            double          CalculateCOCoreMassOnPhase()                                                     { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                            // McCO(TPAGB) = Mc(TPAGB)Same as on phase

            double          CalculateCoreMassAtPhaseEnd()                                                    { return m_CoreMass; }                                                                          // NO-OP
            double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time);
            double          CalculateCoreMassOnPhase()                                                       { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                            // Use class member variables

            double          CalculateEnvelopeMassOnPhase(const double p_Tau)                                 { return m_Mass - std::max(m_COCoreMass, m_HeCoreMass);}                                                                           // NO-OP

            double          CalculateGyrationRadius()                                                        { return 0.1; }                                                                                 // Hurley et al., 2000, after eq 109 for giants. Single number approximation.

            double          CalculateHeCoreMassAtPhaseEnd()                                                  { return m_HeCoreMass; }                                                                        // NO-OP
            double          CalculateHeCoreMassOnPhase()                                                     { return m_CoreMass; }                                                                        // NO-OP

            double          CalculateLambdaDewi();
            double          CalculateLambdaNanjing();

            double          CalculateLuminosityOnPhase(const double p_Time);
            double          CalculateLuminosityOnPhase()                                                     { return CalculateLuminosityOnPhase(m_Age); }                                                   // Use class member variables
            double          CalculateLuminosityAtPhaseEnd()                                                  { return m_Luminosity; }                                                                        // NO-OP

            double          CalculateMcPrime(const double p_Time);

            double          CalculateRadiusAtPhaseEnd()                                                      { return m_Radius; }                                                                            // NO-OP
            double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity)           { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], m_BnCoefficients); }
            double          CalculateRadiusOnPhase()                                                         { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }                                        // Use class member variables
    static  double          CalculateRadiusOnPhase_Static(const double      p_Mass,
                                                          const double      p_Luminosity,
                                                          const double      p_MHeF,
                                                          const DBL_VECTOR &p_BnCoefficients);

            double          CalculateRemnantLuminosity();
            double          CalculateRemnantRadius();

            double          CalculateTauAtPhaseEnd()                                                         { return m_Tau; }                                                                               // NO-OP
            double          CalculateTauOnPhase()                                                            { return m_Tau; }                                                                               // NO-OP

            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) { return m_Temperature; }                                                                       // NO-OP
            double          CalculateTemperatureAtPhaseEnd()                                                 { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }                              // Use class member variables

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                            { CalculateTimescales(m_Mass0, m_Timescales); }                                                 // Use class member variables

            double          ChooseTimestep(const double p_Time);

            ENVELOPE        DetermineEnvelopeType()                                                          { return ENVELOPE::CONVECTIVE; }                                                                // Always CONVECTIVE
            ENVELOPE        DetermineEnvelopeTypeHurley2002()                                                { return ENVELOPE::CONVECTIVE; }                                                                // Always CONVECTIVE

            STELLAR_TYPE    EvolveToNextPhase() { return m_StellarType; }                                                                                                                                    // NO-OP

            bool            IsEndOfPhase()                                                                   { return !ShouldEvolveOnPhase(); }                                                              // Phase ends when envelope loss or going supernova
            bool            IsSupernova()                                                                    { return (utils::Compare(m_COCoreMass, m_GBParams[static_cast<int>(GBP::McSN)]) >= 0 && utils::Compare(m_COCoreMass, m_Mass) < 0); }    // Going supernova if envelope lost and core mass large enough

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
            void            ResolveHeliumFlash() { }                                                                                                                                                         // NO-OP
            STELLAR_TYPE    ResolveRemnantAfterEnvelopeLoss();
            STELLAR_TYPE    ResolveSkippedPhase()                                                            { return m_StellarType; }                                                                       // NO-OP

            bool            ShouldEvolveOnPhase()                                                            { return (utils::Compare(m_COCoreMass, std::min(m_GBParams[static_cast<int>(GBP::McSN)], m_Mass)) < 0); }   // Evolve on TPAGB phase if envelope is not lost and not going supernova
            bool            ShouldSkipPhase()                                                                { return false; }                                                                               // Never skip TPAGB phase

};

#endif // __TPAGB_h__
