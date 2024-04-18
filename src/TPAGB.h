#ifndef __TPAGB_h__
#define __TPAGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
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
            double          CalculateCOCoreMassAtPhaseEnd() const                                                   { return (utils::Compare(m_COCoreMass, m_GBParams[static_cast<int>(GBP::McSN)]) >= 0 && utils::Compare(m_COCoreMass, m_Mass) < 0) ? m_COCoreMass : m_Mass; }
            double          CalculateCOCoreMassOnPhase() const                                                      { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                    // McCO(TPAGB) = Mc(TPAGB)Same as on phase

            double          CalculateCoreMassAtPhaseEnd() const                                                     { return m_CoreMass; }                                                                  // NO-OP
            double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const;
            double          CalculateCoreMassOnPhase() const                                                        { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                    // Use class member variables

            double          CalculateHeCoreMassAtPhaseEnd() const                                                   { return m_HeCoreMass; }                                                                // NO-OP
            double          CalculateHeCoreMassOnPhase() const                                                      { return m_CoreMass; }                                                                  // NO-OP

            double          CalculateLambdaDewi() const;
            double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
            double          CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const;

            double          CalculateLuminosityOnPhase(const double p_Time) const;
            double          CalculateLuminosityOnPhase() const                                                      { return CalculateLuminosityOnPhase(m_Age); }                                           // Use class member variables
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return m_Luminosity; }                                                                // NO-OP

            double          CalculateMcPrime(const double p_Time) const;

            double          CalculateRadiusAtPhaseEnd() const                                                       { return m_Radius; }                                                                    // NO-OP
            double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const            { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], m_BnCoefficients); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }                                // Use class member variables
    static  double          CalculateRadiusOnPhase_Static(const double      p_Mass,
                                                          const double      p_Luminosity,
                                                          const double      p_MHeF,
                                                          const DBL_VECTOR &p_BnCoefficients);

            double          CalculateRemnantLuminosity() const;
            double          CalculateRemnantRadius() const;

            double          CalculateTauAtPhaseEnd() const                                                          { return m_Tau; }                                                                       // NO-OP
            double          CalculateTauOnPhase() const                                                             { return m_Tau; }                                                                       // NO-OP

            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return m_Temperature; }                                                               // NO-OP
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }                      // Use class member variables

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                                         // Use class member variables

            double          ChooseTimestep(const double p_Time) const;

            ENVELOPE        DetermineEnvelopeType() const                                                           { return ENVELOPE::CONVECTIVE; }                                                        // Always CONVECTIVE

            STELLAR_TYPE    EvolveToNextPhase()                                                                     { return m_StellarType; }                                                                                                                                   // NO-OP

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }                                                      // Phase ends when envelope loss or going supernova
            bool            IsSupernova() const                                                                     { return (utils::Compare(m_COCoreMass, m_GBParams[static_cast<int>(GBP::McSN)]) >= 0 && utils::Compare(m_COCoreMass, m_Mass) < 0); } // Going supernova if still has envelope and core mass large enough

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
            void            ResolveHeliumFlash() { }                                                                                                                                                        // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return m_StellarType; }                                                               // NO-OP

            bool            ShouldEvolveOnPhase() const                                                             { return ((utils::Compare(m_COCoreMass, std::min(m_GBParams[static_cast<int>(GBP::McSN)], m_Mass)) < 0) && !ShouldEnvelopeBeExpelledByPulsations()); } // Evolve on TPAGB phase if envelope is not lost and not going supernova
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                                       // Never skip TPAGB phase

};

#endif // __TPAGB_h__
