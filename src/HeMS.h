#ifndef __HeMS_h__
#define __HeMS_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "TPAGB.h"


// JR: todo: revisit this one day - sometimes HeMS is MS, sometimes it is GiantBranch...
// Right now it is GiantBranch - figure out which is more appropriate

class BaseStar;
class TPAGB;

class HeMS: virtual public BaseStar, public TPAGB {

public:

    HeMS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), TPAGB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    HeMS& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically
    static double   CalculateLifetimeOnPhase_Static(const double p_Mass);

    static double   CalculateLuminosityAtZAMS_Static(const double p_Mass);
    static double   CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Tau);
    static double   CalculateLuminosityAtPhaseEnd_Static(const double p_Mass);


    static DBL_DBL  CalculateRadiusAtPhaseEnd_Static(const double p_Mass, const double p_Luminosity);
    static double   CalculateRadiusAtZAMS_Static(const double p_Mass);
    static double   CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Tau);


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                                                                                                                     // Set stellar type
        CalculateTimescales();
        // JR: Age for HeMS is partially calculated before switching -
        // can get here from various places in ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
    }

    void FastForward() {                                                                                                                                                        // Set stellar attributes for stars initialized to this stellar type

        m_Radius                                   = CalculateRadiusAtZAMS_Static(m_MZAMS);
        m_Luminosity                               = CalculateLuminosityAtZAMS_Static(m_MZAMS);
        m_Temperature                              = CalculateTemperatureOnPhase_Static(m_LZAMS, m_RZAMS);

        m_InitialLuminosity                        = m_Luminosity;
        m_InitialRadius                            = m_Radius;
        m_InitialTemperature                       = m_Temperature;
        m_InitialStellarType                       = m_StellarType;
        m_StellarTypePrev                          = m_StellarType;
        m_HeCoreMass                               = m_Mass;
    }

    // member functions - alphabetically
            double          CalculateCOCoreMassAtPhaseEnd() const                                       { return CalculateCOCoreMassOnPhase(); }                                // Same as on phase
            double          CalculateCOCoreMassOnPhase() const                                          { return 0.0; }                                                         // McCO(HeMS) = 0.0

            double          CalculateCoreMassAtPhaseEnd() const                                         { return CalculateHeCoreMassOnPhase(); }                                // Same as on phase
            double          CalculateCoreMassOnPhase() const                                            { return m_COCoreMass; }                                                // Mc(HeMS) = McCOMass

            double          CalculateGyrationRadius() const                                             { return 0.1; }                                                         // JR: todo: Nobody seems sure about this...

            double          CalculateHeCoreMassOnPhase() const                                          { return m_Mass; }                                                      // McHe(HeMS) = Mass
            double          CalculateHeCoreMassAtPhaseEnd() const                                       { return CalculateHeCoreMassOnPhase(); }                                // Same as on phase

            double          CalculateInitialSupernovaMass() const                                       { return GiantBranch::CalculateInitialSupernovaMass(); }                // Use GiantBranch

            double          CalculateLambdaDewi() const                                                 { return 0.5; }
            double          CalculateLambdaNanjingStarTrack(const double p_Mass,
                                                   const double p_Metallicity) const                    { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }                  // Not supported - use BaseStar (0.0 are dummy values)
            double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const                    { return CalculateLuminosityAtPhaseEnd_Static(p_Mass); }
            double          CalculateLuminosityAtPhaseEnd() const                                       { return CalculateLuminosityAtPhaseEnd(m_Mass); }                       // Use class member variables
            double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const   { return CalculateLuminosityOnPhase_Static(p_Mass, p_Tau); }
            double          CalculateLuminosityOnPhase() const                                          { return CalculateLuminosityOnPhase(m_Mass, m_Tau); }                   // Use class member variables

            double          CalculateMassLossRateHurley();
            double          CalculateMassLossRateVink();

            double          CalculateMassTransferRejuvenationFactor() const;

            double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0) const          { return MainSequence::CalculateMomentOfInertia(p_RemnantRadius); }     // Use MainSequence
            double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0) const        { return MainSequence::CalculateMomentOfInertiaAU(p_RemnantRadius); }   // Use MainSequence

            double          CalculatePerturbationMu() const                                             { return 5.0; }                                                         // Hurley et al. 2000, eqs 97 & 98

            double          CalculateRadialExtentConvectiveEnvelope() const                             { return BaseStar::CalculateRadialExtentConvectiveEnvelope(); }         // HeMS stars don't have a convective envelope

            double          CalculateRadiusAtPhaseEnd(const double p_Mass) const                        { return CalculateRadiusAtPhaseEnd_Static(p_Mass); }
            double          CalculateRadiusAtPhaseEnd() const                                           { return CalculateRadiusAtPhaseEnd(m_Mass); }                           // Use class member variables
    static  double          CalculateRadiusAtPhaseEnd_Static(const double p_Mass);
            double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau) const       { return CalculateRadiusOnPhase_Static(p_Mass, p_Tau); }
            double          CalculateRadiusOnPhase() const                                              { return CalculateRadiusOnPhase(m_Mass, m_Tau); }                       // Use class member variables

            double          CalculateTauAtPhaseEnd() const                                              { return 1.0; }
            double          CalculateTauOnPhase() const                                                 { return m_Age / m_Timescales[static_cast<int>(TIMESCALE::tHeMS)]; }

            double          CalculateTemperatureAtPhaseEnd() const                                      { return BaseStar::CalculateTemperatureAtPhaseEnd(); }

            double          CalculateThermalMassLossRate() const                                        { return BaseStar::CalculateThermalMassLossRate(); }                    // Use BaseStar

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                       { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables
    
            double          CalculateZeta(ZETA_PRESCRIPTION p_ZetaPrescription)                         { return OPTIONS->ZetaMainSequence(); }                                 // A HeMS star is treated as any other MS star for Zeta calculation purposes

            double          ChooseTimestep(const double p_Time) const;

            STELLAR_TYPE    EvolveToNextPhase();

            ENVELOPE        DetermineEnvelopeType() const                                               { return ENVELOPE::RADIATIVE; }                                         // Always RADIATIVE

            bool            IsEndOfPhase() const                                                        { return !ShouldEvolveOnPhase(); }
            bool            IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) const;
            bool            IsSupernova() const                                                         { return false; }                                                       // Not here

            void            PerturbLuminosityAndRadius() { }                                                                                                                    // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false);
            void            ResolveHeliumFlash() { }                                                                                                                            // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                       { return m_StellarType; }                                               // NO-OP

            void            SetSNHydrogenContent()                                                      { m_SupernovaDetails.isHydrogenPoor = true; }                           // Always true

            bool            ShouldEvolveOnPhase() const                                                 { return (utils::Compare(m_Tau, 0.0) >= 0 && utils::Compare(m_Tau, 1.0) < 0); } // Evolve on HeMS phase if 0 <= tau < 1.0
            bool            ShouldSkipPhase() const                                                     { return false; }                                                       // Never skip HeMS phase

            void            UpdateInitialMass()                                                         { m_Mass0 = m_Mass; }                                                   // Per Hurley et al. 2000, section 7.1
            void            UpdateAgeAfterMassLoss();                                                                                                                           // Per Hurley et al. 2000, section 7.1

};

#endif // __HeMS_h__
