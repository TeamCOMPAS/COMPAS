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

    HeMS() { m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS; };
    
    HeMS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), TPAGB(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                                                                                                                                         // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                             // Initialise if required
    }

    HeMS* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeMS* clone = new HeMS(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static HeMS* Clone(HeMS& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeMS* clone = new HeMS(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions - alphabetically
    static double   CalculateLifetimeOnPhase_Static(const double p_Mass);

    static double   CalculateLuminosityAtZAMS_Static(const double p_Mass);
    static double   CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Tau);
    static double   CalculateLuminosityAtPhaseEnd_Static(const double p_Mass);

           double   CalculateMassLossRateBelczynski2010();
           double   CalculateMassLossRateMerritt2024();

    static DBL_DBL  CalculateRadiusAtPhaseEnd_Static(const double p_Mass, const double p_Luminosity);
    static double   CalculateRadiusAtZAMS_Static(const double p_Mass);
    static double   CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Tau);
    
           double   CalculateRemnantRadius() const                                                                  { return Radius(); }

    MT_CASE         DetermineMassTransferTypeAsDonor() const                                                        { return MT_CASE::OTHER; }                                                      // Not A, B, C, or NONE


protected:

    void Initialise() {
        CalculateTimescales();
        // JR: Age for HeMS is partially calculated before switching -
        // can get here from various places in ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically
            DBL_DBL         CalculateConvectiveEnvelopeMass() const                                                 { return std::tuple<double, double> (0.0, 0.0); }                           // No convective envelope for naked He stars
            double          CalculateCOCoreMassAtPhaseEnd() const                                                   { return CalculateCOCoreMassOnPhase(); }                                        // Same as on phase
            double          CalculateCOCoreMassOnPhase() const                                                      { return 0.0;  }                                                                // McCO(HeMS) = 0.0

            double          CalculateConvectiveCoreMass() const;
            double          CalculateConvectiveCoreRadius() const;
            double          CalculateCoreMassAtPhaseEnd() const                                                     { return CalculateHeCoreMassOnPhase(); }                                        // Same as on phase /*ILYA*/ To fix, not everything will become CO core
            double          CalculateCoreMassOnPhase() const                                                        { return 0.0; }                                                                 // Mc(HeMS) = 0.0

    static  double          CalculateCoreMass_Luminosity_B_Static()                                                 { return 4.1E4; }
    static  double          CalculateCoreMass_Luminosity_D_Static(const double p_Mass)                              { return 5.5E4 / (1.0 + (0.4 * p_Mass * p_Mass * p_Mass * p_Mass)); }           // pow() is slow - use multiplication
    static  double          CalculateCoreMass_Luminosity_p_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs) { return 5.0; }
    static  double          CalculateCoreMass_Luminosity_q_Static(const double p_Mass, const DBL_VECTOR &p_MassCutoffs) { return 3.0; }

            double          CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const;
            double          CalculateCriticalMassRatioHurleyHjellmingWebbink() const                                { return 0.33; }                                                                // As coded in BSE. Using the inverse owing to how qCrit is defined in COMPAS. See Hurley et al. 2002 sect. 2.6.1 for additional details.

            void            CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams);
            void            CalculateGBParams()                                                                     { CalculateGBParams(m_Mass0, m_GBParams); }                                     // Use class member variables

            double          CalculateHeCoreMassOnPhase() const                                                      { return m_Mass; }                                                              // McHe(HeMS) = Mass
            double          CalculateHeCoreMassAtPhaseEnd() const                                                   { return CalculateHeCoreMassOnPhase(); }                                        // Same as on phase

            // Abundances
            double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                                          { return CalculateHeliumAbundanceCoreOnPhase(); }
            double          CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const;
            double          CalculateHeliumAbundanceCoreOnPhase() const                                             { return CalculateHeliumAbundanceCoreOnPhase(m_Tau); }                          // Use class member variables                                       
            
            double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const                                       { return CalculateHeliumAbundanceSurfaceOnPhase(); }
            double          CalculateHeliumAbundanceSurfaceOnPhase() const                                          { return m_HeliumAbundanceSurface; }                                            // Use class member variables                      

            double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                                        { return CalculateHydrogenAbundanceCoreOnPhase(); } 
            double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const                         { return 0.0; }
            double          CalculateHydrogenAbundanceCoreOnPhase() const                                           { return CalculateHydrogenAbundanceCoreOnPhase(m_Tau); }                        // Use class member variables                                 

            double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const                                     { return CalculateHydrogenAbundanceSurfaceOnPhase(); } 
            double          CalculateHydrogenAbundanceSurfaceOnPhase() const                                        { return m_HydrogenAbundanceSurface; }                                          // Use class member variables

            double          CalculateInitialSupernovaMass() const                                                   { return GiantBranch::CalculateInitialSupernovaMass(); }                        // Use GiantBranch

            double          CalculateLambdaDewi() const                                                             { return 0.5; }
            double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const  { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }                 // Not supported - use BaseStar (0.0 are dummy values)
            double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const                                { return CalculateLuminosityAtPhaseEnd_Static(p_Mass); }
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return CalculateLuminosityAtPhaseEnd(m_Mass); }                               // Use class member variables
            double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const               { return CalculateLuminosityOnPhase_Static(p_Mass, p_Tau); }
            double          CalculateLuminosityOnPhase() const                                                      { return CalculateLuminosityOnPhase(m_Mass, m_Tau); }                           // Use class member variables

            double          CalculateMassLossRateHurley();
            double          CalculateMassLossRateWolfRayetShenar2019() const;

            double          CalculateMassTransferRejuvenationFactor();

            double          CalculateMomentOfInertia() const                                                        { return MainSequence::CalculateMomentOfInertia(); }

            double          CalculatePerturbationMu() const                                                         { return 5.0; }                                                                 // Hurley et al. 2000, eqs 97 & 98

            double          CalculateRadiusAtPhaseEnd(const double p_Mass) const                                    { return CalculateRadiusAtPhaseEnd_Static(p_Mass); }
            double          CalculateRadiusAtPhaseEnd() const                                                       { return CalculateRadiusAtPhaseEnd(m_Mass); }                                   // Use class member variables
    static  double          CalculateRadiusAtPhaseEnd_Static(const double p_Mass);
            double          CalculateRadiusOnPhaseTau(const double p_Mass, const double p_Tau) const                { return CalculateRadiusOnPhase_Static(p_Mass, p_Tau); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhaseTau(m_Mass, m_Tau); }                            // Use class member variables

            double          CalculateTauAtPhaseEnd() const                                                          { return 1.0; }
            double          CalculateTauOnPhase() const                                                             { return m_Age / m_Timescales[static_cast<int>(TIMESCALE::tHeMS)]; }

            double          CalculateTemperatureAtPhaseEnd() const                                                  { return BaseStar::CalculateTemperatureAtPhaseEnd(); }
            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }                 // Same as on phase

            double          CalculateThermalMassLossRate() const                                                    { return BaseStar::CalculateThermalMassLossRate(); }                            // Use BaseStar

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                                 // Use class member variables
    
            double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)                  { return OPTIONS->ZetaMainSequence(); }                                         // A HeMS star is treated as any other MS star for Zeta calculation purposes

            double          ChooseTimestep(const double p_Time) const;

            STELLAR_TYPE    EvolveToNextPhase();

            ENVELOPE        DetermineEnvelopeType() const                                                           { return ENVELOPE::RADIATIVE; }                                                 // Always RADIATIVE

            double          InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, 
                                                   const double p_massTransferEfficiencyBeta)                       { return InterpolateGeEtAlQCrit(); }                                            // The function arguments are irrelavant for He stars, for now
            double          InterpolateGeEtAlQCrit(); 

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }
            bool            IsSupernova() const                                                                     { return false; }                                                               // Not here

            void            PerturbLuminosityAndRadius() { }                                                                                                                                        // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                                // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return m_StellarType; }                                                       // NO-OP

            void            SetSNHydrogenContent()                                                                  { m_SupernovaDetails.isHydrogenPoor = true; }                                   // Always true

            bool            ShouldEnvelopeBeExpelledByPulsations() const                                            { return false; }                                                               // No envelope to lose by pulsations
            bool            ShouldEvolveOnPhase() const                                                             { return (utils::Compare(m_Tau, 0.0) >= 0 && utils::Compare(m_Tau, 1.0) < 0); } // Evolve on HeMS phase if 0 <= tau < 1.0
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                               // Never skip HeMS phase

            void            UpdateInitialMass()                                                                     { m_Mass0 = m_Mass; }                                                           // Per Hurley et al. 2000, section 7.1
            void            UpdateAgeAfterMassLoss();                                                                                                                                               // Per Hurley et al. 2000, section 7.1
};

#endif // __HeMS_h__
