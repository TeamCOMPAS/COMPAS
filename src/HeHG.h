#ifndef __HeHG_h__
#define __HeHG_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeMS.h"


class BaseStar;
class HeMS;

class HeHG: virtual public BaseStar, public HeMS {

public:

    HeHG() { m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP; };
    
    HeHG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeMS(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;                                                                                                                    // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                     // Initialise if required
    }

    HeHG* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeHG* clone = new HeHG(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static HeHG* Clone(HeHG& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeHG* clone = new HeHG(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions
    double  CalculateConvectiveCoreRadius () const                                                          { return std::min(5.0 * CalculateRemnantRadius(), m_Radius); }          // Last paragraph of section 6 of Hurley+ 2000
    static void CalculateGBParams_Static(const double p_Mass0, const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_BnCoefficients, DBL_VECTOR &p_GBParams);
 
    double  CalculateRemnantRadius() const;

protected:

    void Initialise() {
        m_Tau = 0.0;                                                                                      // Start of phase
        CalculateTimescales();                                                                            // Initialise timescales
        // Age for HeHG is calculated before switching -
        // can get here via EvolveOneTimestep() and ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
        
        // Update stellar properties at start of HeHG phase (since core definition changes)
        CalculateGBParams();

        EvolveOnPhase(0.0);
    }


    // member functions - aphabetically
            double          CalculateCOCoreMassAtPhaseEnd() const                                                   { return m_COCoreMass; }  //*ILYA* check                                               // NO-OP
            double          CalculateCOCoreMassOnPhase() const;
    
            double          CalculateConvectiveCoreMass() const                                                     { return m_CoreMass; }

            double          CalculateCoreMassAtBAGB() const                                                         { return m_Mass0; }                                                     // McBAGB = M0 (Hurely et al. 2000, discussion just before eq 89)
            double          CalculateCoreMassAtPhaseEnd() const                                                     { return m_CoreMass; }                                                  // NO-OP
            double          CalculateCoreMassOnPhase() const                                                        { return m_COCoreMass; }                                                // Mc(HeMS) = McCOMass

            double          CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const;
            double          CalculateCriticalMassRatioHurleyHjellmingWebbink() const                                { return 1.28; }                                                        // From BSE. Using the inverse owing to how qCrit is defined in COMPAS. See Hurley et al. 2002 sect. 2.6.1 for additional details.

            void            CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams);
            void            CalculateGBParams()                                                                     { CalculateGBParams(m_Mass0, m_GBParams); }                             // Use class member variables

            double          CalculateHeCoreMassAtPhaseEnd() const                                                   { return CalculateHeCoreMassOnPhase(); }                                // Same as on phase
            double          CalculateHeCoreMassOnPhase() const                                                      { return m_Mass; }                                                      // NO-OP

            double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
            double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const { return CalculateLambdaNanjingStarTrack(0.0, 0.0); }        // 0.0 are dummy values that are not used

            double          CalculateLuminosityOnPhase() const;
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return m_Luminosity; }                                                // NO-OP

            double          CalculateMassTransferRejuvenationFactor()                                               { return 1.0; }

            double          CalculateMomentOfInertia() const                                                        { return GiantBranch::CalculateMomentOfInertia(); }

            double          CalculatePerturbationMu() const;
            double          CalculatePerturbationMuAtPhaseEnd() const                                               { return m_Mu; }                                                        // NO-OP

            double          CalculateRadiusAtPhaseEnd() const                                                       { return m_Radius; }                                                    // NO-OP
            double          CalculateRadiusOnPhase() const;

            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const;
            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                          { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }

            double          CalculateTauAtPhaseEnd() const                                                          { return m_Tau; }                                                       // NO-OP
            double          CalculateTauOnPhase() const                                                             { return 0.0; }

            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return m_Temperature; }                                               // NO-OP
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }      // Use class member variables

            double          CalculateThermalMassLossRate() const                                                    { return GiantBranch::CalculateThermalMassLossRate(); }                 // Skip HeMS

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables
    
            double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)                  { return GiantBranch::CalculateZetaConstantsByEnvelope(p_ZetaPrescription); } // Calculate Zetas as for other giant stars (HeMS stars were an exception)
    
            double          CalculateZetaEquilibrium()                                                              { return -std::numeric_limits<double>::infinity(); }                     // Nuclear timescale MT should be impossible from HG stars that evolve on a thermal timescale

            double          ChooseTimestep(const double p_Time) const;

            ENVELOPE        DetermineEnvelopeType() const;

            STELLAR_TYPE    EvolveToNextPhase();

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }
            bool            IsSupernova() const;
            double          CalculateInitialSupernovaMass() const;
    
            void            PerturbLuminosityAndRadius()                                                            { GiantBranch::PerturbLuminosityAndRadius(); }                          // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                        // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return m_StellarType; }                                               // NO-OP

            bool            ShouldEnvelopeBeExpelledByPulsations() const                                            { return CHeB::ShouldEnvelopeBeExpelledByPulsations(); }                // Envelope of convective star with luminosity to mass ratio beyond threshold should be expelled
            bool            ShouldEvolveOnPhase() const;
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                       // Never skip HeMS phase

            void            UpdateAgeAfterMassLoss()                                                                { GiantBranch::UpdateAgeAfterMassLoss(); }                              // No action for He giants
            void            UpdateInitialMass()                                                                     { GiantBranch::UpdateInitialMass(); }                                   // Skip HeMS
};

#endif // __HeHG_h__
