#ifndef __MainSequence_h__
#define __MainSequence_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"

class BaseStar;

class MainSequence: virtual public BaseStar {

public:

    MainSequence(){};

    MainSequence(const BaseStar& p_BaseStar) : BaseStar(p_BaseStar) {}

    MT_CASE DetermineMassTransferTypeAsDonor() const                                        { return MT_CASE::A; }                                                  // Always case A
    
    const std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> SHIKAUCHI_COEFFICIENTS = InterpolateShikauchiCoefficients(m_Metallicity);                                 // Interpolate Shikauchi coefficients for the given metallicity

protected:
    
    // member variables
    
    double          m_HeliumAbundanceCoreOut      = m_InitialHeliumAbundance;                                                                                       // Helium abundance just outside the core, used for rejuvenation calculations
    double          m_InitialMainSequenceCoreMass = 0.0;                                                                                                            // Initial mass of the mixing core is initialised in MS_gt_07 class

    // member functions - alphabetically
    double          CalculateAlphaL(const double p_Mass) const;
    double          CalculateAlphaR(const double p_Mass) const;

    double          CalculateConvectiveCoreMass() const;
    double          CalculateConvectiveCoreRadius() const;
    DBL_DBL         CalculateConvectiveEnvelopeMass() const;
    double          CalculateBetaL(const double p_Mass) const;
    double          CalculateBetaR(const double p_Mass) const;

    double          CalculateDeltaL(const double p_Mass) const;
    double          CalculateDeltaR(const double p_Mass) const;

    double          CalculateEta(const double p_Mass) const;

    double          CalculateGamma(const double p_Mass) const;

    double          CalculateCOCoreMassAtPhaseEnd() const                                   { return CalculateCOCoreMassOnPhase(); }                                // Same as on phase
    double          CalculateCOCoreMassOnPhase() const                                      { return 0.0; }                                                         // McCO(MS) = 0.0

    double          CalculateCoreMassAtPhaseEnd() const                                     { return (OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::MANDEL) ? MainSequenceCoreMass() : 0.0; }      // Accounts for minimal core mass built up prior to mass loss through mass transfer
    double          CalculateCoreMassOnPhase() const                                        { return 0.0; }                                                         // Mc(MS) = 0.0 (Hurley et al. 2000, just before eq 28)

    double          CalculateHeCoreMassAtPhaseEnd() const                                   { return CalculateCoreMassAtPhaseEnd(); }                               // Same as He core mass
    double          CalculateHeCoreMassOnPhase() const                                      { return 0.0; }                                                         // McHe(MS) = 0.0

    double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                          { return CalculateHeliumAbundanceCoreOnPhase(); }
    double          CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const;                                         
    double          CalculateHeliumAbundanceCoreOnPhase() const                             { return CalculateHeliumAbundanceCoreOnPhase(m_Tau); }                  // Use class member variables                                       
    
    double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const                       { return CalculateHeliumAbundanceSurfaceOnPhase(); }
    double          CalculateHeliumAbundanceSurfaceOnPhase() const                          { return m_InitialHeliumAbundance; }                                    // Use class member variables                      
    
    double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                        { return CalculateHydrogenAbundanceCoreOnPhase(); } 
    double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const;                                                          
    double          CalculateHydrogenAbundanceCoreOnPhase() const                           { return CalculateHydrogenAbundanceCoreOnPhase(m_Tau); }                // Use class member variables                                 
    
    double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const                     { return CalculateHydrogenAbundanceSurfaceOnPhase(); } 
    double          CalculateHydrogenAbundanceSurfaceOnPhase() const                        { return m_InitialHydrogenAbundance; }                                  // Use class member variables
    
    double          CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const                                   { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                      // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const;
    double          CalculateLuminosityOnPhase() const                                      { return CalculateLuminosityOnPhase(m_Age, m_Mass0, m_LZAMS0); }        // Use class member variables
    double          CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore, const double p_Age) const;
    DBL_DBL         CalculateMainSequenceCoreMassShikauchi(const double p_Dt, const double p_MassLossRate);
    double          CalculateInitialMainSequenceCoreMass(const double p_MZAMS) const;
    double          CalculateMomentOfInertia() const                                        { return (0.1 * (m_Mass) * m_Radius * m_Radius); }                      // k2 = 0.1 as defined in Hurley et al. 2000, after eq 109

    double          CalculatePerturbationMu() const                                         { return 5.0; }                                                         // mu(MS) = 5.0 (Hurley et al. 2000, eqs 97 & 98)

    double          CalculateRadialExtentConvectiveEnvelope() const;

    double          CalculateRadiusOnPhaseTau(const double p_Mass, const double p_Tau) const;

    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Time, const double p_RZAMS) const;
    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS) const;
    double          CalculateRadiusAtPhaseEnd() const                                       { return CalculateRadiusAtPhaseEnd(m_Mass, m_RZAMS); }                  // Use class member variables
    double          CalculateRadiusOnPhase() const                                          { return CalculateRadiusOnPhase(m_Mass, m_Age, m_RZAMS0); }             // Use class member variables
    double          CalculateRadiusTransitionToHG(const double p_Mass, const double p_Age, double const p_RZAMS) const;
     
    double          CalculateTauAtPhaseEnd() const                                          { return 1.0; }                                                         // tau = 1.0 at end of MS
    double          CalculateTauOnPhase() const;

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)  { return OPTIONS->ZetaMainSequence(); }

    double          ChooseTimestep(const double p_Time) const;

    void            EvolveOneTimestepPreamble();
    STELLAR_TYPE    EvolveToNextPhase()                                                     { return STELLAR_TYPE::HERTZSPRUNG_GAP; }

    double          InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, 
                                           const double p_massTransferEfficiencyBeta); // RTW do I need a const here?
    
    std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> InterpolateShikauchiCoefficients(const double p_Metallicity) const;
    
    bool            IsEndOfPhase() const                                                    { return !ShouldEvolveOnPhase(); }                                      // Phase ends when age at or after MS timescale

    void            PerturbLuminosityAndRadius() { }                                                                                                                // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);

    bool            ShouldEvolveOnPhase() const                                             { return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]); }    // Evolve on MS phase if age in MS timescale
    
    double          TAMSCoreMass() const;

    void            UpdateInitialMass()                                                     { m_Mass0 = m_Mass; }                                                   // Per Hurley et al. 2000, section 7.1
   
    void            UpdateAfterMerger(double p_Mass, double p_HydrogenMass);
    
    void            UpdateAgeAfterMassLoss();                                                                                                                       // Per Hurley et al. 2000, section 7.1

    void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_MassLossRate);

};

#endif // __MainSequence_h__
