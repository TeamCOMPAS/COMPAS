#ifndef __NS_h__
#define __NS_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"
#include "BH.h"


class BaseStar;
class Remnants;

class NS: virtual public BaseStar, public Remnants {

public:

    NS() {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                             // Set stellar type

        // set NS values based on options provided (so we don't need to do this every timestep)
        // JR: these are good candidates for the "globals/constants" singleton...
        ALFVEN_CONST             = PPOW(2.0 * PI_2 / G_CGS, 1.0 / 7.0);
        NS_MAG_FIELD_LOWER_LIMIT = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField());
        NS_DECAY_MASS_SCALE      = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;
        NS_DECAY_TIME_SCALE      = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR; 
    };
    
    NS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                             // Set stellar type

        // set NS values based on options provided (so we don't need to do this every timestep)
        // JR: these are good candidates for the "globals/constants" singleton...
        ALFVEN_CONST             = PPOW(2.0 * PI_2 / G_CGS, 1.0 / 7.0);
        NS_MAG_FIELD_LOWER_LIMIT = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField());
        NS_DECAY_MASS_SCALE      = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;
        NS_DECAY_TIME_SCALE      = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR; 

        if (p_Initialise) Initialise();                                                                                                                         // Initialise if required
    }

    NS* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        NS* clone = new NS(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static NS* Clone(NS& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        NS* clone = new NS(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member variables

    // static variables that only need to be calculated once
    // be aware that these variables are global (because they're static), and are shared amongst all NS instances
    // JR: these are good candidates for the "globals/constants" singleton...

    static inline double ALFVEN_CONST { 0.0 };                                                                                                                  // Constant for calculating Alfven radius - CGS units (mu0 = 1.0)
    static inline double NS_MAG_FIELD_LOWER_LIMIT { 0.0 };
    static inline double NS_DECAY_MASS_SCALE { 0.0 };
    static inline double NS_DECAY_TIME_SCALE { 0.0 }; 


    // member functions - alphabetically
    static  DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass);

    static double       DeltaJByAccretion_Static(const double p_Mass, const double p_Radius_6, const double p_MagField, const double p_SpinFrequency, const double p_mDot, const double p_Epsilon);

    MT_CASE             DetermineMassTransferTypeAsDonor() const                { return MT_CASE::NONE; }                                                       // Always NONE


protected:
    
    void Initialise() {
        
        // set internal properties to zero to avoid meaningless values
        m_Age        = 0.0;
        m_COCoreMass = 0.0;
        m_HeCoreMass = 0.0;
        m_CoreMass   = 0.0;
        m_Mass0      = 0.0;
        
        EvolveOnPhase(0.0);

        CalculateAndSetPulsarParameters();
    }

    // member variables

    double m_AngularMomentum_CGS;                                                                                                                               // Current angular momentum in CGS - only required in NS class
    double m_MomentOfInertia_CGS;                                                                                                                               // MoI in CGS - only required in NS class


    // member functions - alphabetically
            void            CalculateAndSetPulsarParameters();

            double          CalculateBirthMagneticField();
            double          CalculateBirthSpinPeriod();
    
            double          CalculateCriticalMassRatioHurleyHjellmingWebbink() const { return 0.0; }

    static  double          CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time);
            double          CalculateLuminosityOnPhase() const                  { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age); }                    // Use class member variables

            double          CalculateMassLossRate()                             { return 0.0; }                                                                 // Ensure that NSs don't lose mass in winds
    
    static  double          CalculateMomentOfInertiaCGS_Static(const double p_Mass, const double p_Radius);                                                     // MoI in CGS            
            double          CalculateMomentOfInertiaCGS() const                 { return CalculateMomentOfInertiaCGS_Static(m_Mass * MSOL_TO_G, m_Radius * RSOL_TO_CM); } // MOI in CGS - use member variables
            double          CalculateMomentOfInertia() const                    { return CalculateMomentOfInertiaCGS() / MSOL_TO_G / RSOL_TO_CM / RSOL_TO_CM; } // MoI (default is solar units)

    static  double          CalculateRadiusOnPhaseInKM_Static(const double p_Mass);                                                                             // Radius on phase in km
    static  double          CalculateRadiusOnPhase_Static(const double p_Mass)  { return CalculateRadiusOnPhaseInKM_Static(p_Mass) * KM_TO_RSOL; }              // Radius on phase in Rsol
            double          CalculateRadiusOnPhase() const                      { return CalculateRadiusOnPhase_Static(m_Mass); }                               // Radius on phase in Rsol
    
            double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const { return CalculateRadiusOnPhase(); }                               // not a meaningful calculation for NS, ignore arguments

            double          CalculateSpinDownRate(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius) const;
  
            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) { }                                                              // not a meaningful calculation for NS and BH
            void            CalculateTimescales() { }                                                                                                           // not a meaningful calculation for NS and BH

            double          ChooseTimestep(const double p_Time) const;

            STELLAR_TYPE    EvolveToNextPhase()                                 { return STELLAR_TYPE::BLACK_HOLE; }

            double          ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                                           const double p_CompanionMass     = 0.0,
                                                           const double p_CompanionRadius   = 0.0,
                                                           const double p_CompanionEnvelope = 0.0);
    
            bool            ShouldEvolveOnPhase() const                         { return (m_Mass <= OPTIONS->MaximumNeutronStarMass()); }                       // Evolve as a neutron star unless mass > maximum neutron star mass (e.g. through accretion)
            void            SpinDownIsolatedPulsar(const double p_Stepsize);
            void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                                       const bool   p_RecycledNS,
                                                       const double p_Stepsize,
                                                       const double p_MassGainPerTimeStep,
                                                       const double p_Epsilon);

};

#endif // __NS_h__
