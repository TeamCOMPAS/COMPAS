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

    NS() { m_StellarType = STELLAR_TYPE::NEUTRON_STAR; };
    
    NS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                             // Set stellar type
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


    // member functions - alphabetically
    static  DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass);

    static  double      CalculateRemnantMass_Static(const double p_COCoreMass)  { return 1.17 + (0.09 * p_COCoreMass); }                                        // Hurley et al., eq 92 (JR: Why is this in NS and not Remnants?) *Ilya*

    MT_CASE             DetermineMassTransferTypeAsDonor() const                { return MT_CASE::NONE; }                                                       // Always NONE


protected:
    
    void Initialise() {
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        
        //Set internal properties to zero to avoid meaningless values
        m_Age        = 0.0;
        m_COCoreMass = 0.0;
        m_HeCoreMass = 0.0;
        m_CoreMass   = 0.0;
        m_Mass0      = 0.0;
        
        m_Radius     = CalculateRadiusOnPhase();                                                                                                                // Set the NS radius, in Rsol
        m_Luminosity = CalculateLuminosityOnPhase();                                                                                                            // Set the NS luminosity

        CalculateAndSetPulsarParameters();
    }

    double m_AngularMomentum_CGS;                                                                                                                               // Current angular momentum in CGS - only required in NS class
    double m_MomentOfInertia_CGS;                                                                                                                               // MoI in CGS - only required in NS class


    // member functions - alphabetically
            void            CalculateAndSetPulsarParameters();

            double          CalculateBirthMagneticField();
            double          CalculateBirthSpinPeriod();

    static  double          CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time);
            double          CalculateLuminosityOnPhase() const                  { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age); }                    // Use class member variables

            double          CalculateMassLossRate()                             { return 0.0; }                                                                 // Ensure that NSs don't lose mass in winds
    
            double          CalculateMomentOfInertiaCGS() const;                                                                                                // MoI in CGS
            double          CalculateMomentOfInertia() const                    { return CalculateMomentOfInertiaCGS() / MSOL_TO_G / RSOL_TO_CM / RSOL_TO_CM; } // MoI (default is solar units)

    static  double          CalculateRadiusOnPhaseInKM_Static(const double p_Mass);                                                                             // Radius on phase in km
    static  double          CalculateRadiusOnPhase_Static(const double p_Mass)  { return CalculateRadiusOnPhaseInKM_Static(p_Mass) * KM_TO_RSOL; }              // Radius on phase in Rsol
            double          CalculateRadiusOnPhase() const                      { return CalculateRadiusOnPhase_Static(m_Mass); }                               // Radius on phase in Rsol

            double          CalculateSpinDownRate(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius) const;

            double          ChooseTimestep(const double p_Time) const;

            STELLAR_TYPE    EvolveToNextPhase()                                 { return STELLAR_TYPE::BLACK_HOLE; }

            void            ResolveCommonEnvelopeAccretion(const double p_FinalMass,
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