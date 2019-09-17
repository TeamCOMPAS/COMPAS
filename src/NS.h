#ifndef __NS_h__
#define __NS_h__

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "ONeWD.h"


class BaseStar;
class ONeWD;

class NS: virtual public BaseStar, public ONeWD {

public:

    NS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), ONeWD(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    NS& operator = (const BaseStar &p_BaseStar) {
        static_cast<BaseStar&>(*this) = p_BaseStar;
        Initialise();
        return *this;
    }


    // member functions - alphabetically
    static  DBL_DBL_DBL     CalculateCoreCollapseSNParams_Static(const double p_Mass);

    static  double          CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time);

    static  double          CalculatePulsarBirthSpinPeriod_Static();

    static  double          CalculateRadiusOnPhase_Static(const double p_Mass);

    static  double          CalculateRemnantMass_Static(const double p_CoreMass)    { return 1.17 + (0.09 * p_CoreMass); }                                          // Hurley et al., eq 92


protected:

    void Initialise() {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                                 // Set stellar type
        CalculateTimescales();                                                                                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                                                                                // Set age appropriately
    }


    // member functions - alphabetically
            void            CalculateAndSetPulsarParameters();

            double          CalculateCOCoreMassOnPhase()                            { return m_Mass; }                                                              // return m_Mass

            double          CalculateConvergedMassStepZetaThermal()                 { return 1.0; }                                                                 // For NS & BH  JR: todo: check this - BH seems to be different...

            double          CalculateHeCoreMassOnPhase()                            { return m_Mass; }                                                              // return m_Mass

            double          CalculateLuminosityOnPhase()                            { return CalculateLuminosityOnPhase_Static(m_Mass, m_Age); }                    // Use class member variables

            DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                        const double p_FractionAccreted,
                                                        const double p_AccretorMassRate = 0.0 );

    static  double          CalculateMomentOfInertia_Static(const double p_Mass, const double p_Radius);

    static  double          CalculatePulsarBirthMagneticField_Static();

    static  double          CalculateSpinDownRate_Static(const double p_Omega,
                                                         const double p_MomentOfInteria,
                                                         const double p_MagField,
                                                         const double p_Radius,
                                                         double const p_Alpha);

            ENVELOPE        DetermineEnvelopeType()                                 { return ENVELOPE::REMNANT; }                                                   // Always REMNANT
            ENVELOPE        DetermineEnvelopeTypeHurley2002()                       { return ENVELOPE::REMNANT; }                                                   // Always REMNANT

            void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                                       const double p_Stepsize,
                                                       const double p_MassGainPerTimeStep,
                                                       const double p_Epsilon);

            bool            ShouldEvolveOnPhase()                                   { return true; }                                                                // Always
            bool            ShouldSkipPhase()                                       { return false; }                                                               // Don't skip
};

#endif // __NS_h__
