#ifndef __WhiteDwarfs_h__
#define __WhiteDwarfs_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"


class BaseStar;
class Remnants;

class WhiteDwarfs: virtual public BaseStar, public Remnants {

public:

    WhiteDwarfs(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    WhiteDwarfs& operator = (const BaseStar &p_BaseStar) { static_cast<BaseStar&>(*this) = p_BaseStar; return *this; }


    // member functions
    static  double      CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                         const double p_Time, 
                                                         const double p_Metallicity, 
                                                         const double p_BaryonNumber);

    static  double      CalculateRadiusOnPhase_Static(const double p_Mass);


    void                ResolveShellChange(const double p_AccretedMass);

protected:
    // member variables

            bool                    m_HeShellDetonation;                        // Flag to initialize He-Shell detonation (i.e. as described in Wang. 2018, sect 5 2018RAA....18...49W)
            double                  m_HeShell;                                  // Current WD He-shell size (Msol). Increases through accretion.
            double                  m_HShell;                                   // Current WD H-shell size (Msol). Increases through accretion.
            double                  m_l0Ritter;                                 // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2. Corresponds to L0.
            double                  m_lambdaRitter;                             // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2.
            bool                    m_OffCenterIgnition;                        // Flag for CO WD evolution into ONe WD
            bool                    m_ShouldRejuvenate;                         // Flag for evolution of HeWD back into HeMS
            bool                    m_IsSubChandrasekharTypeIa;                 // Flag for SubCh SN of HeWD
            double                  m_XRitter;                                  // Assumed hydrogen-mass fraction of material being accreted by He WD, as in Ritter 1999, table 2.
            ACCRETION_REGIME        m_AccretionRegime;
            
            // member functions - alphabetically
            double  CalculateAccretionRegime(const bool p_DonorIsHeRich,
                                             const bool p_DonorIsGiant,
                                             const double p_DonorThermalMassLossRate,
                                             const double p_MassLostByDonor);
        
            double      CalculateCOCoreMassOnPhase() const                          { return m_COCoreMass; }                                                // NO-OP

            double      CalculateHeCoreMassOnPhase() const                          { return m_HeCoreMass; }                                                // NO-OP

            double      CalculateEtaH(const double p_MassIntakeRate);

            double      CalculateEtaHe(const double p_MassIntakeRate);

            double      CalculateEtaPTY(const double p_MassIntakeRate);

            double      Calculatel0Ritter() const                                   {return (m_Metallicity > 0.01) ? 1995262.3 : 31622.8; } // Luminosity constant which depends on metallicity in Ritter 1999, eq 10

    virtual DBL_DBL     CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                    const bool   p_IsHeRich)        { return std::make_tuple(0.0, 0.0); }                   // Should never be called
            DBL_DBL     CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                    const double p_AccretorMassRate,
                                                    const bool   p_IsHeRich)        { return CalculateMassAcceptanceRate(p_DonorMassRate, p_IsHeRich); } // Ignore the input accretion rate for WDs

            double      CalculateXRitter() const                                    {return (m_Metallicity > 0.01) ? 0.7 : 0.8 ; }          // Assumed Hydrogen-mass fraction

            double      CalculatelambdaRitter() const                               {return (m_Metallicity > 0.01) ? 8 : 5 ;  }             // Exponent for the assumed core-mass and luminosity relationship in Ritter 1999

            double      CalculateInitialSupernovaMass() const                       { return 0.0; }

            double      CalculateRadiusOnPhase(const double p_Mass) const           { return CalculateRadiusOnPhase_Static(p_Mass); }
            double      CalculateRadiusOnPhase() const                              { return CalculateRadiusOnPhase(m_Mass); }              // Use class member variables

            ENVELOPE    DetermineEnvelopeType() const                               { return ENVELOPE::CONVECTIVE; }                        // Always CONVECTIVE

            bool        IsMassAboveEcsnThreshold() const                            { return (utils::Compare(m_Mass, MECS) > 0); }          // Mass exceeds ECSN threshold mass
            bool        IsMassAboveChandrasekhar() const                            { return (utils::Compare(m_Mass, MCH) > 0); }           // Mass exceeds Chandrasekhar limit 

            STELLAR_TYPE ResolveAIC();  
            STELLAR_TYPE ResolveSNIa();  
            STELLAR_TYPE ResolveHeSD();  
            STELLAR_TYPE ResolveSupernova()                                         { return EvolveToNextPhase(); }                         // SNe for WDs are handled internally to each WD type

            ACCRETION_REGIME WhiteDwarfAccretionRegime() const                      { return m_AccretionRegime; }

};

#endif // __WhiteDwarfs_h__
