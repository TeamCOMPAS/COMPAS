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

    WhiteDwarfs(){};

    WhiteDwarfs(const BaseStar &p_BaseStar) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {}


    // member functions
    static  double  CalculateLuminosityOnPhase_Static(const double p_Mass, 
                                                      const double p_Time, 
                                                      const double p_Metallicity, 
                                                      const double p_BaryonNumber);

    static  double  CalculateRadiusOnPhase_Static(const double p_Mass);

    MT_CASE         DetermineMassTransferTypeAsDonor() const                                { return MT_CASE::OTHER; }                                  // Not A, B, C, or NONE

    ACCRETION_REGIME DetermineAccretionRegime(const double p_DonorThermalMassLossRate, const bool p_HeRich);                                            // Get the current accretion regime. Can also change m_HeShellDetonation and m_OffCenterIgnition flags.
    
    void            ResolveShellChange(const double p_AccretedMass);


protected:
    // member variables

            bool             m_HeShellDetonation;                                                                                                       // Flag to initialize He-Shell detonation (i.e. as described in Wang. 2018, sect 5 2018RAA....18...49W)
            double           m_HeShell;                                                                                                                 // Current WD He-shell size (Msol). Increases through accretion.
            double           m_HShell;                                                                                                                  // Current WD H-shell size (Msol). Increases through accretion.
            double           m_L0Ritter;                                                                                                                // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2. Corresponds to L0.
            double           m_LambdaRitter;                                                                                                            // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2.
            bool             m_OffCenterIgnition;                                                                                                       // Flag for CO WD evolution into ONe WD
            bool             m_ShouldRejuvenate;                                                                                                        // Flag for evolution of HeWD back into HeMS
            bool             m_IsSubChandrasekharTypeIa;                                                                                                // Flag for SubCh SN of HeWD
            double           m_XRitter;                                                                                                                 // Assumed hydrogen-mass fraction of material being accreted by He WD, as in Ritter 1999, table 2.
            ACCRETION_REGIME m_AccretionRegime;
            
            // member functions - alphabetically
            double           CalculateAccretionRegime(const bool   p_DonorIsHeRich,
                                                      const bool   p_DonorIsGiant,
                                                      const double p_DonorThermalMassLossRate,
                                                      const double p_MassLostByDonor);

            double          CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate,
                                               const double p_massTransferEfficiencyBeta)                       { return CalculateCriticalMassRatioHurleyHjellmingWebbink(); }
            double          CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const         { return CalculateCriticalMassRatioHurleyHjellmingWebbink(); }
            double          CalculateCriticalMassRatioGeEtAl(const QCRIT_PRESCRIPTION p_qCritPrescription,
                                                     const double p_massTransferEfficiencyBeta)                 { return CalculateCriticalMassRatioHurleyHjellmingWebbink(); }
            double           CalculateCriticalMassRatioHurleyHjellmingWebbink() const       { return HURLEY_HJELLMING_WEBBINK_QCRIT_WD; }
        
            double           CalculateCOCoreMassOnPhase() const                             { return m_COCoreMass; }                                    // NO-OP

            double           CalculateHeCoreMassOnPhase() const                             { return m_HeCoreMass; }                                    // NO-OP

            double           CalculateHeliumAbundanceCoreOnPhase() const                    { return 0.0; }
            double           CalculateHeliumAbundanceSurfaceOnPhase() const                 { return 0.0; }
            
            double           CalculateHydrogenAbundanceCoreOnPhase() const                  { return 0.0; }
            double           CalculateHydrogenAbundanceSurfaceOnPhase() const               { return 0.0; }

            double           CalculateEtaH(const double p_MassIntakeRate);

            double           CalculateEtaHe(const double p_MassIntakeRate);

            double           CalculateEtaPTY(const double p_MassIntakeRate);

            double           Calculatel0Ritter() const                                      { return (m_Metallicity > 0.01) ? L0_RITTER_HIGH_Z : L0_RITTER_LOW_Z; }

            DBL_DBL          CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                         const bool   p_IsHeRich);          
            DBL_DBL          CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                         const double p_AccretorMassRate,
                                                         const bool   p_IsHeRich)           { return CalculateMassAcceptanceRate(p_DonorMassRate, p_IsHeRich); }

            double           CalculateXRitter() const                                       { return utils::Compare(m_Metallicity, 0.01) > 0 ? 0.7 : 0.8; } // Assumed Hydrogen-mass fraction

            double           CalculateLambdaRitter() const                                  { return utils::Compare(m_Metallicity, 0.01) > 0 ? 8.0 : 5.0; } // Exponent for the assumed core-mass and luminosity relationship in Ritter 1999

            double           CalculateInitialSupernovaMass() const                          { return 0.0; }

            double           CalculateRadiusOnPhase(const double p_Mass) const              { return CalculateRadiusOnPhase_Static(p_Mass); }
            double           CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const { return CalculateRadiusOnPhase(p_Mass); }                // Ignore luminosity argument for WDs
            double           CalculateRadiusOnPhase() const                                 { return CalculateRadiusOnPhase(m_Mass); }                  // Use class member variables

            ENVELOPE         DetermineEnvelopeType() const                                  { return ENVELOPE::CONVECTIVE; }                            // Always CONVECTIVE

            bool             IsMassAboveChandrasekhar() const                               { return (utils::Compare(m_Mass, MCH) > 0); }               // Mass exceeds Chandrasekhar limit 

            STELLAR_TYPE     ResolveAIC();  
            STELLAR_TYPE     ResolveSNIa();  
            STELLAR_TYPE     ResolveHeSD();  
            STELLAR_TYPE     ResolveSupernova()                                             { return EvolveToNextPhase(); }                             // SNe for WDs are handled internally to each WD type

            ACCRETION_REGIME WhiteDwarfAccretionRegime() const                              { return m_AccretionRegime; }

};

#endif // __WhiteDwarfs_h__
