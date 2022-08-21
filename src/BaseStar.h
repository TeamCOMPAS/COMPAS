#ifndef __BaseStar_h__
#define __BaseStar_h__

#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"
#include "vector3d.h"

#include "Options.h"
#include "Log.h"
#include "Errors.h"


class BaseStar {

public:

    BaseStar();

    BaseStar(const unsigned long int p_RandomSeed, 
             const double            p_MZAMS, 
             const double            p_Metallicity, 
             const KickParameters    p_KickParameters,
             const double            p_RotationalVelocity = -1.0); 

    virtual ~BaseStar() {}


    // object identifiers - all classes have these
    OBJECT_ID           ObjectId() const                                                        { return m_ObjectId; }
    OBJECT_TYPE         ObjectType() const                                                      { return m_ObjectType; }
    STELLAR_TYPE        InitialStellarType() const                                              { return m_InitialStellarType; }
    STELLAR_TYPE        StellarType() const                                                     { return m_StellarType; }
    STELLAR_TYPE        StellarTypePrev() const                                                 { return m_StellarTypePrev; }


    // getters - alphabetically
            double              Age() const                                                     { return m_Age; }
            double              AngularMomentum() const                                         { return CalculateGyrationRadius() * m_Radius * RSOL_TO_AU * m_Radius * RSOL_TO_AU * m_Omega; }
            double              BindingEnergy_Fixed() const                                     { return m_BindingEnergies.fixed; }
            double              BindingEnergy_Nanjing() const                                   { return m_BindingEnergies.nanjing; }
            double              BindingEnergy_Loveridge() const                                 { return m_BindingEnergies.loveridge; }
            double              BindingEnergy_LoveridgeWinds() const                            { return m_BindingEnergies.loveridgeWinds; }
            double              BindingEnergy_Kruckow() const                                   { return m_BindingEnergies.kruckow; }
            bool                CHonMS() const                                                  { return m_CHE; }
            double              COCoreMass() const                                              { return m_COCoreMass; }
            double              CoreMass() const                                                { return m_CoreMass; }
            int                 DominantMassLossRate() const                                    { return static_cast<int>(m_DominantMassLossRate); }
            double              Dt() const                                                      { return m_Dt; }
            double              DtPrev() const                                                  { return m_DtPrev; }
            ERROR               Error() const                                                   { return m_Error; }
            bool                ExperiencedAIC() const                                          { return (m_SupernovaDetails.events.past & SN_EVENT::AIC) == SN_EVENT::AIC; }
            bool                ExperiencedCCSN() const                                         { return (m_SupernovaDetails.events.past & SN_EVENT::CCSN) == SN_EVENT::CCSN; }
            bool                ExperiencedECSN() const                                         { return (m_SupernovaDetails.events.past & SN_EVENT::ECSN) == SN_EVENT::ECSN; }
            bool                ExperiencedPISN() const                                         { return (m_SupernovaDetails.events.past & SN_EVENT::PISN) == SN_EVENT::PISN; }
            bool                ExperiencedPPISN() const                                        { return (m_SupernovaDetails.events.past & SN_EVENT::PPISN) == SN_EVENT::PPISN; }
            SN_EVENT            ExperiencedSN_Type() const                                      { return utils::SNEventType(m_SupernovaDetails.events.past); }
            bool                ExperiencedUSSN() const                                         { return (m_SupernovaDetails.events.past & SN_EVENT::USSN) == SN_EVENT::USSN; }
            double              HeCoreMass() const                                              { return m_HeCoreMass; }
            bool                IsAIC() const                                                   { return (m_SupernovaDetails.events.current & SN_EVENT::AIC) == SN_EVENT::AIC; }
            bool                IsCCSN() const                                                  { return (m_SupernovaDetails.events.current & SN_EVENT::CCSN) == SN_EVENT::CCSN; }
    virtual bool                IsDegenerate() const                                            { return false; }   // default is not degenerate - White Dwarfs, NS and BH are degenerate
            bool                IsECSN() const                                                  { return (m_SupernovaDetails.events.current & SN_EVENT::ECSN) == SN_EVENT::ECSN; }
    virtual bool                IsMassRatioUnstable(const double p_AccretorMass, const bool p_AccretorIsDegenerate) const { return false; } // default is stable
            bool                IsOneOf(const STELLAR_TYPE_LIST p_List) const;
            bool                IsPISN() const                                                  { return (m_SupernovaDetails.events.current & SN_EVENT::PISN) == SN_EVENT::PISN; }
            bool                IsPPISN() const                                                 { return (m_SupernovaDetails.events.current & SN_EVENT::PPISN) == SN_EVENT::PPISN; }
            bool                IsUSSN() const                                                  { return (m_SupernovaDetails.events.current & SN_EVENT::USSN) == SN_EVENT::USSN; }
            double              Lambda_Dewi() const                                             { return m_Lambdas.dewi; }
            double              Lambda_Fixed() const                                            { return m_Lambdas.fixed; }
            double              Lambda_Kruckow() const                                          { return m_Lambdas.kruckow; }
            double              Lambda_KruckowBottom() const                                    { return m_Lambdas.kruckowBottom; }
            double              Lambda_KruckowMiddle() const                                    { return m_Lambdas.kruckowMiddle; }
            double              Lambda_KruckowTop() const                                       { return m_Lambdas.kruckowTop; }
            double              Lambda_Loveridge() const                                        { return m_Lambdas.loveridge; }
            double              Lambda_LoveridgeWinds() const                                   { return m_Lambdas.loveridgeWinds; }
            double              Lambda_Nanjing() const                                          { return m_Lambdas.nanjing; }
            bool                LBV_Phase_Flag() const                                          { return m_LBVphaseFlag; }
            double              Luminosity() const                                              { return m_Luminosity; }
            double              Mass() const                                                    { return m_Mass; }
            double              Mass0() const                                                   { return m_Mass0; }
            double              MinimumCoreMass() const                                         { return m_MinimumCoreMass; }
            double              MassPrev() const                                                { return m_MassPrev; }
            STYPE_VECTOR        MassTransferDonorHistory() const                                { return m_MassTransferDonorHistory; }
            std::string         MassTransferDonorHistoryString() const;
            double              Mdot() const                                                    { return m_Mdot; }
            double              Metallicity() const                                             { return m_Metallicity; }
            double              MZAMS() const                                                   { return m_MZAMS; }
            double              Omega() const                                                   { return m_Omega; }
            double              OmegaCHE() const                                                { return m_OmegaCHE; }
            double              OmegaBreak() const                                              { return CalculateOmegaBreak(); }
            double              OmegaPrev() const                                               { return m_OmegaPrev; }
            double              OmegaZAMS() const                                               { return m_OmegaZAMS; }
            COMPAS_VARIABLE     PropertyValue(const T_ANY_PROPERTY p_Property) const;
            double              Pulsar_MagneticField() const                                    { return m_PulsarDetails.magneticField; }
            double              Pulsar_SpinPeriod() const                                       { return m_PulsarDetails.spinPeriod; }
            double              Pulsar_SpinFrequency() const                                    { return m_PulsarDetails.spinFrequency; }
            double              Pulsar_SpinDownRate() const                                     { return m_PulsarDetails.spinDownRate; }
            double              Radius() const                                                  { return m_Radius; }
            double              RadiusPrev() const                                              { return m_RadiusPrev; }
            unsigned long int   RandomSeed() const                                              { return m_RandomSeed; }
            double              RZAMS() const                                                   { return m_RZAMS; }
            double              SN_CoreMassAtCOFormation() const                                { return m_SupernovaDetails.coreMassAtCOFormation; }
            double              SN_COCoreMassAtCOFormation() const                              { return m_SupernovaDetails.COCoreMassAtCOFormation; }
            SupernovaDetailsT   SN_Details() const                                              { return m_SupernovaDetails; }
            double              SN_DrawnKickMagnitude() const                                   { return m_SupernovaDetails.drawnKickMagnitude; }
            double              SN_EccentricAnomaly() const                                     { return m_SupernovaDetails.eccentricAnomaly; }
            double              SN_FallbackFraction() const                                     { return m_SupernovaDetails.fallbackFraction; }
            double              SN_HeCoreMassAtCOFormation() const                              { return m_SupernovaDetails.HeCoreMassAtCOFormation; }
            bool                SN_IsHydrogenPoor() const                                       { return m_SupernovaDetails.isHydrogenPoor; }
            double              SN_KickMagnitude() const                                        { return m_SupernovaDetails.kickMagnitude; }
            double              SN_MeanAnomaly() const                                          { return m_SupernovaDetails.meanAnomaly; }
            double              SN_Phi() const                                                  { return m_SupernovaDetails.phi; }
            double              SN_TotalMassAtCOFormation() const                               { return m_SupernovaDetails.totalMassAtCOFormation; }
            double              SN_TrueAnomaly() const                                          { return m_SupernovaDetails.trueAnomaly; }
            double              SN_Theta() const                                                { return m_SupernovaDetails.theta; }
            SN_EVENT            SN_Type() const                                                 { return utils::SNEventType(m_SupernovaDetails.events.current); }
            double              SN_KickMagnitudeRandom() const                                  { return m_SupernovaDetails.kickMagnitudeRandom; }
            double              Speed() const                                                   { return m_ComponentVelocity.Magnitude(); }
            COMPAS_VARIABLE     StellarPropertyValue(const T_ANY_PROPERTY p_Property) const;
            double              Tau() const                                                     { return m_Tau; }
            double              Temperature() const                                             { return m_Temperature; }
            double              Time() const                                                    { return m_Time; }
            double              Timescale(TIMESCALE p_Timescale) const                          { return m_Timescales[static_cast<int>(p_Timescale)]; }
            double              XExponent() const                                               { return m_XExponent; }


    // setters
            void                SetInitialType(STELLAR_TYPE p_InitialType)                      { m_InitialStellarType = p_InitialType; }                                           // JR Could do some sanity checks here
            void                SetOmega(double p_vRot)                                         { if (p_vRot >= 0.0) m_Omega = p_vRot; };                                           // Do nothing if sanity check fails (JR: I don't really like this, but I think unavoidable - at least for now)

            void                SetSNCurrentEvent(SN_EVENT p_SNEvent)                           { m_SupernovaDetails.events.current |= p_SNEvent; }                                 // Set supernova primary event/state for current timestep
            void                SetSNPastEvent(const SN_EVENT p_SNEvent)                        { m_SupernovaDetails.events.past |= p_SNEvent; }                                    // Set supernova primary event/state for any past timestep
            
            void                UpdateComponentVelocity(const Vector3d p_newVelocity);	

            void                UpdateMassTransferDonorHistory();




    // member functions - alphabetically
            void            ApplyMassTransferRejuvenationFactor()                                               { m_Age *= CalculateMassTransferRejuvenationFactor(); }             // Apply age rejuvenation factor

            void            CalculateBindingEnergies(const double p_CoreMass, const double p_EnvMass, const double p_Radius);

            double          CalculateDynamicalTimescale() const                                                 { return CalculateDynamicalTimescale_Static(m_Mass, m_Radius); }         // Use class member variables

            double          CalculateEddyTurnoverTimescale();

    virtual void            CalculateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams) { }                                                                                      // Default is NO-OP
    virtual void            CalculateGBParams()                                                                 { CalculateGBParams(m_Mass0, m_GBParams); }                         // Use class member variables

    virtual double          CalculateGyrationRadius() const                                                     { return 0.0; }                                                     // Default is 0.0

            double          CalculateInterpolatedQCritGe2015() const;  

            void            CalculateLambdas()                                                                  { CalculateLambdas(m_Mass - m_CoreMass); }                          // Use class member variables
            void            CalculateLambdas(const double p_EnvMass);

    virtual DBL_DBL         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                        const double p_AccretorMassRate = 0.0);

            double          CalculateMassLossValues(const bool p_UpdateMDot = false, const bool p_UpdateMDt = false);                                                               // JR: todo: better name?

    virtual double          CalculateMomentOfInertia(const double p_RemnantRadius = 0.0) const                  { return 0.0; }                                                     // Use inheritance hierarchy
    virtual double          CalculateMomentOfInertiaAU(const double p_RemnantRadius = 0.0) const                { return 0.0; }                                                     // Use inheritance hierarchy
    
            double          CalculateNuclearTimescale() const                                                   { return CalculateNuclearTimescale_Static(m_Mass, m_Luminosity); }  // Use class member variables
    
            double          CalculateOmegaCHE(const double p_MZAMS, const double p_Metallicity) const;

            double          CalculateRadialChange() const                                                       { return (utils::Compare(m_RadiusPrev,0)<=0)? 0 : std::abs(m_Radius - m_RadiusPrev) / m_RadiusPrev; }                    // Return fractional radial change (if previous radius is negative or zero, return 0 to avoid NaN

            double          CalculateRadialExpansionTimescale() const                                           { return CalculateRadialExpansionTimescale_Static(m_StellarType, m_StellarTypePrev, m_Radius, m_RadiusPrev, m_DtPrev); } // Use class member variables
    
            void            CalculateSNAnomalies(const double p_Eccentricity);

            double          CalculateSNKickMagnitude(const double p_RemnantMass, const double p_EjectaMass, const STELLAR_TYPE p_StellarType);

            double          CalculateThermalMassAcceptanceRate(const double p_Radius) const;
            double          CalculateThermalMassAcceptanceRate() const                                          { return CalculateThermalMassAcceptanceRate(m_Radius); }

    virtual double          CalculateThermalMassLossRate() const                                                { return m_Mass / CalculateThermalTimescale(); }                    // Use class member variables - and inheritance hierarchy

    virtual double          CalculateThermalTimescale(const double p_Radius) const;                                                                                                 // Use inheritance hierarchy
    virtual double          CalculateThermalTimescale() const                                                   { return CalculateThermalTimescale(m_Radius); }                     // Use inheritance hierarchy

            double          CalculateTimestep();

    virtual double          CalculateZeta(ZETA_PRESCRIPTION p_ZetaPrescription)                                 { return 0.0; }                                                     // Use inheritance hierarchy

            void            ClearCurrentSNEvent()                                                               { m_SupernovaDetails.events.current = SN_EVENT::NONE; }             // Clear supernova event/state for current timestep

    virtual ENVELOPE        DetermineEnvelopeType() const                                                       { return ENVELOPE::REMNANT; }                                       // Default is REMNANT - but should never be called

            void            IncrementOmega(const double p_OmegaDelta)                                           { m_Omega += p_OmegaDelta; }                                        // Apply delta to current m_Omega

            void            ResolveAccretion(const double p_AccretionMass)                                      { m_Mass = std::max(0.0, m_Mass + p_AccretionMass); }               // Handles donation and accretion - won't let mass go negative

    virtual STELLAR_TYPE    ResolveEnvelopeLoss(bool p_NoCheck = false)                                         { return m_StellarType; }

    virtual void            ResolveMassLoss();

            void            SetStellarTypePrev(const STELLAR_TYPE p_StellarTypePrev)                            { m_StellarTypePrev = p_StellarTypePrev; }

            void            StashSupernovaDetails(const STELLAR_TYPE p_StellarType)                             { LOGGING->StashSSESupernovaDetails(this, p_StellarType); }

    virtual void            UpdateAgeAfterMassLoss() { }                                                                                                                             // Default is NO-OP

            STELLAR_TYPE    UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                              const double p_DeltaMass0,
                                                              const double p_DeltaTime,
                                                              const bool   p_ForceRecalculate);

    virtual void            UpdateInitialMass() { }                                                                                                                                 // Default is NO-OP

    virtual void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                                       const bool   p_RecyclesNS,
                                                       const double p_Stepsize,
                                                       const double p_MassGainPerTimeStep,
                                                       const double p_Epsilon) { }                                                                                                  // Default is NO-OP
    virtual void            UpdateMinimumCoreMass()  {}                                                                                                                 // Only set minimal core mass following Main Sequence mass transfer to MS age fraction of TAMS core mass; default is NO-OP

    
    // printing functions
            bool            PrintDetailedOutput(const int p_Id) const                                           { return OPTIONS->DetailedOutput() ? LOGGING->LogSSEDetailedOutput(this, p_Id) : true; } // Write record to SSE Detailed Output log file
            bool            PrintSupernovaDetails() const                                                       { return LOGGING->LogSSESupernovaDetails(this); }                   // Write record to SSE Supernovae log file
            bool            PrintStashedSupernovaDetails()                                                      { return LOGGING->LogStashedSSESupernovaDetails(this); }            // Write record to SSE Supernovae log file
            bool            PrintSwitchLog() const                                                              { return OPTIONS->SwitchLog() ? LOGGING->LogSSESwitchLog(this) : true; } // Write record to SSE Switchlog log file
            bool            PrintSystemParameters() const                                                       { return LOGGING->LogSSESystemParameters(this); }                   // Write record to SSE System Parameters file

protected:

    OBJECT_ID               m_ObjectId;                                 // Instantiated object's unique object id
    OBJECT_TYPE             m_ObjectType;                               // Instantiated object's object type
    STELLAR_TYPE            m_InitialStellarType;                       // Stellar type at birth, defined in Hurley et al. 2000
    STELLAR_TYPE            m_StellarType;                              // Stellar type defined in Hurley et al. 2000

    ERROR                   m_Error;                                    // Records most recent error encountered for this star

    // member variables - alphabetical in groups

    bool                    m_CHE;                                      // CHE flag - true if the star spent entire MS as a CH star; false if evolved CH->MS

    // Stellar variables
    unsigned long int       m_RandomSeed;                               // Seeds the random number generator for this star

    // Zero Age Main Sequence
    double                  m_LZAMS;                                    // ZAMS Luminosity
    double                  m_MZAMS;                                    // ZAMS Mass
    double                  m_OmegaZAMS;                                // ZAMS Angular Frequency
    double                  m_OmegaCHE;                                 // Minimum angular frequency at which CHE will occur (calculated at ZAMS)
    double                  m_RZAMS;                                    // ZAMS Radius
    double                  m_TZAMS;                                    // ZAMS Temperature

    // Effective Zero Age Main Sequence
    double                  m_LZAMS0;                                   // Effective ZAMS Luminosity
    double                  m_RZAMS0;                                   // Effective ZAMS Radius
    double                  m_TZAMS0;                                   // Effective ZAMS Temperature

    // Current timestep variables
    double                  m_Age;                                      // Current effective age (changes with mass loss/gain)(myrs)
    double                  m_COCoreMass;                               // Current CO core mass (Msol)
    double                  m_CoreMass;                                 // Current core mass (Msol)
    double                  m_Dt;                                       // Current timestep (myrs)
    double                  m_HeCoreMass;                               // Current He core mass (Msol)
    bool                    m_LBVphaseFlag;                             // Flag to know if the star satisfied the conditions, at any point in its evolution, to be considered a Luminous Blue Variable (LBV)
    double                  m_Luminosity;                               // Current luminosity (Lsol)
    double                  m_Mass;                                     // Current mass (Msol)
    double                  m_Mass0;                                    // Current effective initial mass (Msol)
    double                  m_MinimumCoreMass;                          // Minimum core mass at end of main sequence (MS stars have no core in the Hurley prescription)
    double                  m_MinimumLuminosityOnPhase;                 // JR: Only required for CHeB stars, but only needs to be calculated once per star
    double                  m_Mdot;                                     // Current mass loss rate (Msol per ?)
    MASS_LOSS_TYPE          m_DominantMassLossRate;                     // Current dominant mass loss rate
    double                  m_Mu;                                       // Current small envelope parameter mu
    double                  m_Omega;                                    // Current angular frequency (yr-1)
    double                  m_Radius;                                   // Current radius (Rsol)
    double                  m_Tau;                                      // Relative time
    double                  m_Temperature;                              // Current temperature (Tsol)
    double                  m_Time;                                     // Current physical time the star has been evolved(myrs)

    // Previous timestep variables
    double                  m_DtPrev;                                   // Previous timestep
    double                  m_MassPrev;                                 // Previous mass (Msol)
    double                  m_OmegaPrev;                                // Previous angular frequency (yr-1)
    double                  m_RadiusPrev;                               // Previous radius (Rsol)
    STELLAR_TYPE            m_StellarTypePrev;                          // Stellar type at previous timestep

    // Metallicity variables
    double                  m_LogMetallicityRho;                        // logMetallicityXi + 1.0       - called rho in Hurley et al 2000
    double                  m_LogMetallicitySigma;                      // log10(Metallicity)           - called sigma in Hurley et al 2000
    double                  m_LogMetallicityXi;                         // log10(Metallicity / Zsol)    - called xi in Hurley et al 2000
    double                  m_Metallicity;                              // Metallicity

    // Metallicity dependent constants
    double                  m_Alpha1;                                   // alpha1 in Hurly et al. 2000, just after eq 49
    double                  m_Alpha3;                                   // alpha4 in Hurley et al. 2000, just after eq 56
    double                  m_Alpha4;                                   // alpha4 in Hurley et al. 2000, just after eq 57
    double                  m_XExponent;                                // exponent to which R depends on M - 'x' in Hurley et al. 2000, eq 47


    // constants only calculated once
    double                  m_BaryonicMassOfMaximumNeutronStarMass;      // baryonic mass of MaximumNeutronStarMass 

    // JR:
    // I initially implemented the following vectors as unordered_maps.  The code worked
    // quite well, except for one small problem - access times (presumably due to hashing)
    // were prohibitive when accessed hundreds of thousands, and in some cases, millions,
    // of times as we evolve the star.  So I used vectors instead - the code is not as
    // elegant, but performance is better by an order of magnitude

    // Timescales, Giant Branch parameters, mass cutoffs
    DBL_VECTOR              m_GBParams;                                 // Giant Branch Parameters
    DBL_VECTOR              m_MassCutoffs;                              // Mass cutoffs
    DBL_VECTOR              m_Timescales;                               // Timescales

    // Luminosity, Radius, a(n) and b(n) coefficients
    DBL_VECTOR              m_AnCoefficients;                           // a(n) coefficients
    DBL_VECTOR              m_BnCoefficients;                           // b(n) coefficients
    DBL_VECTOR              m_LCoefficients;                            // Luminosity coefficients
    DBL_VECTOR              m_RCoefficients;                            // Radius coefficients

    // Luminosity, Radius and Gamma constants
    // JR:
    // These are calculated in CalculateAnCoefficients()
    // Calculating the a(n) coefficients requires one of the R constants, and the L, R and Gamma
    // constants are calculated using the a(n) coefficients, so it seemed logical to calculate
    // them all in the same function
    DBL_VECTOR              m_GammaConstants;                           // Gamma constants
    DBL_VECTOR              m_LConstants;                               // Luminosity constants
    DBL_VECTOR              m_RConstants;                               // Radius constants

    // Binding energies, Lambdas and Zetas
    BindingEnergiesT        m_BindingEnergies;                          // Binding enery values
    LambdasT                m_Lambdas;                                  // Lambda values

    // Stellar details squirrelled away...
    SupernovaDetailsT       m_SupernovaDetails;                         // Supernova attributes
    PulsarDetailsT          m_PulsarDetails;                            // Pulsar attributes

    // Star vector velocity 
    Vector3d                m_ComponentVelocity;                        // Isolated star velocity vector (binary's center-of-mass velocity for bound binary)

    // Star mass transfer history 
    STYPE_VECTOR            m_MassTransferDonorHistory;                 // List of MT donor stellar types - mostly relevent for binary stars

    // member functions - alphabetically
            void                AgeOneTimestepPreamble(const double p_DeltaTime);

            double              ApplyBlackHoleKicks(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass);

            double              CalculateAlpha1() const;
            double              CalculateAlpha3() const;
            double              CalculateAlpha4() const;

            void                CalculateAnCoefficients(DBL_VECTOR &p_AnCoefficients,
                                                        DBL_VECTOR &p_LConstants,
                                                        DBL_VECTOR &p_RConstants,
                                                        DBL_VECTOR &p_GammaConstants);

    virtual void                CalculateAndSetPulsarParameters() { }                                                                                                                               // NO-OP for most stellar types

            double              CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const;

            void                CalculateBnCoefficients(DBL_VECTOR &p_BnCoefficients);

    virtual double              CalculateCOCoreMassAtPhaseEnd() const                                                   { return m_COCoreMass; }                                                    // Default is NO-OP
    virtual double              CalculateCOCoreMassOnPhase() const                                                      { return m_COCoreMass; }                                                    // Default is NO-OP

    virtual double              CalculateCoreMassAtPhaseEnd() const                                                     { return m_CoreMass; }                                                      // Default is NO-OP
    static  double              CalculateCoreMassGivenLuminosity_Static(const double p_Luminosity, const DBL_VECTOR &p_GBParams);
    virtual double              CalculateCoreMassOnPhase() const                                                        { return m_CoreMass; }                                                      // Default is NO-OP

    static  double              CalculateDynamicalTimescale_Static(const double p_Mass, const double p_Radius);

    virtual double              CalculateEddingtonCriticalRate() const                                                  { return 2.08E-3 / 1.7 * m_Radius * MYR_TO_YEAR; }       // Hurley+, 2002, Eq. (67); should never be called...

            double              CalculateGBRadiusXExponent() const;

    virtual double              CalculateHeCoreMassAtPhaseEnd() const                                                   { return m_HeCoreMass; }                                                    // Default is NO-OP
    virtual double              CalculateHeCoreMassOnPhase() const                                                      { return m_HeCoreMass; }                                                    // Default is NO-OP

    static  double              CalculateHeRateConstant_Static()                                                        { return HE_RATE_CONSTANT; }                                                // Only >= CHeB stars need AHe, but no drama if other stars calculate (retrieve it) - it's only a constant (we could just use the constant inline...)
    static  double              CalculateHHeRateConstant_Static()                                                       { return HHE_RATE_CONSTANT; }                                               // Only TPAGB stars need AHHe, but no drama if other stars calculate (retrieve it) - it's only a constant (we could just use the constant inline...)

    static  double              CalculateInitialEnvelopeMass_Static(const double p_Mass);

    virtual double              CalculateLambdaDewi() const                                                             { SHOW_WARN(ERROR::NO_LAMBDA_DEWI, "Default used: 1.0"); return 1.0; }      // Not supported: show error
            double              CalculateLambdaKruckow(const double p_Radius, const double p_Alpha) const;
            double              CalculateLambdaLoveridgeEnergyFormalism(const double p_EnvMass, const double p_IsMassLoss = false) const;
    virtual double              CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const           { SHOW_WARN(ERROR::NO_LAMBDA_NANJING, "Default used: 1.0"); return 1.0; }   // Not supported: show error
    virtual double              CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const             { SHOW_WARN(ERROR::NO_LAMBDA_NANJING, "Default used: 1.0"); return 1.0; }   // Not supported: show error

            void                CalculateLCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_LCoefficients);

            double              CalculateLifetimeToBAGB(const double p_tHeI, const double p_tHe) const;
            double              CalculateLifetimeToBGB(const double p_Mass) const;

            double              CalculateLogBindingEnergyLoveridge(bool p_IsMassLoss) const;

            double              CalculateLuminosityAtBAGB(double p_Mass) const;
    virtual double              CalculateLuminosityAtPhaseEnd() const                                                   { return m_Luminosity; }                                                    // Default is NO-OP
            double              CalculateLuminosityAtZAMS(const double p_MZAMS);
            double              CalculateLuminosityGivenCoreMass(const double p_CoreMass) const;
    virtual double              CalculateLuminosityOnPhase() const                                                      { return m_Luminosity; }                                                    // Default is NO-OP

            double              CalculateMassAndZInterpolatedLambdaNanjing(const double p_Mass, const double p_Z) const;
            double              CalculateMassInterpolatedLambdaNanjing(const double p_Mass, const int p_Zind) const;
            double              CalculateZInterpolatedLambdaNanjing(const double p_Z, const int p_MassInd) const;

            void                CalculateMassCutoffs(const double p_Metallicity, const double p_LogMetallicityXi, DBL_VECTOR &p_MassCutoffs);

    static  double              CalculateMassLoss_Static(const double p_Mass, const double p_Mdot, const double p_Dt);

            double              CalculateMassLossRate();
    virtual double              CalculateMassLossRateHurley();
            double              CalculateMassLossRateKudritzkiReimers() const;
            double              CalculateMassLossRateLBV(const LBV_PRESCRIPTION p_LBV_prescription);
            double              CalculateMassLossRateLBVHurley(const double p_HD_limit_fac) const;
            double              CalculateMassLossRateLBVBelczynski() const;
            double              CalculateMassLossRateNieuwenhuijzenDeJager() const;
            double              CalculateMassLossRateOB(const double p_Teff);
            double              CalculateMassLossRateVassiliadisWood() const;
    virtual double              CalculateMassLossRateVink();
            double              CalculateMassLossRateWolfRayetZDependent(const double p_Mu) const;
            double              CalculateMassLossRateWolfRayet3() const;                                                                                                                            // JR: Never called - do we need it?
            double              CalculateMassLossRateWolfRayet(const double p_Mu) const;

    virtual double              CalculateMassTransferRejuvenationFactor() const;

            double              CalculateMaximumCoreMass(double p_Mass) const;

    static  double              CalculateNuclearTimescale_Static(const double p_Mass, const double p_Luminosity);

            double              CalculateOmegaBreak() const;

    static  double              CalculateOStarRotationalVelocityAnalyticCDF_Static(const double p_Ve);
    static  double              CalculateOStarRotationalVelocityAnalyticCDFInverse_Static(double p_Ve, void *p_Params);
    static  double              CalculateOStarRotationalVelocity_Static(const double p_Xmin, const double p_Xmax);

            double              CalculatePerturbationB(const double p_Mass) const;
            double              CalculatePerturbationC(double p_Mass) const;
    virtual double              CalculatePerturbationMu() const                                                         { return m_Mu; }                                                            // Default is NO-OP
    virtual double              CalculatePerturbationMuAtPhaseEnd() const                                               { return CalculatePerturbationMuOnPhase(); }                                // Same as on phase
    virtual double              CalculatePerturbationMuOnPhase() const                                                  { return CalculatePerturbationMu(); }
            double              CalculatePerturbationQ(const double p_Radius, const double p_Rc) const;
            double              CalculatePerturbationR(const double p_Mu, const double p_Mass, const double p_Radius, const double p_Rc) const;
            double              CalculatePerturbationS(const double p_Mu, const double p_Mass) const;

    static  double              CalculateRadialExpansionTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                                     const STELLAR_TYPE p_StellarTypePrev,
                                                                     const double       p_Radius,
                                                                     const double       p_RadiusPrev,
                                                                     const double       p_DtPrev);

            virtual double      CalculateRadialExtentConvectiveEnvelope() const                                         { return m_Radius; }                                                        // default for stars with no convective envelope

    virtual double              CalculateRadiusAtPhaseEnd() const                                                       { return m_Radius; }                                                        // Default is NO-OP
            double              CalculateRadiusAtZAMS(const double p_MZAMS) const;
    virtual double              CalculateRadiusOnPhase() const                                                          { return m_Radius; }                                                        // Default is NO-OP
    virtual std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return std::make_tuple(CalculateRadiusOnPhase(), m_StellarType); }

            void                CalculateRCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_RCoefficients);

            double              CalculateRotationalVelocity(double p_MZAMS) const;

    virtual double              CalculateTauOnPhase() const                                                             { return m_Tau; }                                                           // Default is NO-OP
    virtual double              CalculateTauAtPhaseEnd() const                                                          { return m_Tau; }                                                           // Default is NO-OP

    virtual double              CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }
    virtual double              CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }             // Same as on phase
            double              CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) const;
    virtual double              CalculateTemperatureOnPhase() const                                                     { return CalculateTemperatureOnPhase(m_Luminosity, m_Radius); }
    virtual double              CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const     { return CalculateTemperatureOnPhase_Static(p_Luminosity, p_Radius); }
    static  double              CalculateTemperatureOnPhase_Static(const double p_Luminosity, const double p_Radius);

    virtual void                CalculateTimescales()                                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                             // Use class member variables
    virtual void                CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) { }                                                                                              // Default is NO-OP

            double              CalculateZadiabaticHurley2002(const double p_CoreMass) const;
            double              CalculateZadiabaticSPH(const double p_CoreMass) const;
            double              CalculateZadiabatic(ZETA_PRESCRIPTION p_ZetaPrescription);

            double              CalculateZAMSAngularFrequency(const double p_MZAMS, const double p_RZAMS) const;

    virtual double              ChooseTimestep(const double p_Time) const                                               { return m_Dt; }

            double              DrawKickMagnitudeBrayEldridge(const double p_EjectaMass,
                                                              const double p_RemnantMass,
                                                              const double p_Alpha,
                                                              const double p_Beta) const;

            double              DrawKickMagnitudeDistributionFlat(const double p_MaxVK, const double p_Rand) const;
            double              DrawKickMagnitudeDistributionMaxwell(const double p_Sigma, const double p_Rand) const;

            double              DrawRemnantKickMuller(const double p_COCoreMass) const;

            double              DrawRemnantKickMullerMandel(const double p_COCoreMass,
                                                            const double p_Rand,
                                                            const double p_RemnantMass) const;

            double              DrawSNKickMagnitude(const double p_Sigma,
                                                    const double p_COCoreMass,
                                                    const double p_Rand,
                                                    const double p_EjectaMass,
                                                    const double p_RemnantMass) const;

            double              CalculateLambdaNanjing() const;
    virtual void                EvolveOneTimestepPreamble() { };                                                                                                                                    // Default is NO-OP

            STELLAR_TYPE        EvolveOnPhase();

    virtual STELLAR_TYPE        EvolveToNextPhase()                                                                     { return m_StellarType; }

            double              FindLambdaNanjingNearestMassIndex(const double p_Mass) const;

    virtual bool                IsEndOfPhase() const                                                                    { return false; }
    virtual bool                IsSupernova() const                                                                     { return false; }

            double              LimitTimestep(const double p_Dt);

    /*
     * Perturb Luminosity and Radius
     *
     * See Hurley et al. 2000, section 6.3
     *
     * The default is no perturbation - this function does nothing and is called
     * only if the stellar class doesn't define its own perturbation function.
     * See the stellar class perturbation functions for perturbation details specific
     * to the stellar class.
     *
     * Perturbation is disabled by default when DEBUG is enabled - except when
     * DEBUG_PERTURB is defined (see below).  The stellar class perturbation
     * functions are defined away if DEBUG is defined - so this generic Star
     * function is called (and does nothing).
     *
     * If DEBUG_PERTURB is defined then perturbation is not disabled while debbuging.
     * To enable perturbation while DEBUG is enabled, define DEBUG_PERTURB.
     */
    virtual void                PerturbLuminosityAndRadius() { }                                                                                                                                    // NO-OP
    virtual void                PerturbLuminosityAndRadiusAtPhaseEnd()                                                  { PerturbLuminosityAndRadiusOnPhase(); }                                    // Same as on phase
    virtual void                PerturbLuminosityAndRadiusOnPhase()                                                     { PerturbLuminosityAndRadius(); }

            STELLAR_TYPE        ResolveEndOfPhase();
    virtual void                ResolveHeliumFlash() { }
    virtual STELLAR_TYPE        ResolveSkippedPhase()                                                                   { return EvolveToNextPhase(); }                                             // Default is evolve to next phase
    virtual STELLAR_TYPE        ResolveSupernova()                                                                      { return m_StellarType; }                                                   // Default is NO-OP

    virtual void                SetSNHydrogenContent()                                                                  { m_SupernovaDetails.isHydrogenPoor = false; }                              // Default is false

            bool                ShouldBeMasslessRemnant() const                                                         { return (m_Mass <= 0.0 || m_StellarType==STELLAR_TYPE::MASSLESS_REMNANT); }
    virtual bool                ShouldEvolveOnPhase() const                                                             { return true; }
    virtual bool                ShouldSkipPhase() const                                                                 { return false; }                                                           // Default is false

            void                UpdateAttributesAndAgeOneTimestepPreamble(const double p_DeltaMass, const double p_DeltaMass0, const double p_DeltaTime);

};

#endif // __BaseStar_h__
