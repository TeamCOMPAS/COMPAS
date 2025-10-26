// gsl includes
#include <gsl/gsl_roots.h>
#include <gsl/gsl_cdf.h>

#include "Rand.h"
#include "BaseStar.h"
#include "vector3d.h"
#include "BH.h"

// boost includes
#include <boost/math/distributions.hpp>

using std::max;
using std::min;


BaseStar::BaseStar() {

    // initialise member variables

    m_ObjectId           = globalObjectId++;                                                        // unique object id - remains for life of star (even through evolution to other phases)
    m_ObjectPersistence  = OBJECT_PERSISTENCE::PERMANENT;                                           // object persistence - permanent or ephemeral (ephemeral used for ephemeral clones)
    m_InitialStellarType = STELLAR_TYPE::STAR;                                                      // stellar type - changes throughout life of star (through evolution to other phases)
    m_StellarType        = STELLAR_TYPE::STAR;                                                      // stellar type - changes throughout life of star (through evolution to other phases)

    m_Error              = ERROR::NONE;                                                             // clear error flag
}


BaseStar::BaseStar(const unsigned long int p_RandomSeed, 
                   const double            p_MZAMS,
                   const double            p_Metallicity,
                   const KickParameters    p_KickParameters,
                   const double            p_RotationalFrequency) {

    // initialise member variables

    m_ObjectId            = globalObjectId++;                                                       // unique object id - remains for life of star (even through evolution to other phases)
    m_ObjectPersistence   = OBJECT_PERSISTENCE::PERMANENT;                                          // object persistence - permanent or ephemeral (ephemeral used for ephemeral clones)
    m_InitialStellarType  = STELLAR_TYPE::STAR;                                                     // stellar type - changes throughout life of star (through evolution to other phases)
    m_StellarType         = STELLAR_TYPE::STAR;                                                     // stellar type - changes throughout life of star (through evolution to other phases)

    m_Error               = ERROR::NONE;                                                            // clear error flag

    m_CHE                 = false;                                                                  // initially

    m_EvolutionStatus     = EVOLUTION_STATUS::CONTINUE;                                             // initially
    
    // Initialise member variables from input parameters
    // (kick parameters initialised below - see m_SupernovaDetails)
    m_RandomSeed          = p_RandomSeed;
    m_MZAMS               = p_MZAMS;
    m_Metallicity         = p_Metallicity;

    // Initialise metallicity dependent values
    m_Log10Metallicity    = log10(m_Metallicity);


    // Initialise coefficients, parameters and constants

    // initialise m_Timescales vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(TIMESCALE::COUNT); i++) {
        m_Timescales.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);
    }

    // initialise m_GBParams vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(GBP::COUNT); i++) {
        m_GBParams.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);
    }

    // initialise m_MassCutoffs vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(MASS_CUTOFF::COUNT); i++) {
        m_MassCutoffs.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);
    }

    // initialise m_LConstants vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(L_CONSTANTS::COUNT); i++) {
        m_LConstants.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);
    }

    // initialise m_RConstants vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(R_CONSTANTS::COUNT); i++) {
        m_RConstants.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);
    }

    // initialise m_GammaConstants vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(GAMMA_CONSTANTS::COUNT); i++) {
        m_GammaConstants.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);
    }


    // calculate coefficients, constants etc.

    CalculateRCoefficients(LogMetallicityXiHurley(), m_RCoefficients);
    CalculateLCoefficients(LogMetallicityXiHurley(), m_LCoefficients);

    CalculateMassCutoffs(m_Metallicity, LogMetallicityXiHurley(), m_MassCutoffs);

    CalculateAnCoefficients(m_AnCoefficients, m_LConstants, m_RConstants, m_GammaConstants);
    CalculateBnCoefficients(m_BnCoefficients);

    m_XExponent                                = CalculateGBRadiusXExponent();
    m_Alpha1                                   = CalculateAlpha1();
    m_Alpha3                                   = CalculateAlpha3();
    m_Alpha4                                   = CalculateAlpha4();

    // initialise remaining member variables

    // Zero age main sequence parameters
    m_InitialMainSequenceCoreMass              = DEFAULT_INITIAL_DOUBLE_VALUE;                      // initialised in MS_gt_07 class if BRCEK core mass prescription is used
    m_RZAMS                                    = CalculateRadiusAtZAMS(m_MZAMS);
    m_LZAMS                                    = CalculateLuminosityAtZAMS(m_MZAMS);
    m_TZAMS                                    = CalculateTemperatureOnPhase_Static(m_LZAMS, m_RZAMS);

    // Initial abundances
    m_InitialHeliumAbundance                   = CalculateInitialHeliumAbundance();
    m_HeliumAbundanceCore                      = m_InitialHeliumAbundance;
    m_HeliumAbundanceSurface                   = m_InitialHeliumAbundance;

    m_InitialHydrogenAbundance                 = CalculateInitialHydrogenAbundance();
    m_HydrogenAbundanceCore                    = m_InitialHydrogenAbundance;
    m_HydrogenAbundanceSurface                 = m_InitialHydrogenAbundance;

    // Effective initial Zero Age Main Sequence parameters corresponding to Mass0
    m_RZAMS0                                   = m_RZAMS;
    m_LZAMS0                                   = m_LZAMS;

    // Current timestep attributes
    m_Time                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Dt                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Tau                                      = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Age                                      = 0.0;                                               // ensure age = 0.0 at construction (rather than default initial value)
    m_MainSequenceCoreMass                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mass                                     = m_MZAMS;
    m_Mass0                                    = m_MZAMS;
    m_Luminosity                               = m_LZAMS;
    m_Radius                                   = m_RZAMS;
    m_Temperature                              = m_TZAMS;
    m_TotalMassLossRate                        = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_ComponentVelocity						   = Vector3d();
    
    m_OmegaCHE                                 = CalculateOmegaCHE(m_MZAMS, m_Metallicity);
    m_OmegaZAMS                                = p_RotationalFrequency >= 0.0                           // valid rotational frequency passed in?
                                                    ? _2_PI * p_RotationalFrequency                     // yes - convert from cycles/yr to rad/yr and use it
                                                    : CalculateZAMSAngularFrequency(m_MZAMS, m_RZAMS);  // no - calculate it
    m_AngularMomentum                          = CalculateMomentOfInertiaAU() * m_OmegaZAMS;

    m_CoreMass                                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_COCoreMass                               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_HeCoreMass                               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mu                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mdot                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_DominantMassLossRate                     = MASS_LOSS_TYPE::NONE;

    m_MinimumLuminosityOnPhase                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_LBVphaseFlag                             = false;
    m_EnvelopeJustExpelledByPulsations         = false;

    // Previous timestep attributes
    m_StellarTypePrev                          = m_StellarType;
    m_MassPrev                                 = m_MZAMS;
    m_RadiusPrev                               = m_RZAMS;
    m_DtPrev                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    
    // Supernova details

    m_SupernovaDetails.initialKickParameters   = p_KickParameters;

    m_SupernovaDetails.events.current          = SN_EVENT::NONE;
    m_SupernovaDetails.events.past             = SN_EVENT::NONE;

    m_SupernovaDetails.coreMassAtCOFormation   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.coreRadiusAtCOFormation = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.COCoreMassAtCOFormation = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.HeCoreMassAtCOFormation = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.totalMassAtCOFormation  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.totalRadiusAtCOFormation= DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.drawnKickMagnitude      = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.kickMagnitude           = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.isHydrogenPoor          = false;
    m_SupernovaDetails.fallbackFraction        = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.eccentricAnomaly        = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.trueAnomaly             = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.supernovaState          = SN_STATE::NONE;

    m_SupernovaDetails.kickMagnitudeRandom     = p_KickParameters.magnitudeRandom;
    m_SupernovaDetails.theta                   = p_KickParameters.theta;
    m_SupernovaDetails.phi                     = p_KickParameters.phi;
    m_SupernovaDetails.meanAnomaly             = p_KickParameters.meanAnomaly;

    m_SupernovaDetails.rocketKickMagnitude     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.rocketKickPhi           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.rocketKickTheta         = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Calculates the Baryonic mass for which the GravitationalRemnantMass will be equal to the maximumNeutronStarMass (inverse of SolveQuadratic())
    // needed to decide whether to calculate Fryer+2012 for Neutron Star or Black Hole in GiantBranch::CalculateGravitationalRemnantMass()
    // calculate only once for entire simulation of N binaries in the future.
    m_BaryonicMassOfMaximumNeutronStarMass     = (0.075 * OPTIONS->MaximumNeutronStarMass() * OPTIONS->MaximumNeutronStarMass()) + OPTIONS->MaximumNeutronStarMass();

    // Pulsar details
    m_PulsarDetails.magneticField              = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.spinPeriod                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.spinFrequency              = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.spinDownRate               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.birthPeriod                = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.birthSpinDownRate          = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Mass Transfer Donor Type History
    m_MassTransferDonorHistory                 = ST_VECTOR();

}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  CLASS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Determine the value of the requested property of the constituent star (parameter p_Property)
 *
 * The property is a boost variant variable, and is one of the following types:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *      BINARY_PROPERTY         - property of the binary
 *      PROGRAM_OPTION          - program option
 *
 * This function handles properties of type:
 * 
 *    STAR_PROPERTY, STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY
 * 
 * only - anything else will result in an error being thrown and the evolution of the star (or binary)
 * terminated.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  
 *
 *
 * COMPAS_VARIABLE StellarPropertyValue(const T_ANY_PROPERTY p_Property)
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseStar::StellarPropertyValue(const T_ANY_PROPERTY p_Property) const {

    COMPAS_VARIABLE value;

    ANY_STAR_PROPERTY property;

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {

        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : { STAR_PROPERTY      prop = boost::get<STAR_PROPERTY>(p_Property);      property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : { STAR_1_PROPERTY    prop = boost::get<STAR_1_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : { STAR_2_PROPERTY    prop = boost::get<STAR_2_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: { SUPERNOVA_PROPERTY prop = boost::get<SUPERNOVA_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: { COMPANION_PROPERTY prop = boost::get<COMPANION_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;

        default:                                                                                                        // unexpected stellar property type
            // the only ways this can happen are if someone added a stellar type property (into ANY_PROPERTY_TYPE)
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without a
            // warning - this is a code defect, so we flag it as an error and that will result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing property type or find and fix the code defect.

            THROW_ERROR(ERROR::UNEXPECTED_STELLAR_PROPERTY_TYPE);                                                       // throw error
    }

    switch (property) {
        case ANY_STAR_PROPERTY::AGE:                                                value = Age();                                                  break;
        case ANY_STAR_PROPERTY::ANGULAR_MOMENTUM:                                   value = AngularMomentum();                                      break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_FIXED:                               value = CalculateBindingEnergy(OPTIONS->CommonEnvelopeLambda()); break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_NANJING:                             value = CalculateBindingEnergy(CalculateLambdaNanjing());       break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE:                           value = CalculateBindingEnergy(CalculateLambdaLoveridge());     break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS:                     value = CalculateBindingEnergy(CalculateLambdaLoveridge(m_Mass - m_CoreMass, true));  break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_KRUCKOW:                             value = CalculateBindingEnergy(CalculateLambdaKruckow());       break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_CONVECTIVE_ENVELOPE:                 value = CalculateConvectiveEnvelopeBindingEnergy(CalculateConvectiveEnvelopeLambdaPicker(CalculateConvectiveEnvelopeMass()));                   break;
        case ANY_STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE:               value = CHonMS();                                               break;
        case ANY_STAR_PROPERTY::CO_CORE_MASS:                                       value = COCoreMass();                                           break;
        case ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION:           value = SN_COCoreMassAtCOFormation();                           break;
        case ANY_STAR_PROPERTY::CONVECTIVE_ENV_MASS:                                std::tie(value, std::ignore) = CalculateConvectiveEnvelopeMass();  break;
        case ANY_STAR_PROPERTY::CORE_MASS:                                          value = CoreMass();                                             break;
        case ANY_STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION:              value = SN_CoreMassAtCOFormation();                             break;
        case ANY_STAR_PROPERTY::CORE_RADIUS_AT_COMPACT_OBJECT_FORMATION:            value = SN_CoreRadiusAtCOFormation();                           break; 
        case ANY_STAR_PROPERTY::DRAWN_KICK_MAGNITUDE:                               value = SN_DrawnKickMagnitude();                                break;
        case ANY_STAR_PROPERTY::DOMINANT_MASS_LOSS_RATE:                            value = DominantMassLossRate();                                 break;
        case ANY_STAR_PROPERTY::DT:                                                 value = Dt();                                                   break;
        case ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE:                                value = CalculateDynamicalTimescale();                          break;
        case ANY_STAR_PROPERTY::ECCENTRIC_ANOMALY:                                  value = SN_EccentricAnomaly();                                  break;
        case ANY_STAR_PROPERTY::ENV_MASS:                                           value = Mass() - CoreMass();                                    break;
        case ANY_STAR_PROPERTY::ERROR:                                              value = Error();                                                break;
        case ANY_STAR_PROPERTY::EVOL_STATUS:                                        value = EvolutionStatus();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_AIC:                                    value = ExperiencedAIC();                                       break;
        case ANY_STAR_PROPERTY::EXPERIENCED_HeSD:                                   value = ExperiencedHeSD();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_CCSN:                                   value = ExperiencedCCSN();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_ECSN:                                   value = ExperiencedECSN();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_PISN:                                   value = ExperiencedPISN();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_PPISN:                                  value = ExperiencedPPISN();                                     break;
        case ANY_STAR_PROPERTY::EXPERIENCED_SNIA:                                   value = ExperiencedSNIA();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_SN_TYPE:                                value = ExperiencedSN_Type();                                   break;
        case ANY_STAR_PROPERTY::EXPERIENCED_USSN:                                   value = ExperiencedUSSN();                                      break;
        case ANY_STAR_PROPERTY::FALLBACK_FRACTION:                                  value = SN_FallbackFraction();                                  break;
        case ANY_STAR_PROPERTY::MASS_TRANSFER_DONOR_HISTORY:                        value = GetMassTransferDonorHistoryString();                    break;
        case ANY_STAR_PROPERTY::HE_CORE_MASS:                                       value = HeCoreMass();                                           break;
        case ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION:           value = SN_HeCoreMassAtCOFormation();                           break;
        case ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_CORE:                              value = HeliumAbundanceCore();                                  break;
        case ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_SURFACE:                           value = HeliumAbundanceSurface();                               break;
        case ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_CORE:                            value = HydrogenAbundanceCore();                                break;
        case ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_SURFACE:                         value = HydrogenAbundanceSurface();                             break;
        case ANY_STAR_PROPERTY::IS_HYDROGEN_POOR:                                   value = SN_IsHydrogenPoor();                                    break;
        case ANY_STAR_PROPERTY::ID:                                                 value = ObjectId();                                             break;
        case ANY_STAR_PROPERTY::INITIAL_HELIUM_ABUNDANCE:                           value = CalculateInitialHeliumAbundance();                      break;
        case ANY_STAR_PROPERTY::INITIAL_HYDROGEN_ABUNDANCE:                         value = CalculateInitialHydrogenAbundance();                    break;
        case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE:                               value = InitialStellarType();                                   break;
        case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE_NAME:                          value = STELLAR_TYPE_LABEL.at(InitialStellarType());            break;
        case ANY_STAR_PROPERTY::IS_AIC:                                             value = IsAIC();                                                break;
        case ANY_STAR_PROPERTY::IS_CCSN:                                            value = IsCCSN();                                               break;
        case ANY_STAR_PROPERTY::IS_HeSD:                                            value = IsHeSD();                                               break;
        case ANY_STAR_PROPERTY::IS_ECSN:                                            value = IsECSN();                                               break;
        case ANY_STAR_PROPERTY::IS_PISN:                                            value = IsPISN();                                               break;
        case ANY_STAR_PROPERTY::IS_PPISN:                                           value = IsPPISN();                                              break;
        case ANY_STAR_PROPERTY::IS_SNIA:                                            value = IsSNIA();                                               break;
        case ANY_STAR_PROPERTY::IS_USSN:                                            value = IsUSSN();                                               break;
        case ANY_STAR_PROPERTY::KICK_MAGNITUDE:                                     value = SN_KickMagnitude();                                     break;
        case ANY_STAR_PROPERTY::LAMBDA_CONVECTIVE_ENVELOPE:                         value = CalculateConvectiveEnvelopeLambdaPicker(CalculateConvectiveEnvelopeMass());              break;
        case ANY_STAR_PROPERTY::LAMBDA_DEWI:                                        value = CalculateLambdaDewi();                                  break;
        case ANY_STAR_PROPERTY::LAMBDA_FIXED:                                       value = OPTIONS->CommonEnvelopeLambda();                        break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW:                                     value = CalculateLambdaKruckow();                               break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM:                              value = CalculateLambdaKruckow(m_Radius, -1.0);                 break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE:                              value = CalculateLambdaKruckow(m_Radius, -4.0 / 5.0);           break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_TOP:                                 value = CalculateLambdaKruckow(m_Radius, -2.0 / 3.0);           break;
        case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE:                                   value = CalculateLambdaLoveridge(m_Mass - m_CoreMass, false);   break;
        case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE_WINDS:                             value = CalculateLambdaLoveridge(m_Mass - m_CoreMass, true);    break;
        case ANY_STAR_PROPERTY::LAMBDA_NANJING:                                     value = CalculateLambdaNanjing();                               break;
        case ANY_STAR_PROPERTY::LBV_PHASE_FLAG:                                     value = LBV_PhaseFlag();                                        break;
        case ANY_STAR_PROPERTY::LUMINOSITY:                                         value = Luminosity();                                           break;
        case ANY_STAR_PROPERTY::MASS:                                               value = Mass();                                                 break;
        case ANY_STAR_PROPERTY::MASS_0:                                             value = Mass0();                                                break;
        case ANY_STAR_PROPERTY::MDOT:                                               value = Mdot();                                                 break;
        case ANY_STAR_PROPERTY::MEAN_ANOMALY:                                       value = SN_MeanAnomaly();                                       break;
        case ANY_STAR_PROPERTY::METALLICITY:                                        value = Metallicity();                                          break;
        case ANY_STAR_PROPERTY::MOMENT_OF_INERTIA:                                  value = CalculateMomentOfInertia();                             break;
        case ANY_STAR_PROPERTY::MZAMS:                                              value = MZAMS();                                                break;
        case ANY_STAR_PROPERTY::OMEGA:                                              value = Omega() / SECONDS_IN_YEAR;                              break;
        case ANY_STAR_PROPERTY::OMEGA_BREAK:                                        value = OmegaBreak() / SECONDS_IN_YEAR;                         break;
        case ANY_STAR_PROPERTY::OMEGA_ZAMS:                                         value = OmegaZAMS() / SECONDS_IN_YEAR;                          break;
        case ANY_STAR_PROPERTY::PULSAR_MAGNETIC_FIELD:                              value = PulsarMagneticField();                                  break;
        case ANY_STAR_PROPERTY::PULSAR_SPIN_DOWN_RATE:                              value = PulsarSpinDownRate();                                   break;
        case ANY_STAR_PROPERTY::PULSAR_SPIN_FREQUENCY:                              value = PulsarSpinFrequency();                                  break;
        case ANY_STAR_PROPERTY::PULSAR_SPIN_PERIOD:                                 value = PulsarSpinPeriod();                                     break;
        case ANY_STAR_PROPERTY::PULSAR_BIRTH_PERIOD:                                value = PulsarBirthPeriod();                                    break;
        case ANY_STAR_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE:                        value = PulsarBirthSpinDownRate();                              break;
        case ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE:                         value = CalculateRadialExpansionTimescale();                    break;
        case ANY_STAR_PROPERTY::RADIUS:                                             value = Radius();                                               break;
        case ANY_STAR_PROPERTY::RANDOM_SEED:                                        value = RandomSeed();                                           break;
        case ANY_STAR_PROPERTY::ROCKET_KICK_MAGNITUDE:                              value = SN_RocketKickMagnitude();                               break;
        case ANY_STAR_PROPERTY::ROCKET_KICK_PHI:                                    value = SN_RocketKickPhi();                                     break;
        case ANY_STAR_PROPERTY::ROCKET_KICK_THETA:                                  value = SN_RocketKickTheta();                                   break;
        case ANY_STAR_PROPERTY::RZAMS:                                              value = RZAMS();                                                break;
        case ANY_STAR_PROPERTY::SN_TYPE:                                            value = SN_Type();                                              break;
        case ANY_STAR_PROPERTY::SPEED:                                              value = Speed();												break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE:                                       value = StellarType();                                          break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE_NAME:                                  value = STELLAR_TYPE_LABEL.at(StellarType());                   break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE_PREV:                                  value = StellarTypePrev();                                      break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE_PREV_NAME:                             value = STELLAR_TYPE_LABEL.at(StellarTypePrev());               break;
        case ANY_STAR_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER:             value = SN_KickMagnitudeRandom();                               break;
        case ANY_STAR_PROPERTY::SUPERNOVA_PHI:                                      value = SN_Phi();                                               break;
        case ANY_STAR_PROPERTY::SUPERNOVA_THETA:                                    value = SN_Theta();                                             break;
        case ANY_STAR_PROPERTY::TEMPERATURE:                                        value = Temperature() * TSOL;                                   break;
        case ANY_STAR_PROPERTY::THERMAL_TIMESCALE:                                  value = CalculateThermalTimescale();                            break;
        case ANY_STAR_PROPERTY::TIME:                                               value = Time();                                                 break;
        case ANY_STAR_PROPERTY::TIMESCALE_MS:                                       value = Timescale(TIMESCALE::tMS);                              break;
        case ANY_STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION:             value = SN_TotalMassAtCOFormation();                            break;
        case ANY_STAR_PROPERTY::TOTAL_RADIUS_AT_COMPACT_OBJECT_FORMATION:           value = SN_TotalRadiusAtCOFormation();                          break;
        case ANY_STAR_PROPERTY::TRUE_ANOMALY:                                       value = SN_TrueAnomaly();                                       break;
        case ANY_STAR_PROPERTY::TZAMS:                                              value = TZAMS() * TSOL;                                         break;
        case ANY_STAR_PROPERTY::VELOCITY_X:                                         value = VelocityX();											break;
        case ANY_STAR_PROPERTY::VELOCITY_Y:                                         value = VelocityY();											break;
        case ANY_STAR_PROPERTY::VELOCITY_Z:                                         value = VelocityZ();											break;
        case ANY_STAR_PROPERTY::ZETA_HURLEY:                                        value = CalculateZetaAdiabaticHurley2002(m_CoreMass);           break;
        case ANY_STAR_PROPERTY::ZETA_HURLEY_HE:                                     value = CalculateZetaAdiabaticHurley2002(m_HeCoreMass);         break;
        case ANY_STAR_PROPERTY::ZETA_SOBERMAN:                                      value = CalculateZetaAdiabaticSPH(m_CoreMass);                  break;
        case ANY_STAR_PROPERTY::ZETA_SOBERMAN_HE:                                   value = CalculateZetaAdiabaticSPH(m_HeCoreMass);                break;

        default:                                                                                                        // unexpected stellar property
            // the only ways this can happen are if someone added a stellar property (into ANY_STAR_PROPERTY),
            // or allowed users to specify a stellar property (via the logfile definitions file), and it isn't
            // accounted for in this code.  We should not default here, with or without a warning - this is a
            // code defect, so we flag it as an error and that will result in termination of the evolution of
            // the star or binary.
            // The correct fix for this is to add code for the missing property, or prevent it from being 
            // specified in the logfile definitions file.

            THROW_ERROR(ERROR::UNEXPECTED_STELLAR_PROPERTY);                                                            // throw error
    }

    return value;
}


/*
 * Determine the value of the requested property of the star (parameter p_Property)
 *
 * The property is a boost variant variable, and is one of the following types:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *      BINARY_PROPERTY         - property of the binary
 *      PROGRAM_OPTION          - program option
 *
 * This function calls the appropriate helper function to retrieve the value.
 * 
 * This function handles properties of type:
 * 
 *    STAR_PROPERTY, PROGRAM_OPTION
 * 
 * only - anything else will result in an error being thrown and the evolution of the star (or binary)
 * terminated.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  
 *
 *
 * COMPAS_VARIABLE PropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseStar::PropertyValue(const T_ANY_PROPERTY p_Property) const {

    COMPAS_VARIABLE value = 0.0;                                                                                                // default property value

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {                                                          // which property type?

        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY:  value = StellarPropertyValue(p_Property); break;                              // star property
        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION: value = OPTIONS->OptionValue(p_Property); break;                              // program option

        default:                                                                                                                // unexpected property type
            // the only ways this can happen are if someone added a stellar type property (into ANY_PROPERTY_TYPE)
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without a
            // warning - this is a code defect, so we flag it as an error and that will result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing property type or find and fix the code defect.

            THROW_ERROR(ERROR::UNEXPECTED_STELLAR_PROPERTY_TYPE);                                                               // throw error
    }

    return value;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                     COEFFICIENT AND CONSTANT CALCULATIONS ETC.                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate a(n) coefficients
 *
 * a(n) coefficients depend on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster.  This function isn't
 * called too often, but the pattern is the same for others that are called many, many times.

 *
 *
 * void CalculateAnCoefficients(DBL_VECTOR &p_AnCoefficients,
 *                              DBL_VECTOR &p_LConstants,
 *                              DBL_VECTOR &p_RConstants,
 *                              DBL_VECTOR &p_GammaConstants)
 *
 * @param   [IN/OUT]    p_AnCoefficients        a(n) coefficients - calculated here
 * @param   [IN/OUT]    p_LConstants            Luminosity constants - calculated here
 * @param   [IN/OUT]    p_RConstants            Radius constants - calculated here
 * @param   [IN/OUT]    p_GammaConstants        Gamma constants - calculated here
 */
void BaseStar::CalculateAnCoefficients(DBL_VECTOR &p_AnCoefficients,
                                       DBL_VECTOR &p_LConstants,
                                       DBL_VECTOR &p_RConstants,
                                       DBL_VECTOR &p_GammaConstants) {
#define a p_AnCoefficients                                                          // for convenience and readability - undefined at end of function
#define index    coeff.first                                                        // for convenience and readability - undefined at end of function
#define coeff(x) coeff.second[AB_TCoeff::x]                                         // for convenience and readability - undefined at end of function
#define LConstants(x) p_LConstants[static_cast<int>(L_CONSTANTS::x)]                // for convenience and readability - undefined at end of function
#define RConstants(x) p_RConstants[static_cast<int>(R_CONSTANTS::x)]                // for convenience and readability - undefined at end of function
#define GammaConstants(x) p_GammaConstants[static_cast<int>(GAMMA_CONSTANTS::x)]    // for convenience and readability - undefined at end of function

    double Z     = m_Metallicity;
    double xi    = LogMetallicityXiHurley();
    double sigma = LogMetallicitySigma();

    // pow() is slow - use multiplication
    // do these calculations once only - and esp. outside the loop
    double xi_2  = xi * xi;
    double xi_3  = xi * xi_2;
    double xi_4  = xi_2 * xi_2;

    // calculate initial values for a(n) coefficients
    a.push_back(0.0);           // this is a dummy entry - so our index is the same as that in Hurley et al. 2000 (we just ignore the zeroeth entry)
    for (auto coeff: A_COEFF) {
        a.push_back(coeff(ALPHA) + (coeff(BETA) * xi) + (coeff(GAMMA) * xi_2) + (coeff(ETA) * xi_3) + (coeff(MU) * xi_4));
    }

    // Special cases - see Hurley et al. 2000

    a[11] *= a[14];
    a[12] *= a[14];
    a[17] = PPOW(10.0, max((0.097 - (0.1072 * (sigma + 3.0))), max(0.097, min(0.1461, (0.1461 + (0.1237 * (sigma + 2.0)))))));
    a[18] *= a[20];
    a[19] *= a[20];
    a[29] = PPOW(a[29], (a[32]));
    a[33] = min(1.4, 1.5135 + (0.3769 * xi));
    a[33] = max(0.6355 - (0.4192 * xi), max(1.25, a[33]));
    a[42] = min(1.25, max(1.1, a[42]));
    a[44] = min(1.3, max(0.45, a[44]));
    a[49] = max(a[49], 0.145);
    a[50] = min(a[50], (0.306 + (0.053 * xi)));
    a[51] = min(a[51], (0.3625 + (0.062 * xi)));
    a[52] = max(a[52], 0.9);
    a[52] = (utils::Compare(Z, 0.01) > 0) ? min(a[52], 1.0) : a[52];
    a[53] = max(a[53], 1.0);
    a[53] = (utils::Compare(Z, 0.01) > 0) ? min(a[53], 1.1) : a[53];
    a[57] = min(1.4, a[57]);
    a[57] = max((0.6355 - (0.4192 * xi)), max(1.25, a[57]));
    a[62] = max(0.065, a[62]);
    a[63] = (utils::Compare(Z, 0.004) < 0) ? min(0.055, a[63]) : a[63];
    a[66] = max(a[66], min(1.6, -0.308 - (1.046 * xi)));
    a[66] = max(0.8, min(0.8 - (2.0 * xi), a[66]));
    a[68] = max(0.9, min(a[68], 1.0));

    // Need bAlphaR - calculate it now
    RConstants(B_ALPHA_R) = (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61]));                            // Hurley et al. 2000, eq 21a (wrong in the arxiv version - says = a59*M**(a61))

    // Continue special cases

    a[64] = (utils::Compare(a[68], a[66]) > 0) ? RConstants(B_ALPHA_R) : max(0.091, min(0.121, a[64]));
    a[68] = min(a[68], a[66]);
    a[72] = (utils::Compare(Z, 0.01) > 0) ? max(a[72], 0.95) : a[72];
    a[74] = max(1.4, min(a[74], 1.6));
    a[75] = max(1.0, min(a[75], 1.27));
    a[75] = max(a[75], 0.6355 - (0.4192 * xi));
    a[76] = max(a[76], -0.1015564 - (0.2161264 * xi) - (0.05182516 * xi_2));
    a[77] = max((-0.3868776 - (0.5457078 * xi) - (0.1463472 * xi_2)), min(0.0, a[77]));
    a[78] = max(0.0, min(a[78], 7.454 + (9.046 * xi)));
    a[79] = min(a[79], max(2.0, -13.3 - (18.6 * xi)));
    a[80] = max(0.0585542, a[80]);
    a[81] = min(1.5, max(0.4, a[81]));

    LConstants(B_ALPHA_L)   = (a[45] + (a[46] * PPOW(2.0, a[48]))) / (PPOW(2.0, 0.4) + (a[47] * PPOW(2.0, 1.9)));   // Hurley et al. 2000, eq 19a
    LConstants(B_BETA_L)    = max(0.0, (a[54] - (a[55] * PPOW(a[57], a[56]))));                                     // Hurley et al. 2000, eq 20
    LConstants(B_DELTA_L)   = min((a[34] / PPOW(a[33], a[35])), (a[36] / PPOW(a[33], a[37])));                      // Hurley et al. 2000, eq 16

    RConstants(C_ALPHA_R)   = (a[58] * PPOW(a[67], a[60])) / (a[59] + PPOW(a[67], a[61]));                          // Hurley et al. 2000, eq 21a (wrong in the arxiv version)
    RConstants(B_BETA_R)    = (a[69] * 8.0 * M_SQRT2) / (a[70] + PPOW(2.0, a[71]));                                 // Hurley et al. 2000, eq 22a
    RConstants(C_BETA_R)    = (a[69] * 16384.0) / (a[70] + PPOW(16.0, a[71]));                                      // Hurley et al. 2000, eq 22a
    RConstants(B_DELTA_R)   = (a[38] + a[39] * 8.0 * M_SQRT2) / (a[40] * 8.0 + PPOW(2.0, a[41])) - 1.0;             // Hurley et al. 2000, eq 17

    GammaConstants(B_GAMMA) = max(0.0, a[76] + (a[77] * PPOW((1.0 - a[78]), a[79])));                               // Hurley et al. 2000, eq 23 and discussion immediately following - max() confirmed in BSE Fortran code
    GammaConstants(C_GAMMA) = (utils::Compare(a[75], 1.0) <= 0) ? GammaConstants(B_GAMMA) : a[80];                  // Hurley et al. 2000, eq 23 and discussion immediately following - <= 1.0 confirmed in BSE Fortran code

#undef GammaConstants
#undef RConstants
#undef LConstants
#undef coeff
#undef index
#undef a
}


/*
 * Calculate b(n) coefficients
 *
 * b(n) coefficients depend on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster.  This function isn't
 * called too often, but the pattern is the same for others that are called many, many times.
 *
 *
 * void CalculateBnCoefficients(DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN/OUT]    p_BnCoefficients        b(n) coefficients - calculated here
 */
void BaseStar::CalculateBnCoefficients(DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define index    coeff.first                                            // for convenience and readability - undefined at end of function
#define coeff(x) coeff.second[AB_TCoeff::x]                             // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function


    double Z     = m_Metallicity;
    double xi    = LogMetallicityXiHurley();
    double sigma = LogMetallicitySigma();
    double rho   = LogMetallicityRho();

    // pow() is slow - use multiplication
    // do these calculations once only - and esp. outside the loop
    double xi_2  = xi * xi;
    double xi_3  = xi * xi_2;
    double xi_4  = xi_2 * xi_2;
    double xi_5  = xi * xi_4;

    double rho_2 = rho * rho;
    double rho_3 = rho * rho_2;

    b.push_back(0.0);           // this is a dummy entry for b(n) coefficients - so our index is the same as that in Hurley et al. 2000
    for (auto coeff: B_COEFF) {
        b.push_back(coeff(ALPHA) + (coeff(BETA) * xi) + (coeff(GAMMA) * xi_2) + (coeff(ETA) * xi_3) + (coeff(MU) * xi_4));
    }

    // Special Cases - see Hurley et al. 2000

    b[1] = min(0.54, b[1]);
    b[2] = PPOW(10.0, (-4.6739 - (0.9394 * sigma)));
    b[2] = min(max(b[2], (-0.04167 + (55.67 * Z))), (0.4771 - (9329.21 * PPOW(Z, 2.94))));
    b[3] = PPOW(10.0, max(-0.1451, (-2.2794 - (1.5175 * sigma) - (0.254 * sigma * sigma))));
    b[3] = (utils::Compare(Z, 0.004) > 0) ? max(b[3], 0.7307 + (14265.1 * PPOW(Z, 3.395))) : b[3];
    b[4] += 0.1231572 * xi_5;
    b[6] += 0.01640687 * xi_5;
    b[11] = b[11] * b[11];
    b[13] = b[13] * b[13];
    b[14] = PPOW(b[14], b[15]);
    b[16] = PPOW(b[16], b[15]);
    b[17] = (utils::Compare(xi, -1.0) > 0) ? 1.0 - (0.3880523 * PPOW((xi + 1.0), 2.862149)) : 1.0;
    b[24] = PPOW(b[24], b[28]);
    b[26] = 5.0 - (0.09138012 * PPOW(Z, -0.3671407));
    b[27] = PPOW(b[27], (2.0 * b[28]));
    b[31] = PPOW(b[31], b[33]);
    b[34] = PPOW(b[34], b[33]);
    b[36] = b[36] * b[36] * b[36] * b[36];
    b[37] = 4.0 * b[37];
    b[38] = b[38] * b[38] * b[38] * b[38];
    b[40] = max(b[40], 1.0);
    b[41] = PPOW(b[41], b[42]);
    b[44] = b[44] * b[44] * b[44] * b[44] * b[44];
    b[45] = utils::Compare(rho, 0.0) <= 0 ? 1.0 : 1.0 - ((2.47162 * rho) - (5.401682 * rho_2) + (3.247361 * rho_3));
    b[46] = -1.0 * b[46] * log10(massCutoffs(MHeF) / massCutoffs(MFGB));
    b[47] = (1.127733 * rho) + (0.2344416 * rho_2) - (0.3793726 * rho_3);
    b[51] -= 0.1343798 * xi_5;
    b[53] += 0.4426929 * xi_5;
    b[55] = min((0.99164 - (743.123 * PPOW(Z, 2.83))), b[55]);
    b[56] += 0.1140142 * xi_5;
    b[57] -= 0.01308728 * xi_5;

#undef massCutoffs
#undef coeff
#undef index
#undef b
}


/*
 * Calculate all alpha-like metallicity dependent luminosity coefficients
 *
 * Luminosity coefficients depend on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster.  This function isn't
 * called too often, but the pattern is the same for others that are called many, many times.
 *
 *
 * void CalculateLCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_LCoefficients)
 *
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - xi in Hurley et al. 2000
 * @param   [IN/OUT]    p_LCoefficients         Luminosity coefficients - calculated here
 */
void BaseStar::CalculateLCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_LCoefficients) {
#define index    coeff.first                // for convenience and readability - undefined at end of function
#define coeff(x) coeff.second[LR_TCoeff::x] // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication
    // do these calculations once only - and esp. outside the loop
    double xi   = p_LogMetallicityXi;
    double xi_2 = xi * xi;
    double xi_3 = xi * xi_2;
    double xi_4 = xi_2 * xi_2;

    // iterate over Luminosity coefficients constants L_COEFF (see constants.h)
    // These are from table 1 in Tout et al. 1996
    // Each row (indexed by 'index') defines the coefficients of the 5 terms (coefficients 'a', 'b', 'c', 'd' & 'e')
    for(auto coeff: L_COEFF) {
        p_LCoefficients.push_back(coeff(a) + (coeff(b) * xi) + (coeff(c) * xi_2) + (coeff(d) * xi_3) + (coeff(e) * xi_4));
    }

#undef coeff
#undef index
}


/*
 * Calculate all alpha-like metallicity dependent radius coefficients
 *
 * Radius coefficients depend on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster.  This function isn't
 * called too often, but the pattern is the same for others that are called many, many times.
 *
 *
 * void CalculateRCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_RCoefficients)
 *
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - xi in Hurley et al. 2000
 * @param   [IN/OUT]    p_LCoefficients         Radius coefficients - calculated here
 */
void BaseStar::CalculateRCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_RCoefficients) {
#define index    coeff.first                // for convenience and readability - undefined at end of function
#define coeff(x) coeff.second[LR_TCoeff::x] // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication
    // do these calculations once only - and esp. outside the loop
    double xi   = p_LogMetallicityXi;
    double xi_2 = xi * xi;
    double xi_3 = xi * xi_2;
    double xi_4 = xi_2 * xi_2;

    // iterate over Radius coefficients constants R_COEFF (see constants.h)
    // These are from table 2 in Tout et al. 1996
    // Each row (indexed by 'index') defines the coefficients of the 5 terms (coefficients 'a', 'b', 'c', 'd' & 'e')
    for(auto coeff: R_COEFF) {
        p_RCoefficients.push_back(coeff(a) + (coeff(b) * xi) + (coeff(c) * xi_2) + (coeff(d) * xi_3) + (coeff(e) * xi_4));
    }

#undef coeff
#undef index
}


/*
 * Calculate the constant alpha1
 *
 * Hurley et al, 2000, just after eq 49
 *
 * Alpha1 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateAlpha1()
 *
 * @return                                      Metallicity dependent constant alpha1
 */
double BaseStar::CalculateAlpha1() const {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double LHeI_MHeF = (b[11] + (b[12] * PPOW(massCutoffs(MHeF), 3.8))) / (b[13] + (massCutoffs(MHeF) * massCutoffs(MHeF)));
    return ((b[9] * PPOW(massCutoffs(MHeF), b[10])) - LHeI_MHeF) / LHeI_MHeF;

#undef massCutoffs
#undef b
}


/*
 * Calculate the constant alpha3
 *
 * Hurley et al. 2000, just after eq 56
 *
 * Alpha3 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateAlpha3()
 *
 * @return                                      Metallicity dependent constant alpha3
 */
double BaseStar::CalculateAlpha3() const {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double LBAGB = (b[31] + (b[32] * PPOW(massCutoffs(MHeF), (b[33] + 1.8)))) / (b[34] + PPOW(massCutoffs(MHeF), b[33]));
    return ((b[29] * PPOW(massCutoffs(MHeF), b[30])) - LBAGB) / LBAGB;

#undef massCutoffs
#undef b
}


/*
 * Calculate the constant alpha4
 *
 * Hurley et al. 2000, just after eq 57
 *
 * Alpha4 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateAlpha4()
 *
 * @return                                      Metallicity dependent constant alpha4
 */
double BaseStar::CalculateAlpha4() const {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double MHeF      = massCutoffs(MHeF);
    double MHeF_5    = MHeF * MHeF * MHeF * MHeF * MHeF;    // pow() is slow - use multiplication
    double tBGB_MHeF = CalculateLifetimeToBGB(MHeF);        // tBGB for mass M = MHeF
    double tHe_MHeF  = tBGB_MHeF * (b[41] * PPOW(MHeF, b[42]) + b[43] * MHeF_5) / (b[44] + MHeF_5);
    
    return ((tHe_MHeF - b[39]) / b[39]);

#undef massCutoffs
#undef b
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate mass cutoffs:
 *
 *   MHook: the metallicity dependent mass above which a hook appears on the MS
 *   MHeF : the metallicity dependent maximum initial mass for which He ignites degenerately in the He Flash
 *   MFGB : the metallicity dependent maximum mass at which He ignites degenerately on the First Giant Branch (FGB)
 *
 * Mass cutoffs depend on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster.  This function isn't
 * called too often, but the pattern is the same for others that are called many, many times.
 *
 *
 * void CalculateMassCutoffs(const double p_Metallicity, const double p_LogMetallicityXi, DBL_VECTOR &p_MassCutoffs)
 *
 * @param   [IN]        p_Metallicity           Metallicity Z (Z = 0.02 = Zsol)
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - xi in Hurley et al. 2000
 * @param   [IN/OUT]    p_MassCutoffs           Mass cutoffs - calculated here
 */
void BaseStar::CalculateMassCutoffs(const double p_Metallicity, const double p_LogMetallicityXi, DBL_VECTOR &p_MassCutoffs) {
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double xi_2 = p_LogMetallicityXi * p_LogMetallicityXi;                          // pow() is slow - use multiplication

    massCutoffs(MHook) = 1.0185 + (0.16015 * p_LogMetallicityXi) + (0.0892 * xi_2); // MHook - Hurley et al. 2000, eq 1
    massCutoffs(MHeF)  = 1.995 + (0.25 * p_LogMetallicityXi) + (0.087 * xi_2);      // MHeF - Hurley et al. 2000, eq 2

    double top         = 13.048 * PPOW((p_Metallicity / ZSOL_HURLEY), 0.06);
    double bottom      = 1.0 + (0.0012 * PPOW((ZSOL_HURLEY / p_Metallicity), 1.27));
    massCutoffs(MFGB)  = top / bottom;                                              // MFGB - Hurley et al. 2000, eq 3

    massCutoffs(MCHE)  = 100.0;                                                     // MCHE - Mandel/Butler - CHE calculation

#undef massCutoffs
}


/*
 * Calculate the parameter x for the Giant Branch
 *
 * X depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 * Hybrid of b5 and b7 from Hurley et al. 2000
 * Hurley et al. 2000, eq 47
 *
 *
 * double CalculateGBRadiusXExponent()
 *
 * @return                                      'x' exponent to which Radius depends on Mass (at constant Luminosity)- 'x' in Hurley et al. 2000, eq 47
 */
double BaseStar::CalculateGBRadiusXExponent() const {

    // pow()is slow - use multiplication
    double xi   = LogMetallicityXiHurley();
    double xi_2 = xi * xi;
    double xi_3 = xi_2 * xi;
    double xi_4 = xi_2 * xi_2;

    return 0.30406 + (0.0805 * xi) + (0.0897 * xi_2) + (0.0878 * xi_3) + (0.0222 * xi_4);   // Hurley et al. 2000, eq 47
}


/*
 * Calculate the perturbation parameter b
 *
 * Hurley et al. 2000, eq 103
 *
 *
 * double CalculatePerturbationB(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Perturbation parameter b
 */
double BaseStar::CalculatePerturbationB(const double p_Mass) const {
    return 0.002 * max(1.0, (2.5 / p_Mass));
}


/*
 * Calculate the perturbation parameter c
 *
 * Hurley et al. 2000, eq 104
 *
 *
 * double CalculatePerturbationC(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Perturbation parameter c
 */
double BaseStar::CalculatePerturbationC(double p_Mass) const {
    return 0.006 * max(1.0, (2.5 / p_Mass));
}


/*
 * Calculate the perturbation parameter s
 *
 * Hurley et al. 2000, eq 101
 *
 *
 * double CalculatePerturbationS(const double p_Mass)
 *
 * @param   [IN]    p_Mu                        Perturbation parameter mu
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Perturbation parameter s
 */
double BaseStar::CalculatePerturbationS(const double p_Mu, const double p_Mass) const {

    double b      = CalculatePerturbationB(p_Mass);
    double b_3    = b * b * b;                      // pow() is slow - use multiplication
    double mu_b_3 = p_Mu * p_Mu * p_Mu / b_3;       // calculate once, use many times...

    return ((1.0 + b_3) * mu_b_3) / (1.0 + mu_b_3);
}


/*
 * Calculate the perturbation parameter q
 *
 * Hurley et al. 2000, eq 105
 *
 *
 * double CalculatePerturbationQ(const double p_Radius, const double p_Rc)
 *
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @param   [IN]    p_Rc                        Radius that the remnant would have if the star immediately lost its envelope (in Rsol)
 * @return                                      Perturbation parameter q
 */
double BaseStar::CalculatePerturbationQ(const double p_Radius, const double p_Rc) const {
    return log(p_Radius / p_Rc); // really is natural log
}


/*
 * Calculate the perturbation parameter r
 *
 * Hurley et al. 2000, eq 102
 *
 *
 * double CalculatePerturbationR(const double p_Mu, const double p_Mass, const double p_Radius, const double p_Rc)
 *
 * @param   [IN]    p_Mu                        Perturbation parameter mu
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @param   [IN]    p_Rc                        Radius that the remnant would have if the star immediately lost its envelope (in Rsol)
 * @return                                      Perturbation parameter r
 */
double BaseStar::CalculatePerturbationR(const double p_Mu, const double p_Mass, const double p_Radius, const double p_Rc) const {

    double r = 0.0;

    if (utils::Compare(p_Mu, 0.0) > 0 && utils::Compare(p_Radius, p_Rc) > 0) {  // only if mu > 0 and radius is larger than core radius, otherwise r = 0 and perturbed radius = core radius

        double c      = CalculatePerturbationC(p_Mass);
        double c_3    = c * c * c;                                              // pow() is slow - use multiplication
        double mu_c_3 = p_Mu * p_Mu * p_Mu / c_3;                               // calculate once

        double q        = CalculatePerturbationQ(p_Radius, p_Rc);
        double exponent = min((0.1 / q), (-14.0 / log10(p_Mu)));                // Hurley et al. 2000 is just 0.1 / q, but the Hurley sse code does this (`rpertf()` in `zfuncs.f`) - no explanation.

        r = ((1.0 + c_3) * mu_c_3 * PPOW((p_Mu), exponent)) / ((1.0 + mu_c_3));
    }

    return r;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                LAMBDA CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Proposed fit for the common envelope lambda parameter
 * Kruckow et al. 2016 (arXiv:1610.04417), fig 1
 *
 * Spectrum fit to the region bounded by the upper and lower limits as shown in Kruckow+ 2016.
 * Fit as presented in Vigna-Gomez et al. 2018 (arXiv:1805.07974)
 *
 *
 * double CalculateLambdaKruckow(const double p_Radius, const double p_Alpha)
 *
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @param   [IN]    p_Alpha                     Power
 * @return                                      Common envelope lambda parameter
 */
double BaseStar::CalculateLambdaKruckow(const double p_Radius, const double p_Alpha) const {

	double alpha = max(-2.0 / 3.0, min(-1.0, p_Alpha));             // clamp alpha to [-1.0, -2/3]

	return 1600.0 * PPOW(0.00125, -alpha) * PPOW(p_Radius, alpha);
}


/* 
 * Wrapper function to return Nanjing lambda based on options
 * 
 * 
 * double BaseStar::CalculateLambdaNanjing()
 * 
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateLambdaNanjing() const {

    double lambda = 0.0;                                                                                        // return value

    double mass = m_MZAMS;
    if (OPTIONS->CommonEnvelopeLambdaNanjingUseRejuvenatedMass()) mass = m_Mass0;                               // use rejuvenated mass to calculate lambda instead of true birth mass
    
    if (OPTIONS->CommonEnvelopeLambdaNanjingEnhanced()) {                                                       // if using enhanced Nanjing lambdas
        STELLAR_POPULATION pop = utils::Compare(m_Metallicity, LAMBDA_NANJING_ZLIMIT) < 0 ? STELLAR_POPULATION::POPULATION_II : STELLAR_POPULATION::POPULATION_I;
        if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMass()) {
            if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMetallicity()) {
                lambda = BaseStar::CalculateMassAndZInterpolatedLambdaNanjing(mass, m_Metallicity);
            }
            else {
                lambda = BaseStar::CalculateMassInterpolatedLambdaNanjing(mass, pop);
            }
        }
        else {
            int massIndex = BaseStar::FindLambdaNanjingNearestMassIndex(mass);                                  // do not interpolate in mass, so need to use nearest mass bin
            if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMetallicity()) {
                lambda = BaseStar::CalculateZInterpolatedLambdaNanjing(m_Metallicity, massIndex);
            }
            else {
                lambda = BaseStar::CalculateLambdaNanjingEnhanced(massIndex, pop);
            }
        }
    }
    else { 
        lambda = CalculateLambdaNanjingStarTrack(mass, m_Metallicity);
    }

    return lambda;
}


/* 
 * Calculate mass- and metallicity-interpolated Nanjing lambda
 * 
 * 
 * double BaseStar::CalculateMassAndZInterpolatedLambdaNanjing(const double p_Mass, const double p_Z)
 * 
 * @param   [IN]    p_Mass                      Mass / Msun to evaluate lambda with
 * @param   [IN]    p_Z                         Metallicity
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateMassAndZInterpolatedLambdaNanjing(const double p_Mass, const double p_Z) const {

    double lambda = 0.0;                                                                                        // return value

    if (utils::Compare(m_Metallicity, LAMBDA_NANJING_POPII_Z) < 0) {
        lambda = BaseStar::CalculateMassInterpolatedLambdaNanjing(p_Mass, STELLAR_POPULATION::POPULATION_II);   // use lambda for pop. II metallicity
    }
    else if (utils::Compare(m_Metallicity, LAMBDA_NANJING_POPI_Z) > 0) {
        lambda = BaseStar::CalculateMassInterpolatedLambdaNanjing(p_Mass, STELLAR_POPULATION::POPULATION_I);    // use lambda for pop. I metallicity
    }
    else {                                                                                                      // linear interpolation in logZ between pop. I and pop. II metallicities
        double lambdaLow = BaseStar::CalculateMassInterpolatedLambdaNanjing(p_Mass, STELLAR_POPULATION::POPULATION_II);
        double lambdaUp  = BaseStar::CalculateMassInterpolatedLambdaNanjing(p_Mass, STELLAR_POPULATION::POPULATION_I);
        lambda           = lambdaLow + (m_Log10Metallicity - LAMBDA_NANJING_POPII_LOGZ) / (LAMBDA_NANJING_POPI_LOGZ - LAMBDA_NANJING_POPII_LOGZ) * (lambdaUp - lambdaLow);
    }

    return lambda;
}


/* 
 * Interpolate Nanjing lambda in mass for a given metallicity
 * 
 * 
 * double BaseStar::CalculateMassInterpolatedLambdaNanjing(const double p_Mass, const int p_StellarPop)
 * 
 * @param   [IN]    p_StellarPop                Stellar population (POP_I or POP_II)
 * @param   [IN]    p_Mass                      Mass / Msun to evaluate lambda with
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateMassInterpolatedLambdaNanjing(const double p_Mass, const STELLAR_POPULATION p_StellarPop) const {

    double lambda = 0.0;                                                                                        // return value

    INT_VECTOR ind = utils::BinarySearch(NANJING_MASSES, p_Mass);
    int low        = ind[0];
    int up         = ind[1];

    if ( (low < 0) && (up >= 0) ) {                                                                             // mass below range calculated by Xu & Li (2010)
        lambda = CalculateLambdaNanjingEnhanced(0, p_StellarPop);                                               // use lambda for minimum mass
    }
    else if ( (low >= 0) && (up < 0) ) {                                                                        // mass above range calculated by Xu & Li (2010)
        lambda = CalculateLambdaNanjingEnhanced(NANJING_MASSES.size() - 1, p_StellarPop);                       // use lambda for maximum mass
    }
    else if (low == up) {                                                                                       // mass is exactly equal to the mass of a model evolved by Xu & Li (2010)
        lambda = CalculateLambdaNanjingEnhanced(low, p_StellarPop);
    }
    else {                                                                                                      // linear interpolation between upper and lower mass bins
        double lambdaLow = CalculateLambdaNanjingEnhanced(low, p_StellarPop);
        double lambdaUp  = CalculateLambdaNanjingEnhanced(up, p_StellarPop);
        lambda           = lambdaLow + (p_Mass - NANJING_MASSES[low]) / (NANJING_MASSES[up] - NANJING_MASSES[low]) * (lambdaUp - lambdaLow);
    }

    return lambda;
}


/* 
 * Interpolate Nanjing lambda in metallicity for a given mass
 * 
 * 
 * double BaseStar::alculateZInterpolatedLambdaNanjing(const double p_Z, const int p_MassInd)
 * 
 * @param   [IN]    p_Z                         Metallicity
 * @param   [IN]    p_MassIndex                 Index specifying donor mass (see NANJING_MASSES in constants.h)
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateZInterpolatedLambdaNanjing(const double p_Z, const int p_MassIndex) const {

    double lambda = 0.0;                                                                                        // return value
    
    if (utils::Compare(m_Metallicity, LAMBDA_NANJING_POPII_Z) < 0) {
        lambda = CalculateLambdaNanjingEnhanced(p_MassIndex, STELLAR_POPULATION::POPULATION_II);                // use lambda for pop. II metallicity
    }
    else if (utils::Compare(m_Metallicity, LAMBDA_NANJING_POPI_Z) > 0) {
        lambda = CalculateLambdaNanjingEnhanced(p_MassIndex, STELLAR_POPULATION::POPULATION_I);                 // use lambda for pop. I metallicity
    }
    else {                                                                                                      // linear interpolation in logZ between pop. I and pop. II metallicities
        double lambdaLow = CalculateLambdaNanjingEnhanced(p_MassIndex, STELLAR_POPULATION::POPULATION_II);
        double lambdaUp  = CalculateLambdaNanjingEnhanced(p_MassIndex, STELLAR_POPULATION::POPULATION_I);
        lambda           = lambdaLow + (m_Log10Metallicity - LAMBDA_NANJING_POPII_LOGZ) / (LAMBDA_NANJING_POPI_LOGZ - LAMBDA_NANJING_POPII_LOGZ) * (lambdaUp - lambdaLow);
    }

    return lambda;
}


/* 
 * Returns index in NANJING_MASSES corresponding to nearest mass model computed by Xu & Li (2010)
 * 
 * 
 * double BaseStar::FindLambdaNanjingNearestMassIndex(const double p_Mass)
 * 
 * @param   [IN]    p_Mass                      Mass
 * @return                                      Index in NANJING_MASSES
 */ 
double BaseStar::FindLambdaNanjingNearestMassIndex(const double p_Mass) const {

    if (p_Mass < NANJING_MASSES_MIDPOINTS[0]) return 0.0;                                                       // if M < midpoint of lowermost bin, use lambda for the 1 Msun model
    
    if (p_Mass >= NANJING_MASSES_MIDPOINTS.back()) return NANJING_MASSES.size() - 1.0;                          // if M >= midpoint of uppermost bin, use lambda for the 100 Msun model

    return utils::BinarySearch(NANJING_MASSES_MIDPOINTS, p_Mass)[1];                                            // search for upper and lower mass bin edges
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 ZETA CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate zeta, the adiabatic donor radial response to mass loss
 *
 * double BaseStar::CalculateZetaAdiabatic() 
 *
 * @return                                      Adiabatic exponent zeta = dlnR/dlnM
 */
double BaseStar::CalculateZetaAdiabatic() { 
                                                                                
    double zetaStar = 0.0;                                                              // return value

    switch (OPTIONS->StellarZetaPrescription()) {

        case ZETA_PRESCRIPTION::SOBERMAN: 
        case ZETA_PRESCRIPTION::HURLEY:   
        case ZETA_PRESCRIPTION::ARBITRARY:
            zetaStar = CalculateZetaConstantsByEnvelope(OPTIONS->StellarZetaPrescription());
            break;

        default:                                                                        // unknown prescription
            // the only way this can happen is if someone added a ZETA_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ZETA_PRESCRIPTION);                              // throw error
    }

    return zetaStar;
}


/*
 * Calculate the Adiabatic Exponent per Hurley et al. 2002
 *
 *
 * double CalculateZetaAdiabaticHurley2002(const double p_CoreMass) const
 *
 * @param   [IN]    p_CoreMass                  Core mass of the star (Msol)
 * @return                                      Adiabatic exponent
 */
double BaseStar::CalculateZetaAdiabaticHurley2002(const double p_CoreMass) const{
    
    if (utils::Compare(p_CoreMass, m_Mass) >= 0) return 0.0;                        // if the object is all core, the calculation is meaningless

    double m = p_CoreMass / m_Mass;
    double x = BaseStar::CalculateGBRadiusXExponent();                              // x from Hurley et al 2000, Eq. 47 - Depends on composition

    return -x + (2.0 * m * m * m * m * m);
}


/*
 * Calculate the Adiabatic Exponent per Soberman, Phinney, vdHeuvel 1997
 *
 *
 * double CalculateZetaAdiabaticSPH(const double p_CoreMass) const
 *
 * @param   [IN]    p_CoreMass                  Core mass of the star (Msol)
 * @return                                      Adiabatic exponent
 */
double BaseStar::CalculateZetaAdiabaticSPH(const double p_CoreMass) const {
    
    if (utils::Compare(p_CoreMass, m_Mass) >= 0) return 0.0;                        // if the object is all core, the calculation is meaningless (and would result in division by zero)

    double m           = p_CoreMass / m_Mass;                                       // eq (57) Soberman, Phinney, vdHeuvel (1997)
    double oneMinusM   = 1.0 - m;
    double oneMinusM_6 = oneMinusM * oneMinusM * oneMinusM * oneMinusM * oneMinusM * oneMinusM;

    return ((2.0 / 3.0) * m / oneMinusM) - ((1.0 / 3.0) * (oneMinusM / (1.0 + (m + m)))) - (0.03 * m) + (0.2 * m / (1.0 + (1.0 / oneMinusM_6))); // eq (61) Soberman, Phinney, vdHeuvel (1997)
}


/*
 * Calculate the critical mass ratio for unstable mass transfer
 *
 * double BaseStar::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, const double p_massTransferEfficiencyBeta)
 *
 * @param   [IN]    p_AccretorIsDegenerate       Whether or not the accretor is a degenerate star
 * @param   [IN]    p_massTransferEfficiencyBeta Mass transfer accretion efficiency
 * @return                                       Critical mass ratio
 */
double BaseStar::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, const double p_massTransferEfficiencyBeta) {
    
        double qCrit = 0.0;                                                                 // return value

        switch (OPTIONS->QCritPrescription()) {

            case QCRIT_PRESCRIPTION::NONE:                                                  // a safe placeholder; really, we should not be calling this function if there is the prescription is "NONE"
                qCrit = 0.0;
                break;
                
            case QCRIT_PRESCRIPTION::GE: 
            case QCRIT_PRESCRIPTION::GE_IC:
                qCrit = CalculateCriticalMassRatioGeEtAl(OPTIONS->QCritPrescription(), p_massTransferEfficiencyBeta);   
                break;

            case QCRIT_PRESCRIPTION::CLAEYS:
                qCrit = CalculateCriticalMassRatioClaeys14(p_AccretorIsDegenerate);
                break;

            case QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK:
                qCrit = CalculateCriticalMassRatioHurleyHjellmingWebbink();
                break;
        
            default:                                                                        // unknown prescription
                // the only way this can happen is if someone added a QCRIT_PRESCRIPTION
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_QCRIT_PRESCRIPTION);                             // throw error
        }

        return qCrit;
}

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate luminosity at ZAMS (in Lsol)
 * Tout et al. 1996, eq 1
 *
 *
 * double CalculateLuminosityAtZAMS(const double p_MZAMS)
 *
 * @param   [IN]    p_MZAMS                     Zero age main sequence mass in Msol
 * @return                                      Luminosity in Lsol (LZAMS)
 */
double BaseStar::CalculateLuminosityAtZAMS(const double p_MZAMS) {
#define coeff(x) m_LCoefficients[static_cast<int>(L_Coeff::x)]   // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication where it makes sense
    // sqrt() is much faster than pow()
    double m_0_5 = std::sqrt(p_MZAMS);
    double m_2   = p_MZAMS * p_MZAMS;
    double m_3   = m_2 * p_MZAMS;
    double m_5   = m_3 * m_2;
    double m_5_5 = m_5 * m_0_5;
    double m_7   = m_5 * m_2;
    double m_8   = m_7 * p_MZAMS;
    double m_9_5 = m_8 * p_MZAMS * m_0_5;
    double m_11  = m_8 * m_3;

    double top    = (coeff(ALPHA) * m_5_5) + (coeff(BETA) * m_11);
    double bottom = (coeff(GAMMA) + m_3) + (coeff(DELTA) * m_5) + (coeff(EPSILON) * m_7) + (coeff(ZETA) * m_8) + (coeff(ETA) * m_9_5);

    return top / bottom;

#undef coeff
}


/*
 * Calculate luminosity at the base of the Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 56
 *
 *
 * double CalculateLuminosityAtBAGB(double p_Mass)
 *
 * @param   [IN]    p_Mass                      (Effective) mass in Msol
 * @return                                      Luminosity at BAGB in Lsol
 */
double BaseStar::CalculateLuminosityAtBAGB(double p_Mass) const {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    return (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0)
            ? (b[29] * PPOW(p_Mass, b[30])) / (1.0 + (m_Alpha3 * exp(15.0 * (p_Mass - massCutoffs(MHeF)))))
            : (b[31] + (b[32] * PPOW(p_Mass, (b[33] + 1.8)))) / (b[34] + PPOW(p_Mass, b[33]));

#undef massCutoffs
#undef b
}


/*
 * Calculate luminosity for a given core mass, used for AGB stars
 *
 * Hurley et al. 2000, eq 37
 *
 *
 * double CalculateLuminosityGivenCoreMass(const double p_CoreMass)
 *
 * @param   [IN]    p_CoreMass                  Core mass in Msol
 * @return                                      Luminosity in Lsol
 */
double BaseStar::CalculateLuminosityGivenCoreMass(const double p_CoreMass) const {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function
    return min((gbParams(B) * PPOW(p_CoreMass, gbParams(q))), (gbParams(D) * PPOW(p_CoreMass, gbParams(p))));
#undef gbParams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate radius at ZAMS in units of Rsol
 * Tout et al. 1996, eq 2
 *
 *
 * double CalculateRadiusAtZAMS(const double p_MZAMS)
 *
 * @param   [IN]    p_MZAMS                     Zero age main sequence mass in Msol
 * @return                                      Radius in units of Rsol (RZAMS)
 * Uses class member m_RCoefficients as radius coefficients
 */
double BaseStar::CalculateRadiusAtZAMS(const double p_MZAMS) const {
#define coeff(x) m_RCoefficients[static_cast<int>(R_Coeff::x)]  // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication where it makes sense
    // std::sqrt() is much faster than pow()
    double m_0_5  = std::sqrt(p_MZAMS);
    double m_2    = p_MZAMS * p_MZAMS;
    double m_2_5  = m_2 * m_0_5;
    double m_6    = m_2 * m_2 * m_2;
    double m_6_5  = m_6 * m_0_5;
    double m_8    = m_6 * m_2;
    double m_8_5  = m_8 * m_0_5;
    double m_11   = m_8 * m_2 * p_MZAMS;
    double m_18_5 = m_6 * m_6 * m_6 * m_0_5;
    double m_19   = m_11 * m_8;
    double m_19_5 = m_19 * m_0_5;

    double top    = (coeff(THETA) * m_2_5 ) + (coeff(IOTA) * m_6_5) + (coeff(KAPPA) * m_11) + (coeff(LAMBDA) * m_19) + (coeff(MU) * m_19_5);
    double bottom =  coeff(NU) + (coeff(XI) * m_2) + (coeff(OMICRON) * m_8_5 ) + m_18_5 + (coeff(PI) * m_19_5);

    return top / bottom;

#undef coeff
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the maximum core mass
 *
 * Hurley et al. 2000, eq 89
 *
 *
 * double CalculateMaximumCoreMass(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Maximum core mass in Msol (McMax)
 */
double BaseStar::CalculateMaximumCoreMass(const double p_Mass) const {
    return min(((1.45 * p_Mass) - 0.31), p_Mass);
}


/*
 * Calculate the initial convective envelope mass
 *
 * Hurley et al. 2000, just after eq 111
 *
 *
 * double CalculateInitialEnvelopeMass_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      ZAMS envelope mass - Menv in Hurley et al. 2000
 */
double BaseStar::CalculateInitialEnvelopeMass_Static(const double p_Mass) {

    double envMass = 0.0;

    if (utils::Compare(p_Mass, 0.35) < 0) {                 // star is fully convective so Menv = M
        envMass = p_Mass;
    }
    else if (utils::Compare(p_Mass, 1.25) < 0) {
        double brackets = (1.25 - p_Mass) / 0.9;            // pow() is slow - use multiplication
        envMass         = 0.35 * brackets * brackets;       // Hurley et al. 2000, just after eq 111
    }

    return envMass;
}


/*
 * Calculate mass loss rate enhancement for rapidly rotating stars
 *
 * Langer 1998 (https://ui.adsabs.harvard.edu/abs/1998A%26A...329..551L/abstract) eq 3
 * 
 * The exponent originally comes from Bjorkman & Cassinelli 1993 (https://ui.adsabs.harvard.edu/abs/1993ApJ...409..429B/abstract),
 * based on a fit to data from Friend & Abbott 1986 (https://ui.adsabs.harvard.edu/abs/1986ApJ...311..701F/abstract) 
 *
 * double CalculateMassLossRateEnhancementRotation()
 *
 * @return                                      Mass loss enhancement factor for rapidly rotating stars
 */
double BaseStar::CalculateMassLossRateEnhancementRotation() {
    return OPTIONS->EnableRotationallyEnhancedMassLoss() ? PPOW((1.0 - Omega() / OmegaBreak()), -0.43) : 1.0;   // default is no enhancement
}


/*
 * Calculate the mass loss rate on the AGB based on the Mira pulsation period (P0)
 *
 * Hurley et al. 2000, just after eq 106 (from Vassiliadis and Wood 1993)
 *
 *
 * double CalculateMassLossRateVassiliadisWood()
 *
 * @return                                      Mass loss rate on AGB in Msol per year
 */
double BaseStar::CalculateMassLossRateVassiliadisWood() const {

    double logP0      = min(3.3, (-2.07 - (0.9 * log10(m_Mass)) + (1.94 * log10(m_Radius))));
    double P0         = PPOW(10.0, (logP0));        // in their fortran code, Hurley et al. take P0 to be min(p0, 2000.0): implemented here as a minimum power
    double logMdot_VW = -11.4 + (0.0125 * (P0 - 100.0 * max((m_Mass - 2.5), 0.0)));
    double Mdot_VW    = PPOW(10.0, (logMdot_VW));

    return min(Mdot_VW, (1.36E-9 * m_Luminosity));
}


/*
 * Calculate mass loss rate on the GB and beyond
 *
 * Hurley et al. 2000, eq 106 (based on a prescription taken from Kudritzki and Reimers 1978)
 *
 *
 * double CalculateMassLossRateKudritzkiReimers()
 *
 * @return                                      Kudritzki and Reimers mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateKudritzkiReimers() const {
    return 4.0E-13 * (MASS_LOSS_ETA * m_Luminosity * m_Radius / m_Mass);    // shouldn't be eta squared like in paper!
}


/*
 * Calculate the mass-loss rate for massive stars (L > 4000 L_sol) using the
 * Nieuwenhuijzen & de Jager 1990 prescription, modified by a metallicity
 * dependent factor (Kudritzki et al 1989).
 *
 * Hurley et al. 2000, just after eq 106
 *
 *
 * double CalculateMassLossRateNieuwenhuijzenDeJager()
 *
 * @return                                      Nieuwenhuijzen & de Jager mass-loss rate for massive stars (in Msol yr^-1)
 */
double BaseStar::CalculateMassLossRateNieuwenhuijzenDeJager() const {
    
    double rate = 0.0;
    
    if (utils::Compare(m_Luminosity, NJ_MINIMUM_LUMINOSITY) > 0) {          // check for minimum luminosity
        double smoothTaper = min(1.0, (m_Luminosity - 4000.0) / 500.0);     // smooth taper between no mass loss and mass loss
        rate = std::sqrt((m_Metallicity / ZSOL_HURLEY)) * smoothTaper * 9.6E-15 * PPOW(m_Radius, 0.81) * PPOW(m_Luminosity, 1.24) * PPOW(m_Mass, 0.16);
    } else {
        rate = 0.0;
    }

    return rate;
}


/*
 * Calculate the opacity for this star (e.g., to determine the Eddington luminosity)
 *
 * See text surrounding Equation 6 in Bjorklund et al. 2022 (https://arxiv.org/abs/2203.08218)
 *
 * double CalculateOpacity_Static(const double p_HeliumAbundanceSurface)
 * 
 * @return                                      Opacity in SI units (m^2/kg)
 *
 */
double BaseStar::CalculateOpacity_Static(const double p_HeliumAbundanceSurface) {

    const double iHe   = 2.0;                                           // Helium ionisation state - For O stars, doubly ionised helium
    double YHe         = p_HeliumAbundanceSurface;                      // Star's surface helium abundance 

    double kappa_e_cgs = 0.4 * (1.0 + iHe * YHe) / (1.0 + 4.0 * YHe);   // cgs units cm^2/g
    double kappa_e_SI  = kappa_e_cgs * OPACITY_CGS_TO_SI;               // Convert to SI units - m^2/kg

    return kappa_e_SI;
}


/*
 * Calculate the Eddington Luminosity L_edd for this star
 *
 * See e.g., above Equation 6 in Bjorklund et al. 2022 (https://arxiv.org/abs/2203.08218)
 * 
 * double CalculateEddingtonLuminosity_Static(const double p_Mass, const double p_HeliumAbundanceSurface)
 * 
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_HeliumAbundanceSurface    Helium abundance
 * @return                                      Eddington luminosity in solar luminosities
 *
 */
double BaseStar::CalculateEddingtonLuminosity_Static(const double p_Mass, const double p_HeliumAbundanceSurface) {

    double kappa_SI = CalculateOpacity_Static(p_HeliumAbundanceSurface);    // Determine opacity
    
    double top      = 4.0 * M_PI * G * C * p_Mass * MSOL_TO_KG;
    double bot      = kappa_SI;
    double L_Edd    = top / bot;

    return L_Edd;
}


/*
 * Calculate the Eddington factor (L/L_Edd) as required by CalculateMassLossRateBjorklund
 * see text surrounding Equation 6 in https://arxiv.org/abs/2203.08218
 * 
 * 
 * double CalculateMassLossRateBjorklundEddingtonFactor()
 *
 * @return                                      Eddington factor
 */
double BaseStar::CalculateMassLossRateBjorklundEddingtonFactor() const {

    const double YHe = 0.1;                                                 // assumed constant by Bjorklund et al.
        
    double Ledd      = CalculateEddingtonLuminosity_Static(m_Mass, YHe);    // W

    double LoverLedd = (m_Luminosity * LSOLW) / Ledd;                       // Dimensionless

    return LoverLedd;
}


/*
 * Calculate the mass loss rate for massive OB stars according to the prescription from Bjorklund et al. 2022
 * See Equation 7 and surrounding text in https://arxiv.org/abs/2203.08218
 * 
 * This prescription is calibrated to the following ranges:
 * 10^4.5 < L / Lsol < 10^6
 * 15,000 < Teff / K < 50,000
 * 15 < M / Msol < 80
 * Z Zsol, Z_LMC = 0.5 * Zsol and Z_SMC = 0.2 * Zsol, with Zsol = 0.014 
 * 
 *
 * double CalculateMassLossRateOBBjorklund2022()
 *
 * @return                                      Bjorklund mass-loss rate for massive stars (in Msol yr^-1)
 */
double BaseStar::CalculateMassLossRateOBBjorklund2022() const {

    double Gamma   = CalculateMassLossRateBjorklundEddingtonFactor();

    double logZ    = log10(m_Metallicity / 0.014);
    double logL    = log10(m_Luminosity / 1.0E6);
    double Teff    = m_Temperature * TSOL;                                  // convert effective temperature to Kelvin
    double logTeff = log10(Teff / 45000.0);           

    double Meff    = m_Mass * (1.0 - Gamma);
    double logMeff = log10(Meff / 45.0);

    // Constants, q depends on logTeff
    const double constC = -5.52;
    const double m      = 2.39;
    const double n      = -1.48;
    const double p      = 2.12;
    double q            = 0.75 - (1.87 * logTeff);

    // Equation 7 in Bjorklund et al. 2022
    double logMdot = constC + (m * logL) + (n * logMeff) + (p * logTeff) + (q * logZ);

    return PPOW(10.0, logMdot);
}


/*
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
 *
 * Sets class member variable m_LBVphaseFlag if necessary
 * 
 *  
 * double CalculateMassLossRateLBV(const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription)
 *
 * @param   [IN]    p_LBVprescription           Which LBV prescription to use
 * @return                                      LBV-like mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateLBV(const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription) {

    double rate = 0.0;                                                                                                          // default return value

    double HDlimitfactor = m_Radius * std::sqrt(m_Luminosity) * 1.0E-5;                                                         // calculate factor by which the star is above the HD limit
    if ((utils::Compare(m_Luminosity, LBV_LUMINOSITY_LIMIT_STARTRACK) > 0) && (utils::Compare(HDlimitfactor, 1.0) > 0)) {       // check if luminous blue variable
		m_LBVphaseFlag         = true;                                                                                          // mark the star as LBV
        m_DominantMassLossRate = MASS_LOSS_TYPE::LBV;                                                                           // set the dominant mass loss rate
        
        switch (p_LBVprescription) {                                                                                            // decide which LBV prescription to use

            case LBV_MASS_LOSS_PRESCRIPTION::ZERO:
                rate = 0.0;
                break;
            case LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD:
            case LBV_MASS_LOSS_PRESCRIPTION::HURLEY:
                rate = CalculateMassLossRateLBVHurley(HDlimitfactor);
                break;
            case LBV_MASS_LOSS_PRESCRIPTION::BELCZYNSKI:
                rate = CalculateMassLossRateLBVBelczynski();
                break;

            default:                                                                                                            // unknown prescription
                // the only ways this can happen are if someone added an LBV_MASS_LOSS_PRESCRIPTION
                // and it isn't accounted for in this code, or if there is a defect in the code that causes
                // this function to be called with a bad parameter.  We should not default here, with or without
                // a warning.
                // We are here because the user chose a prescription this code doesn't account for, or as a result
                // of a code defect, and either of those should be flagged as an error and result in termination of
                // the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option, or find and fix the code defect.

                THROW_ERROR(ERROR::UNKNOWN_LBV_MASS_LOSS_PRESCRIPTION);                                                         // throw error
        }
    } else {
        rate = 0.0;                                                                                                             // no winds if it isn't an LBV star!
    }

    return rate;
}


/*
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
 *
 * Hurley+ 2000 Section 7.1 a few equation after Eq. 106 (Equation not labelled)
 * 
 *
 * double CalculateMassLossRateLBVHurley(const double p_HDlimitfactor)
 *
 * @param   [IN]    p_HDlimitfactor             Factor by which star is above Humphreys-Davidson limit
 * @return                                      LBV-like mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateLBVHurley(const double p_HDlimitfactor) const {
    double v = p_HDlimitfactor - 1.0;
    return 0.1 * v * v * v * ((m_Luminosity / 6.0E5) - 1.0);
}


/*
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
 *
 * Belczynski et al. 2010, eq 8
 * 
 *
 * double CalculateMassLossRateLBVBelczynski()
 *
 * @return                                      LBV-like mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateLBVBelczynski() const {
    return OPTIONS->LuminousBlueVariableFactor() * 1.0E-4;
}


/*
 * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0).
 *
 * Hurley et al. 2000, just after eq 106 (taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998)
 *
 * Note that the reduction of this formula is imposed in order to match the observed number of black holes in binaries (Hurley et al 2000)
 *
 *
 * double CalculateMassLossRateWolfRayet(const double p_Mu)
 *
 * @param   [IN]    p_Mu                        Small envelope parameter (see Hurley et al. 2000, eq 97 & 98)
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayet(const double p_Mu) const {
    // In Hurley's fortran code there is a parameter here 'hewind' which by default is 1.0 -
    // can be set to zero to disable this particular part of winds. We instead opt for all winds on or off.
    double rate = 0.0;
    if (utils::Compare(p_Mu, 1.0) < 0) {
        rate = PPOW(m_Luminosity, 1.5) * (1.0 - p_Mu) * 1.0E-13;
    }
    return rate;
}


/*
 * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0).
 *
 * Belczynski et al. 2010, eq 9 (taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998)
 *
 * Note that the reduction of this formula is imposed in order to match the observed number of black holes in binaries (Hurley et al 2000)
 *
 *
 * double CalculateMassLossRateWolfRayetZDependent(const double p_Mu)
 *
 * @param   [IN]    p_Mu                        Small envelope parameter (see Hurley et al. 2000, eq 97 & 98)
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayetZDependent(const double p_Mu) const {
    // I think StarTrack may still do something different here,
    // there are references to Hamann & Koesterke 1998 and Vink and de Koter 2005
    // TW - Haven't seen StarTrack but I think H&K gives the original equation and V&dK gives the Z dependence
    double rate = 0.0;
    if (utils::Compare(p_Mu, 1.0) < 0) {
        rate = 1.0E-13 * PPOW(m_Luminosity, 1.5) * PPOW(m_Metallicity / ZSOL_ANDERS, 0.86) * (1.0 - p_Mu);
    }
    return rate;
}


/*
 * Calculate mass loss rate for massive OB stars using the Vink et al 2001 prescription
 *
 * Vink et al. 2001, eqs 24 & 25
 * Belczynski et al. 2010, eqs 6 & 7
 *
 *
 * double CalculateMassLossRateOBVink2001()
 *
 * @return                                      Mass loss rate for hot OB stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateOBVink2001() const {

    double rate = 0.0;                                                                                          // default return value

    double teff = m_Temperature * TSOL;  

    if (utils::Compare(teff, VINK_MASS_LOSS_MINIMUM_TEMP) >= 0 && utils::Compare(teff, VINK_MASS_LOSS_BISTABILITY_TEMP) <= 0) {
        double v = 1.3;                                                                                                // v_inf/v_esc
        v        = v * PPOW(m_Metallicity / ZSOL_ANDERS, OPTIONS->ScaleTerminalWindVelocityWithMetallicityPower());    // Scale Vinf with metallicity

        double logMdotOB = -6.688 +
                           (2.210 * log10(m_Luminosity / 1.0E5)) -
                           (1.339 * log10(m_Mass / 30.0)) -
                           (1.601 * log10(v / 2.0)) +
                           (0.85  * LogMetallicityXiAnders()) +
                           (1.07  * log10(teff / 20000.0));

        rate = PPOW(10.0, logMdotOB);

    }
    else if (utils::Compare(teff, VINK_MASS_LOSS_BISTABILITY_TEMP) > 0) {
        SHOW_WARN_IF(utils::Compare(teff, VINK_MASS_LOSS_MAXIMUM_TEMP) > 0, ERROR::HIGH_TEFF_WINDS);            // show warning if winds being used outside comfort zone

        double v = 2.6;                                                                                                // v_inf/v_esc
        v        = v * PPOW(m_Metallicity / ZSOL_ANDERS, OPTIONS->ScaleTerminalWindVelocityWithMetallicityPower());    // Scale Vinf with metallicity

        double logMdotOB = -6.697 +
                           (2.194 * log10(m_Luminosity / 1.0E5)) -
                           (1.313 * log10(m_Mass / 30.0)) -
                           (1.226 * log10(v / 2.0)) +
                           (0.85  * LogMetallicityXiAnders()) +
                           (0.933 * log10(teff / 40000.0)) -
                           (10.92 * log10(teff / 40000.0) * log10(teff/40000.0));

        rate = PPOW(10.0, logMdotOB);

    }
    else {
        SHOW_WARN(ERROR::LOW_TEFF_WINDS, "Mass Loss Rate = 0.0");                                               // too cold to use winds - show warning.
    }

    return rate;
}


/*
 * Calculate mass loss rate for massive OB stars using the Vink+Sander 2021 update
 * https://arxiv.org/pdf/2103.12736.pdf
 * features two bi-stability jumps, at T1 and T2
 * offset = {"cold":-5.99,"inter":-6.688,"hot":-6.697}
 *
 * 
 * double CalculateMassLossRateOBVinkSander2021()
 *
 * @return                                            Mass loss rate for hot OB stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateOBVinkSander2021() const {

    double rate = 0.0;                                                                                          // default return value

    const double zExp2001 = 0.85;
    const double zExp     = 0.42;

    double teff    = m_Temperature * TSOL;  
    double Gamma   = EDDINGTON_PARAMETER_FACTOR * m_Luminosity / m_Mass;
    double charrho = -14.94 + (3.1857 * Gamma) + (zExp * LogMetallicityXiAnders());
    double T2      = ( 61.2 + (2.59 * charrho) ) * 1000.0;                                                      // typically around 25000.0, higher jump first as in Vink python recipe
    double T1      = ( 100.0 + (6.0 * charrho) ) * 1000.0;                                                      // typically around 20000.0, has similar behavior when fixed

    double logL5  = log10(m_Luminosity / 1.0E5);
    double logM30 = log10(m_Mass / 30.0);
    double logT40 = log10(teff / 40000.0);
    double logT20 = log10(teff / 20000.0);

    if (utils::Compare(teff, VINK_MASS_LOSS_MINIMUM_TEMP) >= 0 && utils::Compare(teff, T1) <= 0) {

        double V         = 0.7;                                                                                 // v_inf/v_esc
        double logMdotOB = -5.99 +
                           (2.210 * logL5) -
                           (1.339 * logM30) -
                           (1.601 * log10(V / 2.0)) +
                           (zExp2001 * LogMetallicityXiAnders()) +
                           (1.07  * logT20);

        rate = PPOW(10.0, logMdotOB);
    }
    else if (utils::Compare(teff, T1) > 0 && utils::Compare(teff, T2) <= 0) {
        SHOW_WARN_IF(utils::Compare(teff, VINK_MASS_LOSS_MAXIMUM_TEMP) > 0, ERROR::HIGH_TEFF_WINDS);            // show warning if winds being used outside comfort zone

        double V         = 1.3;                                                                                 // v_inf/v_esc
        double logMdotOB = -6.688 +
                           (2.210 * logL5) -
                           (1.339 * logM30) -
                           (1.601 * log10(V / 2.0)) +
                           (zExp2001  * LogMetallicityXiAnders()) +
                           (1.07  * logT20);

        rate = PPOW(10.0, logMdotOB);
    }
    else if (utils::Compare(teff, T2) > 0) {
        SHOW_WARN_IF(utils::Compare(teff, VINK_MASS_LOSS_MAXIMUM_TEMP) > 0, ERROR::HIGH_TEFF_WINDS);            // show warning if winds being used outside comfort zone

        double V         = 2.6;                                                                                 // v_inf/v_esc
        double logMdotOB = -6.697 +
                           (2.194 * logL5) -
                           (1.313 * logM30) -
                           (1.226 * log10(V / 2.0)) +
                           (zExp  * LogMetallicityXiAnders()) +
                           (0.933 * logT40) -
                           (10.92 * logT40 * logT40);

        rate = PPOW(10.0, logMdotOB);
    }
    else {
        SHOW_WARN(ERROR::LOW_TEFF_WINDS, "Mass Loss Rate = 0.0");                                               // too cold to use winds - show warning.
    }

    return rate;
}


/*
 * Calculate mass loss rate for massive OB stars using the Krticka+ 2018 prescription
 *
 * https://arxiv.org/pdf/1712.03321.pdf
 *
 * 
 * double CalculateMassLossRateOBKrticka2018()
 *
 * @return                                      Mass loss rate for hot OB stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateOBKrticka2018() const {
    
    double logMdot = -5.70 + 0.50 * LogMetallicityXiAsplund() + (1.61 - 0.12 * LogMetallicityXiAsplund()) * log10(m_Luminosity / 1.0E6);

    return PPOW(10.0, logMdot);
}


/*
 * Calculate mass loss rate for RSG stars using the Beasor+2020 prescription
 *
 * https://arxiv.org/pdf/2001.07222.pdf eq 4.
 * 
 * fit corrected slightly in Decin 2023, eq E.1 
 * https://arxiv.org/pdf/2303.09385.pdf
 * 
 * corrected again by Beasor+2023, https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.2460B/abstract
 *
 * 
 * double CalculateMassLossRateRSGBeasor2020()
 *
 * @return                                      Mass loss rate for RSG stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateRSGBeasor2020() const {

    double logMdot = (-21.5 - 0.15 * m_MZAMS) + (3.6 * log10(m_Luminosity));        // further correction by Beasor+

    return PPOW(10.0, logMdot);
}


/*
 * Calculate mass loss rate for RSG stars using the Decin2023 prescription
 * 
 *  https://arxiv.org/pdf/2303.09385.pdf eq 6.
 *
 * 
 * double CalculateMassLossRateRSGDecin2023()
 *
 * @return                                      Mass loss rate for RSG stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateRSGDecin2023() const {
    return PPOW(10.0, -20.63 - 0.16 * m_MZAMS + 3.47 * log10(m_Luminosity));
}


/*
 * Calculate mass loss rate for RSG stars using the Yang 2023 prescription
 *  Third order polynomial in log Luminosity.
 *  https://arxiv.org/pdf/2303.09385.pdf eq 6.
 *
 * 
 * double CalculateMassLossRateRSGYang2023()
 *
 * @return                                      Mass loss rate for RSG stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateRSGYang2023() const {

    double logL    = log10(m_Luminosity);
    double logL_2  = logL * logL;
    double logMdot = 0.45 * logL_2 * logL - 5.26 * logL_2 + 20.93 * logL - 34.56;

    return PPOW(10.0, logMdot);
}


/*
 * Calculate mass loss rate for RSG stars using the Kee + 2021 prescription
 *
 * https://arxiv.org/pdf/2101.03070.pdf eqs 5, 13, 14, 25. 
 *
 * 
 * double CalculateMassLossRateRSGKee2021()
 *
 * @return                                      Mass loss rate for RSG stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateRSGKee2021() const {

    const double vturb    = 1.5E4;                                                                              // turbulent velocity, m/s, for a typical RSG
    const double k_b      = 1.38E-23;                                                                           // Boltzmann Constant in J K^-1
    const double sigma    = 5.67E-8;                                                                            // Stefan Boltzmann constant W m^-2 K^-4
    const double m_h      = 1.67E-27;                                                                           // mass of hydrogen in Kg
    const double kappa    = 0.01 * OPACITY_CGS_TO_SI;                                                           // Given after Eq. 16 

    double teff           = TSOL * m_Temperature;                                                               // in K
    
    double R_SI           = sqrt((m_Luminosity * LSOLW) / (4.0 * M_PI * sigma * PPOW(teff, 4.0)));
    double M_SI           = m_Mass * MSOL_TO_KG;
    double cs             = sqrt(k_b * teff / m_h);
    double gamma          = (kappa * m_Luminosity * LSOLW) / (4.0 * M_PI * G * C * M_SI);
    double vesc           = sqrt(2.0 * G * (M_SI) / (R_SI));                                                    // m/s, not vesc,eff

    double Rpmod          = G * (M_SI) * (1.0 - gamma) / (2.0 * ((cs * cs) + (vturb * vturb)));                 // modified parker radius, in m
    double rho            = (4.0 / 3.0) * (Rpmod / (kappa * (R_SI) * (R_SI))) * 
                            (exp(-(2.0 * Rpmod / (R_SI)) + (3.0 / 2.0))) / (1.0 - exp(-2.0 * Rpmod / (R_SI)));

    double MdotAnalytical = 4.0 * M_PI * rho * sqrt(cs * cs + vturb * vturb) * Rpmod * Rpmod;                   // in kg/s
    double factor         = PPOW(((vturb / 17000.0) / (vesc / 60000.0)), 1.30);                                 // non-isothermal correction factor

    return factor * MdotAnalytical * SECONDS_IN_YEAR / MSOL_TO_KG; 
}   


/*
 *  Calculate mass loss rate for RSG stars using the Vink and Sabhahit 2023 prescription
 *  A kinked function of L and M
 *  https://arxiv.org/pdf/2309.08657.pdf eqs 1 and 2
 *
 * 
 * double CalculateMassLossRateRSGVinkSabhahit2023()
 *
 * @return                                      Mass loss rate for RSG stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateRSGVinkSabhahit2023() const {

    const double logLkink = 4.6;

    double logL = log10(m_Luminosity);
    double logM = log10(m_Mass);

    double logMdot;
    if (utils::Compare(logL, logLkink) < 0) {
        logMdot = -8.0 + 0.7 * logL - 0.7 * logM;
    }
    else {
        logMdot = -24.0 + 4.77 * logL - 3.99 * logM;
    }

    return PPOW(10.0, logMdot);
}


/*
 * Calculate mass loss rate for very massive (>100 Msol) OB stars using the Bestenlehner 2020 prescription
 *
 * https://arxiv.org/pdf/2002.05168.pdf
 *
 * 
 * double CalculateMassLossRateVMSBestenlehner2020()
 *
 * @return                                      Mass loss rate for hot OB stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateVMSBestenlehner2020() const {

    const double alpha       = 0.39;                                        // CAK force multiplier
    const double logMdotZero = -4.78;                                       // from substituting LogMdotTrans and Gamma_e trans into eq 12. 

    double gamma   = EDDINGTON_PARAMETER_FACTOR * m_Luminosity / m_Mass;    // Eddington Parameter, not metallicity specific as in the publication
    double logMdot = logMdotZero + ((1.0 / alpha) + 0.5) * log10(gamma) - (((1.0 - alpha) / alpha) + 2.0) * log10(1.0 - gamma);

    return PPOW(10.0, logMdot);
}


/*
 * Calculate the mass loss rate for very massive OB stars using a fit to the Vink 2011 mass loss rate
 *
 * https://arxiv.org/pdf/1105.0556.pdf
 *
 * 
 * double CalculateMassLossRateVMSVink2011()
 *
 * @return                                      Mass loss rate for very massive stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateVMSVink2011() const {

    double rate = CalculateMassLossRateOBVink2001();

    double Gamma    = EDDINGTON_PARAMETER_FACTOR * m_Luminosity / m_Mass;                                       // Eddington Parameter, independent of surface composition

    if (utils::Compare(Gamma, 0.5) > 0) {                                                                       // apply this correction to high gamma only
        double logMdotdiff = 0.04468 + (0.3091 * Gamma) + (0.2434 * Gamma * Gamma);
        rate = PPOW(10.0, (logMdotdiff + log10(rate)));
    }

    return rate;
}


/*
 * Calculate mass loss rate for very massive stars using the Sabhahit 2023 prescription
 *
 * https://arxiv.org/pdf/2306.11785.pdf
 *
 * 
 * double CalculateMassLossRateVMSSabhahit2023()
 *
 * @return                                      Mass loss rate in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateVMSSabhahit2023() {

    double gamma       = EDDINGTON_PARAMETER_FACTOR * m_Luminosity / m_Mass;                                    // Eddington Parameter, independent of surface composition
    double Mswitch     = PPOW(m_Metallicity, -1.574) * 0.0615 + 18.10;                                          // obtained from a powerlaw fit to table 2, given teff=45kK
    double Lswitch     = PPOW(10, (-1.91 * m_Log10Metallicity + 2.36));                                         // loglinear fits to table 2 
    double Mdotswitch  = PPOW(10, (-1.86 * m_Log10Metallicity - 8.90));
    double gammaswitch = EDDINGTON_PARAMETER_FACTOR * Lswitch / Mswitch;

    double Mdot; 
    if (utils::Compare(gamma, gammaswitch) > 0) {
        Mdot = Mdotswitch * PPOW((m_Luminosity / Lswitch) , 4.77) * PPOW((m_Mass/Mswitch) , -3.99);
    }
    else {
        Mdot = CalculateMassLossRateOB(OPTIONS->OBMassLossPrescription());                                      // not in the VMS regime according to Sabhahit+ 2023, fall back to default OB mass loss prescription
    }

    return Mdot;
}


/*
 * Calculate mass loss for main sequence stars. 
 * Switches prescription based on program options. 
 *
 * 
 * double CalculateMassLossRateOB(const OB_MASS_LOSS_PRESCRIPTION p_OB_MassLossPrescription)
 *
 * @param   [IN]    p_OBMassLossPrescription    OB Mass loss prescription to use
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateOB(const OB_MASS_LOSS_PRESCRIPTION p_OB_MassLossPrescription) {

    double rate = 0.0;                                                                                          // default return value                                                      

    m_DominantMassLossRate = MASS_LOSS_TYPE::OB;                                                                // set dominant mass loss rate
    
    switch (p_OB_MassLossPrescription) {                                                                        // decide which prescription to use
        case OB_MASS_LOSS_PRESCRIPTION::ZERO         : rate = 0.0; break;
        case OB_MASS_LOSS_PRESCRIPTION::VINK2001     : rate = CalculateMassLossRateOBVink2001(); break;
        case OB_MASS_LOSS_PRESCRIPTION::VINK2021     : rate = CalculateMassLossRateOBVinkSander2021(); break;
        case OB_MASS_LOSS_PRESCRIPTION::BJORKLUND2022: rate = CalculateMassLossRateOBBjorklund2022(); break;
        case OB_MASS_LOSS_PRESCRIPTION::KRTICKA2018  : rate = CalculateMassLossRateOBKrticka2018(); break;

        default:                                                                                                // unknown prescription
            // the only ways this can happen are if someone added an OB_MASS_LOSS_PRESCRIPTION
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without
            // a warning.
            // We are here because the user chose a prescription this code doesn't account for, or as a result
            // of a code defect, and either of those should be flagged as an error and result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option, or find and fix the code defect.

            THROW_ERROR(ERROR::UNKNOWN_OB_MASS_LOSS_PRESCRIPTION);                                              // throw error
    }

    return rate;
}


/*
 * Calculate mass loss for RSG stars (Red Supergiant). 
 * Switches prescription based on program options. 
 * 
 * 
 * double CalculateMassLossRateRSG(const RSG_MASS_LOSS_PRESCRIPTION p_RSG_MassLossPrescription)
 *
 * @param   [IN]    p_RSG_MassLossPrescription  RSG Mass loss prescription to use
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateRSG(const RSG_MASS_LOSS_PRESCRIPTION p_RSG_MassLossPrescription) {

    double rate = 0.0;                                                                                          // default return value                                                      

    switch (p_RSG_MassLossPrescription) {                                                                       // decide which prescription to use
        case RSG_MASS_LOSS_PRESCRIPTION::ZERO            : rate = 0.0; break;
        case RSG_MASS_LOSS_PRESCRIPTION::VINKSABHAHIT2023: rate = CalculateMassLossRateRSGVinkSabhahit2023(); break;            
        case RSG_MASS_LOSS_PRESCRIPTION::BEASOR2020      : rate = CalculateMassLossRateRSGBeasor2020(); break;
        case RSG_MASS_LOSS_PRESCRIPTION::DECIN2023       : rate = CalculateMassLossRateRSGDecin2023(); break;
        case RSG_MASS_LOSS_PRESCRIPTION::YANG2023        : rate = CalculateMassLossRateRSGYang2023(); break;            
        case RSG_MASS_LOSS_PRESCRIPTION::KEE2021         : rate = CalculateMassLossRateRSGKee2021(); break;
        case RSG_MASS_LOSS_PRESCRIPTION::NJ90            : rate = CalculateMassLossRateNieuwenhuijzenDeJager(); break;

        default:                                                                                                // unknown prescription
            // the only ways this can happen are if someone added an RSG_MASS_LOSS_PRESCRIPTION
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without
            // a warning.
            // We are here because the user chose a prescription this code doesn't account for, or as a result
            // of a code defect, and either of those should be flagged as an error and result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option, or find and fix the code defect.

            THROW_ERROR(ERROR::UNKNOWN_RSG_MASS_LOSS_PRESCRIPTION);                                             // throw error
    }

    return rate;
}


/*
 * Calculate mass loss for very massive MS stars, >100Msol. 
 * Switches prescription based on program options. 
 *
 * 
 * double CalculateMassLossRateVMS(const VMS_MASS_LOSS_PRESCRIPTION p_VMS_MassLossPrescription)
 *
 * @param   [IN]    p_VMS_MassLossPrescription  VMS Mass loss prescription to use
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateVMS(const VMS_MASS_LOSS_PRESCRIPTION p_VMS_MassLossPrescription) {

    double rate = 0.0;                                                      

    switch (p_VMS_MassLossPrescription) {                                                                       // decide which prescription to use
        case VMS_MASS_LOSS_PRESCRIPTION::ZERO            : rate = 0.0; break;
        case VMS_MASS_LOSS_PRESCRIPTION::BESTENLEHNER2020: rate = CalculateMassLossRateVMSBestenlehner2020(); break;
        case VMS_MASS_LOSS_PRESCRIPTION::VINK2011        : rate = CalculateMassLossRateVMSVink2011(); break;
        case VMS_MASS_LOSS_PRESCRIPTION::SABHAHIT2023    : rate = CalculateMassLossRateVMSSabhahit2023(); break;

        default:                                                                                                // unknown prescription
            // the only ways this can happen are if someone added a VMS_MASS_LOSS_PRESCRIPTION
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without
            // a warning.
            // We are here because the user chose a prescription this code doesn't account for, or as a result
            // of a code defect, and either of those should be flagged as an error and result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option, or find and fix the code defect.

            THROW_ERROR(ERROR::UNKNOWN_VMS_MASS_LOSS_PRESCRIPTION);                                             // throw error
    }

    return rate;
}


/*
 * Calculate the mass-loss rate for Wolf-Rayet stars according to the
 * prescription of Sander & Vink 2020 (https://arxiv.org/abs/2009.01849)
 * 
 * Use the luminosity prescription given by Equation 13 (see section 3.4.1)
 * 
 * 
 * double CalculateMassLossRateWolfRayetSanderVink2020(const double p_Mu)
 *
 * @param   [IN]    p_Mu                        Small envelope parameter (see Hurley et al. 2000, eq 97 & 98)
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayetSanderVink2020(const double p_Mu) const {

    double Mdot = 0.0;                                                                                      // default return value

    if (utils::Compare(p_Mu, 1.0) < 0) {

        double logL = log10(m_Luminosity);
        double logZ = LogMetallicityXiAnders();

        // Calculate alpha, L0 and Mdot10
        double alpha     = 0.32 * logZ + 1.4;                                                               // Equation 18 in Sander & Vink 2020
        double logL0     = -0.87 * logZ + 5.06;                                                             // Equation 19 in Sander & Vink 2020
        double logMdot10 = -0.75 * logZ - 4.06;                                                             // Equation 20 in Sander & Vink 2020

        if (utils::Compare(logL0, logL) <= 0) {                                                             // No mass loss for L < L0
            // Equation 13 in Sander & Vink 2020
            double logMdot = alpha * log10(logL - logL0) + 0.75 * (logL - logL0 - 1.0) + logMdot10;
            Mdot           = PPOW(10.0, logMdot);
        }
    }

    return Mdot;
}


/*
 * Calculate the correction to the mass-loss rates for Wolf-Rayet stars 
 * as a function of effective temperature, according to the
 * prescription of Sander et al. 2023 (https://arxiv.org/abs/2301.01785)
 * 
 * Use the correction given in Eq. 18, with the effective temperature
 * (what they refer to as T_\star in Eq. 1) as T_eff,crit
 * 
 * 
 * double CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(const double p_Mdot)
 *
 * @param   [IN]    p_Mdot                      Uncorrected mass-loss rate (in Msol yr^{-1})
 * @return                                      Corrected mass-loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(const double p_Mdot) const {

    if (p_Mdot <= 0.0) return 0.0;                                  // nothing to adjust
    
    const double teffRef = 141.0E3;                                 // reference effective temperature in Kelvin
    const double teffMin = 100.0E3;                                 // minimum effective temperature in Kelvin to apply correction

    double teff                = m_Temperature * TSOL;              // get effective temperature in Kelvin
    double logMdotUncorrected  = log10(p_Mdot);                     // uncorrected mass-loss rate
    double logMdotCorrected    = 0.0;

    // Only apply to sufficiently hot stars
    if (utils::Compare(teff, teffMin) > 0) {
        logMdotCorrected = logMdotUncorrected - 6.0 * log10(teff / teffRef);
    }
    else {
        logMdotCorrected = logMdotUncorrected;
    }
    
    return PPOW(10.0, logMdotCorrected);
}


/*
 * Calculate the mass-loss rate for helium stars according to the
 * prescription of Vink 2017 (https://ui.adsabs.harvard.edu/abs/2017A%26A...607L...8V/abstract)
 * 
 * See their Eq. 1
 * 
 * 
 * double CalculateMassLossRateHeliumStarVink2017()
 *
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateHeliumStarVink2017() const {

    double logMdot = -13.3 + (1.36 * log10(m_Luminosity)) + (0.61 * LogMetallicityXiAnders());    // Vink 2017 Eq. 1.

    return PPOW(10.0, logMdot);
}


/*
 * Calculate the mass-loss rate for Wolf--Rayet stars according to the
 * prescription of Shenar et al. 2019 (https://ui.adsabs.harvard.edu/abs/2019A%26A...627A.151S/abstract)
 * 
 * See their Eq. 6 and Table 5
 * 
 * We use the fitting coefficients for hydrogen rich WR stars (e.g., WNh)
 * The C4 (X_He) term is = 0 and is omitted
 * 
 * 
 * double CalculateMassLossRateWolfRayetShenar2019()
 *
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayetShenar2019() const {

    double teff = m_Temperature * TSOL;

    // For H-rich WR stars (X_H > 0.4)
    const double C1 = -6.78;
    const double C2 =  0.66;
    const double C3 = -0.12;
    const double C5 =  0.74;

    double logMdot = C1 + (C2 * log10(m_Luminosity)) + (C3 * log10(teff)) + (C5 * m_Log10Metallicity); 

    return PPOW(10.0, logMdot);
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * According to Hurley et al. 2000
 *
 * 
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double BaseStar::CalculateMassLossRateHurley() {
    return CalculateMassLossRateNieuwenhuijzenDeJager();
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star at the current evolutionary phase
 * Follows the implementation in StarTrack
 *
 * 
 * double CalculateMassLossRateBelczynski2010()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double BaseStar::CalculateMassLossRateBelczynski2010() {
    m_DominantMassLossRate = MASS_LOSS_TYPE::NONE;                                                                  // reset dominant mass loss rate

    double LBVRate = CalculateMassLossRateLBV(OPTIONS->LBVMassLossPrescription());                                  // start with LBV winds (can be, and is often, 0.0)
    double otherWindsRate = 0.0;

    if (m_DominantMassLossRate != MASS_LOSS_TYPE::LBV || 
        OPTIONS->LBVMassLossPrescription() == LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD ) {                            // check whether we should add other winds to the LBV winds (always for HURLEY_ADD prescription, only if not in LBV regime for others)

        double teff = m_Temperature * TSOL;                                                                         // change to Kelvin so it can be compared with values as stated in Vink prescription
        if (utils::Compare(teff, VINK_MASS_LOSS_MINIMUM_TEMP) < 0) {                                                // cool stars, add Hurley et al 2000 winds
            otherWindsRate = CalculateMassLossRateHurley() * OPTIONS->CoolWindMassLossMultiplier();                 // apply cool wind mass loss multiplier
        }
        else  {                                                                                                     // hot stars, add Vink et al. 2001 winds (ignoring bistability jump)
            
            otherWindsRate = CalculateMassLossRateOBVink2001();
            m_DominantMassLossRate = MASS_LOSS_TYPE::OB;
            
            // If user wants to transition between OB and WR mass loss rates
            if (OPTIONS->ScaleMassLossWithSurfaceHeliumAbundance()) {
                otherWindsRate     = EnhanceWindsWithWolfRayetContribution(otherWindsRate, BaseStar::CalculateMassLossRateWolfRayetZDependent(0.0), false);
            }
        }

        if (utils::Compare(LBVRate, otherWindsRate) > 0) {                                                          // which is dominant?
            m_DominantMassLossRate = MASS_LOSS_TYPE::LBV;                                                           // set LBV dominant again in case Hurley or OB overwrote it
        }
    }

    // BSE and StarTrack have some multiplier they apply here
    return LBVRate + otherWindsRate;
}


/*
 * Calculate the mass loss rate according to the updated framework.
 *
 * The structure is similar to the CalculateMassLossRateBelczynski2010() wrapper (previous default).
 * Mass loss rates are divided into several classes: RSG winds, cool star winds,
 * very massive star (VMS) winds, OB star winds.
 * Furthermore, LBV winds are computed separately, and, if non-zero, either replace other mass loss
 * or are added to other wind mass loss if LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD is used.
 * 
 *
 * double CalculateMassLossRateMerritt2025()
 * 
 * @return                  Mass loss rate in Msol per year
 */
double BaseStar::CalculateMassLossRateMerritt2025() {

    m_DominantMassLossRate = MASS_LOSS_TYPE::NONE;

    double LBVRate         = CalculateMassLossRateLBV(OPTIONS->LBVMassLossPrescription());                          // start with LBV winds (can be, and is often, 0.0)
    double otherWindsRate  = 0.0;
    double teff            = TSOL * m_Temperature;

    // calculate other winds rate
    if (m_DominantMassLossRate != MASS_LOSS_TYPE::LBV || 
        OPTIONS->LBVMassLossPrescription() == LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD ) {                            // check whether we should add other winds to the LBV winds (always for HURLEY_ADD prescription, only if not in LBV regime for others)

        if ((utils::Compare(teff, RSG_MAXIMUM_TEMP) < 0) &&                                                         // teff < max temp for RSG winds?
            (utils::Compare(m_MZAMS, MASSIVE_THRESHOLD) >= 0) &&                                                    // ZAMS mass at or above massive threshold?
            (IsOneOf(GIANTS) || m_StellarType == STELLAR_TYPE::HERTZSPRUNG_GAP)) {                                  // must be core helium burning giant(CHeB, FGB, EAGB, TPAGB), or HG
            otherWindsRate         = CalculateMassLossRateRSG(OPTIONS->RSGMassLossPrescription());                  // yes - use RSG mass loss rate
            m_DominantMassLossRate = MASS_LOSS_TYPE::RSG;                                                           // set dominant mass loss rate
        }                                                                      
        else if (utils::Compare(teff, VINK_MASS_LOSS_MINIMUM_TEMP) < 0) {                                           // cool stars, add Hurley et al 2000 winds (NJ90)
            otherWindsRate = CalculateMassLossRateHurley() * OPTIONS->CoolWindMassLossMultiplier();                 // apply cool wind mass loss multiplier
        }
        else if (utils::Compare(m_Mass, VMS_MASS_THRESHOLD) >= 0) {                                                 // mass at or above VMS winds threshold?
            otherWindsRate         = CalculateMassLossRateVMS(OPTIONS->VMSMassLossPrescription());                  // yes - use VMS mass loss rate
            m_DominantMassLossRate = MASS_LOSS_TYPE::VMS;                                                           // set dominant mass loss rate
            // If user wants to transition between OB/VMS and WR mass loss rates
            if (OPTIONS->ScaleMassLossWithSurfaceHeliumAbundance()) {
                otherWindsRate     = EnhanceWindsWithWolfRayetContribution(otherWindsRate, 0.0, true);
            }
        }
        else {                                                                                                      // otherwise...
            otherWindsRate         = CalculateMassLossRateOB(OPTIONS->OBMassLossPrescription());                    // use OB mass loss rate
            m_DominantMassLossRate = MASS_LOSS_TYPE::OB;                                                            // set dominant mass loss rate
            // If user wants to transition between OB and WR mass loss rates
            if (OPTIONS->ScaleMassLossWithSurfaceHeliumAbundance()) {
                otherWindsRate     = EnhanceWindsWithWolfRayetContribution(otherWindsRate, 0.0, true);
            }
        }

        if (utils::Compare(LBVRate, otherWindsRate) > 0) {                                                          // which is dominant?
            m_DominantMassLossRate = MASS_LOSS_TYPE::LBV;                                                           // set LBV dominant again in case Hurley or OB overwrote it
        }
    }

    return LBVRate + otherWindsRate;
}


/*
 * EnhanceWindsWithWolfRayetContribution
 *
 * Enhance the winds with a contribution due to WR finds (see CalculateMassLossFractionWR)
 *
 * double EnhanceWindsWithWolfRayetContribution (double p_OtherWindsRate, double p_WolfRayetRate, bool p_RecalculateWolfRayetRate)
 *
 * @param       p_OtherWindsRate                Wind mass loss rate due to OB or VMS winds to avoid recompution
 * @param       p_WolfRayetRate                 Wind mass loss rate due to WR winds if already computed
 * @param       p_RecalculateWolfRayetRate      If true, recompute the WR winds by cloning the star as a HeMS star
 * @return                                      Total wind mass loss rate
 */
double BaseStar::EnhanceWindsWithWolfRayetContribution (double p_OtherWindsRate, double p_WolfRayetRate, bool p_RecalculateWolfRayetRate) {
    
    double MdotWR = p_WolfRayetRate;
    double fractionWR = CalculateMassLossFractionWR(m_HeliumAbundanceSurface);
    double totalWindsRate = 0.0;
    
    if (p_RecalculateWolfRayetRate && fractionWR > 0.0 ) {
        MdotWR        = BaseStar::CalculateMassLossRateWolfRayetZDependent(0.0);
        // *Jeff* This used to clone the star as a HeMS star and query its CalculateMassLossRateMerritt2025(); eventually, let's switch to a static function to calculate the luminosity of the WR star
    }
    
    // Combine each of these prescriptions according to the OB wind fraction
    totalWindsRate = ((1.0 - fractionWR) * p_OtherWindsRate) + (fractionWR * MdotWR);
    
    if ( fractionWR * MdotWR > (1.0 - fractionWR) * p_OtherWindsRate) {
        m_DominantMassLossRate =MASS_LOSS_TYPE::WR;
    }
    
    return totalWindsRate;
}

/*
 * CalculateMassLossFractionWR
 *
 * @brief
 * Calculate the fraction of mass loss attributable to WR mass loss, per Yoon et al. 2006
 *
 * The model described in Yoon et al. 2006 (also Szecsi et al. 2015) uses OB mass loss while the
 * He surface abundance is below 0.55, WR mass loss when the surface He abundance is above 0.7,
 * and linearly interpolate when the He surface abundance is between those limits.
 *
 * This function calculates the fraction of mass loss attributable to OB mass loss, based on
 * the He surface abundance and the abundance limits described in Yoon et al. 2006.  The value
 * returned will be 1.0 if 100% of the mass loss is attributable to OB mass lass, 0.0 if 100% of
 * the mass loss is attributable to WR mass loss, and in the range (0.0, 1.0) if the mass loss is
 * a mix of OB and WR.
 *
 *
 * double CalculateMassLossFractionWR(const double p_HeAbundanceSurface) const
 *
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Fraction of mass loss attributable to WR mass loss
 */
double BaseStar::CalculateMassLossFractionWR(const double p_HeAbundanceSurface) const {

    constexpr double limOB = 0.55;                                          // per Yoon et al. 2006
    constexpr double limWR = 0.70;                                          // per Yoon et al. 2006

    return std::min(1.0, std::max (0.0, (p_HeAbundanceSurface - limOB) / (limWR - limOB)));
}

/*
 * Calculate mass loss rate
 *
 * Calls relevant mass loss function based on mass loss prescription given in program options (OPTIONS->massLossPrescription)
 *
 *
 * double CalculateMassLossRate()
 *
 * @return                                      Mass loss rate
 */
double BaseStar::CalculateMassLossRate() {

    double mDot = 0.0;                                                                                          // default return value

    if (OPTIONS->MassLossPrescription() != MASS_LOSS_PRESCRIPTION::ZERO) {                                      // mass loss enabled?
                                                                                                                // yes
        double LBVRate;
        double otherWindsRate;

        switch (OPTIONS->MassLossPrescription()) {                                                              // which prescription?

            case MASS_LOSS_PRESCRIPTION::ZERO:
                mDot = 0.0;
                break;

            case MASS_LOSS_PRESCRIPTION::HURLEY:
                LBVRate        = CalculateMassLossRateLBV(LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD);
                otherWindsRate = CalculateMassLossRateHurley();
                if (utils::Compare(LBVRate, otherWindsRate) > 0) {
                    m_DominantMassLossRate = MASS_LOSS_TYPE::LBV;
                }
                mDot = LBVRate + otherWindsRate;
                break;

            case MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010:
                mDot = CalculateMassLossRateBelczynski2010();
                break;

            case MASS_LOSS_PRESCRIPTION::MERRITT2025:
                mDot = CalculateMassLossRateMerritt2025();
                break;

            default:                                                                                                // unknown prescription
                // the only way this can happen is if someone added a MASS_LOSS_PRESCRIPTION
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_MASS_LOSS_PRESCRIPTION);                                                 // throw error
        }

        mDot *= OPTIONS->OverallWindMassLossMultiplier();                                                           // apply overall wind mass loss multiplier
    }
    
    mDot = min(mDot, MAXIMUM_WIND_MASS_LOSS_RATE);                                                                  // cap winds at a maximum mass loss rate (typically 0.1 solar masses per year) to avoid convergence issues
    
    UpdateTotalMassLossRate(-mDot);                                                                                 // update total mass loss rate
    
    return mDot;
}


/*
 * Calculate values for mDot and mass assuming mass loss is applied
 *
 * Class member variable m_Mdot is updated directly by this function if required (see parameters)
 * Class member variable m_Mass is not updated directly by this function - the calculated mass is returned as the functional return
 *
 * - calculates mass loss
 * - calculates new mass loss rate (mDot) to match (possibly limited) mass loss
 * - calculates new mass (mass) based on (possibly limited) mass loss
 * - returns existing value for mass if mass loss not being used (program option)
 *
 *
 * double CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot)
 *
 * @param   [IN]    p_Dt                        time step (Myr)
 * @param   [IN]    p_UpdateMDot                flag to indicate whether the class member variable m_Mdot should be updated (default is false)
 * @return                                      calculated mass (mSol)
 */
double BaseStar::CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot) {

    double mass = m_Mass;

    if (OPTIONS->MassLossPrescription() != MASS_LOSS_PRESCRIPTION::ZERO) {      // mass loss enabled?
                                                                                // yes
        double mDot     = CalculateMassLossRate();                              // calculate mass loss rate
        double massLoss = max(0.0, mDot * p_Dt * 1.0E6);                        // calculate mass loss; mass loss rate given in Msol per year, times are in Myr so need to multiply by 10^6
        if (p_UpdateMDot) m_Mdot = mDot;                                        // update class member variable if necessary

        if (OPTIONS->CheckPhotonTiringLimit()) {
            double lim = m_Luminosity / (G_SOLAR_YEAR * m_Mass / m_Radius);     // calculate the photon tiring limit in Msol yr^-1 using Owocki & Gayley 1997, equation slightly clearer in Owocki+2004 Eq. 20
            massLoss   = std::min(massLoss, lim);                               // limit mass loss to the photon tiring limit
            if (p_UpdateMDot) m_Mdot = massLoss / p_Dt / 1.0E6;                 // update class member variable if necessary
        }

        mass -= massLoss;                                                       // new mass based on mass loss
    }

    return mass;
}


/*
 * Resolve mass loss
 *
 * - calculates mass loss rate
 * - calculates (and limits) mass loss
 * - resets mass loss rate (m_Mdot) to match (possibly limited) mass loss
 * - calculates and sets new mass (m_Mass) based on (possibly limited) mass loss
 * - applies mass rejuvenation factor and calculates new age
 * - updates angular momentum of mass-losing star
 *
 *
 * void ResolveMassLoss(double p_Dt))
 *
 * @param   [IN]    p_Dt                    time step
 */
void BaseStar::ResolveMassLoss(double p_Dt) {

    if (OPTIONS->MassLossPrescription() != MASS_LOSS_PRESCRIPTION::ZERO) {                          // mass loss enabled?
                                                                                                    // yes
        double mass = CalculateMassLossValues(p_Dt, true);                                          // calculate new values assuming mass loss applied

        double angularMomentumChange = (2.0 / 3.0) * (mass - m_Mass) * m_Radius * RSOL_TO_AU * m_Radius * RSOL_TO_AU * Omega();
                
        // JR: this is here to keep attributes in sync BSE vs SSE
        // Supernovae are caught in UpdateAttributesAndAgeOneTimestep()
        // Don't resolve envelope loss here (JR: we're not going to switch anyway... need to revisit this)
        STELLAR_TYPE st = UpdateAttributesAndAgeOneTimestep(mass - m_Mass, 0.0, 0.0, false, false); // recalculate stellar attributes
        if (st != m_StellarType) {                                                                  // should switch?
            SHOW_WARN(ERROR::SWITCH_NOT_TAKEN);                                                     // show warning if we think we should switch again...
        }

        UpdateInitialMass();                                                                        // update effective initial mass (MS, HG & HeMS)
        UpdateAgeAfterMassLoss();                                                                   // update age (MS, HG & HeMS)
        ApplyMassTransferRejuvenationFactor();                                                      // apply age rejuvenation factor
        SetAngularMomentum(m_AngularMomentum + angularMomentumChange);
    }
}


/*
 * Calculate core mass for a given luminosity using the Mc - L relation
 *
 * Hurley et al. 2000, eqs 37 & 38
 *
 *
 * double BaseStar::CalculateCoreMassGivenLuminosity_Static(const double p_Luminosity, const DBL_VECTOR &p_GBParams)
 *
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @param   [IN]    p_GBParams                  Giant Branch parameters
 * @return                                      Core mass in Msol
 */
double BaseStar::CalculateCoreMassGivenLuminosity_Static(const double p_Luminosity, const DBL_VECTOR &p_GBParams) {
#define gbParams(x) p_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    return (utils::Compare(p_Luminosity, gbParams(Lx)) > 0)
            ? PPOW((p_Luminosity / gbParams(B)), (1.0 / gbParams(q)))
            : PPOW((p_Luminosity / gbParams(D)), (1.0 / gbParams(p)));

#undef gbParams
}


/*
 * Calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the accretion efficiency parameter
 *
 *
 * The maximum acceptance rate of the accretor star during mass transfer is based on stellar type: this function
 * is for main sequence stars, or stars that have evolved off of the main sequence but are not yet remnants.
 *
 * Mass transfer is assumed Eddington limited for BHs and NSs.  The formalism of Nomoto/Claeys is used for WDs.
 *
 * For non compact objects:
 *
 *    1) Kelvin-Helmholtz (thermal) timescale if THERMAL (thermally limited) mass transfer efficiency
 *    2) Choose a fraction of the mass rate that will be effectively accreted for FIXED fraction mass transfer (as in StarTrack)
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_AccretorMassRate          Thermal mass transfer rate of the accretor (this star)
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
DBL_DBL BaseStar::CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate) {

    double acceptanceRate   = 0.0;                                                                      // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                                      // accretion fraction - default  = 0.0

    switch (OPTIONS->MassTransferAccretionEfficiencyPrescription()) {

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED:                                   // thermally limited mass transfer
            acceptanceRate   = min(OPTIONS->MassTransferCParameter() * p_AccretorMassRate, p_DonorMassRate);
            fractionAccreted = acceptanceRate / p_DonorMassRate;
            break;

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::HAMSTARS:                                            // thermally limited mass transfer, following Lau+, 2024
            // the mass transfer C parameter is fit using the data from Figure (1) of Lau+, 2024
            acceptanceRate   = min(PPOW(10.0, (4.0 / PPOW((Mass() + 0.2), 0.3) - 0.6)) * p_AccretorMassRate, p_DonorMassRate);
            fractionAccreted = acceptanceRate / p_DonorMassRate;
            break;

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::FIXED_FRACTION:                                      // fixed fraction of mass accreted, as in StarTrack
            fractionAccreted = OPTIONS->MassTransferFractionAccreted();
            acceptanceRate   = min(p_DonorMassRate, fractionAccreted * p_DonorMassRate);
            break;

        default:                                                                                        // unknown prescription
            // the only way this can happen is if someone added an MT_ACCRETION_EFFICIENCY_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION);                           // throw error
    }

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/*
 * Calculate thermal mass acceptance rate
 *
 *
 * double CalculateThermalMassAcceptanceRate(const double p_Radius)
 *
 * @param   [IN]    p_Radius                    Radius of the accretor (Rsol) [typically called with Roche Lobe radius]
 * @return                                      Thermal mass acceptance rate
 */
double BaseStar::CalculateThermalMassAcceptanceRate(const double p_Radius) {
    double acceptanceRate = 0.0;                                                                    // thermal mass acceptance rate

    switch (OPTIONS->MassTransferThermallyLimitedVariation()) {
        case MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE:
            acceptanceRate = (m_Mass - m_CoreMass) / CalculateThermalTimescale(p_Radius);           // uses provided accretor radius (should be Roche lobe radius in practice)
            break;
        case MT_THERMALLY_LIMITED_VARIATION::C_FACTOR:
            acceptanceRate = CalculateThermalMassLossRate();
            break;

        default:                                                                                    // unknown prescription
            // the only way this can happen is if someone added an MT_THERMALLY_LIMITED_VARIATION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_THERMALLY_LIMITED_VARIATION);                             // throw error
    }

    return acceptanceRate;
}


/*
 * Calculate the mass accreted by a Neutron Star given mass and radius of companion
 *
 * JR: todo: flesh-out this documentation
 * JR: is this just for NS?  Maybe a change of name... Or move it to NS.cpp? *Ilya*
 *
 *
 * double CalculateMassAccretedForCO(const double p_Mass,
 *                                   const double p_CompanionMass,
 *                                   const double p_CompanionRadius,
 *                                   const double p_CompanionEnvelope)
 *
 * @param   [IN]    p_Mass                      The mass of the accreting star (Msol)
 * @param   [IN]    p_CompanionMass             The mass of the companion star (Msol)
 * @param   [IN]    p_CompanionRadius           The radius of the companion star (Rsol)
 * @param   [IN]    p_CompanionEnvelope         Envelope of the companion pre-CE
 * @return                                      Mass accreted by the Neutron Star (Msol)
 */
double BaseStar::CalculateMassAccretedForCO(const double p_Mass,
                                            const double p_CompanionMass,
                                            const double p_CompanionRadius,
                                            const double p_CompanionEnvelope) const {

     double deltaMass;

     switch (OPTIONS->CommonEnvelopeMassAccretionPrescription()) {                                              // which prescription?

        case CE_ACCRETION_PRESCRIPTION::ZERO:                                                                   // ZERO
            deltaMass = 0.0;
            break;

        case CE_ACCRETION_PRESCRIPTION::CONSTANT:                                                               // CONSTANT
            deltaMass = OPTIONS->CommonEnvelopeMassAccretionConstant();                                         // use program option
            break;

        case CE_ACCRETION_PRESCRIPTION::UNIFORM:                                                                // UNIFORM
            deltaMass = RAND->Random(OPTIONS->CommonEnvelopeMassAccretionMin(), OPTIONS->CommonEnvelopeMassAccretionMax()); // uniform random distribution - Oslowski+ (2011)
            break;

        case CE_ACCRETION_PRESCRIPTION::MACLEOD: {                                                              // MACLEOD
                                                                                                                // linear regression estimated from Macleod+ (2015)
            double mm = -1.0714285714285712E-05;                                                                // gradient of the linear fit for gradient
            double cm =  0.00012057142857142856;                                                                // intercept of the linear fit for gradient
            double mc =  0.01588571428571428;                                                                   // gradient of the linear fit for intercept
            double cc = -0.15462857142857137;                                                                   // intercept of the linear fir for intercept
            double m  = mm * p_CompanionMass + cm;                                                              // gradient of linear fit for mass
            double c  = mc * p_CompanionMass + cc;                                                              // intercept of linear fit for mass

            // calculate mass to accrete and clamp to minimum and maximum from program options
            deltaMass = std::min(OPTIONS->CommonEnvelopeMassAccretionMax(), std::max(OPTIONS->CommonEnvelopeMassAccretionMin(), m * p_CompanionRadius + c));
            } break;

        case CE_ACCRETION_PRESCRIPTION::CHEVALIER:                                                              // CHEVALIER
                                                                                                                // Model 2 from van Son et al. 2020
            deltaMass = (p_Mass * p_CompanionMass) / (2.0 * (p_Mass + p_CompanionMass)) ;                       // Hoyle-Lyttleton accretion rate times inspiral time
            break;

        default:                                                                                                // unknown prescription
            // the only way this can happen is if someone added a CE_ACCRETION_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_CE_ACCRETION_PRESCRIPTION);                                              // throw error
    }

    deltaMass = std::min(p_CompanionEnvelope, deltaMass);                                                       // clamp the mass accretion to be no more than the envelope of the companion pre CE

    return deltaMass;
}


/*
 * Construct string representing the Mass Transfer Donor History vector
 *
 * This is so that a string is passed to the output, not a vector of stellar types.
 *
 * std::string BaseStar::GetMassTransferDonorHistoryString() 
 *
 * @return                              string of dash-separated stellar type numbers
 */
std::string BaseStar::GetMassTransferDonorHistoryString() const {

    ST_VECTOR   mtHistVec = m_MassTransferDonorHistory;      
    std::string mtHistStr = "";

    if (mtHistVec.empty()) {                                                            // this star was never a donor for MT
        mtHistStr = "NA";
    }
    else {                                                                              // this star was a donor, return the stellar type string
        for (size_t ii = 0; ii < mtHistVec.size(); ii++) {
            mtHistStr += std::to_string(static_cast<int>(mtHistVec[ii])) + "-";         // create string of stellar type followed by dash
        }
        mtHistStr.pop_back();                                                           // remove final dash
    }
    return mtHistStr;
}


/*
 * Add new MT event to event history - only for donor stars
 *
 * void BaseStar::UpdateMassTransferDonorHistory()
 *
 */
void BaseStar::UpdateMassTransferDonorHistory() {

    if (m_MassTransferDonorHistory.empty()) {                                           // no history?
        m_MassTransferDonorHistory.push_back(m_StellarType);                            // yes - first event
    }
    else if (!utils::IsOneOf(m_StellarType, { m_MassTransferDonorHistory.back() })) {   // first MT as current stellar type?
        m_MassTransferDonorHistory.push_back(m_StellarType);                            // yes - new event
    }
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              TEMPERATURE CALCULATIONS                             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the effective temperature of the star on the phase
 *
 *
 * double CalculateTemperatureOnPhase_Static(const double p_Luminosity, const double p_Radius)
 *
 * @param   [IN]    p_Luminosity                Luminosity of the star (Lsol)
 * @param   [IN]    p_Radius                    Radius of the star (Rsol)
 * @return                                      Effective temperature of the star (Tsol)
 */
double BaseStar::CalculateTemperatureOnPhase_Static(const double p_Luminosity, const double p_Radius) {
    return std::sqrt(std::sqrt(p_Luminosity)) / std::sqrt(p_Radius);   // sqrt() is much faster than pow()
}


/*
 * Calculate the effective temperature of the star in Kelvin, given the luminosity of the
 * star (in Lsol) and the radius of the star (in Rsol)
 *
 *
 * double CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius)
 *
 * @param   [IN]    p_Luminosity                Luminosity of the star (Lsol)
 * @param   [IN]    p_Radius                    Radius of the star (Rsol)
 * @return                                      Effective temperature of the star (Kelvin)
 */
double BaseStar::CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) const {
    return CalculateTemperatureOnPhase(p_Luminosity, p_Radius) * TSOL;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                      ROTATIONAL / FREQUENCY CALCULATIONS ETC.                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the analytic cumulative distribution function (CDF) for the
 * equatorial rotational velocity of single O stars.
 *
 * From Equation 1-4 of Ramirez-Agudelo et al 2013 https://arxiv.org/abs/1309.2929
 * Modelled as a mixture of a gamma component and a normal component.
 *
 * Uses the Boost library
 *
 *
 * double CalculateOStarRotationalVelocityAnalyticCDF_Static(const double p_Ve)
 * 
 * @param   [IN]    p_vE                        Rotational velocity (in km s^-1) at which to calculate CDF
 * @return                                      The CDF for the rotational velocity
 */
double BaseStar::CalculateOStarRotationalVelocityAnalyticCDF_Static(const double p_Ve) {

    constexpr double alpha  = 4.82;
    constexpr double beta   = 1.0 / 25.0;
    constexpr double mu     = 205.0;
    constexpr double sigma  = 190.0;
    constexpr double iGamma = 0.43;

    boost::math::inverse_gamma_distribution<> gammaComponent(alpha, beta); // (shape, scale) = (alpha, beta)
    boost::math::normal_distribution<> normalComponent(mu, sigma);
    
    // Compute CDF at zero rotational velocity -- the CDF should be relative to this quantity
    double CDFzero = (iGamma * boost::math::cdf(gammaComponent, 0.0)) + ((1.0 - iGamma) * boost::math::cdf(normalComponent, 0.0));
    
    double CDFunnormalised = (iGamma * boost::math::cdf(gammaComponent, p_Ve)) + ((1.0 - iGamma) * boost::math::cdf(normalComponent, p_Ve));
    
    return ((CDFunnormalised-CDFzero) / (1.0 - CDFzero));

}


/*
 * Calculate rotational velocity from the analytic cumulative distribution function (CDF)
 * for the equatorial rotational velocity of single O stars.
 *
 * Uses inverse sampling and root finding
 *
 * Ramirez-Agudelo et al. 2013 https://arxiv.org/abs/1309.2929
 *
 *
 * double CalculateOStarRotationalVelocity
 *
 * @return                                      Rotational velocity in km s^-1
 */
double BaseStar::CalculateOStarRotationalVelocity() {

    double desiredCDF            = RAND->Random();                                                      // Random desired CDF

    const boost::uintmax_t maxit = ADAPTIVE_RV_MAX_ITERATIONS;                                          // Limit to maximum iterations.
    boost::uintmax_t it          = maxit;                                                               // Initially our chosen max iterations, but updated with actual.

    // find root
    // we use an iterative algorithm to find the root here:
    //    - if the root finder throws an exception, we stop and return a negative value for the root (indicating no root found)
    //    - if the root finder reaches the maximum number of (internal) iterations, we stop and return a negative value for the root (indicating no root found)
    //    - if the root finder returns a solution, we check that func(solution) = 0.0 +/ ROOT_ABS_TOLERANCE
    //       - if the solution is acceptable, we stop and return the solution
    //       - if the solution is not acceptable, we reduce the search step size and try again
    //       - if we reach the maximum number of search step reduction iterations, or the search step factor reduces to 1.0 (so search step size = 0.0),
    //         we stop and return a negative value for the root (indicating no root found)
   
    double guess      = 100.0;                                                                          // guess at 100 km s^-1 (arbitrary initial guess)

    double factorFrac = ADAPTIVE_RV_SEARCH_FACTOR_FRAC;                                                 // search step size factor fractional part
    double factor     = 1.0 + factorFrac;                                                               // factor to determine search step size (size = guess * factor)
    
    std::pair<double, double> root(-1.0, -1.0);                                                         // initialise root - default return
    std::size_t tries = 0;                                                                              // number of tries
    bool done         = false;                                                                          // finished (found root or exceed maximum tries)?
    ERROR error       = ERROR::NONE;
    OStarRotationVelocityFunctor<double> func = OStarRotationVelocityFunctor<double>(desiredCDF);
    while (!done) {                                                                                     // while no error and acceptable root found

        bool isRising = true;                                                                           //guess for direction of search; CDF increases monotonically

        // run the root finder
        // regardless of any exceptions or errors, display any problems as a warning, then
        // check if the root returned is within tolerance - so even if the root finder
        // bumped up against the maximum iterations, or couldn't bracket the root, use
        // whatever value it ended with and check if it's good enough for us - not finding
        // an acceptable root should be the exception rather than the rule, so this strategy
        // shouldn't cause undue performance issues.
        try {
            error = ERROR::NONE;
            root  = boost::math::tools::bracket_and_solve_root(func, guess, factor, isRising, utils::BracketTolerance, it); // find root
            // root finder returned without raising an exception
            if (error != ERROR::NONE) { SHOW_WARN(error); }                                             // root finder encountered an error
            else if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_RV_ITERATIONS); }                         // too many root finder iterations
        }
        catch(std::exception& e) {                                                                      // catch generic boost root finding error
            // root finder exception
            // could be too many iterations, or unable to bracket root - it may not
            // be a hard error - so no matter what the reason is that we are here,
            // we'll just emit a warning and keep trying
            if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_RV_ITERATIONS); }                              // too many root finder iterations
            else             { SHOW_WARN(ERROR::ROOT_FINDER_FAILED, e.what()); }                        // some other problem - show it as a warning
        }

        // we have a solution from the root finder - it may not be an acceptable solution
        // so we check if it is within our preferred tolerance
        if (fabs(func(root.first + (root.second - root.first) / 2.0)) <= ROOT_ABS_TOLERANCE) {          // solution within tolerance?
            done = true;                                                                                // yes - we're done
        }
        else if (fabs(func(root.first)) <= ROOT_ABS_TOLERANCE) {                                        // solution within tolerance at endpoint 1?
            root.second=root.first;
            done = true;                                                                                // yes - we're done
        }
        else if (fabs(func(root.second)) <= ROOT_ABS_TOLERANCE) {                                       // solution within tolerance at endpoint 2?
            root.first=root.second;
            done = true;                                                                                // yes - we're done
        }
        else {                                                                                          // no - try again
            // we don't have an acceptable solution - reduce search step size and try again
            factorFrac /= 2.0;                                                                          // reduce fractional part of factor
            factor      = 1.0 + factorFrac;                                                             // new search step size
            tries++;                                                                                    // increment number of tries
            if (tries > ADAPTIVE_RV_MAX_TRIES || fabs(factor - 1.0) <= ROOT_ABS_TOLERANCE) {            // too many tries, or step size 0.0?
                // we've tried as much as we can - fail here with -ve return value
                root.first  = -1.0;                                                                     // yes - set error return
                root.second = -1.0;
                SHOW_WARN(ERROR::TOO_MANY_RV_TRIES);                                                    // show warning
                done = true;                                                                            // we're done
            }
        }
    }
    
    return root.first + (root.second - root.first) / 2.0;                                               // Midway between brackets is our result, if necessary we could return the result as an interval here.
}


/*
 * Calculate the initial rotational velocity (in km s^-1 ) of a star with ZAMS mass MZAMS
 *
 * Distribution used is determined by program option "rotationalVelocityDistribution"
 *
 *
 * double CalculateRotationalVelocity(double p_MZAMS)
 *
 * @param   [IN]    p_MZAMS                     Zero age main sequence mass in Msol
 * @return                                      Initial equatorial rotational velocity in km s^-1 - vRot in Hurley et al. 2000
 */
double BaseStar::CalculateRotationalVelocity(double p_MZAMS) {

    double vRot = 0.0;

    switch (OPTIONS->RotationalVelocityDistribution()) {                                            // which prescription?

        case ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO:                                                // ZERO
            vRot = 0.0;
            break;

        case ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY:                                              // HURLEY

            // Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
            vRot = (330.0 * PPOW(p_MZAMS, 3.3)) / (15.0 + PPOW(p_MZAMS, 3.45));
            break;

        case ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES:                                           // VLTFLAMES

            // Rotational velocity based on VLT-FLAMES survey.
            // For O-stars (taken to be above 16 Msol), use results
            // of Ramirez-Agudelo et al. (2013) https://arxiv.org/abs/1309.2929 (single stars)
            // and Ramirez-Agudelo et al. (2015) https://arxiv.org/abs/1507.02286 (spectroscopic binaries)
            // For B-stars (taken to be between 2 and 16 Msol) use results
            // of Dufton et al. (2013) https://arxiv.org/abs/1212.2424
            // For lower mass stars, default back to  Hurley et al. 2000 distribution for now

            if (utils::Compare(p_MZAMS, 16.0) >= 0) {
                vRot = CalculateOStarRotationalVelocity();
                vRot = max(vRot, 0.0);                                                              // Set to no rotation if no positive solution found; warning already raised
            }
            else if (utils::Compare(p_MZAMS, 2.0) >= 0) {
                vRot = utils::InverseSampleFromTabulatedCDF(RAND->Random(), BStarRotationalVelocityCDFTable);
            }
            else {
                // Don't know what better to use for low mass stars so for now
                // default to Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
                vRot = (330.0 * PPOW(p_MZAMS, 3.3)) / (15.0 + PPOW(p_MZAMS, 3.45));
            }
            break;

        default:                                                                                    // unknown prescription
            // the only way this can happen is if someone added a ROTATIONAL_VELOCITY_DISTRIBUTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_VROT_PRESCRIPTION);                                          // throw error
    }

    return vRot;
}


/*
 * Calculate the initial angular frequency (in yr^-1) of a star with
 * ZAMS mass and radius MZAMS and RZAMS respectively
 *
 * Hurley et al. 2000, eq 108
 *
 *
 * double CalculateZAMSAngularFrequency(const double p_MZAMS, const double p_RZAMS)
 *
 * @param   [IN]    p_MZAMS                     Zero age main sequence mass in Msol
 * @param   [IN]    p_RZAMS                     Zero age main sequence radius in Rsol
 * @return                                      Initial angular frequency in yr^-1 - omega in Hurley et al. 2000
 */
double BaseStar::CalculateZAMSAngularFrequency(const double p_MZAMS, const double p_RZAMS) {
    double vRot = CalculateRotationalVelocity(p_MZAMS);
    return utils::Compare(vRot, 0.0) == 0 ? 0.0 : 45.35 * vRot / p_RZAMS;               // Hurley et al. 2000, eq 108
}


/*
 * Calculate the break up angular velocity of a star in rad/yr units, where [G] = 4*pi^2 AU^3 yr^-2 Msol^-1
 *
 *
 * double CalculateOmegaBreak() const
 *
 * @return                                      Break up angular velocity (rad yr^-1)
 */
double BaseStar::CalculateOmegaBreak() const {
    constexpr double RSOL_TO_AU_3 = RSOL_TO_AU * RSOL_TO_AU * RSOL_TO_AU;
	return _2_PI * std::sqrt(m_Mass / (RSOL_TO_AU_3 * m_Radius * m_Radius * m_Radius));
}


/*
 * Calculate the minimum rotational frequency (in yr^-1) at which CHE will occur
 * for a star with ZAMS mass MZAMS
 *
 * Mandel's fit from Butler 2018
 *
 *
 * double CalculateOmegaCHE(const double p_MZAMS, const double p_Metallicity)
 *
 * @param   [IN]        p_MZAMS                 Zero age main sequence mass in Msol
 * @param   [IN]        p_Metallicity           Metallicity of the star
 * @return                                      Minimum angular frequency in rad*yr^-1
 */
double BaseStar::CalculateOmegaCHE(const double p_MZAMS, const double p_Metallicity) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double mRatio = p_MZAMS;                                                                        // in MSol, so ratio is just p_MZAMS

    // calculate omegaCHE(M, Z = 0.004)
    double omegaZ004 = 0.0;
    if (utils::Compare(p_MZAMS, massCutoffs(MCHE)) <= 0) {
        for (std::size_t i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(mRatio, i) / PPOW(mRatio, 0.4);
        }
    }
    else {
        for (std::size_t i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(100.0, i) / PPOW(mRatio, 0.4);
        }
    }

    // calculate omegaCHE(M, Z)
    return (1.0 / ((0.09 * log(p_Metallicity / 0.004)) + 1.0) * omegaZ004) * SECONDS_IN_YEAR;       // in rads/yr

#undef massCutoffs
}


/*
 * Calculate the Dynamical tides contribution to the l=2, (n,m) = [(1,0), (1,2), (2,2), (3,2)] imaginary components of the 
 * potential tidal Love number
 *
 * Gravity Waves, Core boundary:
 * Zahn, 1977, Eq. (5.5) , with the value of E_2 coming from Kushnir et al., 2017, by comparing Eq. (8) to Eq. (1)
 * 
 * Gravity Waves, Envelope Boundary:
 * Ahuir et al, 2021, Eq. (131)
 * 
 * Inertial Waves, Convective Envelope:
 * Ogilvie, 2013, Eq. (B3)
 * 
 * DBL_DBL_DBL_DBL CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number, Dynamical tides only (unitless)
 */
DBL_DBL_DBL_DBL BaseStar::CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {
    
    double coreMass = CalculateConvectiveCoreMass();

    double envMass, envMassMax;
    std::tie(envMass, envMassMax) = CalculateConvectiveEnvelopeMass();
    
    double radIntershellMass = m_Mass - coreMass - envMass;                                             // refers to the combined mass of non-convective layers

    double radiusAU              = m_Radius * RSOL_TO_AU;
    double coreRadiusAU          = CalculateConvectiveCoreRadius() * RSOL_TO_AU;
    double convectiveEnvRadiusAU = CalculateRadialExtentConvectiveEnvelope() * RSOL_TO_AU;
    double radiusIntershellAU    = radiusAU - convectiveEnvRadiusAU;                                    // Outer radial coordinate of radiative intershell

    // There should be no Dynamical tides if the entire star is convective, i.e. if there are no convective-radiative boundaries. 
    // If so, return 0.0 for all dynamical components of ImKnm.
    // This condition should be true for low-mass MS stars (<= 0.35 Msol) at ZAMS.
    if (utils::Compare(radIntershellMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || utils::Compare(radiusIntershellAU, coreRadiusAU) <= 0) {
        return std::make_tuple(0.0, 0.0, 0.0, 0.0);                           
    }
 
    double R_3            = radiusAU * radiusAU * radiusAU;
    double R3OverG_M      = (R_3 / G_AU_Msol_yr / m_Mass);
    double sqrtR3OverG_M  = std::sqrt(R3OverG_M);

    double k10GravityCore = 0.0;                                                                        // gravity Wave dissipation, core boundary
    double k12GravityCore = 0.0;
    double k22GravityCore = 0.0;
    double k32GravityCore = 0.0;

    double k10GravityEnv  = 0.0;                                                                        // gravity Wave dissipation, envelope boundary
    double k12GravityEnv  = 0.0;
    double k22GravityEnv  = 0.0;
    double k32GravityEnv  = 0.0;

    double k22InertialEnv = 0.0;                                                                        // inertial Wave dissipation, envelope
    
    double omegaSpin      = Omega();
    double twoOmegaSpin   = omegaSpin + omegaSpin;

    double w10 = p_Omega;
    double w12 = ((p_Omega) - twoOmegaSpin);
    double w22 = ((p_Omega + p_Omega) - twoOmegaSpin);
    double w32 = ((p_Omega + p_Omega + p_Omega) - twoOmegaSpin);
        
    // Assume that GW dissipation from core boundary is only efficient if the radiative region extends to the surface, i.e. there is no convective envelope.
    if (utils::Compare(coreRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0 && utils::Compare(coreMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0 && utils::Compare(convectiveEnvRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0 && utils::Compare(envMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0) {                   
        constexpr double beta2Dynamical         = 1.0;
        constexpr double rhoFactorDynamcial     = 0.1;
        double coreRadiusOverRadius   = coreRadiusAU / radiusAU;
        double coreRadiusOverRadius_3 = coreRadiusOverRadius * coreRadiusOverRadius * coreRadiusOverRadius;
        double coreRadiusOverRadius_9 = coreRadiusOverRadius_3 * coreRadiusOverRadius_3 * coreRadiusOverRadius_3;
        double massOverCoreMass       = m_Mass / coreMass;
        double E2Dynamical            = (2.0 / 3.0) * coreRadiusOverRadius_9 * massOverCoreMass * std::cbrt(massOverCoreMass) * beta2Dynamical * rhoFactorDynamcial;

        // (l=2, n=1, m=0), Gravity Wave dissipation from core boundary
        double s10     = w10 * sqrtR3OverG_M;
        double s10_4_3 = s10 * std::cbrt(s10);
        double s10_8_3 = s10_4_3 * s10_4_3;
        k10GravityCore = E2Dynamical *  std::copysign(s10_8_3, w10);
        if (std::isnan(k10GravityCore)) k10GravityCore = 0.0;

        // (l=2, n=1, m=2), Gravity Wave dissipation from core boundary
        double s12     = w12 * sqrtR3OverG_M;
        double s12_4_3 = s12 * std::cbrt(s12);
        double s12_8_3 = s12_4_3 * s12_4_3;
        k12GravityCore = E2Dynamical * std::copysign(s12_8_3, w12);
        if (std::isnan(k12GravityCore)) k12GravityCore = 0.0;

        // (l=2, n=2, m=2), Gravity Wave dissipation from core boundary
        double s22     = w22 * sqrtR3OverG_M;
        double s22_4_3 = s22 * std::cbrt(s22);
        double s22_8_3 = s22_4_3 * s22_4_3;
        k22GravityCore = E2Dynamical * std::copysign(s22_8_3, w22);
        if (std::isnan(k22GravityCore)) k22GravityCore = 0.0;

        // (l=2, n=3, m=2), Gravity Wave dissipation from core boundary
        double s32     = w32 * sqrtR3OverG_M;
        double s32_4_3 = s32 * std::cbrt(s32);
        double s32_8_3 = s32_4_3 * s32_4_3;
        k32GravityCore = E2Dynamical * std::copysign(s32_8_3, w32);
        if (std::isnan(k32GravityCore)) k32GravityCore = 0.0;    
    }

    double rint_3 = radiusIntershellAU * radiusIntershellAU * radiusIntershellAU;
    double rc_3   = coreRadiusAU * coreRadiusAU * coreRadiusAU;
    double gamma  = (envMass / (R_3 - rint_3)) / (radIntershellMass / (rint_3 - rc_3));

    // There is no GW or IW dissipation from the envelope boundary if no convective envelope
    if ((utils::Compare(convectiveEnvRadiusAU / radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0) || (utils::Compare(envMass / m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0)) {    

        constexpr double dynPrefactor     = 3.207452512782476;                                                        // 3^(11/3) * Gamma(1/3)^2 / 40 PI
        constexpr double m_l_factor_22    = 0.183440402716368;                                                        // m * (l(l+1))^{-4/3}
        double cbrtdNdlnr       = std::cbrt(G_AU_Msol_yr * radIntershellMass / radiusIntershellAU / (radiusAU - radiusIntershellAU) / (radiusAU - radiusIntershellAU));
        
        double alpha            = radiusIntershellAU / radiusAU;
        double oneMinusAlpha    = 1.0 - alpha;
        double beta             = radIntershellMass / m_Mass;

        double alpha_2          = alpha * alpha;
        double alpha_3          = alpha_2 * alpha;
        double alpha_5          = alpha_3 * alpha_2;
        double alpha_11         = alpha_5 * alpha_5 * alpha;
        double oneMinusAlpha_2  = oneMinusAlpha * oneMinusAlpha;
        double oneMinusAlpha_3  = 1.0 - alpha_3;
        double beta_2           = beta * beta;

    
        double oneMinusGamma    = 1.0 - gamma;
        double oneMinusGamma_2  = oneMinusGamma * oneMinusGamma;
        double alpha_2_3Minus_1 = (alpha * 2.0 / 3.0) - 1.0;

        // Assume GW dissipation from the envelope boundary only acts if the radiative zone extends to the core, i.e. if there is no convective core.
        if (utils::Compare(coreRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0) {                              
            double Epsilon       = alpha_11 * envMass / m_Mass * oneMinusGamma_2 * alpha_2_3Minus_1 * alpha_2_3Minus_1 / beta_2 / oneMinusAlpha_3 / oneMinusAlpha_2;

            // (l=2, n=1, m=0), Gravity Wave dissipation from envelope boundary is always 0.0 since m * (l(l+1))^{-4/3} = 0

            // (l=2, n=1, m=2), Gravity Wave dissipation from envelope boundary
            double w12_4_3       = w12 * std::cbrt(w12);
            double w12_8_3       = w12_4_3 * w12_4_3;
            k12GravityEnv        = dynPrefactor * m_l_factor_22 * std::copysign(w12_8_3, w12) * R3OverG_M * Epsilon / cbrtdNdlnr;
            if (std::isnan(k12GravityEnv)) k12GravityEnv = 0.0;  

            // (l=2, n=2, m=2), Gravity Wave dissipation from envelope boundary
            double w22_4_3       = w22 * std::cbrt(w22);
            double w22_8_3       = w22_4_3 * w22_4_3;
            k22GravityEnv        = dynPrefactor * m_l_factor_22 * std::copysign(w22_8_3, w22)* R3OverG_M * Epsilon / cbrtdNdlnr;
            if (std::isnan(k22GravityEnv)) k22GravityEnv = 0.0;  

            // (l=2, n=3, m=2), Gravity Wave dissipation from envelope boundary
            double w32_4_3       = w32 * std::cbrt(w32);
            double w32_8_3       = w32_4_3 * w32_4_3;
            k32GravityEnv        = dynPrefactor * m_l_factor_22 * std::copysign(w32_8_3, w32) * R3OverG_M * Epsilon / cbrtdNdlnr;
            if (std::isnan(k32GravityEnv)) k32GravityEnv = 0.0;  
        }

        // (l=2, n=2, m=2), Inertial Wave dissipation, convective envelope
        // IW dissipation is only efficient for highly spinning stars, as in Esseldeurs, et al., 2024 
        if (utils::Compare(twoOmegaSpin, p_Omega) >= 0) {                                                                            
            double epsilonIW_2       = omegaSpin * omegaSpin * R3OverG_M;
            double oneMinusAlpha_4 = oneMinusAlpha_2 * oneMinusAlpha_2;
            double bracket1          = 1.0 + (2.0 * alpha) + (3.0 * alpha_2) + (3.0 * alpha_3 / 2.0);
            double bracket2          = 1.0 + (oneMinusGamma / gamma) * alpha_3;
            double bracket3          = 1.0 + (3.0 * gamma / 2.0) + (5.0 * alpha_3 / (2.0 * gamma) * (1.0 + (gamma / 2.0) - (3.0* gamma * gamma / 2.0))) - (9.0 / 4.0 * oneMinusGamma * alpha_5);
            k22InertialEnv           = (100.0 * M_PI / 63.0) * epsilonIW_2 * (alpha_5 / (1.0 - alpha_5)) * oneMinusGamma_2 * oneMinusAlpha_4 * bracket1 * bracket1 * bracket2 / bracket3 / bracket3;
            k22InertialEnv           = std::copysign(k22InertialEnv, w22);
            if (std::isnan(k22InertialEnv)) k22InertialEnv = 0.0;  
        }
    }

    // return ImKnmDynamical
    return std::make_tuple(k10GravityCore + k10GravityEnv, k12GravityCore + k12GravityEnv, k22GravityCore + k22GravityEnv + k22InertialEnv, k32GravityCore + k32GravityEnv);
}


/*
 * Calculate the Equilibrium tides contribution to the l=2, (n,m) = [(1,0), (1,2), (2,2), (3,2)] imaginary components of the 
 * potential tidal Love number
 * 
 * Barker (2020), Eqs. (20) to (27), (l=2, n=2, m=2 mode only).
 *
 * DBL_DBL_DBL_DBL CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number, Equilibrium tides only (unitless)
 */
DBL_DBL_DBL_DBL BaseStar::CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {

    // Viscous dissipation
    // No contribution from convective core; only convective envelope.

    double rEnvAU = CalculateRadialExtentConvectiveEnvelope() * RSOL_TO_AU;
    double envMass, envMassMax;
    std::tie(envMass, envMassMax) = CalculateConvectiveEnvelopeMass();
    
    double rOutAU = m_Radius * RSOL_TO_AU;                                                      // outer boundary of convective envelope
    double rInAU  = (rOutAU - rEnvAU);                                                          // inner boundary of convective envelope
    
    if (utils::Compare(rEnvAU/rOutAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || utils::Compare(envMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || std::isnan(envMass)) return std::make_tuple(0.0, 0.0, 0.0, 0.0);           // skip calculations if there is no convective envelope (to avoid Imk22 = NaN)

    double rOut_2  = rOutAU * rOutAU;
    double rOut_3  = rOut_2 * rOutAU;
    double rOut_5  = rOut_2 * rOut_3;
    double rOut_7  = rOut_2 * rOut_5;
    double rOut_9  = rOut_2 * rOut_7;

    double rIn_2  = rInAU * rInAU;
    double rIn_3  = rIn_2 * rInAU;
    double rIn_5  = rIn_2 * rIn_3;
    double rIn_7  = rIn_2 * rIn_5;
    double rIn_9  = rIn_2 * rIn_7;

    double omegaSpin      = Omega();
    double twoOmegaSpin   = omegaSpin + omegaSpin;

    double rhoConv        = envMass / (4.0 * M_PI * (rOut_3 - rIn_3) / 3.0);
    double lConv          = rEnvAU;                                                              // set length scale to height of convective envelope
    double tConv          = CalculateEddyTurnoverTimescale();
    double vConv          = lConv / tConv;
    double omegaConv      = 1.0 / tConv;                                                         // absent factor of 2*PI, following Barker (2020)
    double vl             = vConv * lConv;
    double M_2            = m_Mass * m_Mass;

    double vl_5           = 5.0 * vl;
    double vl25OverRoot20 = vl * (25.0 / std::sqrt(20.0));
    double vlOver2        = 0.5 * vl;

    double w10 = p_Omega;
    double w12 = ((p_Omega) - (twoOmegaSpin));
    double w22 = ((p_Omega + p_Omega) - (twoOmegaSpin));
    double w32 = ((p_Omega + p_Omega + p_Omega) - (twoOmegaSpin));

    double k2_prefactor   = (224.0 * M_PI / 15.0) * (rOut_9 - rIn_9) * rhoConv / G_AU_Msol_yr / M_2 / rOut_5;

    // (l=2, n=1, m=0), Viscous dissipation, convective envelope
    double omega_t_10            = std::abs(w10);                                               
    double omega_tOverOmega_c_10 = omega_t_10 / omegaConv;
    double nuTidal10             = vl_5;
    if (utils::Compare(omega_tOverOmega_c_10, 5.0) > 0) {             
        nuTidal10 = vl25OverRoot20 / omega_tOverOmega_c_10 / omega_tOverOmega_c_10;
    }
    else if (utils::Compare(omega_tOverOmega_c_10, 0.01) > 0) {
        nuTidal10 = vlOver2 / std::sqrt(omega_tOverOmega_c_10);    
    }
    double k10Equilibrium    = k2_prefactor * nuTidal10 * omega_t_10;
    if (std::isnan(k10Equilibrium)) k10Equilibrium = 0.0;
    if (w10 < 0.0) k10Equilibrium = -std::abs(k10Equilibrium);


    // (l=2, n=1, m=2), Viscous dissipation, convective envelope
    double omega_t_12            = std::abs(w12);                                               
    double omega_tOverOmega_c_12 = omega_t_12 / omegaConv;
    double nuTidal12             = vl_5;
    if (utils::Compare(omega_tOverOmega_c_12, 5.0) > 0) {             
        nuTidal12 = vl25OverRoot20 / omega_tOverOmega_c_12 / omega_tOverOmega_c_12;
    }
    else if (utils::Compare(omega_tOverOmega_c_12, 0.01) > 0) {
        nuTidal12 = vlOver2 / std::sqrt(omega_tOverOmega_c_12);    
    }
    double k12Equilibrium    = k2_prefactor * nuTidal12 * omega_t_12;
    if (std::isnan(k12Equilibrium)) k12Equilibrium = 0.0;
    if (w12 < 0) k12Equilibrium = -std::abs(k12Equilibrium);


    // (l=2, n=2, m=2), Viscous dissipation, convective envelope
    double omega_t_22            = std::abs(w22);                                               
    double omega_tOverOmega_c_22 = omega_t_22 / omegaConv;
    double nuTidal22             = vl_5;
    if (utils::Compare(omega_tOverOmega_c_22, 5.0) > 0) {             
        nuTidal22 = vl25OverRoot20 / omega_tOverOmega_c_22 / omega_tOverOmega_c_22;
    }
    else if (utils::Compare(omega_tOverOmega_c_22, 0.01) > 0) {
        nuTidal22 = vlOver2 / std::sqrt(omega_tOverOmega_c_22);    
    }
    double k22Equilibrium    = k2_prefactor * nuTidal22 * omega_t_22;
    if (std::isnan(k22Equilibrium)) k22Equilibrium = 0.0;
    if (w22 < 0.0) k22Equilibrium = -std::abs(k22Equilibrium);


    // (l=2, n=3, m=2), Viscous dissipation, convective envelope
    double omega_t_32              = std::abs(w32);                                               
    double omega_t_over_omega_c_32 = omega_t_32 / omegaConv;
    double nuTidal32               = vl_5;
    if (utils::Compare(omega_t_over_omega_c_32, 5.0) > 0) {             
        nuTidal32 = vl25OverRoot20 / omega_t_over_omega_c_32 / omega_t_over_omega_c_32;
    }
    else if (utils::Compare(omega_t_over_omega_c_32, 0.01) > 0) {
        nuTidal32 = vlOver2 / std::sqrt(omega_t_over_omega_c_32);    
    }
    double k32Equilibrium    = k2_prefactor * nuTidal32 * omega_t_32;
    if (std::isnan(k32Equilibrium)) k32Equilibrium = 0.0;
    if (w32 < 0.0) k32Equilibrium = -std::abs(k32Equilibrium);

    // return ImKnmEquilibrium
    return std::make_tuple(k10Equilibrium, k12Equilibrium, k22Equilibrium, k32Equilibrium);
}


/*
 * Calculate the l=2, (n,m) = [(1,0), (1,2), (2,2), (3,2)] imaginary components of the potential tidal Love number 
 * by combining Equilibrium and Dynamical tidal contributions.
 *
 * DBL_DBL_DBL_DBL CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number (unitless)
 */
DBL_DBL_DBL_DBL BaseStar::CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {
    
    double Imk10Dynamical, Imk12Dynamical, Imk22Dynamical, Imk32Dynamical;
    std::tie(Imk10Dynamical, Imk12Dynamical, Imk22Dynamical, Imk32Dynamical) = CalculateImKnmDynamical(p_Omega, p_SemiMajorAxis, p_M2);

    double Imk10Equilibrium, Imk12Equilibrium, Imk22Equilibrium, Imk32Equilibrium;
    std::tie(Imk10Equilibrium, Imk12Equilibrium, Imk22Equilibrium, Imk32Equilibrium) = CalculateImKnmEquilibrium(p_Omega, p_SemiMajorAxis, p_M2);
    
    // return combined ImKnm terms;
    return std::make_tuple(Imk10Dynamical + Imk10Equilibrium, Imk12Dynamical + Imk12Equilibrium, Imk22Dynamical + Imk22Equilibrium, Imk32Dynamical + Imk32Equilibrium);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate lifetime to the Base of the Giant Branch (end of the Hertzsprung Gap)
 * For high mass stars, t_BGB = t_HeI.
 *
 * Hurley et al. 2000, eq 4 (plotted in Hurley et al. 2000, fig 5)
 *
 *
 * double CalculateLifetimeToBGB(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Lifetime to the Base of the Giant Branch in Myr
 */
double BaseStar::CalculateLifetimeToBGB(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function
    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    double m_2   = p_Mass * p_Mass;
    double m_4   = m_2 * m_2;
    double m_5_5 = m_4 * p_Mass * std::sqrt(p_Mass);
    double m_7   = m_4 * m_2 * p_Mass;

    return (a[1] + (a[2] * m_4) + (a[3] * m_5_5) + m_7) / ((a[4] * m_2) + (a[5] * m_7));

#undef a
}


/*
 * Calculate lifetime to the Base of the Asymptotic Giant Branch
 * tBAGB = tHeI + tHe
 *
 *
 * double CalculateLifetimeToBAGB(const double p_tHeI, const double p_tHe)
 *
 * @param   [IN]    p_tHeI                      Time to helium ignition
 * @param   [IN]    p_tHe                       Time to helium burning
 * @return                                      Lifetime to Base of the Asymptotic Giant Branch in Myr
 */
double BaseStar::CalculateLifetimeToBAGB(const double p_tHeI, const double p_tHe) const {
    return p_tHeI + p_tHe;
}


/*
 * Calculate dynamical timescale
 *
 * Kalogera & Webbink 1996, eq 1
 *
 *
 * double CalculateDynamicalTimescale_Static(const double p_Mass, const double p_Radius)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @return                                      Dynamical timescale in Myr
 */
double BaseStar::CalculateDynamicalTimescale_Static(const double p_Mass, const double p_Radius) {
    return 5.0 * 1.0E-5 * p_Radius * std::sqrt(p_Radius) * YEAR_TO_MYR / std::sqrt(p_Mass);   // sqrt() is much faster than pow()
}


/*
 * Calculate thermal timescale
 *
 * pre-factor from Kalogera & Webbink 1996 (https://arxiv.org/abs/astro-ph/9508072), equation 2, 
 * combined with p_Mass * p_EnvMass case from equation 61 from https://arxiv.org/abs/astro-ph/0201220 for k in {2,3,4,5,6,8,9}
 * [note that equation 61 of BSE (https://arxiv.org/abs/astro-ph/0201220) approximates this with a value a factor of 3 smaller]
 * 
 * 
 * double CalculateThermalTimescale(const double p_Radius) const
 *
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @return                                      Thermal timescale in Myr
 *
 * The p_Radius parameter is to accommodate the call (of this function) in BaseBinaryStar::CalculateMassTransfer()
*/
double BaseStar::CalculateThermalTimescale(const double p_Radius) const {   
    return 31.4 * m_Mass * (m_Mass == m_CoreMass ? m_Mass : m_Mass - m_CoreMass) / (p_Radius * m_Luminosity); // G*Msol^2/(Lsol*Rsol) ~ 31.4 Myr (~ 30 Myr in Kalogera & Webbink)
}


/*
 * Calculate radial expansion timescale
 *
 *
 * double CalculateRadialExpansionTimescale_Static(const STELLAR_TYPE p_StellarType,
 *                                                 const STELLAR_TYPE p_StellarTypePrev,
 *                                                 const double       p_Radius,
 *                                                 const double       p_RadiusPrev,
 *                                                 const double       p_DtPrev)
 *
 * @param   [IN]    p_StellarType               Current stellar type of star
 * @param   [IN]    p_StellarTypePrev           Previous stellar type of star
 * @param   [IN]    p_Radius                    Current radius of star in Rsol
 * @param   [IN]    p_RadiusPrev                Previous radius of star in Rsol
 * @param   [IN]    p_DtPrev                    Previous timestep in Myr
 * @return                                      Radial expansion timescale in Myr
 *                                              Returns -1.0 if radial expansion timescale can't be calculated
 *                                              (i.e. stellar type has changed or radius has not changed)
 */
double BaseStar::CalculateRadialExpansionTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                          const STELLAR_TYPE p_StellarTypePrev,
                                                          const double       p_Radius,
                                                          const double       p_RadiusPrev,
                                                          const double       p_DtPrev) {

    return (p_StellarTypePrev == p_StellarType && utils::Compare(p_RadiusPrev, p_Radius) != 0)
            ? (p_DtPrev * p_RadiusPrev) / fabs(p_Radius - p_RadiusPrev)
            : -1.0;
}


/*
 * Calculate the radial expansion timescale in the mass transfer regime
 * We do not use CalculateRadialExpansionTimescale(), since in the process of mass transfer the previous radius
 * is determined by binary evolution, not nuclear timescale evolution
 *
 *
 * double CalculateRadialExpansionTimescaleDuringMassTransfer()
 *
 * @return                                      Radial expansion timescale
 */
double BaseStar::CalculateRadialExpansionTimescaleDuringMassTransfer() {
    
    // We create and age it slightly to determine how the radius will change.
    // To be sure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    BaseStar *clone = Clone(OBJECT_PERSISTENCE::EPHEMERAL, false);                              // do not re-initialise the clone

    double timestep = std::max(1000.0 * NUCLEAR_MINIMUM_TIMESTEP, m_Age / 1.0E6);
    clone->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, timestep, true, false);
    double radiusAfterAging = clone->Radius();
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    return timestep * m_Radius / fabs(m_Radius - radiusAfterAging);
}


/*
 * Calculate mass change timescale
 *
 *
 * double CalculateMassChangeTimescale_Static(const STELLAR_TYPE p_StellarType,
 *                                            const STELLAR_TYPE p_StellarTypePrev,
 *                                            const double       p_Mass,
 *                                            const double       p_MassPrev,
 *                                            const double       p_DtPrev)
 *
 * @param   [IN]    p_StellarType               Current stellar type of star
 * @param   [IN]    p_StellarTypePrev           Previous stellar type of star
 * @param   [IN]    p_Mass                      Current mass of star in Msol
 * @param   [IN]    p_MassPrev                  Previous radius of star in Msol
 * @param   [IN]    p_DtPrev                    Previous timestep in Myr
 * @return                                      Mass change timescale in Myr
 *                                              Returns -1.0 if mass change timescale can't be calculated
 *                                              (i.e. stellar type has changed or mass has not changed)
 */
double BaseStar::CalculateMassChangeTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                     const STELLAR_TYPE p_StellarTypePrev,
                                                     const double       p_Mass,
                                                     const double       p_MassPrev,
                                                     const double       p_DtPrev) {

    return p_StellarTypePrev == p_StellarType && utils::Compare(p_MassPrev, p_Mass) != 0
            ? (p_DtPrev * p_MassPrev) / fabs(p_Mass - p_MassPrev)
            : -1.0;
}



/*
 * Calculate the eddy turnover timescale
 * Hurley+2002, sec. 2.3, particularly eq. 31 of subsec. 2.3.1
 *
 *
 * double CalculateEddyTurnoverTimescale()
 *
 * @return                                      eddy turnover timescale (yr)
 */
double BaseStar::CalculateEddyTurnoverTimescale() const {

	double rEnv	= CalculateRadialExtentConvectiveEnvelope();
    double mEnv, mEnvmax;
    std::tie(mEnv, mEnvmax) = CalculateConvectiveEnvelopeMass();
    return 0.4311 * cbrt((mEnv * rEnv * (m_Radius - (0.5 * rEnv))) / (3.0 * m_Luminosity));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                SUPERNOVA FUNCTIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Draw a kick magnitude in km s^-1 from a Maxwellian distribution of the form:
 *
 *
 * double DrawKickMagnitudeDistributionMaxwell(const double p_Sigma, const double p_Rand)
 *
 * @param   [IN]    p_Sigma                     Distribution scale parameter - affects the spread of the distribution
 * @param   [IN]    p_Rand                      Random number between 0 and 1 used for drawing from the inverse CDF of the Maxwellian
 * @return                                      Drawn kick magnitude (km s^-1)
 */
double BaseStar::DrawKickMagnitudeDistributionMaxwell(const double p_Sigma, const double p_Rand) const {
    return p_Sigma * std::sqrt(gsl_cdf_chisq_Pinv(p_Rand, 3)); // a Maxwellian is a chi distribution with three degrees of freedom
}


/*
 * Draw a kick magnitude in km s^-1 from a uniform distribution between 0 and parameter p_MaxVK
 *
 *
 * double DrawKickMagnitudeDistributionFlat(const double p_MaxVK, const double p_Rand)
 *
 * @param   [IN]    p_MaxVK                     Maximum kick magnitude in km s^-1 to draw
 * @param   [IN]    p_Rand                      Random number between 0 and 1 used for drawing from the distribution
 * @return                                      Drawn kick magnitude (km s^-1)
 */
double BaseStar::DrawKickMagnitudeDistributionFlat(const double p_MaxVK, const double p_Rand) const {
    return p_Rand * p_MaxVK;
}


/*
 * Draw a kick magnitude in km s^-1 per Bray & Eldridge 2016, 2018
 *
 * See:
 *    https://arxiv.org/abs/1605.09529
 *    https://arxiv.org/abs/1804.04414
 *
 *
 * double DrawKickMagnitudeBrayEldridge(const double p_EjectaMass,
 *                                     const double p_RemnantMass,
 *                                     const double p_Alpha,
 *                                     const double p_Beta)
 *
 * @param   [IN]    p_EjectaMass                Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @param   [IN]    p_RemnantMass               Mass of the remnant (Msol)
 * @param   [IN]    p_Alpha                     Fitting coefficient (see Bray & Eldridge 2016, 2018)
 * @param   [IN]    p_Beta                      Fitting coefficient (see Bray & Eldridge 2016, 2018)
 * @return                                      Drawn kick magnitude (km s^-1)
 */
double BaseStar::DrawKickMagnitudeBrayEldridge(const double p_EjectaMass,
                                               const double p_RemnantMass,
                                               const double p_Alpha,
                                               const double p_Beta) const {

    return p_Alpha * (p_EjectaMass / p_RemnantMass) + p_Beta;
}


/*
 * Draw kick magnitude per Muller et al. 2016 as presented in eq. B5 of Vigna-Gomez et al. 2018 (arXiv:1805.07974)
 *
 * BHs do not get natal kicks
 *
 * double DrawRemnantKickMuller(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                Carbon Oxygen core mass of exploding star (Msol)
 * @return                                      Drawn kick magnitude (km s^-1)
 */
double BaseStar::DrawRemnantKickMuller(const double p_COCoreMass) const {

    double	remnantKick     = 0.0;	            // units km/s
	double	lowerRegimeKick = 35.0;		        // following Vigna-Gomez et al. 2018 using 35 km/s as the peak of a low-kick Maxwellian (e.g. USSN, ECSN)

	     if (utils::Compare(p_COCoreMass, 1.372) <  0) remnantKick = 0.0;
	else if (utils::Compare(p_COCoreMass, 1.49 ) <  0) remnantKick = lowerRegimeKick + (1000.0 * (p_COCoreMass - 1.372));
    else if (utils::Compare(p_COCoreMass, 1.65 ) <  0) remnantKick = 90.0 + (650.0 * (p_COCoreMass - 1.49));
	else if (utils::Compare(p_COCoreMass, 2.4  ) <  0) remnantKick = 100.0 + (175.0  * (p_COCoreMass - 1.65));
    else if (utils::Compare(p_COCoreMass, 3.2  ) <  0) remnantKick = 200.0 + (550.0 * (p_COCoreMass - 2.4));
    else if (utils::Compare(p_COCoreMass, 3.6  ) <  0) remnantKick = 80.0 + (120.0  * (p_COCoreMass - 3.2));
    else if (utils::Compare(p_COCoreMass, 4.05 ) <  0) remnantKick = 0.0;                                                // going to be a Black Hole
    else if (utils::Compare(p_COCoreMass, 4.6  ) <  0) remnantKick = 350.0 + (50.0  * (p_COCoreMass - 4.05));
    else if (utils::Compare(p_COCoreMass, 5.7  ) <  0) remnantKick = 0.0;                                                // going to be a Black Hole
    else if (utils::Compare(p_COCoreMass, 6.0  ) <  0) remnantKick = 275.0 - (300.0  * (p_COCoreMass - 5.7));
    else if (utils::Compare(p_COCoreMass, 6.0  ) >= 0) remnantKick = 0.0;                                                // going to be a Black Hole

    return remnantKick;
}


/*
 * Draw kick magnitude per Mandel and Mueller, 2020
 *
 * double DrawRemnantKickMullerMandel(const double p_COCoreMass, 
 *                                    const double p_Rand,
 *                                    const double p_RemnantMass)
 * 
 * @param   [IN]    p_COCoreMass                Carbon Oxygen core mass of exploding star (Msol)
 * @param   [IN]    p_Rand                      Random number between 0 and 1 used for drawing from the distribution
 * @param   [IN]    p_RemnantMass               Mass of the remnant (Msol)
 * @return                                      Drawn kick magnitude (km s^-1)
 */
double BaseStar::DrawRemnantKickMullerMandel(const double p_COCoreMass, 
                                             const double p_Rand,
                                             const double p_RemnantMass) const {					
	double remnantKick;
	double muKick      = 0.0;
    double sigmaKick   = 0.0;

	if (utils::Compare(p_RemnantMass, OPTIONS->MaximumNeutronStarMass()) <  0) {
		muKick = max(OPTIONS->MullerMandelKickMultiplierNS() * (p_COCoreMass - p_RemnantMass) / p_RemnantMass, 0.0);
        sigmaKick = OPTIONS->MullerMandelSigmaKickNS();
	}
	else {
		muKick = max(OPTIONS->MullerMandelKickMultiplierBH() * (p_COCoreMass - p_RemnantMass) / p_RemnantMass, 0.0);
        sigmaKick = OPTIONS->MullerMandelSigmaKickBH();
	}

    double quantile0 = gsl_cdf_gaussian_P(-1.0, sigmaKick);  //quantile of -1 in the Gaussian CDF; the goal is to draw from the cut-off Gaussian since the kick must exceed 0
    double rand = quantile0 + p_Rand * (1.0 - quantile0);
    remnantKick = muKick * (1.0 + gsl_cdf_gaussian_Pinv(rand, sigmaKick));
    
    // Mandel * Mueller 2020 call for USSN kicks to be treated in the same way as CCSN kicks; however, if this override flag is set, set the USSN kick to be equal to the user-provided magnitude
    if (utils::SNEventType(m_SupernovaDetails.events.current) == SN_EVENT::USSN && OPTIONS->USSNKicksOverrideMandelMuller() ) {
        remnantKick = OPTIONS->KickMagnitudeDistributionSigmaForUSSN();
    }
    
	return remnantKick;
}


/*
 * Draw a kick magnitude from the user-specified distribution
 *
 *
 * double DrawSNKickMagnitude(const double p_Sigma,
 *                           const double p_COCoreMass,
 *                           const double p_Rand,
 *                           const double p_EjectaMass,
 *                           const double p_RemnantMass)
 *
 * @param   [IN]    p_Sigma                     Distribution scale parameter - affects the spread of the distribution
 * @param   [IN]    p_COCoreMass                Carbon Oxygen core mass of exploding star (Msol)
 * @param   [IN]    p_Rand                      Random number between 0 and 1 used for drawing from the distribution
 * @param   [IN]    p_EjectaMass                Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @param   [IN]    p_RemnantMass               Mass of the remnant (Msol)
 * @return                                      Drawn kick magnitude (km s^-1)
 */
double BaseStar::DrawSNKickMagnitude(const double p_Sigma,
                                     const double p_COCoreMass,
                                     const double p_Rand,
                                     const double p_EjectaMass,
                                     const double p_RemnantMass) {
	double kickMagnitude;

    switch (OPTIONS->KickMagnitudeDistribution()) {                                             // which distribution

        case KICK_MAGNITUDE_DISTRIBUTION::ZERO:                                                 // ZERO
            kickMagnitude = 0.0;
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN:                                           // MAXWELLIAN, MAXWELL
            kickMagnitude = DrawKickMagnitudeDistributionMaxwell(p_Sigma, p_Rand);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::FLAT:                                                 // FLAT
            kickMagnitude = DrawKickMagnitudeDistributionFlat(OPTIONS->KickMagnitudeDistributionMaximum(), p_Rand);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::FIXED:                                                // FIXED
            kickMagnitude = p_Sigma;
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::BRAYELDRIDGE:                                         // BRAY ELDRIDGE
            kickMagnitude = DrawKickMagnitudeBrayEldridge(p_EjectaMass, p_RemnantMass, BRAY_ELDRIDGE_CONSTANT_VALUES.at(BRAY_ELDRIDGE_CONSTANT::ALPHA), BRAY_ELDRIDGE_CONSTANT_VALUES.at(BRAY_ELDRIDGE_CONSTANT::BETA));
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::MULLER2016:                                           // MULLER2016
            kickMagnitude = DrawRemnantKickMuller(p_COCoreMass);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::MULLER2016MAXWELLIAN: {                               // MULLER2016-MAXWELLIAN

            double mullerSigma = DrawRemnantKickMuller(p_COCoreMass) / std::sqrt(3.0);

            kickMagnitude = DrawKickMagnitudeDistributionMaxwell(mullerSigma, p_Rand);
            } break;

        case  KICK_MAGNITUDE_DISTRIBUTION::MULLERMANDEL:                                        // MULLERMANDEL
            kickMagnitude = DrawRemnantKickMullerMandel(p_COCoreMass, p_Rand, p_RemnantMass);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::LOGNORMAL: {                                          // LOGNORMAL
            // Only draw Disberg & Mandel (2025) kicks for CCSN or PPISN (if they receive a kick)
            // use Maxwellians with sigma set in CalculateSNKickMagnitude() for other SN types
            SN_EVENT thisSNevent = utils::SNEventType(m_SupernovaDetails.events.current);       // current SN event
            if (thisSNevent == SN_EVENT::CCSN || (thisSNevent == SN_EVENT::PPISN && OPTIONS->NatalKickForPPISN())) {
                // maximum kick of DISBERG_MANDEL_MAX_KICK = 1000 km/s, following Disberg & Mandel (2025)
                // normalise the draw so that the kick is uniformly drawn from the inverse CDF below this maximum
                double cdfMax = gsl_cdf_lognormal_P(DISBERG_MANDEL_MAX_KICK, DISBERG_MANDEL_MU, DISBERG_MANDEL_SIGMA);
                kickMagnitude = gsl_cdf_lognormal_Pinv(p_Rand*cdfMax, DISBERG_MANDEL_MU, DISBERG_MANDEL_SIGMA);
            }
            else
                kickMagnitude = DrawKickMagnitudeDistributionMaxwell(p_Sigma, p_Rand);
            } break;
            
        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added a KICK_MAGNITUDE_DISTRIBUTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION);                  // throw error
    }

    return kickMagnitude / OPTIONS->KickScalingFactor();
}


/*
 * Calculate supernova kick magnitude
 * Based on the current supernova event type and user-specified kick magnitude distributions
 *
 *
 * double BaseStar::CalculateSNKickMagnitude(const double p_RemnantMass, const double p_EjectaMass, const STELLAR_TYPE p_StellarType)
 *
 * @param   [IN]    p_RemnantMass               The mass of the remnant (Msol)
 * @param   [IN]    p_EjectaMass                Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @param   [IN]    p_StellarType		        Expected remnant type
 * @return                                      Kick magnitude
 */
double BaseStar::CalculateSNKickMagnitude(const double p_RemnantMass, const double p_EjectaMass, const STELLAR_TYPE p_StellarType) {

    ERROR error = ERROR::NONE;

    SN_EVENT thisSNevent = utils::SNEventType(m_SupernovaDetails.events.current);                   // current SN event
	double   vK          = 0.0;                                                                     // kick magnitude

    if (!m_SupernovaDetails.initialKickParameters.magnitudeSpecified ||                             // user did not supply kick magnitude, or
         m_SupernovaDetails.initialKickParameters.magnitudeRandomSpecified) {                       // ... wants to draw magnitude using supplied random number

        double sigma;
        switch (thisSNevent) {                                                                      // what type of supernova event happening now?

		    case SN_EVENT::ECSN:                                                                    // ECSN may have a separate kick prescription
			    sigma = OPTIONS->KickMagnitudeDistributionSigmaForECSN();
                break;

		    case SN_EVENT::USSN:                                                                    // USSN may have a separate kick prescription
			    sigma = OPTIONS->KickMagnitudeDistributionSigmaForUSSN();
                break;

		    case SN_EVENT::AIC:                                                                     // AIC have 0 kick 
		    case SN_EVENT::SNIA:                                                                    // SNIA have 0 kick 
		    case SN_EVENT::HeSD:                                                                    // HeSD have 0 kick 
			    sigma = 0.0;
                break;

            case SN_EVENT::CCSN:
		    case SN_EVENT::PPISN:                                                                   // draw a random kick magnitude from the user selected distribution - sigma based on whether compact object is a NS or BH

                if(thisSNevent == SN_EVENT::PPISN && !OPTIONS->NatalKickForPPISN()) {
                    sigma = 0.0;
                }
                else {
                    if (p_StellarType == STELLAR_TYPE::NEUTRON_STAR)                                // neutron star?
                        sigma = OPTIONS->KickMagnitudeDistributionSigmaCCSN_NS();                   // yes
                    else if (p_StellarType == STELLAR_TYPE::BLACK_HOLE)                             // no - black hole?
                        sigma = OPTIONS->KickMagnitudeDistributionSigmaCCSN_BH();                   // yes
                    else                                                                            // unexpected stellar type - shouldn't happen                                        
                        error = ERROR::UNEXPECTED_STELLAR_TYPE;                                     // set error value
                }
                break;

            case SN_EVENT::PISN:                                                                    // not expected here
                error = ERROR::UNEXPECTED_SN_EVENT;                                                 // set error value
                break;

            case SN_EVENT::NONE:                                                                    // no supernova event - shouldn't be here...
                error = ERROR::EXPECTED_SN_EVENT;                                                   // set error value
                break;

            default:                                                                                // unknown SN event type
                // the only ways this can happen are if someone added an SN_EVENT
                // and it isn't accounted for in this code, or if there is a defect in the code that causes
                // m_SupernovaDetails.events.current to be invalid.  We should not default here, with or without
                // a warning.
                // We are here because an SN event occurred that this code doesn't account for, or as a result
                // of a code defect, and either of those should be flagged as an error and result in termination of
                // the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option, or find and fix the code defect.

                error = ERROR::UNKNOWN_SN_EVENT;                                                    // set error value
	    }
    
	    if (error == ERROR::NONE) {                                                                 // check for errors
                                                                                                    // no errors - draw kick magnitude
            vK = DrawSNKickMagnitude(sigma, m_SupernovaDetails.COCoreMassAtCOFormation, m_SupernovaDetails.kickMagnitudeRandom, p_EjectaMass, p_RemnantMass);
        }
    }
    else {                                                                                          // user supplied kick parameters and wants to use supplied kick magnitude, so ...
        vK = m_SupernovaDetails.initialKickParameters.magnitude;                                    // ... use it 
    }

	if (error == ERROR::NONE) {                                                                     // check for errors
                                                                                                    // no errors
        m_SupernovaDetails.drawnKickMagnitude = vK;                                                 // drawn kick magnitude

        if (thisSNevent == SN_EVENT::CCSN && utils::IsOneOf(p_StellarType, { STELLAR_TYPE::BLACK_HOLE })) { // core-collapse supernova event this timestep, and remnant is black hole?
            vK = BH::ReweightSupernovaKickByMass_Static(vK, m_SupernovaDetails.fallbackFraction, m_Mass);   // yes - re-weight kick by mass of remnant according to user specified black hole kicks option, if relevant (default is no reweighting)
        }
        else {                                                                                      // otherwise
            m_SupernovaDetails.fallbackFraction = 0.0;                                              // set fallback fraction to zero
        }
        m_SupernovaDetails.kickMagnitude = vK;                                                      // updated kick magnitude
    }
    else {                                                                                          // error occurred
        THROW_ERROR(error);                                                                         // throw error
    }

    return vK;
}


/*
 * Calculate eccentric anomaly and true anomaly - uses kepler's equation
 *
 * Modifies class member variables m_SupernovaDetails.eccentricAnomaly and m_SupernovaDetails.trueAnomaly
 *
 *
 * void CalculateSNAnomalies(const double p_Eccentricity)
 *
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary
 */
void BaseStar::CalculateSNAnomalies(const double p_Eccentricity) {

    ERROR error = ERROR::NONE;

    std::tie(error, m_SupernovaDetails.eccentricAnomaly, m_SupernovaDetails.trueAnomaly) = utils::SolveKeplersEquation(m_SupernovaDetails.meanAnomaly, p_Eccentricity);

         if (error == ERROR::NO_CONVERGENCE) { THROW_ERROR(error, "Solving Kepler's equation"); }       // no convergence - throw error   
    else if (error == ERROR::OUT_OF_BOUNDS ) { THROW_ERROR(error, "Eccentric anomaly"); }               // eccentric anomaly out of bounds - throw error

    return;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              ENERGY RELATED FUNCTIONS                             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 *	Calculate the absolute value of the binding energy of the envelope of the star
 *
 *
 * double CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda)
 *
 * @param   [IN]    p_CoreMass                  Core mass of the star (Msol)
 * @param   [IN]    p_EnvMass                   Envelope mass of the star (Msol)
 * @param   [IN]    p_Radius                    Radius of the star (Rsol)
 * @param   [IN]    p_Lambda                    Dimensionless parameter defining the binding energy
 * @return                                      Binding energy (erg)
 */
double BaseStar::CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const {

    double bindingEnergy = 0.0;                                                         // default

	if (p_Radius <= 0.0) {                                                              // positive radius?
        THROW_ERROR(ERROR::RADIUS_NOT_POSITIVE, "Binding energy = 0.0");
	}
	else if (p_Lambda <= 0.0) {                                                         // positive lambda?
        THROW_ERROR(ERROR::LAMBDA_NOT_POSITIVE, "Binding energy = 0.0");
	}
	else {                                                                              // calculate binding energy
        // convert to CGS where necessary
        double radius    = p_Radius * RSOL_TO_CM;
        double coreMass  = p_CoreMass * MSOL_TO_G;
        double envMass   = p_EnvMass * MSOL_TO_G;

        double totalMass = coreMass + envMass;                                          // total mass

		bindingEnergy    = G_CGS * totalMass * envMass / (p_Lambda * radius);           // erg
	}

	return bindingEnergy;
}


/*
 * Calculate convective envelope binding energy for the two-stage Hirai & Mandel (2022) common envelope formalism
 *
 *
 * double CalculateConvectiveEnvelopeBindingEnergy(const double p_TotalMass, const double p_ConvectiveEnvelopeMass, const double p_Radius, const double p_Lambda) const
 *
 * @param   [IN]    p_TotalMass                 Total mass of the star (Msol)
 * @param   [IN]    p_ConvectiveEnvelopeMass    Mass of the convective outer envelope  (Msol)
 * @param   [IN]    p_Radius                    Radius of the star (Rsol)
 * @param   [IN]    p_Lambda                    Lambda parameter for the convective envelope
 * @return                                      Binding energy (erg)
 */
double BaseStar::CalculateConvectiveEnvelopeBindingEnergy(const double p_TotalMass, const double p_ConvectiveEnvelopeMass, const double p_Radius, const double p_Lambda) const {
    return CalculateBindingEnergy(p_TotalMass - p_ConvectiveEnvelopeMass, p_ConvectiveEnvelopeMass, p_Radius, p_Lambda);
}


/*
 * Approximates the lambda binding energy parameter of the outer convective envelope
 *
 * This is needed for the Hirai & Mandel (2022) two-stage CE formalism.
 * Follows the fits of Picker, Hirai, Mandel (2024), arXiv:2402.13180 for lambda_He
 *
 *
 * double BaseStar::CalculateConvectiveEnvelopeLambdaPicker(const DBL_DBL p_convectiveEnvelopeMass) const
 *
 * @param   [IN]    p_convectiveEnvelopeMass    Tuple: Mass of the outer convective envelope shell and Maximum mass of the outer convective envelope shell at the onset of carbon burning
 * @return                                      Lambda binding energy parameter for the outer convective envelope
 */
double BaseStar::CalculateConvectiveEnvelopeLambdaPicker(const DBL_DBL p_convectiveEnvelopeMass) const {
    
    double envMass, envMassMax;
    std::tie(envMass, envMassMax) = p_convectiveEnvelopeMass;
    double m2         = 0.0023 * m_Log10Metallicity * m_Log10Metallicity + 0.0088 * m_Log10Metallicity + 0.013;         // Eq. (12) and Table 1 of Picker, Hirai, Mandel (2024)
    double b1         = m2 * m_Mass - 0.23;                                                                             // Eq. (11) of Picker+ (2024)
    double logLambda  = envMass / envMassMax > 0.3
                            ? 0.42 * envMass / envMassMax + b1
                            : 0.3 * 0.42 + b1;
    
    return exp(logLambda);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Determines if the star is one of a list of stellar types passed
 *
 *
 * bool IsOneOf(const STELLAR_TYPE_LIST p_List)
 *
 * @param   [IN]    p_List                      List of stellar types
 * @return                                      Boolean - true if star is in list, false if not
 */
bool BaseStar::IsOneOf(const STELLAR_TYPE_LIST p_List) const {
    for (auto elem: p_List) {
        if (m_StellarType == elem) return true;
    }
	return false;
}


/*
 * Calculate next timestep for stellar evolution
 *
 * Timestep based on stellar type, age, etc.
 *
 *
 * double CalculateTimestep()
 *
 * @return                                      Timestep
 */
double BaseStar::CalculateTimestep() {
    
    double radialExpansionTimescale = CalculateRadialExpansionTimescale();
    double massChangeTimescale      = CalculateMassChangeTimescale();
    double dt                       = 0.0;

    if (massChangeTimescale > 0.0)                                                                           // non-positive means it could not be computed (e.g., just after stellar type change)
        dt = OPTIONS->MassChangeFraction() * massChangeTimescale;

    if (radialExpansionTimescale > 0.0)                                                                      // non-positive means it could not be computed (e.g., just after stellar type change)
        dt = dt <= 0.0 ? OPTIONS->RadialChangeFraction() * radialExpansionTimescale : min(dt, OPTIONS->RadialChangeFraction() * radialExpansionTimescale);
    
    // the GBParams and Timescale calculations need to be done
    // before the timestep calculation - since the binary code
    // calls this functiom, the GBParams and Timescale functions
    // are called here
    CalculateGBParams();                                                                                    // calculate giant branch parameters
    CalculateTimescales();                                                                                  // calculate timescales

    dt = dt <= 0.0 ? ChooseTimestep(m_Age) : min(dt, ChooseTimestep(m_Age));
    
    // there is a chance that mass loss from winds is much faster than previously estimated if, say, LBV winds have turned on
    // we therefore precompute the mass loss rate to avoid taking an overly long timestep, despite the extra computational costs
    double massChangeWinds = m_Mass - CalculateMassLossValues(dt, false);
    if(utils::Compare(massChangeWinds, 0.0) != 0)
        dt = min(dt, OPTIONS->MassChangeFraction() * (dt * m_Mass / fabs(massChangeWinds)));
            
    dt = max(round(dt / TIMESTEP_QUANTUM) * TIMESTEP_QUANTUM, NUCLEAR_MINIMUM_TIMESTEP);
        
    return dt;
}


/*
 * Set parameters required before ageing one timestep - modify star attributes
 *
 *
 * void AgeOneTimestepPreamble(const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaTime                 Timestep in Myr
 */
void BaseStar::AgeOneTimestepPreamble(const double p_DeltaTime) {

    if (p_DeltaTime > 0.0) {                        // if dt > 0    (don't use utils::Compare() here)
        m_Time += p_DeltaTime;                      // advance physical simulation time
        m_Age  += p_DeltaTime;                      // advance age of star
        m_Dt    = p_DeltaTime;                      // set timestep
    }
    EvolveOneTimestepPreamble();                    // per stellar type
}


/*
 * Set parameters required before updating stellar attributes (via evolution) - modify star attributes
 *
 * Will apply mass changes to m_Mass and/or m_Mass0.  Note discussion in documentation for
 * UpdateAttributes() in Star.cpp - changing m_Mass0 is a bit of a kludge and should be fixed.
 *
 * If the timestep (p_DeltaTime) is > 0 then m_Mass and m_Radius will be saved as m_MassPrev and m_RadiusPrev respectively.
 *
 *
 * void UpdateAttributesAndAgeOneTimestepPreamble(const double p_DeltaMass, const double p_DeltaMass0, const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_DeltaTime                 Timestep in Myr
 */
void BaseStar::UpdateAttributesAndAgeOneTimestepPreamble(const double p_DeltaMass, const double p_DeltaMass0, const double p_DeltaTime) {

    if (utils::Compare(p_DeltaMass,  0.0) != 0) { m_Mass  = max(0.0, m_Mass  + p_DeltaMass);  }     // update mass as required (only change if delta != 0 with tolerance) and prevent -ve
    if (utils::Compare(p_DeltaMass0, 0.0) != 0) { m_Mass0 = max(0.0, m_Mass0 + p_DeltaMass0); }     // update mass0 as required (only change if delta != 0 with tolerance) and prevent -ve

    // record some current values before they are (possibly) changed by evolution
    if (p_DeltaTime > 0.0) {                                                                        // don't use utils::Compare() here
            m_StellarTypePrev = m_StellarType;
            m_MassPrev        = m_Mass;
            m_RadiusPrev      = m_Radius;
    }
    
    // the GBParams and Timescale calculations need to be done before taking the timestep - since
    // the binary code ultimately calls this via UpdateAttributesAndAgeOneTimestep(), the GBParams
    // and Timescale functions are called here.
    //
    // JR: todo: we should revisit where and how often we recalculate GBParams and Timescales.  The
    // problem is that there are multiple entry points into the calculate/take timestep code that it
    // isn't always obvious where we need to do this...  A project for another time.

    CalculateGBParams();                                                                            // calculate giant branch parameters
    if (p_DeltaTime > 0.0) CalculateTimescales();                                                   // calculate timescales if necessary
}


/*
 * Apply mass changes if required, age the star one timestep, advance the simulation time, and update the
 * attributes of the star.
 *
 * The star's attributes (Age, Radius, Luminosity etc.) are calculated and updated as required.
 *
 * Free parameters in the update process are the star's mass (m_Mass), initial mass (m_Mass0), the star's age
 * (m_Age) and the simulation time attribute (m_Time):
 *
 *    - if required, the star's mass is changed by the amount passed as the p_DeltaMass parameter before other
 *      attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to the
 *      star's mass before the attributes of the star are calculated.
 *
 *    - if required, the star's effective initial mass is changed by the amount passed as the p_DeltaMass0 parameter before
 *      other attributes are updated.  The p_DeltaMass0 parameter may be zero, in which case no change is made to
 *      the star's mass before the attributes of the star are calculated.  The Mass0 attribute in Hurley et al. 2000 is overloaded by the introduction
 *      of mass loss (see section 7.1).
 *
 *    - if required, the star is aged by the amount passed as the p_DeltaTime parameter, and the simulation time is
 *      advanced by the same amount, before other attributes are updated.  The p_deltaTime parameter may be zero,
 *      in which case no change is made to the star's age or the physical time attribute.
 *
 *
 * Before updating attributes we check whether the star:
 *    - is a massless remnant - we don't update attributes of massless remnants
 *    - has become a supernova - if so we resolve the supernova and do not update the attributes
 *    - should skip this phase for this timestep (checked after applying p_DeltaMass, p_DeltaMass0 and p_DeltaTime)
 *
 * If none of the above are true the star evolves on phase for the specified timestep (which may be 0, in which case
 * the star's attributes other than age are re-calculated), then the need to evolve the star off phase is checked.
 *
 * If p_DeltaMass, p_DeltaMass0 and p_DeltaTime are all passed as zero the checks for massless remnant and supernova
 * are performed (and consequential changes made), but no other changes to the star's attributes are made - unless
 * the p_ForceRecalculate parameter is set true.
 *
 * The functional return is the stellar type to which the star should evolve.  The returned stellar type is just the
 * stellar type of the star upon entry if it should remain on phase.  The star's stellar type is not changed here.
 *
 *
 * STELLAR_TYPE UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
 *                                                const double p_DeltaMass0,
 *                                                const double p_DeltaTime,
 *                                                const bool   p_ForceRecalculate,
 *                                                const bool   p_ResolveEnvelopeLoss)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_DeltaTime                 The timestep to evolve in Myr
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the three deltas are 0.0
 *                                              (optional, default = false)
 * @param   [IN]    p_ResolveEnvelopeLoss       Specifies whether envelope loss should be resolved here
 *                                              (optional, default = true)
 *                                              JR: this is a bit of a kludge to resolve problems introduced by modifying stellar attributes in
 *                                                  anticipation of switching stellar type, but using those attributes before the actual switch
 *                                                  to the new stellar type - we need to resolve those situations in the code.
 *                                                  The bottom line is that this parameter is a temporary fix while I develop the permanent fix.
 * @return                                      Stellar type to which star should evolve
 */
STELLAR_TYPE BaseStar::UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                         const double p_DeltaMass0,
                                                         const double p_DeltaTime,
                                                         const bool   p_ForceRecalculate,
                                                         const bool   p_ResolveEnvelopeLoss) {
    STELLAR_TYPE stellarType = m_StellarType;                                                   // default is no change

    if (ShouldBeMasslessRemnant()) {                                                            // do not update the star if it lost all of its mass
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;
    }
    else {
        stellarType = ResolveSupernova();                                                       // handle supernova
        if (stellarType == m_StellarType) {                                                     // still on phase?
            
            UpdateAttributesAndAgeOneTimestepPreamble(p_DeltaMass, p_DeltaMass0, p_DeltaTime);  // apply mass changes and save current values if required

            if (p_ForceRecalculate                     ||                                       // force recalculate?
                utils::Compare(p_DeltaMass,  0.0) != 0 ||                                       // mass change? or...
                utils::Compare(p_DeltaMass0, 0.0) != 0 ||                                       // mass0 change? or...
                               p_DeltaTime         > 0) {                                       // age/time advance? (don't use utils::Compare() here)
                                                                                                // yes - update attributes
                AgeOneTimestepPreamble(p_DeltaTime);                                            // advance dt, age, simulation time if necessary (don't use utils::Compare() here)

                if (ShouldSkipPhase()) stellarType = ResolveSkippedPhase();                     // skip phase if required
                else {                                                                          // not skipped - execute phase
                    stellarType = EvolveOnPhase(p_DeltaTime);                                   // evolve on phase
                    if (stellarType == m_StellarType) {                                         // need to switch to new stellar type?
                        stellarType = ResolveEndOfPhase(p_ResolveEnvelopeLoss);                 // no - check for need to move off phase
                    }   
                }
            }
        }
    }

    return stellarType;                                                                         // stellar type to which star should evolve
}


/*
 * Evolve the star on it's current phase - take one timestep on the current phase
 *
 *
 * STELLAR_TYPE EvolveOnPhase(const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaTime                 Timestep in Myr
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::EvolveOnPhase(const double p_DeltaTime) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (ShouldEvolveOnPhase()) {                                                    // evolve timestep on phase
        
        UpdateMainSequenceCoreMass(p_DeltaTime, -m_Mdot);                           // update core mass, relevant for MS stars

        m_Tau        = CalculateTauOnPhase();

        m_COCoreMass = CalculateCOCoreMassOnPhase();
        m_CoreMass   = CalculateCoreMassOnPhase();
        m_HeCoreMass = CalculateHeCoreMassOnPhase();

        m_Luminosity = CalculateLuminosityOnPhase();

        // Calculate abundances
        m_HeliumAbundanceCore      = CalculateHeliumAbundanceCoreOnPhase();
        m_HeliumAbundanceSurface   = CalculateHeliumAbundanceSurfaceOnPhase();
        m_HydrogenAbundanceCore    = CalculateHydrogenAbundanceCoreOnPhase();
        m_HydrogenAbundanceSurface = CalculateHydrogenAbundanceSurfaceOnPhase();  
        
        std::tie(m_Radius, stellarType) = CalculateRadiusAndStellarTypeOnPhase();   // radius and possibly new stellar type

        m_Mu = CalculatePerturbationMuOnPhase();

        PerturbLuminosityAndRadiusOnPhase();

        m_Temperature = CalculateTemperatureOnPhase();

        if (p_DeltaTime > 0.0) {
            STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();                   // resolve envelope loss if it occurs - possibly new stellar type
            if (thisStellarType != m_StellarType) {                                 // thisStellarType overrides stellarType (from CalculateRadiusAndStellarTypeOnPhase())
                stellarType = thisStellarType;
            }
        }
    }

    return stellarType;
}


/*
 * Evolve the star onto the next phase if necessary - take one timestep at the end of the current phase
 *
 *
 * STELLAR_TYPE ResolveEndOfPhase(const bool p_ResolveEnvelopeLoss)
 *
 * @param   [IN]    p_ResolveEnvelopeLoss       Specifies whether envelope loss should be resolved here
 *                                              (optional, default = true)
 *                                              JR: this is a bit of a kludge to resolve problems introduced by modifying stellar attributes in
 *                                                  anticipation of switching stellar type, but using those attributes before the actual switch
 *                                                  to the new stellar type - we need to resolve those situations in the code.
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::ResolveEndOfPhase(const bool p_ResolveEnvelopeLoss) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (IsEndOfPhase()) {                                                       // end of phase

        if (p_ResolveEnvelopeLoss) stellarType = ResolveEnvelopeLoss();         // if required, resolve envelope loss if it occurs

        if (stellarType == m_StellarType) {                                     // staying on phase?
            
            m_Tau         = CalculateTauAtPhaseEnd();

            m_COCoreMass  = CalculateCOCoreMassAtPhaseEnd();
            m_CoreMass    = CalculateCoreMassAtPhaseEnd();
            m_HeCoreMass  = CalculateHeCoreMassAtPhaseEnd();

            m_Luminosity  = CalculateLuminosityAtPhaseEnd();

            m_Radius      = CalculateRadiusAtPhaseEnd();

            m_Mu          = CalculatePerturbationMuAtPhaseEnd();

            PerturbLuminosityAndRadiusAtPhaseEnd();

            m_Temperature = CalculateTemperatureAtPhaseEnd();

            stellarType   = EvolveToNextPhase();                                // determine the stellar type to which the star should evolve
        }
    }

    return stellarType;
}
