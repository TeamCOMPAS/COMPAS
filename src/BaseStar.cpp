// gsl includes
#include <gsl/gsl_roots.h>

// boos includes
#include <boost/math/distributions.hpp>

#include "Rand.h"
#include "BaseStar.h"

using std::max;
using std::min;


BaseStar::BaseStar() {

    // initialise member variables

    m_ObjectId           = globalObjectId++;                           // unique object id - remains for life of star (even through evolution to other phases)
    m_ObjectType         = OBJECT_TYPE::BASE_STAR;                     // object type - remains for life of star (even through evolution to other phases)
    m_InitialStellarType = STELLAR_TYPE::STAR;                         // stellar type - changes throughout life of star (through evolution to other phases)
    m_StellarType        = STELLAR_TYPE::STAR;                         // stellar type - changes throughout life of star (through evolution to other phases)

    m_Error              = ERROR::NOT_INITIALISED;                     // clear error flag
}


BaseStar::BaseStar(const unsigned long int p_RandomSeed, 
                   const double            p_MZAMS, 
                   const double            p_Metallicity, 
                   const KickParameters    p_KickParameters,
                   const double            p_LBVfactor, 
                   const double            p_WolfRayetFactor) {

    // initialise member variables

    m_ObjectId            = globalObjectId++;                           // unique object id - remains for life of star (even through evolution to other phases)
    m_ObjectType          = OBJECT_TYPE::BASE_STAR;                     // object type - remains for life of star (even through evolution to other phases)
    m_InitialStellarType  = STELLAR_TYPE::STAR;                         // stellar type - changes throughout life of star (through evolution to other phases)
    m_StellarType         = STELLAR_TYPE::STAR;                         // stellar type - changes throughout life of star (through evolution to other phases)

    m_Error               = ERROR::NONE;                                // clear error flag

    m_CHE                 = false;                                      // initially
    
    // Initialise member variables from input parameters
    // (kick parameters initialised below - see m_SupernovaDetails)
    m_RandomSeed          = p_RandomSeed;
    m_MZAMS               = p_MZAMS;
    m_Metallicity         = std::min(std::max(p_Metallicity, 0.0), 1.0);
    m_LBVfactor           = p_LBVfactor;
    m_WolfRayetFactor     = p_WolfRayetFactor;


    // Initialise metallicity dependent values
    m_LogMetallicityXi    = log10(m_Metallicity / ZSOL);
    m_LogMetallicitySigma = log10(m_Metallicity);
    m_LogMetallicityRho   = m_LogMetallicityXi + 1.0;


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

    CalculateRCoefficients(m_LogMetallicityXi, m_RCoefficients);
    CalculateLCoefficients(m_LogMetallicityXi, m_LCoefficients);

    CalculateMassCutoffs(m_Metallicity, m_LogMetallicityXi, m_MassCutoffs);

    CalculateAnCoefficients(m_AnCoefficients, m_LConstants, m_RConstants, m_GammaConstants);
    CalculateBnCoefficients(m_BnCoefficients);

    m_XExponent                                = CalculateGBRadiusXExponent();
    m_Alpha1                                   = CalculateAlpha1();
    m_Alpha3                                   = CalculateAlpha3();
    m_Alpha4                                   = CalculateAlpha4();

    // timescales
    m_DynamicalTimescale                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_ThermalTimescale                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_NuclearTimescale                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_RadialExpansionTimescale                 = DEFAULT_INITIAL_DOUBLE_VALUE;

    // initialise remaining member variables

    // Zero age main sequence parameters
    m_RZAMS                                    = CalculateRadiusAtZAMS(m_MZAMS);
    m_LZAMS                                    = CalculateLuminosityAtZAMS(m_MZAMS);
    m_TZAMS                                    = CalculateTemperatureOnPhase_Static(m_LZAMS, m_RZAMS);

    m_OmegaCHE                                 = CalculateOmegaCHE(m_MZAMS, m_Metallicity);
    m_OmegaZAMS                                = CalculateZAMSAngularFrequency(m_MZAMS, m_RZAMS);

    // Effective initial Zero Age Main Sequence parameters corresponding to Mass0
    m_RZAMS0                                   = m_RZAMS;
    m_LZAMS0                                   = m_LZAMS;

    // Current timestep attributes
    m_Time                                     = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Dt                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Tau                                      = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Age                                      = 0.0;           // ensure age = 0.0 at construction (rather than default initial value)
    m_Mass                                     = m_MZAMS;
    m_Mass0                                    = m_MZAMS;
    m_Luminosity                               = m_LZAMS;
    m_Radius                                   = m_RZAMS;
    m_Temperature                              = m_TZAMS;
    m_EnvMass                                  = CalculateInitialEnvelopeMass_Static(m_Mass);

    m_CoreMass                                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_COCoreMass                               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_HeCoreMass                               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mu                                       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_CoreRadius                               = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Mdot                                     = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_Omega                                    = m_OmegaZAMS;
    m_OmegaBreak                               = CalculateOmegaBreak();
    m_AngularMomentum                          = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_MomentOfInertia                          = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_MinimumLuminosityOnPhase                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_LBVphaseFlag                             = false;

    // Previous timestep attributes
    m_StellarTypePrev                          = m_StellarType;
    m_MassPrev                                 = m_MZAMS;
    m_RadiusPrev                               = m_RZAMS;
    m_DtPrev                                   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_OmegaPrev                                = m_OmegaZAMS;

    // Lambdas
	m_Lambdas.dewi                             = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.fixed                            = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.kruckow                          = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.kruckowBottom                    = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.kruckowMiddle                    = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.kruckowTop                       = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.loveridge                        = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.loveridgeWinds                   = DEFAULT_INITIAL_DOUBLE_VALUE;
	m_Lambdas.nanjing                          = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Zetas
    m_Zetas.hurley                             = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Zetas.hurleyHe                           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Zetas.nuclear                            = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Zetas.soberman                           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Zetas.sobermanHe                         = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_Zetas.thermal                            = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Binding energies
    m_BindingEnergies.fixed                    = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BindingEnergies.nanjing                  = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BindingEnergies.loveridge                = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BindingEnergies.loveridgeWinds           = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_BindingEnergies.kruckow                  = DEFAULT_INITIAL_DOUBLE_VALUE;

    // Supernova detais

    m_SupernovaDetails.initialKickParameters   = p_KickParameters;

    m_SupernovaDetails.events.current          = SN_EVENT::NONE;
    m_SupernovaDetails.events.past             = SN_EVENT::NONE;

    m_SupernovaDetails.coreMassAtCOFormation   = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.COCoreMassAtCOFormation = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.HeCoreMassAtCOFormation = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.totalMassAtCOFormation  = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.drawnKickVelocity       = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.kickVelocity            = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.hydrogenContent         = HYDROGEN_CONTENT::RICH;
    m_SupernovaDetails.fallbackFraction        = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.eccentricAnomaly        = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_SupernovaDetails.trueAnomaly             = DEFAULT_INITIAL_DOUBLE_VALUE;

    m_SupernovaDetails.supernovaState          = SN_STATE::NONE;

    if (p_KickParameters.supplied) {
        m_SupernovaDetails.kickVelocityRandom  = p_KickParameters.useVelocityRandom ? p_KickParameters.velocityRandom : DEFAULT_INITIAL_DOUBLE_VALUE;
        m_SupernovaDetails.theta               = p_KickParameters.theta;
        m_SupernovaDetails.phi                 = p_KickParameters.phi;
        m_SupernovaDetails.meanAnomaly         = p_KickParameters.meanAnomaly;
    }
    else {
        m_SupernovaDetails.kickVelocityRandom  = RAND->Random();
        std::tie(m_SupernovaDetails.theta, m_SupernovaDetails.phi) = DrawKickDirection();
        m_SupernovaDetails.meanAnomaly         = RAND->Random(0.0, _2_PI);
    }

    // Calculates the Baryonic mass for which the GravitationalRemnantMass will be equal to the maximumNeutronStarMass (inverse of SolveQuadratic())
    // needed to decide whether to calculate Fryer+2012 for Neutron Star or Black Hole in GiantBranch::CalculateGravitationalRemnantMass()
    m_BaryonicMassOfMaximumNeutronStarMass      = (0.075 * OPTIONS->MaximumNeutronStarMass() * OPTIONS->MaximumNeutronStarMass()) + OPTIONS->MaximumNeutronStarMass();
    // calculate only once for entire simulation of N binaries in the future.

    // Pulsar details
    m_PulsarDetails.magneticField              = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.spinPeriod                 = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.spinFrequency              = DEFAULT_INITIAL_DOUBLE_VALUE;
    m_PulsarDetails.spinDownRate               = DEFAULT_INITIAL_DOUBLE_VALUE;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  CLASS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Determine the value of the requested property of the constitusnt star (parameter p_Property)
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
 *    STAR_PROPERTY, STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY only.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  The type of the
 * functional return is a tuple: std::tuple<bool, COMPAS_VARIABLE_TYPE>.  This type
 * is COMPAS_VARIABLE by typedef.
 *
 * The bool returned indicates whether the property value was retrieved ok: true = yes, fales = no
 * The COMPAS_VARIABLE_TYPE variable returned is a boost variant variable, the value of which is the
 * value of the underlying primitive variable.
 *
 *
 * COMPAS_VARIABLE StellarPropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseStar::StellarPropertyValue(const T_ANY_PROPERTY p_Property) const {

    bool ok = true;                                                                                                     // default is no error

    COMPAS_VARIABLE_TYPE value;

    ANY_STAR_PROPERTY property;

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {

        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : { STAR_PROPERTY      prop = boost::get<STAR_PROPERTY>(p_Property);      property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : { STAR_1_PROPERTY    prop = boost::get<STAR_1_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : { STAR_2_PROPERTY    prop = boost::get<STAR_2_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: { SUPERNOVA_PROPERTY prop = boost::get<SUPERNOVA_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: { COMPANION_PROPERTY prop = boost::get<COMPANION_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;

        default:                                                                                                        // unknown property type
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            SHOW_WARN(ERROR::UNKNOWN_PROPERTY_TYPE);                                                                    // show warning
    }

    if (ok) {
        switch (property) {
            case ANY_STAR_PROPERTY::AGE:                                                value = Age();                                                  break;
            case ANY_STAR_PROPERTY::ANGULAR_MOMENTUM:                                   value = AngularMomentum();                                      break;
            case ANY_STAR_PROPERTY::BINDING_ENERGY_FIXED:                               value = BindingEnergy_Fixed();                                  break;
            case ANY_STAR_PROPERTY::BINDING_ENERGY_NANJING:                             value = BindingEnergy_Nanjing();                                break;
            case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE:                           value = BindingEnergy_Loveridge();                              break;
            case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS:                     value = BindingEnergy_LoveridgeWinds();                         break;
            case ANY_STAR_PROPERTY::BINDING_ENERGY_KRUCKOW:                             value = BindingEnergy_Kruckow();                                break;
            case ANY_STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE:               value = CHonMS();                                               break;
            case ANY_STAR_PROPERTY::CO_CORE_MASS:                                       value = COCoreMass();                                           break;
            case ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION:           value = SN_COCoreMassAtCOFormation();                           break;
            case ANY_STAR_PROPERTY::CORE_MASS:                                          value = CoreMass();                                             break;
            case ANY_STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION:              value = SN_CoreMassAtCOFormation();                             break;
            case ANY_STAR_PROPERTY::DRAWN_KICK_VELOCITY:                                value = SN_DrawnKickVelocity();                                 break;
            case ANY_STAR_PROPERTY::DT:                                                 value = Dt();                                                   break;
            case ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE:                                value = DynamicalTimescale();                                   break;
            case ANY_STAR_PROPERTY::ECCENTRIC_ANOMALY:                                  value = SN_EccentricAnomaly();                                  break;
            case ANY_STAR_PROPERTY::ENV_MASS:                                           value = EnvMass();                                              break;
            case ANY_STAR_PROPERTY::ERROR:                                              value = Error();                                                break;
            case ANY_STAR_PROPERTY::EXPERIENCED_CCSN:                                   value = ExperiencedCCSN();                                      break;
            case ANY_STAR_PROPERTY::EXPERIENCED_ECSN:                                   value = ExperiencedECSN();                                      break;
            case ANY_STAR_PROPERTY::EXPERIENCED_PISN:                                   value = ExperiencedPISN();                                      break;
            case ANY_STAR_PROPERTY::EXPERIENCED_PPISN:                                  value = ExperiencedPPISN();                                     break;
            case ANY_STAR_PROPERTY::EXPERIENCED_SN_TYPE:                                value = ExperiencedSN_Type();                                   break;
            case ANY_STAR_PROPERTY::EXPERIENCED_USSN:                                   value = ExperiencedUSSN();                                      break;
            case ANY_STAR_PROPERTY::FALLBACK_FRACTION:                                  value = SN_FallbackFraction();                                  break;
            case ANY_STAR_PROPERTY::HE_CORE_MASS:                                       value = HeCoreMass();                                           break;
            case ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION:           value = SN_HeCoreMassAtCOFormation();                           break;
            case ANY_STAR_PROPERTY::HYDROGEN_POOR:                                      value = SN_HydrogenContent() == HYDROGEN_CONTENT::POOR;         break;
            case ANY_STAR_PROPERTY::HYDROGEN_RICH:                                      value = SN_HydrogenContent() == HYDROGEN_CONTENT::RICH;         break;
            case ANY_STAR_PROPERTY::ID:                                                 value = ObjectId();                                             break;
            case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE:                               value = InitialStellarType();                                   break;
            case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE_NAME:                          value = STELLAR_TYPE_LABEL.at(InitialStellarType());            break;
            case ANY_STAR_PROPERTY::IS_CCSN:                                            value = IsCCSN();                                               break;
            case ANY_STAR_PROPERTY::IS_ECSN:                                            value = IsECSN();                                               break;
            case ANY_STAR_PROPERTY::IS_PISN:                                            value = IsPISN();                                               break;
            case ANY_STAR_PROPERTY::IS_PPISN:                                           value = IsPPISN();                                              break;
            case ANY_STAR_PROPERTY::IS_USSN:                                            value = IsUSSN();                                               break;
            case ANY_STAR_PROPERTY::KICK_VELOCITY:                                      value = SN_KickVelocity();                                      break;
            case ANY_STAR_PROPERTY::LAMBDA_DEWI:                                        value = Lambda_Dewi();                                          break;
            case ANY_STAR_PROPERTY::LAMBDA_FIXED:                                       value = Lambda_Fixed();                                         break;
            case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW:                                     value = Lambda_Kruckow();                                       break;
            case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM:                              value = Lambda_Kruckow();                                       break;
            case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE:                              value = Lambda_Kruckow();                                       break;
            case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_TOP:                                 value = Lambda_Kruckow();                                       break;
            case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE:                                   value = Lambda_Loveridge();                                     break;
            case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE_WINDS:                             value = Lambda_LoveridgeWinds();                                break;
            case ANY_STAR_PROPERTY::LAMBDA_NANJING:                                     value = Lambda_Nanjing();                                       break;
            case ANY_STAR_PROPERTY::LBV_PHASE_FLAG:                                     value = LBV_Phase_Flag();                                       break;
            case ANY_STAR_PROPERTY::LUMINOSITY:                                         value = Luminosity();                                           break;
            case ANY_STAR_PROPERTY::MASS:                                               value = Mass();                                                 break;
            case ANY_STAR_PROPERTY::MASS_0:                                             value = Mass0();                                                break;
            case ANY_STAR_PROPERTY::MDOT:                                               value = Mdot();                                                 break;
            case ANY_STAR_PROPERTY::MEAN_ANOMALY:                                       value = SN_MeanAnomaly();                                       break;
            case ANY_STAR_PROPERTY::METALLICITY:                                        value = Metallicity();                                          break;
            case ANY_STAR_PROPERTY::MZAMS:                                              value = MZAMS();                                                break;
            case ANY_STAR_PROPERTY::NUCLEAR_TIMESCALE:                                  value = NuclearTimescale();                                     break;
            case ANY_STAR_PROPERTY::OMEGA:                                              value = Omega();                                                break;
            case ANY_STAR_PROPERTY::OMEGA_BREAK:                                        value = OmegaBreak();                                           break;
            case ANY_STAR_PROPERTY::OMEGA_ZAMS:                                         value = OmegaZAMS();                                            break;
            case ANY_STAR_PROPERTY::PULSAR_MAGNETIC_FIELD:                              value = Pulsar_MagneticField();                                 break;
            case ANY_STAR_PROPERTY::PULSAR_SPIN_DOWN_RATE:                              value = Pulsar_SpinDownRate();                                  break;
            case ANY_STAR_PROPERTY::PULSAR_SPIN_FREQUENCY:                              value = Pulsar_SpinFrequency();                                 break;
            case ANY_STAR_PROPERTY::PULSAR_SPIN_PERIOD:                                 value = Pulsar_SpinPeriod();                                    break;
            case ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE:                         value = RadialExpansionTimescale();                             break;
            case ANY_STAR_PROPERTY::RADIUS:                                             value = Radius();                                               break;
            case ANY_STAR_PROPERTY::RANDOM_SEED:                                        value = RandomSeed();                                           break;
            case ANY_STAR_PROPERTY::RECYCLED_NEUTRON_STAR:                              value = ExperiencedRecycledNS();                                break;
            case ANY_STAR_PROPERTY::RLOF_ONTO_NS:                                       value = ExperiencedRLOFOntoNS();                                break;
            case ANY_STAR_PROPERTY::RUNAWAY:                                            value = ExperiencedRunaway();                                   break;
            case ANY_STAR_PROPERTY::RZAMS:                                              value = RZAMS();                                                break;
            case ANY_STAR_PROPERTY::SN_TYPE:                                            value = SN_Type();                                              break;
            case ANY_STAR_PROPERTY::STELLAR_TYPE:                                       value = StellarType();                                          break;
            case ANY_STAR_PROPERTY::STELLAR_TYPE_NAME:                                  value = STELLAR_TYPE_LABEL.at(StellarType());                   break;
            case ANY_STAR_PROPERTY::STELLAR_TYPE_PREV:                                  value = StellarTypePrev();                                      break;
            case ANY_STAR_PROPERTY::STELLAR_TYPE_PREV_NAME:                             value = STELLAR_TYPE_LABEL.at(StellarTypePrev());               break;
            case ANY_STAR_PROPERTY::SUPERNOVA_KICK_VELOCITY_MAGNITUDE_RANDOM_NUMBER:    value = SN_KickVelocityRandom();                                break;
            case ANY_STAR_PROPERTY::SUPERNOVA_PHI:                                      value = SN_Phi();                                               break;
            case ANY_STAR_PROPERTY::SUPERNOVA_THETA:                                    value = SN_Theta();                                             break;
            case ANY_STAR_PROPERTY::TEMPERATURE:                                        value = Temperature();                                          break;
            case ANY_STAR_PROPERTY::THERMAL_TIMESCALE:                                  value = ThermalTimescale();                                     break;
            case ANY_STAR_PROPERTY::TIME:                                               value = Time();                                                 break;
            case ANY_STAR_PROPERTY::TIMESCALE_MS:                                       value = Timescale(TIMESCALE::tMS);                              break;
            case ANY_STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION:             value = SN_TotalMassAtCOFormation();                            break;
            case ANY_STAR_PROPERTY::TRUE_ANOMALY:                                       value = SN_TrueAnomaly();                                       break;
            case ANY_STAR_PROPERTY::ZETA_HURLEY:                                        value = Zeta_Hurley();                                          break;
            case ANY_STAR_PROPERTY::ZETA_HURLEY_HE:                                     value = Zeta_HurleyHe();                                        break;
            case ANY_STAR_PROPERTY::ZETA_NUCLEAR:                                       value = Zeta_Nuclear();                                         break;
            case ANY_STAR_PROPERTY::ZETA_SOBERMAN:                                      value = Zeta_Soberman();                                        break;
            case ANY_STAR_PROPERTY::ZETA_SOBERMAN_HE:                                   value = Zeta_SobermanHe();                                      break;
            case ANY_STAR_PROPERTY::ZETA_THERMAL:                                       value = Zeta_Thermal();                                         break;

        default:                                                                                                        // unknown property
            ok    = false;                                                                                              // that's not ok...
            value = "UNKNOWN";                                                                                          // default value
            SHOW_WARN(ERROR::UNKNOWN_STELLAR_PROPERTY);                                                                 // show warning
        }
    }

    return std::make_tuple(ok, value);
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
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  The type of the
 * functional return is a tuple: std::tuple<bool, COMPAS_VARIABLE_TYPE>.  This type
 * is COMPAS_VARIABLE by typedef.
 *
 * The bool returned indicates whether the property value was retireved ok: true = yes, fales = no
 * The COMPAS_VARIABLE_TYPE variable returned is a boost variant variable, the value of which is the
 * value of the underlying primitive variable.
 *
 *
 * COMPAS_VARIABLE PropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseStar::PropertyValue(const T_ANY_PROPERTY p_Property) const {

    bool ok = true;                                                                                                             // default is no error

    COMPAS_VARIABLE_TYPE value = 0.0;                                                                                           // default property value

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {                                                          // which property type?

        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY:                                                                                // star property
            std::tie(ok, value) = StellarPropertyValue(p_Property);
            break;

        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION:                                                                               // program option
            std::tie(ok, value) = OPTIONS->OptionValue(p_Property);
            break;

        default:                                                                                                                // unknown property type
            ok    = false;                                                                                                      // that's not ok...
            value = "UNKNOWN";                                                                                                  // default value
            SHOW_WARN(ERROR::UNKNOWN_PROPERTY_TYPE  );                                                                          // show warning
    }

    return std::make_tuple(ok, value);
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
    double xi    = m_LogMetallicityXi;
    double sigma = m_LogMetallicitySigma;

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
    a[17] = pow(10.0, max((0.097 - (0.1072 * (sigma + 3.0))), max(0.097, min(0.1461, (0.1461 + (0.1237 * (sigma + 2.0)))))));
    a[18] *= a[20];
    a[19] *= a[20];
    a[29] = pow(a[29], (a[32]));
    a[33] = min(1.4, 1.5135 + (0.3769 * xi));
    a[42] = min(1.25, max(1.1, a[42]));
    a[44] = min(1.3, max(0.45, a[44]));
    a[49] = max(a[49], 0.145);
    a[50] = min(a[50], (0.306 + (0.053 * xi)));
    a[51] = min(a[51], (0.3625 + (0.062 * xi)));
    a[52] = (utils::Compare(Z, 0.01) > 0) ? min(a[52], 1.0) : max(a[52], 0.9);
    a[53] = (utils::Compare(Z, 0.01) > 0) ? min(a[53], 1.1) : max(a[53], 1.0);
    a[57] = min(1.4, a[57]);
    a[57] = max((0.6355 - (0.4192 * xi)), max(1.25, a[57]));
    a[62] = max(0.065, a[62]);
    a[63] = (utils::Compare(Z, 0.004) < 0) ? min(0.055, a[63]) : a[63];
    a[66] = max(a[66], min(1.6, -0.308 - (1.046 * xi)));
    a[66] = max(0.8, min(0.8 - (2.0 * xi), a[66]));
    a[68] = max(0.9, min(a[68], 1.0));

    // Need bAlpaR - calculate it now
    RConstants(B_ALPHA_R) = (a[58] * pow(a[66], a[60])) / (a[59] + pow(a[66], a[61]));                          // Hurley et al. 2000, eq 21a (wrong in the arxiv version - says = a59*M**(a61))

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

    LConstants(B_ALPHA_L)   = (a[45] + (a[46] * pow(2.0, a[48]))) / (pow(2.0, 0.4) + (a[47] * pow(2.0, 1.9)));  // Hurley et al. 2000, eq 19a
    LConstants(B_BETA_L)    = max(0.0, (a[54] - (a[55] * pow(a[57], a[56]))));                                  // Hurley et al. 2000, eq 20
    LConstants(B_DELTA_L)   = min((a[34] / pow(a[33], a[35])), (a[36] / pow(a[33], a[37])));                    // Hurley et al. 2000, eq 16

    RConstants(C_ALPHA_R)   = (a[58] * pow(a[67], a[60])) / (a[59] + pow(a[67], a[61]));                        // Hurley et al. 2000, eq 21a (wrong in the arxiv version)
    RConstants(B_BETA_R)    = (a[69] * 8.0 * M_SQRT2) / (a[70] + pow(2.0, a[71]));                              // Hurley et al. 2000, eq 22a
    RConstants(C_BETA_R)    = (a[69] * 16384.0) / (a[70] + pow(16.0, a[71]));                                   // Hurley et al. 2000, eq 22a
    RConstants(B_DELTA_R)   = (a[38] + (a[39] * 8.0 * M_SQRT2)) / ((a[40] * 8.0) + pow(2.0, a[41]));            // Hurley et al. 2000, eq 17

    GammaConstants(B_GAMMA) = a[76] + (a[77] * pow((1.0 - a[78]), a[79]));                                      // Hurley et al. 2000, eq 23
    GammaConstants(C_GAMMA) = (utils::Compare(a[75], 1.0) == 0) ? GammaConstants(B_GAMMA) : a[80];              // Hurley et al. 2000, eq 23

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
    double xi    = m_LogMetallicityXi;
    double sigma = m_LogMetallicitySigma;
    double rho   = m_LogMetallicityRho;

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
    b[2] = pow(10.0, (-4.6739 - (0.9394 * sigma)));
    b[2] = min(max(b[2], (-0.04167 + (55.67 * Z))), (0.4771 - (9329.21 * pow(Z, 2.94))));
    b[3] = max(-0.1451, (-2.2794 - (1.5175 * sigma) - (0.254 * sigma * sigma)));
    b[3] = (utils::Compare(Z, 0.004) > 0) ? max(b[3], 0.7307 + (14265.1 * pow(Z, 3.395))) : pow(10.0, b[3]);
    b[4] += 0.1231572 * xi_5;
    b[6] += 0.01640687 * xi_5;
    b[11] = b[11] * b[11];
    b[13] = b[13] * b[13];
    b[14] = pow(b[14], b[15]);
    b[16] = pow(b[16], b[15]);
    b[17] = (utils::Compare(xi, -1.0) > 0) ? 1.0 - (0.3880523 * pow((xi + 1.0), 2.862149)) : 1.0;
    b[24] = pow(b[24], b[28]);
    b[26] = 5.0 - (0.09138012 * pow(Z, -0.3671407));
    b[27] = pow(b[27], (2.0 * b[28]));
    b[31] = pow(b[31], b[33]);
    b[34] = pow(b[34], b[33]);
    b[36] = b[36] * b[36] * b[36] * b[36];
    b[37] = 4.0 * b[37];
    b[38] = b[38] * b[38] * b[38] * b[38];
    b[40] = max(b[40], 1.0);
    b[41] = pow(b[41], b[42]);
    b[44] = b[44] * b[44] * b[44] * b[44] * b[44];
    b[45] = utils::Compare(rho, 0.0) <= 0 ? 1.0 : 1.0 - ((2.47162 * rho) - (5.401682 * rho_2) + (3.247361 * rho_3));
    b[46] = -1.0 * b[46] * log10(massCutoffs(MHeF) / massCutoffs(MFGB));
    b[47] = (1.127733 * rho) + (0.2344416 * rho_2) - (0.3793726 * rho_3);
    b[51] -= 0.1343798 * xi_5;
    b[53] += 0.4426929 * xi_5;
    b[55] = min((0.99164 - (743.123 * pow(Z, 2.83))), b[55]);
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
double BaseStar::CalculateAlpha1() {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double LHeI_MHeF = (b[11] + (b[12] * pow(massCutoffs(MHeF), 3.8))) / (b[13] + (massCutoffs(MHeF) * massCutoffs(MHeF)));
    return ((b[9] * pow(massCutoffs(MHeF), b[10])) - LHeI_MHeF) / LHeI_MHeF;

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
double BaseStar::CalculateAlpha3() {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double LBAGB = (b[31] + (b[32] * pow(massCutoffs(MHeF), (b[33] + 1.8)))) / (b[34] + pow(massCutoffs(MHeF), b[33]));
    return ((b[29] * pow(massCutoffs(MHeF), b[30])) - LBAGB) / LBAGB;

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
double BaseStar::CalculateAlpha4() {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double MHeF      = massCutoffs(MHeF);
    double MHeF_5    = MHeF * MHeF * MHeF * MHeF * MHeF;    // pow() is slow - use multiplication
    double tBGB_MHeF = CalculateLifetimeToBGB(MHeF);        // tBGB for mass M = MHeF

    return tBGB_MHeF * (((b[41] * pow(MHeF, b[42])) + (b[43] * MHeF_5)) / (b[44] + MHeF_5));

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

    double top         = 13.048 * pow((p_Metallicity / ZSOL), 0.06);
    double bottom      = 1.0 + (0.0012 * pow((ZSOL / p_Metallicity), 1.27));
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
 * @return                                      'x' exponent to which Radius depends on Mass (at conatant Luminosity)- 'x' in Hurley et al. 2000, eq 47
 */
double BaseStar::CalculateGBRadiusXExponent() {

    // pow()is slow - use multiplication
    double xi   = m_LogMetallicityXi;
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
double BaseStar::CalculatePerturbationB(const double p_Mass) {
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
double BaseStar::CalculatePerturbationC(double p_Mass) {
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
double BaseStar::CalculatePerturbationS(const double p_Mu, const double p_Mass) {

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
double BaseStar::CalculatePerturbationQ(const double p_Radius, const double p_Rc) {
    return log(p_Radius / p_Rc); // really is natural log
}


/*
 * Calculate the perturbation parameter r
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
double BaseStar::CalculatePerturbationR(const double p_Mu, const double p_Mass, const double p_Radius, const double p_Rc) {

    double r = 0.0;

    if(utils::Compare(p_Mu, 0.0) >= 0) {                            // only if mu >= 0

        double c      = CalculatePerturbationC(p_Mass);
        double c_3    = c * c * c;                                  // pow() is slow - use multiplication
        double mu_c_3 = p_Mu * p_Mu * p_Mu / c_3;                   // calculate once

        double q        = CalculatePerturbationQ(p_Radius, p_Rc);
        double exponent = min((0.1 / q), (-14.0 / log10(p_Mu)));    // JR: todo: Hurley et al. 2000 is just 0.1 / q ?

        r = ((1.0 + c_3) * mu_c_3 * pow((p_Mu), exponent)) / ((1.0 + mu_c_3));
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
 * Made by Alejandro Vigna-Gomez
 *
 *
 * double CalculateLambdaKruckow(const double p_Radius, const double p_Alpha)
 *
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @param   [IN]    p_Alpha                     Power
 * @return                                      Common envelope lambda parameter
 */
double BaseStar::CalculateLambdaKruckow(const double p_Radius, const double p_Alpha) {

	double alpha = max(-2.0 / 3.0, min(-1.0, p_Alpha));             // clamp alpha to [-1.0, -2/3]

	return 1600.0 * pow(0.00125, -alpha) * pow(p_Radius, alpha);
}


/*
 * Calculate the binding energy of the envelope
 * Loveridge et al. 2011
 *
 * This function computes log[BE/erg] as a function of log[Z], Mzams, M, log[R/Ro] and GB.
 * Electronic tables, program and further information in: http://astro.ru.nl/~sluys/index.php?title=BE
 *
 *
 * double CalculateLogBindingEnergyLoveridge(bool p_IsMassLoss)
 *
 * @param   [IN]    p_IsMassLoss                Boolean indicating whether mass-loss correction should be applied
 * @return                                      log binding energy in ergs
 */
double BaseStar::CalculateLogBindingEnergyLoveridge(bool p_IsMassLoss) {

    // find closest metallicity covered by Loveridge et al. 2011
    // (see LOVERIDGE_METALLICITY and LOVERIDGE_METALLICITYValue)

    int lMetallicity = 0;
    double minDiff   = std::numeric_limits<double>::max();

    // initialise m_MassCutoffs vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(LOVERIDGE_METALLICITY::COUNT); i++) {
        double thisDiff = std::abs(m_Metallicity - std::get<1>(LOVERIDGE_METALLICITY_VALUE[i]));
        if (utils::Compare(thisDiff, minDiff) < 0) {
            lMetallicity = i;
            minDiff      = thisDiff;
        }
    }

    // Determine the evolutionary stage of the star (see LOVERIDGE_GROUP)

    LOVERIDGE_GROUP lGroup;

    if (utils::Compare(m_Mass, LOVERIDGE_LM_HM_CUTOFFS[lMetallicity]) > 0) {                // mass > low mass / high mass cutoff?
        lGroup = LOVERIDGE_GROUP::HM;                                                       // yes, group is HM - High Mass
    }
    else {                                                                                  // no - low mass
        if (utils::Compare(m_COCoreMass, 0.0) > 0) {                                        // CO core exists?
            lGroup = LOVERIDGE_GROUP::LMA;                                                  // yes, group is LMA - Low mass on the AGB
        }
        else {                                                                              // no - low mass star on RGB

            // calculate early / late cutoff for low mass RGB stars
            constexpr double deltaM   = 1.0E-5;                                             // JR: todo: what is this for?  Should it be in constants.h?
                      double cutOff   = 0.0;
                      int    exponent = 0;
            for (auto const& aCoefficient: LOVERIDGE_LM1_LM2_CUTOFFS[lMetallicity]) {
                cutOff += aCoefficient * utils::intPow(log10(m_Mass + deltaM), exponent++);
            }

            // set evolutionary stage based on cutoff
            lGroup = utils::Compare(log10(m_Radius), cutOff) > 0 ? LOVERIDGE_GROUP::LMR2 : LOVERIDGE_GROUP::LMR1;
        }
    }

    // calculate log10(binding energy)
    constexpr double deltaR           = 1E-5;                                               // JR: todo: what is this for?  Should it be in constants.h?
              double logBindingEnergy = 0.0;
    for (auto const& lCoefficients: LOVERIDGE_COEFFICIENTS[lMetallicity][static_cast<int>(lGroup)]) {
        logBindingEnergy += lCoefficients.alpha_mr * utils::intPow(log10(m_Mass), lCoefficients.m) * utils::intPow(log10(m_Radius + deltaR), lCoefficients.r);
    }

    double MZAMS_Mass = (m_MZAMS - m_Mass) / m_MZAMS;
    logBindingEnergy *= p_IsMassLoss ? 1.0 + (0.25 * MZAMS_Mass * MZAMS_Mass) : 1.0;        // apply mass-loss correction factor (lambda)

    constexpr double logBE0 = 33.29866;                                                     // JR: todo: what is this for?  Should it be in constants.h?
    logBindingEnergy += logBE0;

	return logBindingEnergy;
}


/*
 * Calculata lambda parameter from the so-called energy formalism of CE (Webbink 1984).
 *
 * Binding energy from detailed models (Loveridge et al. 2011) is given in [E]=ergs, so use cgs
 *
 *
 * double CalculateLambdaLoveridgeEnergyFormalism(const double p_EnvMass, const double p_IsMassLoss)
 *
 * @param   [IN]    p_EnvMass                   Envelope mass (Msol)
 * @param   [IN]    p_IsMassLoss                Boolean indicating whether mass-loss correction should be applied
 * @return                                      Common envelope lambda parameter
 */
double BaseStar::CalculateLambdaLoveridgeEnergyFormalism(const double p_EnvMass, const double p_IsMassLoss) {   // JR: todo: p_IsMassLoss is ignored - do we need it?  I've set the default to false

    double bindingEnergy = pow(10.0, CalculateLogBindingEnergyLoveridge(false));
    return bindingEnergy > 0.0 ? (G_CGS * m_Mass * MSOL_TO_G * p_EnvMass * MSOL_TO_G) / (m_Radius * RSOL_TO_AU * AU_TO_CM * bindingEnergy) : pow(10.0, -20);   // JR: todo: fix this
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 ZETA CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the radius-time exponent zeta
 * Derived by the nuclear evolution of the star according to SSE, calculated using a fake time step and recomputing the stars radius.
 *
 *
 * double CalculateZetaNuclear(const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaTime                 Timestep (Myr)
 * @return                                      Radius-Time exponent Zeta for the nuclear timescale
 */
double BaseStar::CalculateZetaNuclear(const double p_DeltaTime) {

    BaseStar* starCopy = new BaseStar(*this);	                                                                                // copy of star - about to be updated for fake mass loss

    double radiusBeforeTimeStep = starCopy->Radius();                                                                           // radius before timestep
    double ageBeforeTimeStep    = starCopy->Age();

    SHOW_ERROR_IF(utils::Compare(radiusBeforeTimeStep, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "Before fake timestep");     // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(ageBeforeTimeStep,    0.0) <  0, ERROR::AGE_NEGATIVE_ONCE,        "Before fake timestep");     // show error if age < 0

    double logRadiusBeforeTimeStep = utils::Compare(radiusBeforeTimeStep, 0.0) > 0 ? log(radiusBeforeTimeStep) : 0.0;

    starCopy->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);                                                           // allow star to respond to previous mass loss changes      JR: todo: is this really necessary?

    // With the current way of doing things, this star has already lost a bunch of mass due to the mass loss caller,
    // need to update radius of star to respond to this, before we evolve it for a short amount of time assuming no mass loss.
    // When Alejandro moves the massLossCaller function, we will need to do things differently, and this will need fixing.
    // JR: todo: what does this actually mean?

    ageBeforeTimeStep = starCopy->Age();                                                                                        // age before timestep      JR: todo: why don't we get the (possibly) new radius here too?
    starCopy->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, p_DeltaTime, false);                                                  // age the star 'p_DeltaTime' Myr

    double radiusAfterTimeStep = starCopy->Radius();
    double ageAfterTimeStep    = starCopy->Age();

    delete starCopy; starCopy = nullptr;

    SHOW_ERROR_IF(utils::Compare(radiusAfterTimeStep, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "After fake timestep");       // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(ageAfterTimeStep,    0.0) <  0, ERROR::AGE_NEGATIVE_ONCE,        "After fake timestep");       // show error if age < 0

    double logRadiusAfterTimeStep = utils::Compare(radiusAfterTimeStep, 0.0) > 0 ? log(radiusAfterTimeStep) : 0.0;

    return (logRadiusAfterTimeStep - logRadiusBeforeTimeStep) / (ageAfterTimeStep - ageBeforeTimeStep);

}


/*
 * Calculate the radius response of a star to mass loss on the nuclear timescale as characterised
 * by the mass-radius exponent Zeta nuclear
 *
 * This calculation is good for all stars except Black Holes
 *
 *
 * double CalculateConvergedTimestepZetaNuclear()
 *
 * @return                                      Mass-radius exponent Zeta nuclear (= dlnR/dt)
 */
double BaseStar::CalculateConvergedTimestepZetaNuclear() {

    double timestep  = ZETA_NUCLEAR_TIMESTEP;                                                                                   // timestep
    double tolerance = 1.0 - ZETA_NUCLEAR_TOLERANCE;                                                                            // max change allowed in zeta nuclear between iterations - signals convergence

    double thisZetaNuclear;                                                                                                     // value of mass-radius exponent calculated at current iteration
    double lastZetaNuclear = 0.0;                                                                                               // value of mass-radius exponent calculated at previous iteration

    int    i    = 0;
    bool   done = false;
    while (!done) {                                                                                                             // iterate until convergence or max iterations reached

        thisZetaNuclear = CalculateZetaNuclear(timestep);                                                                       // calculate zeta nuclear for this percentage mass loss

        if (i >= ZETA_NUCLEAR_ITERATIONS || std::abs(lastZetaNuclear / thisZetaNuclear) > tolerance) {                          // max iterations or convergence?   (Don't use Compare() here)
            done = true;                                                                                                        // yes - we're done
        }
        else {                                                                                                                  // no - go again
            lastZetaNuclear = thisZetaNuclear;                                                                                  // remember this zeta nuclear value for comparison next iteration
            timestep       /= ZETA_NUCLEAR_ITERATIONS;                                                                          // reduce timestep   JR: todo: is this factor right?
            i++;                                                                                                                // increment number of iterations
        }
    }

    SHOW_ERROR_IF(std::abs(lastZetaNuclear / thisZetaNuclear) <= tolerance, ERROR::NO_CONVERGENCE, "Zeta Nuclear calculation"); // show error if no convergence

    return thisZetaNuclear;
}


/*
 * Calculate the radius-mass exponent Zeta, assuming the star has had time to recover its thermal equilibrium.
 * Zeta is calculated using a fake mass change step and recomputing the star's attributes.
 *
 *
 * double CalculateZetaThermal(double p_PercentageMassChange)
 *
 * @param   [IN]    p_PercentageMassChange      Percentage of mass the star should artificially lose/gain in order to calcluate Zeta (thermal)
 *                                              (sign of p_PercentageMassChange determine if loss or gain - -ve is loss, +ve is gain)
 * @return                                      Radius-Mass exponent Zeta for the thermal timescale
 */
double BaseStar::CalculateZetaThermal(double p_PercentageMassChange) {

    BaseStar* starCopy = new BaseStar(*this);	                                                                                // copy of star - about to be updated for fake mass loss

    starCopy->UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);                                                           // allow star to respond to previous mass loss changes      JR: todo: is this really necessary?
    starCopy->UpdateAttributesAndAgeOneTimestep(-(starCopy->Mass() * FAKE_MASS_LOSS_PERCENTAGE / 100.0), 0.0, 0.0, false);      // apply fake mass loss                                     JR: todo: why do we do this...?

    // record properties of the star before fake mass change
    double radiusBeforeMassLoss = starCopy->Radius();                                                                           // radius before fake mass change       JR: todo: didn't we just do fake mass loss...?
    double massBeforeMassLoss   = starCopy->MassPrev();                                                                         // mass before fake mass change - due to order of updating radius bug     JR: todo: check this

    SHOW_ERROR_IF(utils::Compare(radiusBeforeMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "Before fake mass change");  // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(massBeforeMassLoss,   0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE,   "Before fake mass change");  // show error if mass <= 0

    double logRadiusBeforeMassLoss = utils::Compare(radiusBeforeMassLoss, 0.0) > 0 ? log(radiusBeforeMassLoss) : 0.0;
    double logMassBeforeMassLoss   = utils::Compare(massBeforeMassLoss,   0.0) > 0 ? log(massBeforeMassLoss  ) : 0.0;

    double deltaMass           = starCopy->Mass() * p_PercentageMassChange / 100.0;                                             // fake mass change - sign of p_PercentageMassChange determines whether loss or gain
    double massAfterMassLoss   = starCopy->Mass() + deltaMass;                                                                  // mass after (just) fake mass change
    starCopy->UpdateAttributesAndAgeOneTimestep(starCopy->Mass() * p_PercentageMassChange / 100.0, 0.0, 0.0, false);            // apply fake mass change and recalculate attributes of star
    double radiusAfterMassLoss = starCopy->m_Radius;                                                                            // radius after fake mass change
    delete starCopy; starCopy  = nullptr;

    SHOW_ERROR_IF(utils::Compare(radiusAfterMassLoss, 0.0) <= 0, ERROR::RADIUS_NOT_POSITIVE_ONCE, "After fake mass change");    // show error if radius <= 0
    SHOW_ERROR_IF(utils::Compare(massAfterMassLoss,   0.0) <= 0, ERROR::MASS_NOT_POSITIVE_ONCE,   "After fake mass change");    // show error if mass <= 0

    double logRadiusAfterMassLoss = utils::Compare(radiusAfterMassLoss, 0.0) > 0 ? log(radiusAfterMassLoss) : 0.0;
    double logMassAfterMassLoss   = utils::Compare(massAfterMassLoss,   0.0) > 0 ? log(massAfterMassLoss  ) : 0.0;

    return (logRadiusAfterMassLoss - logRadiusBeforeMassLoss) / (logMassAfterMassLoss - logMassBeforeMassLoss);                 // zeta thermal
}


/*
 * Calculate the radius response of a star to mass loss on the thermal timescale as characterised
 * by the mass-radius exponent Zeta thermal
 *
 * This calculation is good for MS, HG, HeMS, HeWD, COWD, and ONeWD stars
 *
 *
 * double CalculateConvergedMassStepZetaThermal()
 *
 * @return                                      Mass-radius exponent Zeta thermal (= dlnR/dlnM)
 */
double BaseStar::CalculateConvergedMassStepZetaThermal() {

    double percentageMassLoss = -ZETA_THERMAL_PERCENTAGE_MASS_CHANGE;                                                       // percentage mass loss
    double tolerance          = 1.0 - ZETA_THERMAL_TOLERANCE;                                                               // max change allowed in zeta thermal between iterations - signals convergence

    double thisZetaThermal;                                                                                                 // value of mass-radius exponent calculated at current iteration
    double lastZetaThermal = 0.0;                                                                                           // value of mass-radius exponent calculated at previous iteration

    int    i    = 0;
    bool   done = false;
    while (!done) {                                                                                                         // iterate until convergence or max iterations reached

        thisZetaThermal = CalculateZetaThermal(percentageMassLoss);                                                         // calculate zeta thermal for this percentage mass loss

        if (i >= ZETA_THERMAL_ITERATIONS || std::abs(lastZetaThermal / thisZetaThermal) > tolerance) {                      // max iterations or convergence?   (Don't use Compare() here)
            done = true;                                                                                                    // yes - we're done
        }
        else {                                                                                                              // no - go again
            lastZetaThermal     = thisZetaThermal;                                                                          // remember this zeta thermal value for comparison next iteration
            percentageMassLoss /= ZETA_THERMAL_ITERATIONS;                                                                  // reduce mass loss   JR: todo: is this factor right?
            i++;                                                                                                            // increment number of iterations
        }
    }

    SHOW_ERROR_IF(std::abs(lastZetaThermal / thisZetaThermal) <= tolerance, ERROR::NO_CONVERGENCE, "Zeta Thermal calculation");    // show error if no convergence

    return thisZetaThermal;
}


/*
 * Calculate the Adiabatic Exponent per Hurley et al. 2002
 *
 *
 * double CalculateZadiabaticHurley2002(const double p_CoreMass)
 *
 * @param   [IN]    p_CoreMass                  Core mass of the star (Msol)
 * @return                                      Adiabatic exponent
 */
double BaseStar::CalculateZadiabaticHurley2002(const double p_CoreMass) {
    double m = p_CoreMass / m_Mass;
    double x = -0.3;                                    // Depends on composition, should use x from Hurley et al 2000
    return -x + (2.0 * m * m * m * m * m);
}


/*
 * Calculate the Adiabatic Exponent per Soberman, Phinney, vdHeuvel 1997
 *
 *
 * double CalculateZadiabaticSPH(const double p_CoreMass)
 *
 * @param   [IN]    p_CoreMass                  Core mass of the star (Msol)
 * @return                                      Adiabatic exponent
 */
double BaseStar::CalculateZadiabaticSPH(const double p_CoreMass) {
    double m           = p_CoreMass / m_Mass;                                                                                                       // eq (57) Soberman, Phinney, vdHeuvel (1997)
    double oneMinusM   = 1.0 - m;
    double oneMinusM_6 = oneMinusM * oneMinusM * oneMinusM * oneMinusM * oneMinusM * oneMinusM;
    return ((2.0 / 3.0) * m / oneMinusM) - ((1.0 / 3.0) * (oneMinusM / (1.0 + (m + m)))) - (0.03 * m) + (0.2 * m / (1.0 + (1.0 / oneMinusM_6)));    // eq (61) Soberman, Phinney, vdHeuvel (1997)
}


/*
 * Calculate all Lambdas
 *
 * ALEJANDRO - 04/11/2016 - Lambda calculations as tracker for binding energy;
 *
 * lambda is recalculated in the CommonEnvelopeEvent function.                                          // JR: todo: different envelope mass in CEE function?  Esp. for loveridge?
 * Maybe we can give it as an argument later, as it is already being calculated here.                   // JR: todo: different envelope mass in CEE function?  Esp. for loveridge?
 *
 *
 * void CalculateLambdas(const double p_EnvMass)
 *
 * @param   [IN]    p_EnvMass                   Envelope mass of the star (Msol)
 */
void BaseStar::CalculateLambdas(const double p_EnvMass) {

    m_Lambdas.fixed          = OPTIONS->CommonEnvelopeLambda();
	m_Lambdas.nanjing        = CalculateLambdaNanjing();
	m_Lambdas.loveridge      = CalculateLambdaLoveridgeEnergyFormalism(p_EnvMass, false);     // JR: todo: (1) arg 2 (ismassloss) is ignored
	m_Lambdas.loveridgeWinds = CalculateLambdaLoveridgeEnergyFormalism(p_EnvMass, true);      // JR: todo: (1) arg 2 (ismassloss) is ignored
	m_Lambdas.kruckow        = CalculateLambdaKruckow(m_Radius, OPTIONS->CommonEnvelopeSlopeKruckow());
	m_Lambdas.kruckowTop     = CalculateLambdaKruckow(m_Radius, -2.0 / 3.0);
	m_Lambdas.kruckowMiddle  = CalculateLambdaKruckow(m_Radius, -4.0 / 5.0);
	m_Lambdas.kruckowBottom  = CalculateLambdaKruckow(m_Radius, -1.0);
	m_Lambdas.dewi           = CalculateLambdaDewi();
}


/*
 * Calculate all Zetas
 *
 * ALEJANDRO - 16/10/2017 - Zeta calculations as a mass transfer tracker.
 *
 *
 * void CalculateZetas()
 */
void BaseStar::CalculateZetas() {
	m_Zetas.thermal    = CalculateConvergedMassStepZetaThermal();
	m_Zetas.nuclear    = CalculateConvergedTimestepZetaNuclear();
	m_Zetas.soberman   = CalculateZadiabaticSPH(m_CoreMass);
	m_Zetas.sobermanHe = CalculateZadiabaticSPH(m_HeCoreMass);
	m_Zetas.hurley     = CalculateZadiabaticHurley2002(m_CoreMass);
	m_Zetas.hurleyHe   = CalculateZadiabaticHurley2002(m_HeCoreMass);
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
    double m_0_5 = sqrt(p_MZAMS);
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
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at BAGB in Lsol
 */
double BaseStar::CalculateLuminosityAtBAGB(double p_Mass) {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    return (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0)
            ? (b[29] * pow(p_Mass, b[30])) / (1.0 + (m_Alpha3 * exp(15.0 * (p_Mass - massCutoffs(MHeF)))))
            : (b[31] + (b[32] * pow(p_Mass, (b[33] + 1.8)))) / (b[34] + pow(p_Mass, b[33]));

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
double BaseStar::CalculateLuminosityGivenCoreMass(const double p_CoreMass) {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    return min((gbParams(B) * pow(p_CoreMass, gbParams(q))), (gbParams(D) * pow(p_CoreMass, gbParams(p))));

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
double BaseStar::CalculateRadiusAtZAMS(const double p_MZAMS) {
#define coeff(x) m_RCoefficients[static_cast<int>(R_Coeff::x)]  // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication where it makes sense
    // sqrt() is much faster than pow()
    double m_0_5  = sqrt(p_MZAMS);
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
 * double CalculateMaximumCoreMass(double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Maximum core mass in Msol (McMax)
 */
double BaseStar::CalculateMaximumCoreMass(double p_Mass) {
    return min(((1.45 * p_Mass) - 0.31), p_Mass);
}


/*
 * Calculate the core mass at which the AGB phase is terminated in a SN/loss of envelope
 *
 * Hurley et al. 2000, eq 75
 *
 *
 * double Star::CalculateMaximumCoreMassSN()
 *
 * @return                                      Maximum core mass before supernova (McSN)
 */
double BaseStar::CalculateMaximumCoreMassSN() {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]        // for convenience and readability - undefined at end of function

    return max(MECS, (0.773 * gbParams(McBAGB)) - 0.35);        // Mch constant in constants.h

#undef gbParams
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
 * Calculate envelope mass on phase as a function of time
 *
 * Given just after Eq 111 in Hurley et al. 2000
 *
 * This function works for most phases.  Stellar types MS_lte_07, MS_gt_07 and HertzsprungGap
 * have specialised functions.
 *
 *
 * double CalculateEnvelopeMassOnPhase(const double p_Tau)
 *
 * @param   [IN]    p_Tau                       Relative lifetime
 * @return                                      ZAMS envelope mass - Menv in Hurley et al. 2000
 *
 * Parameter p_Tau not required here - but we want same signature for all classes.  (JR: revisit this)
 */
double BaseStar::CalculateEnvelopeMassOnPhase(const double p_Tau) {
    return m_Mass - m_CoreMass; // For most phases with core envelope separation
}


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * JR: This function returns 1.0 in all cases....
 *
 *
 * double CalculateMassTransferRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double BaseStar::CalculateMassTransferRejuvenationFactor() {

    double fRej;
    switch (OPTIONS->MassTransferRejuvenationPrescription()) {                          // which prescription

        case MT_REJUVENATION_PRESCRIPTION::NONE:                                        // NONE - use default Hurley et al. 2000 prescription = 1.0
            fRej = 1.0;
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                   // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf

            if (utils::Compare(m_Mass, m_MassPrev) <= 0) {                              // Rejuvenation factor is unity for mass losing stars
                fRej = 1.0;
            }
            else {
                fRej = 1.0;
            }
            break;

        default:                                                                        // unknown prescription - use default Hurley et al. 2000 prescription = 1.0
            SHOW_WARN(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION, "Using default fRej = 1.0");     // show warning
            fRej = 1.0;
    }

    return fRej;
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
double BaseStar::CalculateMassLossRateVassiliadisWood() {

    double logP0      = min(3.3, (-2.07 - (0.9 * log10(m_Mass)) + (1.94 * log10(m_Radius))));
    double P0         = pow(10.0, (logP0)); // In their fortran code, Hurley et al take P0 to be min(p0, 2000.0), implemented here as a minimum power
    double logMdot_VW = -11.4 + (0.0125 * (P0 - 100.0 * max((m_Mass - 2.5), 0.0)));
    double Mdot_VW    = pow(10.0, (logMdot_VW));

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
double BaseStar::CalculateMassLossRateKudritzkiReimers() {
    return 4.0E-13 * (MASS_LOSS_ETA * m_Luminosity * m_Radius / m_Mass);    // Shouldn't be eta squared like in paper!
}


/*
 * Calculate the mass loss rate for massive stars (L > 4000 L_sol) using the
 * Nieuwenhuijzen & de Jager 1990 prescription, modified by a metallicity
 * dependent factor (Kudritzki et al 1989).
 *
 * Hurley et al. 2000, just after eq 106
 *
 *
 * double CalculateMassLossRateNieuwenhuijzenDeJagerStatic()
 *
 * @return                                      Nieuwenhuijzen & de Jager mass loss rate for massive stars (in Msol yr^-1)
 */
double BaseStar::CalculateMassLossRateNieuwenhuijzenDeJager() {
    double smoothTaper = min(1.0, (m_Luminosity - 4000.0) / 500.0); // Smooth taper between no mass loss and mass loss
    return sqrt((m_Metallicity / ZSOL)) * smoothTaper * 9.6E-15 * pow(m_Radius, 0.81) * pow(m_Luminosity, 1.24) * pow(m_Mass, 0.16);
}


/*
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
 *
 * Ref?
 *
 *
 * double CalculateMassLossRateLBV()
 *
 * @return                                      LBV-like mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateLBV() {
    return 0.1 * pow(((1.0E-5 * m_Radius * sqrt(m_Luminosity)) - 1.0), 3.0) * ((m_Luminosity / 6.0E5) - 1.0);
}


/*
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
 *
 * Belczynski et al. 2010, eq 8
 *
 *
 * double CalculateMassLossRateLBV2(const double p_Flbv)
 *
 * @param   [IN]    p_Flbv                      Multiplicitive constant multiplying base rate of 1E-4 Msol yr^-1 for LBV mass loss
 *                                              (Belczynski et al 2010 sets this as 1.5, Mennekens & Vanbeveren 2014 set this as 10)
 * @return                                      LBV-like mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateLBV2(const double p_Flbv) {
    return p_Flbv * 1.0E-4;
}


/*
 * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0).
 *
 * Hurley et al. 2000, just after eq 106 (taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998)
 *
 * Note that the reduction of this formula is imposed in order to match the observed number of black holes in binaries (Hurley et al 2000)
 *
 *
 * double CalculateMassLossRateWolfRayetLike(const double p_Mu)
 *
 * @param   [IN]    p_Mu                        Small envelope parameter (see Hurley et al. 2000, eq 97 & 98)
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayetLike(const double p_Mu) {
    // In the fortran code there is a parameter here hewind which by default is 1.0 -
    // can be set to zero to disable this particular part of winds. We instead opt for all winds on or off.
    return pow(m_Luminosity, 1.5) * (1.0 - p_Mu) * 1.0E-13;
}


/*
 * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0).
 *
 * Belczynski et al. 2010, eq 9 (taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998)
 *
 * Note that the reduction of this formula is imposed in order to match the observed number of black holes in binaries (Hurley et al 2000)
 *
 *
 * double CalculateMassLossRateWolfRayet2(const double p_Mu)
 *
 * @param   [IN]    p_Mu                        Small envelope parameter (see Hurley et al. 2000, eq 97 & 98)
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayet2(const double p_Mu) {
    // I think StarTrack may still do something different here,
    // there are references to Hamann & Koesterke 1998 and Vink and de Koter 2005

    return m_WolfRayetFactor * 1.0E-13 * pow(m_Luminosity, 1.5) * pow(m_Metallicity / ZSOL, 0.86) * (1.0 - p_Mu);
}


/*
 * Calculate the Wolf-Rayet like mass loss rate independent of WR star composition as given by Nugis & Lamers 2000
 *
 * Belczynski et al. 2010, eq 10.  We do not use this equation by default.
 *
 *
 * double CalculateMassLossRateWolfRayet3()
 *
 * @return                                      Mass loss rate (in Msol yr^{-1})
 */
double BaseStar::CalculateMassLossRateWolfRayet3() {
    return exp(-5.73 + (0.88 * log(m_Mass)));
}


/*
 * Calculate mass loss rate for massive OB stars using the Vink et al 2001 prescription
 *
 * Vink et al. 2001, eqs 24 & 25
 * Belczynski et al. 2010, eqs 6 & 7
 *
 *
 * double CalculateMassLossRateOB(const double p_Teff)
 *
 * @param   [IN]    p_Teff                      Effective temperature in K
 * @return                                      Mass loss rate for hot OB stars in Msol yr^-1
 */
double BaseStar::CalculateMassLossRateOB(const double p_Teff) {

    double rate;

    if (utils::Compare(p_Teff, 12500.0) >= 0 && utils::Compare(p_Teff, 25000.0) <= 0) {
        double V         = 1.3;                                                             // v_inf/v_esc

        double logMdotOB = -6.688                             +
                           (2.210 * log10(m_Luminosity / 1.0E5)) -
                           (1.339 * log10(m_Mass / 30.0))        -
                           (1.601 * log10(V / 2.0))              +
                           (0.85  * log10(m_Metallicity / ZSOL)) +
                           (1.07  * log10(p_Teff / 20000.0));

        rate = pow(10.0, logMdotOB);
    }
    else if (utils::Compare(p_Teff, 25000.0) > 0) {
        SHOW_WARN_IF(utils::Compare(p_Teff, 50000.0) > 0, ERROR::HIGH_TEFF_WINDS);          // show warning if winds being used outside comfort zone

        double V         = 2.6;                                                             // v_inf/v_esc

        double logMdotOB = -6.697 +
                           (2.194 * log10(m_Luminosity / 1.0E5)) -
                           (1.313 * log10(m_Mass / 30.0))        -
                           (1.226 * log10(V / 2.0))              +
                           (0.85  * log10(m_Metallicity / ZSOL)) +
                           (0.933 * log10(p_Teff / 40000.0))     -
                           (10.92 * log10(p_Teff / 40000.0) * log10(p_Teff/40000.0));

        rate = pow(10.0, logMdotOB);
    }
    else{
        SHOW_WARN(ERROR::LOW_TEFF_WINDS, "Mass Loss Rate = 0.0");                           // too cold to use winds - show warning
        rate = 0.0;
    }

    return rate;
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * According to Hurley et al. 2000
 *
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double BaseStar::CalculateMassLossRateHurley() {
    return (utils::Compare(m_Luminosity, 4.0E3) > 0) ? CalculateMassLossRateNieuwenhuijzenDeJager() : 0.0;          // JR: todo: make this a constant
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star at the current evolutionary phase
 * According to Vink - based on implementation in StarTrack
 *
 * double CalculateMassLossRateVink()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double BaseStar::CalculateMassLossRateVink() {
    double rate;

    double tmp = m_Radius * sqrt(m_Luminosity) * 1.0E-5;
    if ((utils::Compare(m_Luminosity, LBV_LUMINOSITY_LIMIT_STARTRACK) > 0) && (utils::Compare(tmp, 1.0) > 0)) {     // luminous blue variable
		m_LBVphaseFlag = true;                                                                                      // ... is true

        rate = CalculateMassLossRateLBV2(m_LBVfactor);                                                              // calculate mass loss rate
    }
    else {
        double teff = m_Temperature * TSOL;                                                                         // change to Kelvin so it can be compared with values as stated in Vink prescription

        if (utils::Compare(teff, 12500.0) < 0) {                                                                    // cool stars, use Hurley et al 2000 winds  JR: todo: make this a constant
            rate = CalculateMassLossRateHurley();
        }
        else  {                                                                                                     // hot stars, use Vink et al. 2001 winds (ignoring bistability jump)
            rate = CalculateMassLossRateOB(teff);
        }
    }

    // BSE and StarTrack have some mulptilier they apply here

    return rate;
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

    double mDot = 0.0;
    if (OPTIONS->UseMassLoss()) {

        switch (OPTIONS->MassLossPrescription()) {                                  // which prescription?

            case MASS_LOSS_PRESCRIPTION::HURLEY:                                    // HURLEY
                mDot = CalculateMassLossRateHurley();
                break;

            case MASS_LOSS_PRESCRIPTION::VINK:                                      // VINK
                mDot = CalculateMassLossRateVink();
                break;

            default:                                                                // unknown mass loss prescription
                SHOW_WARN(ERROR::UNKNOWN_MASS_LOSS_PRESCRIPTION, "Using HURLEY");   // show warning
                mDot = CalculateMassLossRateHurley();                               // use HURLEY
        }
    }

    return mDot;
}


/*
 * Calculate mass loss given mass loss rate - uses current timestep (m_Dt)
 * Returned mass loss is limited to 1% of current mass
 *
 *
 * double CalculateMassLoss_Static(const double p_Mass, const double p_Mdot, const double p_Dt)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Mdot                      Mass loss rate
 * @param   [IN]    p_Dt                        Timestep
 * @return                                      Mass loss
 */
double BaseStar::CalculateMassLoss_Static(const double p_Mass, const double p_Mdot, const double p_Dt) {
    return max(0.0, min(p_Mdot * p_Dt * 1.0E6, p_Mass * MAXIMUM_MASS_LOSS_FRACTION));   // Mass loss rate given in Msol per year, times are in Myr so need to multiply by 10^6
}


/*
 * Calculate values for dt, mDot and mass assuming mass loss is applied
 *
 * Class member variables m_Mdot and m_Dt are updated directly by this function if required (see paramaters)
 * Class member variables m_Mass is not updated directly by this function - the calculated mass is returned as the functional return
 *
 * - calculates (and limits) mass loss
 * - calculate new timestep (dt) and mass loss rate (mDot) to match (possibly limited) mass loss
 * - calculates new mass (mass) based on (possibly limited) mass loss
 *
 * Returns existing value for mass if mass loss not being used (program option)
 *
 *
 * double CalculateMassLossValues()                                                             // JR: todo: pick a better name for this...
 * @return                                      calculated mass (mSol)
 */
double BaseStar::CalculateMassLossValues(const bool p_UpdateMDot, const bool p_UpdateMDt) {

    double dt   = m_Dt;
    double mDot = m_Mdot;
    double mass = m_Mass;

    if (OPTIONS->UseMassLoss()) {                                           // only if using mass loss (program option)

        mDot = CalculateMassLossRate();                                     // calculate mass loss rate
        double massLoss = CalculateMassLoss_Static(mass, mDot, dt);         // calculate mass loss - limited to (mass * MAXIMUM_MASS_LOSS_FRACTION)

        // could do this without the test - we know the mass loss may already
        // have been limited.  This way is probably marginally faster
        if (utils::Compare(massLoss, (mass * MAXIMUM_MASS_LOSS_FRACTION)) < 0) {
            mass -= massLoss;                                               // new mass based on mass loss
        }
        else {
            dt    = massLoss / (mDot * 1.0E6);                              // new timestep to match limited mass loss
            mDot  = massLoss / (dt * 1.0E6);                                // new mass loss rate to match limited mass loss
            mass -= massLoss;                                               // new mass based on limited mass loss

            if (p_UpdateMDt) m_Dt = dt;                                     // update class member variable if necessary
        }

        if (p_UpdateMDot) m_Mdot = mDot;                                    // update class member variable if necessary
    }

    return mass;
}


/*
 * Resolve mass loss
 *
 * - calculates mass loss rate
 * - calculates (and limits) mass loss
 * - resets timestep (m_Dt) and mass loss rate (m_Mdot) to match (possibly limited) mass loss
 * - calculates and sets new mass (m_Mass) based on (possibly limited) mass loss
 * - applies mass rejuvenation factor and calculates new age
 *
 *
 * double ResolveMassLoss()
 */
void BaseStar::ResolveMassLoss() {

    if (OPTIONS->UseMassLoss()) {
        m_Mass = CalculateMassLossValues(true, true);                           // calculate new values assuming mass loss applied

        UpdateInitialMass();                                                    // update initial mass (MS, HG & HeMS)  JR: todo: fix this kludge one day - mass0 is overloaded, and isn't always "initial mass"
        UpdateAgeAfterMassLoss();                                               // update age (MS, HG & HeMS)
        ApplyMassTransferRejuvenationFactor();                                  // apply age rejuvenation factor
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
            ? pow((p_Luminosity / gbParams(B)), (1.0 / gbParams(q)))
            : pow((p_Luminosity / gbParams(D)), (1.0 / gbParams(p)));

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
 *    1) Kelvin-Helmholtz (thermal) timescale if THERMALLY_LIMITED_MASS_TRANSFER
 *    2) Choose a fraction of the mass rate that will be effectively accreted for FIXED_FRACTION_MASS_TRANSFER (as in StarTrack)
 *    3) Disk vs impact accretion for CENTRIFUGALLY_LIMITED_MASS_TRANSFER
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_FractionAccreted, const double p_AccretorMassRate)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_FractionAccreted          Accretion efficiency parameter - may be returned unchanged
 * @param   [IN]    p_AccretorMassRate          Thermal mass loss rate of the accretor (this star)
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
DBL_DBL BaseStar::CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_FractionAccreted, const double p_AccretorMassRate) {

    double acceptanceRate   = 0.0;                                                          // acceptance mass rate - default = 0.0
    double fractionAccreted = p_FractionAccreted;                                           // accretion fraction - default is unchanged

    switch (OPTIONS->MassTransferAccretionEfficiencyPrescription()) {

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED:                       // thermally limited mass transfer:

            if(utils::Compare(p_AccretorMassRate, p_DonorMassRate) < 0) {                   // the accretor cannot accrete at the rate the donor is providing mass

                acceptanceRate = min(OPTIONS->MassTransferCParameter() * p_AccretorMassRate, p_DonorMassRate);

                switch (OPTIONS->MassTransferThermallyLimitedVariation()) {

                    case MT_THERMALLY_LIMITED_VARIATION::C_FACTOR:
                        fractionAccreted = min(1.0, OPTIONS->MassTransferCParameter() * (p_AccretorMassRate / p_DonorMassRate));
                        break;

                    case MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE:
                        fractionAccreted = min(1.0, OPTIONS->MassTransferCParameter() * (p_AccretorMassRate / p_DonorMassRate));
                        break;

                    default:                                                                // unknown thermally limited variation - shouldn't happen
                        m_Error = ERROR::UNKNOWN_MT_THERMALLY_LIMITED_VARIATION;            // set error value
                        SHOW_WARN(m_Error);                                                 // warn that an error occurred
                }
            }
            else {                                                                          // accretor can accrete more than the donor is providing
                acceptanceRate   = p_DonorMassRate;
                fractionAccreted = 1.0;
            }
            break;

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::FIXED_FRACTION:                          // fixed fraction of mass accreted, as in StarTrack
            acceptanceRate = min(p_DonorMassRate, p_FractionAccreted * p_DonorMassRate);
            break;

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::CENTRIFUGALLY_LIMITED:                   // centrifugally limited mass transfer
            acceptanceRate = 0.0;
            break;

        default:                                                                            // unknown mass transfer accretion efficiency prescription - shouldn't happen
            m_Error = ERROR::UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION;                  // set error value
            SHOW_WARN(m_Error);                                                             // warn that an error occurred
    }

    return std::make_tuple(acceptanceRate, fractionAccreted);
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
    return sqrt(sqrt(p_Luminosity)) / sqrt(p_Radius);   // sqrt() is much faster than pow()
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
double BaseStar::CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) {
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
 * @param   [IN]    p_vE                        Rotational velocity (in km s^-1) at which to calculate CDF
 * @return                                      The CDF for the rotational velocity
 */
double BaseStar::CalculateOStarRotationalVelocityAnalyticCDF_Static(const double p_Ve) {

    double alpha  = 4.82;
    double beta   = 1.0 / 25.0;
    double mu     = 205.0;
    double sigma  = 190.0;
    double iGamma = 0.43;

    boost::math::inverse_gamma_distribution<> gammaComponent(alpha, beta); // (shape, scale) = (alpha, beta)
    boost::math::normal_distribution<> normalComponent(mu, sigma);

	return (iGamma * boost::math::cdf(gammaComponent, p_Ve)) + ((1.0 - iGamma) * boost::math::cdf(normalComponent, p_Ve));
}


/*
 * Calculate the inverse of the analytic cumulative distribution function (CDF) for the
 * equatorial rotational velocity of single O stars.
 *
 * (i.e. calculate the inverse of CalculateOStarRotationalVelocityAnalyticCDF_Static())
 *
 *
 * double CalculateOStarRotationalVelocityAnalyticCDFInverse_Static(const double p_Ve, const void *p_Params)
 * @param   [IN]    p_vE                        Rotational velocity (in km s^-1) - value of the kick vk which we want to find
 * @param   [IN]    p_Params                    Pointer to RotationalVelocityParams structure containing y, the CDF draw U(0,1)
 * @return                                      Inverse CDF
 *                                              Should be zero when p_Ve = vk, the value of the kick to draw
 */
double BaseStar::CalculateOStarRotationalVelocityAnalyticCDFInverse_Static(double p_Ve, void* p_Params) {
    RotationalVelocityParams* params = (RotationalVelocityParams*) p_Params;
    return CalculateOStarRotationalVelocityAnalyticCDF_Static(p_Ve) - params->u;
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
 * double CalculateOStarRotationalVelocity_Static(double p_U, double p_Xmin, double p_Xmax)
 *
 * @param   [IN]    p_Xmin                      Minimum value for root
 * @param   [IN]    p_Xmax                      Maximum value for root
 * @return                                      Rotational velocity in km s^-1
 */
double BaseStar::CalculateOStarRotationalVelocity_Static(const double p_Xmin, const double p_Xmax) {

    double xMin = p_Xmin;
    double xMax = p_Xmax;

    double result = xMin;

    double maximumInverse = CalculateOStarRotationalVelocityAnalyticCDF_Static(xMax);
    double minimumInverse = CalculateOStarRotationalVelocityAnalyticCDF_Static(xMin);

    double rand = RAND->Random();

    while(utils::Compare(rand, maximumInverse) > 0) {
        xMax          *= 2.0;
        maximumInverse = CalculateOStarRotationalVelocityAnalyticCDF_Static(xMax);
    }

    if(utils::Compare(rand, minimumInverse) >= 0) {

        const gsl_root_fsolver_type *T;
        gsl_root_fsolver            *s;
        gsl_function                 F;

    	RotationalVelocityParams     params = {rand};

	    F.function = &CalculateOStarRotationalVelocityAnalyticCDFInverse_Static;
	    F.params   = &params;

	    // gsl_root_fsolver_brent
	    // gsl_root_fsolver_bisection
	    T = gsl_root_fsolver_brent;
	    s = gsl_root_fsolver_alloc(T);

	    gsl_root_fsolver_set(s, &F, xMin, xMax);

	    int status  = GSL_CONTINUE;
        int iter    = 0;
        int maxIter = 100;

    	while (status == GSL_CONTINUE && iter < maxIter) {
        	iter++;
        	status = gsl_root_fsolver_iterate(s);
        	result = gsl_root_fsolver_root(s);
        	xMin   = gsl_root_fsolver_x_lower(s);
        	xMax   = gsl_root_fsolver_x_upper(s);
        	status = gsl_root_test_interval(xMin, xMax, 0, 0.001);
        }

    	gsl_root_fsolver_free(s);   // de-allocate memory for root solver
    }

    return result;
}


/*
 * Calculate the inital rotational velocity (in km s^-1 ) of a star with ZAMS mass MZAMS
 *
 * Distribution used is determined by program option "rotationalVelocityDistribution"
 *
 *
 * double CalculateRotationalVelocity(double p_MZAMS)
 *
 * @param   [IN]    p_MZAMS                     Zero age main sequence mass in Msol
 * @return                                      Initial equatorial rotational velocity in km s^-1 - vRot in Hurley at al. 2000
 */
double BaseStar::CalculateRotationalVelocity(double p_MZAMS) {

    double vRot = 0.0;

    switch (OPTIONS->RotationalVelocityDistribution()) {                                // which prescription?

        case ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO: break;                             // ZERO

        case ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY:                                  // HURLEY

            // Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
            vRot = (330.0 * pow(p_MZAMS, 3.3)) / (15.0 + pow(p_MZAMS, 3.45));
            break;

        case ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES:                               // VLTFLAMES

            // Rotational velocity based on VLT-FLAMES survey.
            // For O-stars use results of Ramirez-Agudelo et al. (2013) https://arxiv.org/abs/1309.2929 (single stars)
            // and Ramirez-Agudelo et al. (2015) https://arxiv.org/abs/1507.02286 (spectroscopic binaries)
            // For B-stars use results of Dufton et al. (2013) https://arxiv.org/abs/1212.2424
            // For lower mass stars, I don't know what updated results there are so default back to
            // Hurley et al. 2000 distribution for now

            if (utils::Compare(p_MZAMS, 16.0) >= 0) {
                vRot = CalculateOStarRotationalVelocity_Static(0.0, 800.0);
            }
            else if (utils::Compare(p_MZAMS, 2.0) >= 0) {
                vRot = utils::InverseSampleFromTabulatedCDF(RAND->Random(), BStarRotationalVelocityCDFTable);
            }
            else {
                // Don't know what better to use for low mass stars so for now
                // default to Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
                vRot = (330.0 * pow(p_MZAMS, 3.3)) / (15.0 + pow(p_MZAMS, 3.45));
            }
            break;

        default:                                                                        // unknown rorational velocity prescription
            SHOW_WARN(ERROR::UNKNOWN_VROT_PRESCRIPTION, "Using default vRot = 0.0");     // show warning
    }
    return vRot;
}


/*
 * Calculate the initial angular frequency (in yr^-1) of a star with
 * ZAMS mass and radius MZAMS and RZAMS respectively
 *
 * Hurley at al. 2000, eq 108
 *
 *
 * double CalculateRotationalAngularFrequency(const double p_MZAMS, const double p_RZAMS)
 *
 * @param   [IN]    p_MZAMS                     Zero age main sequence mass in Msol
 * @param   [IN]    p_RZAMS                     Zero age main sequence radius in Rsol
 * @return                                      Initial angular frequency in yr^-1 - omega in Hurley et al. 2000
 */
double BaseStar::CalculateZAMSAngularFrequency(const double p_MZAMS, const double p_RZAMS) {
    double vRot = CalculateRotationalVelocity(p_MZAMS);
    return utils::Compare(vRot, 0.0) == 0 ? 0.0 : 45.35 * vRot / p_RZAMS;    // Hurley at al. 2000, eq 108       JR: todo: added check for vRot = 0
}


/*
 * Calculate the break up frequency of star in yr^-1 units, where [G] = 4*pi^2 AU^3 yr^-2 Msol^-1
 *
 *
 * double CalculateOmegaBreak()
 *
 * @return                                      Break up frequency (yr^-1)
 */
double BaseStar::CalculateOmegaBreak() {
    constexpr double _2_PI_2      = _2_PI * _2_PI;
    constexpr double RSOL_TO_AU_3 = RSOL_TO_AU * RSOL_TO_AU * RSOL_TO_AU;

	return sqrt(_2_PI_2 * m_Mass / (RSOL_TO_AU_3 * m_Radius * m_Radius * m_Radius));
}


/*
 * Calculate the minimum rotational frequency (in yr^-1) at which CHE will occur
 * for a star with ZAMS mass MZAMS
 *
 * Mendel's fit from Butler 2018
 *
 *
 * double CalculateOmegaCHE(const double p_MZAMS, const double p_Metallicity)
 *
 * @param   [IN]        p_MZAMS                 Zero age main sequence mass in Msol
 * @param   [IN]        p_Metallicity           Metallicity of the star
 * @return                                      Initial angular frequency in rad*s^-1
 */
double BaseStar::CalculateOmegaCHE(const double p_MZAMS, const double p_Metallicity) {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double mRatio = p_MZAMS;                                                                        // in MSol, so ratio is just p_MZAMS

    // calculate omegaCHE(M, Z = 0.004)
    double omegaZ004 = 0.0;
    if (utils::Compare(p_MZAMS, massCutoffs(MCHE)) <= 0) {
        for (std::size_t i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(mRatio, i) / pow(mRatio, 0.4);
        }
    }
    else {
        for (std::size_t i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(100.0, i) / pow(mRatio, 0.4);
        }
    }

    // calculate omegaCHE(M, Z)
    return (1.0 / ((0.09 * log(p_Metallicity / 0.004)) + 1.0) * omegaZ004) * SECONDS_IN_YEAR;   // in rads/yr

#undef massCutoffs
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate lifetime to the Base of the Giant Branch (end of the Hertzsprung Gap)
 * For high mass stars, t_BGB = t_HeI. (JR: there is no check here...)
 *
 * Hurley et al. 2000, eq 4 (plotted in Hurley et al. 2000, fig 5)
 *
 *
 * double CalculateLifetimeToBGB(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Lifetime to the Base of the Giant Branch in Myr
 */
double BaseStar::CalculateLifetimeToBGB(const double p_Mass) {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    double m_2   = p_Mass * p_Mass;
    double m_4   = m_2 * m_2;
    double m_5_5 = m_4 * p_Mass * sqrt(p_Mass);
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
double BaseStar::CalculateLifetimeToBAGB(const double p_tHeI, const double p_tHe) {
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
    return 5.0 * 1.0E-5 * p_Radius * sqrt(p_Radius) * YEAR_TO_MYR / sqrt(p_Mass);   // sqrt() is much faster than pow()
}


/*
 * Calculate nuclear timescale
 *
 * Kalogera & Webbink 1996, eq 3
 *
 *
 * double CalculateNuclearTimescale_Static(const double p_Mass, const double p_Luminosity)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @return                                      Dynamical timescale in Myr
 */
double BaseStar::CalculateNuclearTimescale_Static(const double p_Mass, const double p_Luminosity) {
    return 1.0E10 * p_Mass * YEAR_TO_MYR / p_Luminosity;
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
 * @return                                      Radial expansion timescale in yr
 *                                              Returns -1.0 if radial expansion timescale can't be calculated
 *                                              (i.e. stellar type has changed or radius has not changed)
 */
double BaseStar::CalculateRadialExpansionTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                          const STELLAR_TYPE p_StellarTypePrev,
                                                          const double       p_Radius,
                                                          const double       p_RadiusPrev,
                                                          const double       p_DtPrev) {

	return p_StellarTypePrev == p_StellarType && utils::Compare(p_RadiusPrev, p_Radius) != 0
            ? (p_DtPrev * MYR_TO_YEAR * p_RadiusPrev) / (p_Radius - p_RadiusPrev)
            : -1.0;
}


/*
 * Calculate all timescales
 *
 * Calculates:
 *     - dynamical timescale
 *     - thermal timescale
 *     - nuclear timesclae
 *     - radial expansion timescale
 *
 * and sets class member variables appropriately
 *
 *
 * void CalculateAllTimescales()
 */
void BaseStar::CalculateAllTimescales() {

    m_DynamicalTimescale	   = CalculateDynamicalTimescale();
    m_ThermalTimescale 		   = CalculateThermalTimescale();
    m_NuclearTimescale 		   = CalculateNuclearTimescale();
    m_RadialExpansionTimescale = CalculateRadialExpansionTimescale();
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
double BaseStar::CalculateEddyTurnoverTimescale() {

	double rEnv	= CalculateRadialExtentConvectiveEnvelope();

	return 0.4311 * pow((m_Mass * rEnv * (m_Radius - (0.5 * rEnv))) / (3.0 * m_Luminosity), 1.0 / 3.0);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                SUPERNOVA FUNCTIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the kick given to a black hole based on users chosen assumptions about black hole kicks,
 * fallback and the magnitude of a kick drawn from a distribution
 *
 * Current options are:
 *
 *    FULL    : Black holes receive the same kicks as neutron stars
 *    REDUCED : Black holes receive the same momentum kick as a neutron star, but downweighted by the black hole mass
 *    ZERO    : Black holes receive zero natal kick
 *    FALLBACK: Black holes receive a kick downweighted by the amount of mass falling back onto them
 *
 *
 *  double ApplyBlackHoleKicks(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass)
 *
 * @param   [IN]    p_vK                        Kick velocity that would otherwise be applied to a neutron star
 * @param   [IN]    p_FallbackFraction          Fraction of mass that falls back onto the proto-compact object
 * @param   [IN]    p_BlackHoleMass             Mass of remnant (in Msol)
 * @return                                      Kick velocity
 */
 double BaseStar::ApplyBlackHoleKicks(const double p_vK, const double p_FallbackFraction, const double p_BlackHoleMass) {

    double vK;

    switch (OPTIONS->BlackHoleKicksOption()) {                      // which BH kicks option specified?

        case BLACK_HOLE_KICK_OPTION::FULL:                          // BH receives full kick - no adjustment necessary
            vK = p_vK;
            break;

        case BLACK_HOLE_KICK_OPTION::REDUCED:                       // Kick is reduced by the ratio of the black hole mass to neutron star mass i.e. v_bh = ns/bh  *v_ns
            vK = p_vK * NEUTRON_STAR_MASS / p_BlackHoleMass;
            break;

        case BLACK_HOLE_KICK_OPTION::ZERO:
            vK = 0.0;                                               // BH Kicks are set to zero regardless of BH mass or kick velocity drawn.
            break;

        case BLACK_HOLE_KICK_OPTION::FALLBACK:                      // Using the so-called 'fallback' prescription for BH kicks
            vK = p_vK * (1.0 - p_FallbackFraction);
            break;

        default:                                                    // unknown BH kick option - shouldn't happen
            vK = p_vK;                                              // return vK unchanged
            m_Error = ERROR::UNKNOWN_BH_KICK_OPTION;                // set error value
            SHOW_WARN(m_Error);                                     // warn that an error occurred
    }

    return vK;
}


/*
 * Inverse sampling from the Maxwell CDF for MCMC and importance sampling
 * Generates a random sample from the distribution
 *
 * From https://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution,
 *
 * the Maxwell CDF is:
 *
 * erf(x/(sqrt(2.0)*sigma)) - (sqrt(2.0/pi) * x * exp(-(x*x)/(2.0*sigma*sigma))/sigma)
 *
 * where erf is the error function which has the form:
 *
 *   erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
 *
 * We use the gsl error function https://www.gnu.org/software/gsl/manual/html_node/Error-Function.html
 *
 *
 * double InverseSampleFromMaxwellCDF_Static(const double p_X, const double p_Sigma)
 *
 * @param   [IN]    p_X                         The value of X for which the distribution function is evaluated
 * @param   [IN]    p_Sigma                     Distribution scale parameter - affects the spread of the distribution
 * @return                                      Maxwell CDF value
 */
double BaseStar::InverseSampleFromMaxwellCDF_Static(const double p_X, const double p_Sigma) {
    return gsl_sf_erf(p_X / (M_SQRT2 * p_Sigma)) - (SQRT_M_2_PI * p_X * exp(-(p_X * p_X) / (2.0 * p_Sigma * p_Sigma)) / p_Sigma);
}


/*
 * Calculate the inverse of the Maxwell CDF
 *
 *
 * double CalculateInverseMaxwellCDF_Static(const double p_X, void* p_Paramsa)
 *
 * @param   [IN]    p_X                         Value of the kick vk which we want to find
 * @param   [IN]    p_Params                    Structure containing:
 *                                                 y    : the CDF draw U(0,1), and
 *                                                 sigma: which sets the characteristic kick of the Maxwellian distribution
 * @return                                      Inverse Maxwell CDF value (Should be zero when x = vk, the value of the kick to draw)
 */
double BaseStar::CalculateInverseMaxwellCDF_Static(const double p_X, void* p_Params) {
    KickVelocityParams* params = (KickVelocityParams*)p_Params;
    return InverseSampleFromMaxwellCDF_Static(p_X, params->sigma) -params->y;
}


/*
 * Draw a kick velocity in km s^-1 from a Maxwellian distribution of the form:
 *
 *    sqrt(2/pi) * (x * x * exp(-(x * x)/(2.0 * a * a)) / (a*a*a))
 *
 * Roughly factor of 10 slower than old method, 10^-4 s vs 10^-5 s for old method.
 *
 * Implement option to only do it this way if using mcmc or importance sampling,
 * otherwise use old method     JR: todo: what is the old method?
 *
 *   Inverse sampling done using root finding in GSL, adapted from the examples here:
 *   https://www.gnu.org/software/gsl/doc/html/roots.html
 *
 *
 * double DrawKickVelocityDistributionMaxwell(const double p_Sigma, const double p_Rand)
 *
 * @param   [IN]    p_Sigma                     Distribution scale parameter - affects the spread of the distribution
 * @param   [IN]    p_Rand                      Random number between 0 and 1 used for drawing from the inverse CDF of the Maxwellian
 * @return                                      Drawn kick velocity (km s^-1)
 */
double BaseStar::DrawKickVelocityDistributionMaxwell(const double p_Sigma, const double p_Rand) {

    double result = 0.0;

    bool useOldMethod = false;
    if (useOldMethod) {                                     // Old method - JR: todo: is this !OPTIONS->useMCMC ?

        double rms = 0.0;
        for (int i=0; i<3; i++) {                           // generate each component of a 3 dimensional velocity vector
            double vel = RAND->RandomGaussian(p_Sigma);
            rms       += vel * vel;                         // sum the squares (because we need the magnitude
        }

        result = sqrt(rms);                                 // magnitude
    }
    else {

    	double xMin = 0.0;
    	double xMax = 5.0 * p_Sigma;                        // options.kickmax or whatever it is called.        JR: todo: find out...

        double maximumInverse = InverseSampleFromMaxwellCDF_Static(xMin, p_Sigma);

        double rand = p_Rand;
        while (utils::Compare(rand, maximumInverse) > 0) {
            xMax *= 2.0;
            maximumInverse = InverseSampleFromMaxwellCDF_Static(xMax, p_Sigma);
        }
        rand = min(rand, maximumInverse);

        const gsl_root_fsolver_type *T;
        gsl_root_fsolver            *s;
        gsl_function                 F;

        KickVelocityParams           params = {rand, p_Sigma};

        F.function = &CalculateInverseMaxwellCDF_Static;
        F.params   = &params;

        // gsl_root_fsolver_brent
        // gsl_root_fsolver_bisection
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);

        gsl_root_fsolver_set (s, &F, xMin, xMax);

	    int status  = GSL_CONTINUE;
        int iter    = 0;
        int maxIter = 100;

    	while (status == GSL_CONTINUE && iter < maxIter) {
            iter++;
            status = gsl_root_fsolver_iterate (s);
            result = gsl_root_fsolver_root (s);
            xMin   = gsl_root_fsolver_x_lower (s);
            xMax   = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval (xMin, xMax, 0, 0.001);
        }

        gsl_root_fsolver_free (s);                          // de-allocate memory for root solver
    }

    return result;
}


/*
 * Draw a kick velocity in km s^-1 from a uniform distribution between 0 and parameter p_MaxVK
 *
 *
 * double DrawKickVelocityDistributionFlat(const double p_MaxVK, const double p_Rand)
 *
 * @param   [IN]    p_MaxVK                     Maximum kick velocity in km s^-1 to draw
 * @param   [IN]    p_Rand                      Random number between 0 and 1 used for drawing from the distribution
 * @return                                      Drawn kick velocity (km s^-1)
 */
double BaseStar::DrawKickVelocityDistributionFlat(const double p_MaxVK, const double p_Rand) {

    bool useOldMethod = false;          //  JR: todo: program option?

    double rand = useOldMethod ? RAND->Random() : p_Rand;

    return rand * p_MaxVK;
}


/*
 * Draw a kick velocity in km s^-1 per Bray & Eldridge 2016, 2018
 *
 * See:
 *    https://arxiv.org/abs/1605.09529
 *    https://arxiv.org/abs/1804.04414
 *
 *
 * double DrawKickVelocityBrayEldridge(const double p_EjectaMass,
 *                                     const double p_RemnantMass,
 *                                     const double p_Alpha,
 *                                     const double p_Beta)
 *
 * @param   [IN]    p_EjectaMass                Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @param   [IN]    p_RemnantMass               Mass of the remnant (Msol)
 * @param   [IN]    p_Alpha                     Fitting coefficient (see Bray & Eldridge 2016, 2018)
 * @param   [IN]    p_Beta                      Fitting coefficient (see Bray & Eldridge 2016, 2018)
 * @return                                      Drawn kick velocity (km s^-1)
 */
double BaseStar::DrawKickVelocityBrayEldridge(const double p_EjectaMass,
                                              const double p_RemnantMass,
                                              const double p_Alpha,
                                              const double p_Beta) {

    return p_Alpha * (p_EjectaMass / p_RemnantMass) + p_Beta;
}


/*
 * Draw kick velocity per Muller et al. 2016
 *
 * For BH assume complete fallback so no ejecta hence no rocket effect
 * If BH doesnt assumes complete fallback, e.g. neutrino mass loss, a Blauuw Kick should be calculated
 *
 * ALEJANDRO - 17/03/2017 - It seems from Bernhard's notes he gives us V_kick = V_kick_3D.
 * Checking with some of the highest kicks he gives, we get V_kick>1000 km s^-1
 *
 *
 * double DrawRemnantKickMuller(const double p_COCoreMass)
 *
 * @param   [IN]    p_COCoreMass                Carbon Oxygen core mass of exploding star (Msol)
 * @return                                      Drawn kick velocity (km s^-1)
 */
double BaseStar::DrawRemnantKickMuller(const double p_COCoreMass) {

    double	remnantKick = 0.0;	                // units km/s
	double	lowerRegimeKick = 70.0;		        // Bernhard proposes to use 10 km s^-1 to replicate ECSN. He quotes Bray & Eldridge 2016 on this.

	     if (utils::Compare(p_COCoreMass, 1.44) <  0) remnantKick = 0.0;
	else if (utils::Compare(p_COCoreMass, 1.49) <  0) remnantKick = lowerRegimeKick + (2000.0 * (p_COCoreMass - 1.372));
    else if (utils::Compare(p_COCoreMass, 1.65) <  0) remnantKick = 180.0 + (1300.0 * (p_COCoreMass - 1.49));
	else if (utils::Compare(p_COCoreMass, 2.4 ) <  0) remnantKick = 250.0 + (350.0  * (p_COCoreMass - 1.65));
    else if (utils::Compare(p_COCoreMass, 3.2 ) <  0) remnantKick = 400.0 + (1100.0 * (p_COCoreMass - 2.4));
    else if (utils::Compare(p_COCoreMass, 3.6 ) <  0) remnantKick = 160.0 + (240.0  * (p_COCoreMass - 3.2));
    else if (utils::Compare(p_COCoreMass, 4.05) <  0) remnantKick = 0.0;
    else if (utils::Compare(p_COCoreMass, 4.6 ) <  0) remnantKick = 700.0 + (100.0  * (p_COCoreMass - 4.05));
    else if (utils::Compare(p_COCoreMass, 5.7 ) <  0) remnantKick = 0.0;
    else if (utils::Compare(p_COCoreMass, 6.0 ) <  0) remnantKick = 550.0 - (600.0  * (p_COCoreMass - 5.7));
    else if (utils::Compare(p_COCoreMass, 6.0 ) >= 0) remnantKick = 0.0;

    return remnantKick;
}


/*
 * Draw a kick velocity from the user-specified distribution
 *
 *
 * double DrawSNKickVelocity(const double p_Sigma,
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
 * @return                                      Drawn kick velocity (km s^-1)
 */
double BaseStar::DrawSNKickVelocity(const double p_Sigma,
                                    const double p_COCoreMass,
                                    const double p_Rand,
                                    const double p_EjectaMass,
                                    const double p_RemnantMass) {
	double kickVelocity;

    switch (OPTIONS->KickVelocityDistribution()) {                                              // which distribution

        case KICK_VELOCITY_DISTRIBUTION::MAXWELLIAN:
        case KICK_VELOCITY_DISTRIBUTION::MAXWELL:
            kickVelocity = DrawKickVelocityDistributionMaxwell(p_Sigma, p_Rand);                // MAXWELLIAN, MAXWELL
            break;

        case KICK_VELOCITY_DISTRIBUTION::FLAT:                                                  // FLAT
            kickVelocity = DrawKickVelocityDistributionFlat(OPTIONS->KickVelocityDistributionMaximum(), p_Rand);
            break;

        case KICK_VELOCITY_DISTRIBUTION::ZERO:                                                  // ZERO
            kickVelocity = 0.0;
            break;

        case KICK_VELOCITY_DISTRIBUTION::FIXED:                                                 // FIXED
            kickVelocity = p_Sigma;
            break;

        case KICK_VELOCITY_DISTRIBUTION::BRAYELDRIDGE:                                          // BRAY ELDRIDGE
            kickVelocity = DrawKickVelocityBrayEldridge(p_EjectaMass, p_RemnantMass, BRAY_ELDRIDGE_CONSTANT_VALUES.at(BRAY_ELDRIDGE_CONSTANT::ALPHA), BRAY_ELDRIDGE_CONSTANT_VALUES.at(BRAY_ELDRIDGE_CONSTANT::BETA));
            break;

        case KICK_VELOCITY_DISTRIBUTION::MULLER2016:                                            // MULLER2016
            kickVelocity = DrawRemnantKickMuller(p_COCoreMass);
            break;

        case KICK_VELOCITY_DISTRIBUTION::MULLER2016MAXWELLIAN: {                                // MULLER2016-MAXWELLIAN

            double mullerSigma = DrawRemnantKickMuller(p_COCoreMass) / sqrt(3.0);

            kickVelocity = DrawKickVelocityDistributionMaxwell(mullerSigma, p_Rand);
            } break;

        default:                                                                                // unknown distribution
            SHOW_WARN(ERROR::UNKNOWN_KICK_VELOCITY_DISTRIBUTION, "Using default: MAXWELL");     // show warning
            kickVelocity = DrawKickVelocityDistributionMaxwell(p_Sigma, p_Rand);
    }

    return kickVelocity / OPTIONS->KickScalingFactor();

}


/*
 * Calculate supernova kick velocity
 * Based on the current supernova event type and user-specified kick velocity distributions
 *
 *
 * double BaseStar::CalculateSNKickVelocity(const double p_RemnantMass, const double p_EjectaMass)
 *
 * @param   [IN]    p_RemnantMass               The mass of the remnant (Msol)
 * @param   [IN]    p_EjectaMass                Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @return                                      Kick velocity
 */
double BaseStar::CalculateSNKickVelocity(const double p_RemnantMass, const double p_EjectaMass) {
    ERROR error = ERROR::NONE;
	double vK;

    if (!m_SupernovaDetails.initialKickParameters.supplied ||                                       // user did not supply kick parameters, or
        (m_SupernovaDetails.initialKickParameters.supplied &&                                       // user did supply kick parameters but ...
         m_SupernovaDetails.initialKickParameters.useVelocityRandom)) {                             // ... wants to draw velocity using supplied random number

        double sigma;
        switch (utils::SNEventType(m_SupernovaDetails.events.current)) {                            // what type of supernova event happening now?

		    case SN_EVENT::ECSN:                                                                    // ALEJANDRO - 04/05/2017 - Allow for ECSN to have kicks different than zero. Still, should be low kicks. Default set to zero.  (JR: todo: check default = 30.0?)
			    sigma = OPTIONS->KickVelocityDistributionSigmaForECSN();
                break;

		    case SN_EVENT::USSN:                                                                    // ALEJANDRO - 25/08/2017 - Allow for USSN to have a separate kick.
			    sigma = OPTIONS->KickVelocityDistributionSigmaForUSSN();
                break;

		    case SN_EVENT::CCSN:                                                                    // draw a random kick velocity from the user selected distribution - sigma based on whether compact object is a NS or BH

                switch (m_StellarType) {                                                            // which stellar type?
                    case STELLAR_TYPE::NEUTRON_STAR:
                        sigma = OPTIONS->KickVelocityDistributionSigmaCCSN_NS();
                        break;

                    case STELLAR_TYPE::BLACK_HOLE:
                        sigma = OPTIONS->KickVelocityDistributionSigmaCCSN_BH();
                        break;

                    default:                                                                        // unknown stellar type - shouldn't happen
                        error = ERROR::UNKNOWN_STELLAR_TYPE;
                }

                break;

            case SN_EVENT::PISN:                                                                    // not expected here
            case SN_EVENT::PPISN:                                                                   // not expected here
                error = ERROR::UNEXPECTED_SN_EVENT;
                break;

            case SN_EVENT::NONE:                                                                    // no supernova event - shouldn't be here...
                error = ERROR::EXPECTED_SN_EVENT;
                break;

		    default:                                                                                // unknown supernova event - shouldn't happen
                error = ERROR::UNKNOWN_SN_EVENT;
	    }
    
	    if (error == ERROR::NONE) {                                                                 // check for errors
                                                                                                    // no errors - draw kick velocity
            vK = DrawSNKickVelocity(sigma, 
                                    m_SupernovaDetails.COCoreMassAtCOFormation, 
                                    m_SupernovaDetails.kickVelocityRandom,
                                    p_EjectaMass, 
                                    p_RemnantMass);
        }
    }
    else {                                                                                          // user supplied kick parameters and wants tu use supplied kick velocity, so ...
        vK = m_SupernovaDetails.initialKickParameters.velocity;                                     // ... use it 
    }


	if (error == ERROR::NONE) {                                                                     // check for errors

        m_SupernovaDetails.drawnKickVelocity = vK;                                                  // drawn kick velocity

        if (utils::SNEventType(m_SupernovaDetails.events.current) == SN_EVENT::CCSN) {              // core-collapse supernova event this timestep?
            vK = ApplyBlackHoleKicks(vK, m_SupernovaDetails.fallbackFraction, m_Mass);              // re-weight kicks by mass of remnant according to user specified black hole kicks option
        }
        else {                                                                                      // otherwise
            m_SupernovaDetails.fallbackFraction = 0.0;                                              // set fallback fraction to zero
        }
        m_SupernovaDetails.kickVelocity = vK;                                                       // updated kick velocity
    }
    else {                                                                                          // error occurred
        vK = 0.0;                                                                                   // set kick velocity to zero
        m_Error = error;                                                                            // set error value
        SHOW_WARN(m_Error);                                                                         // warn that an error occurred
    }

    return vK;
}


 /*
  * Draw the angular components of the supernova kick theta and phi according to user specified options.

    Parameters
    -----------
    options : programOptions
        User specified options
    r : gsl_rng
        Used for the random number generation
    theta : double
        Polar angle for kick (pointer)
    phi : double
        Azimuthal angle for kick (pointer)
    kickDirectionPower : double
        Power law for the POWER angle distribution

    Returns
    --------
    updates values of theta and phi
    */
DBL_DBL BaseStar::DrawKickDirection() {

    double theta;                                                                                               // theta, angle out of the plane
    double phi;                                                                                                 // phi, angle in the plane

    double delta = 1.0 * DEGREE;                                                                                // small angle () in radians - could be set by user in options

    double rand = RAND->Random();                                                                               // do this here to be consistent with legacy code - allows comparison tests (won't work for long - soon there will be too many changes to the code...)

    switch (OPTIONS->KickDirectionDistribution()) {                                                             // which kick direction distribution?

        case KICK_DIRECTION_DISTRIBUTION::ISOTROPIC:                                                            // ISOTROPIC: Draw theta and phi isotropically
            theta = acos(1.0 - (2.0 * RAND->Random())) - (4*std::atan(1.0) / 2.0);//M_PI_2;
            phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
            break;

        case KICK_DIRECTION_DISTRIBUTION::POWERLAW: {                                                           // POWERLAW: Draw phi uniform in [0,2pi], theta according to a powerlaw
                                                                                                                // (power law power = 0 = isotropic, +infinity = kick along pole, -infinity = kick in plane)
            // Choose magnitude of power law distribution -- if using a negative power law that blows up at 0,
            // need a lower cutoff (currently set at 1E-6), check it doesn't affect things too much
            // JR: todo: shsould these be in constants.h?
            double magnitude_of_cos_theta = utils::InverseSampleFromPowerLaw(OPTIONS->KickDirectionPower(), 1.0, 1E-6);

            if (OPTIONS->KickDirectionPower() < 0.0) magnitude_of_cos_theta = 1.0 - magnitude_of_cos_theta;     // don't use utils::Compare() here

            double actual_cos_theta = magnitude_of_cos_theta;

            if (RAND->Random() < 0.5) actual_cos_theta = -magnitude_of_cos_theta;                               // By taking the magnitude of cos theta we lost the information about whether it was up or down, so we put that back in randomly here

            actual_cos_theta = min(1.0, max(-1.0, actual_cos_theta));                                           // clamp to [-1.0, 1.0]

            theta = acos(actual_cos_theta);                                                                     // set the kick angle out of the plane theta
            phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
            } break;

        case KICK_DIRECTION_DISTRIBUTION::INPLANE:                                                              // INPLANE: Force the kick to be in the plane theta = 0
            theta = 0.0;                                                                                        // force the kick to be in the plane - theta = 0
            phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
            break;

        case KICK_DIRECTION_DISTRIBUTION::PERPENDICULAR:                                                        // PERPENDICULAR: Force kick to be along spin axis
            theta = M_PI_2;                                                                                     // pi/2 UP
            if (rand >= 0.5) theta = -theta;                                                                    // switch to DOWN
            phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
            break;

        case KICK_DIRECTION_DISTRIBUTION::POLES:                                                                // POLES: Direct the kick in a small cone around the poles
            if (rand < 0.5) theta = M_PI_2 - fabs(RAND->RandomGaussian(delta));                                 // UP - slightly less than or equal to pi/2
            else            theta = fabs(RAND->RandomGaussian(delta)) - M_PI_2;                                 // DOWN - slightly more than or equal to -pi/2

            phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
            break;

        case KICK_DIRECTION_DISTRIBUTION::WEDGE:                                                                // WEDGE: Direct kick into a wedge around the horizon (theta = 0)
            theta = RAND->RandomGaussian(delta);                                                                // Gaussian around 0 with a deviation delta
            phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
            break;

        default:                                                                                                // unknown kick direction distribution - use ISOTROPIC
            theta = acos(1.0 - (2.0 * RAND->Random())) - M_PI_2;
            phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
            SHOW_WARN(ERROR::UNKNOWN_KICK_DIRECTION_DISTRIBUTION, " Using default");                            // warn that an error occurred
    }

    return std::make_tuple(theta, phi);
}


/*
 * Solve Kepler's Equation using root finding techniques. Here we use Newton-Raphson.
 *
 * For a definition of all the anomalies see here:
 *
 *    https://en.wikipedia.org/wiki/Mean_anomaly
 *    https://en.wikipedia.org/wiki/True_anomaly
 *    https://en.wikipedia.org/wiki/Eccentric_anomaly
 *
 *
 * DBL_DBL SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity)
 *
 * @param   [IN]    p_MeanAnomaly               The mean anomaly
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary
 * @return                                      Tuple containing the eccentric anomaly and the true anomaly
 */
DBL_DBL BaseStar::SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity) {

    double e = p_Eccentricity;
    double M = p_MeanAnomaly;
    double E = p_MeanAnomaly;                                                                                                       // inital guess at E is M - correct for e = 0

    double kepler = E - (e * sin(E)) - M;                                                                                           // let f(E) = 0.  Equation (92) in my "A simple toy model" document

    int iteration = 0;
    while (std::abs(kepler) >= NEWTON_RAPHSON_EPSILON && iteration++ < MAX_KEPLER_ITERATIONS) {                                     // repeat the approximation until E is within the specified error of the true value, or max iterations exceeded
        double keplerPrime = 1.0 - (e * cos(E));                                                                                    // derivative of f(E), f'(E).  Equation (94) in my "A simple toy model" document
        E = E - kepler / keplerPrime;
        kepler = E - (e * sin(E)) - M;                                                                                              // let f(E) = 0.  Equation (92) in my "A simple toy model" document
    }

    if (iteration >= MAX_KEPLER_ITERATIONS) SHOW_ERROR(ERROR::NO_CONVERGENCE, "Solving Kepler's equation");                         // show error

    double nu = 2.0 * atan((sqrt((1.0 + e) / (1.0 - e))) * tan(0.5*E));                                                             // convert eccentric anomaly into true anomaly.  Equation (96) in my "A simple toy model" document

         if (utils::Compare(E, M_PI) >= 0 && utils::Compare(E, _2_PI) <= 0) nu += _2_PI;                                            // add 2PI if necessary
    else if (utils::Compare(E, 0.0)  <  0 || utils::Compare(E, _2_PI) >  0) SHOW_WARN(ERROR::OUT_OF_BOUNDS, "Eccentric anomaly");   // out of bounds - show warning

    return std::make_tuple(E, nu);
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
    std::tie(m_SupernovaDetails.eccentricAnomaly, m_SupernovaDetails.trueAnomaly) = SolveKeplersEquation(m_SupernovaDetails.meanAnomaly, p_Eccentricity);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              ENERGY RELATED FUNCTIONS                             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 *	Calculate the absolute value of the binding energy core to the envelope of the star
 *	ALEJANDRO - 08/03/2017
 *
 *
 * double CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda)
 *
 * @param   [IN]    p_CoreMass                  Core mass of the star (Msol)
 * @param   [IN]    p_EnvMass                   Envelope mass of the star (Msol)
 * @param   [IN]    p_Radius                    Radius of the star (Rsol)
 * @return                                      Binding energy (ergs)
 */
double BaseStar::CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda) {

    double bindingEnergy = 0.0;                                                         // default

	if (utils::Compare(p_Radius, 0.0) <= 0) {                                           // positive radius?
        SHOW_WARN(ERROR::RADIUS_NOT_POSITIVE, "Binding energy = 0.0");                  // warn radius not positive
	}
	else if (utils::Compare(p_Lambda, 0.0) <= 0) {                                      // positive lambda?
        // Not necesarily zero as sometimes lambda is made 0, or maybe weird values for certain parameters of the fit. Not sure about the latter.
//        SHOW_WARN(ERROR::LAMBDA_NOT_POSITIVE, "Binding energy = 0.0");                  // warn lambda not positive
	}
	else {                                                                              // calculate binding energy
        // convert to CGS where necessary
        double radius    = p_Radius * RSOL_TO_CM;
        double coreMass  = p_CoreMass * MSOL_TO_G;
        double envMass   = p_EnvMass * MSOL_TO_G;

        double totalMass = coreMass + envMass;                                          // total mass

		bindingEnergy    = G_CGS * totalMass * envMass / (p_Lambda * radius);           // ergs
	}

	return bindingEnergy;
}


/*
 * Calculate all binding energies
 *
 *
 * void CalculateBindingEnergies(const double p_CoreMass, const double p_EnvMass, const double p_Radius)
 *
 * @param   [IN]    p_CoreMass                  Core mass of the star (Msol)
 * @param   [IN]    p_EnvMass                   Envelope mass of the star (Msol)
 * @param   [IN]    p_Radius                    Radius of the star (Rsol)
 */
void BaseStar::CalculateBindingEnergies(const double p_CoreMass, const double p_EnvMass, const double p_Radius) {
    m_BindingEnergies.fixed          = CalculateBindingEnergy(p_CoreMass, p_EnvMass, p_Radius, m_Lambdas.fixed);
	m_BindingEnergies.nanjing        = CalculateBindingEnergy(p_CoreMass, p_EnvMass, p_Radius, m_Lambdas.nanjing);
	m_BindingEnergies.loveridge      = CalculateBindingEnergy(p_CoreMass, p_EnvMass, p_Radius, m_Lambdas.loveridge);
	m_BindingEnergies.loveridgeWinds = CalculateBindingEnergy(p_CoreMass, p_EnvMass, p_Radius, m_Lambdas.loveridgeWinds);
	m_BindingEnergies.kruckow        = CalculateBindingEnergy(p_CoreMass, p_EnvMass, p_Radius, m_Lambdas.kruckow);
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
 * Limit timestep to 1% mass change
 *
 * May modify m_Mdot  -  JR: todo: revisit this
 *
 *
 * double LimitTimestep()
 *
 * @return                                      Timestep
 */
double BaseStar::LimitTimestep(const double p_Dt) {

    double dt = p_Dt;

    // COEN NEIJSSEL 16/11/2016: cap timestep according to max 1% change in mass star.
    if (OPTIONS->UseMassLoss()) {
        double mDot     = CalculateMassLossRate();                                          // First, calculate mass loss rate
        double massLoss = CalculateMassLoss_Static(m_Mass, mDot, dt);                       // Next, calculate mass loss - limited to (mass * MAXIMUM_MASS_LOSS_FRACTION)

        if (utils::Compare(massLoss, 0.0) > 0) {                                            // No change if no mass loss
            double dtWind = massLoss / (mDot * 1.0E6);                                      // Calculate timestep to match (possibly limited) mass loss
                   dt     = min(dt, dtWind);                                                // choose dt
                   m_Mdot = massLoss / (dt * 1.0E6);                                        // Reset mass loss rate to match (possibly limited) mass loss
        }
    }

    return dt;
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

    // the GBParams and Timescale calculations need to be done
    // before the timestep calculation - since the binary code
    // calls this functiom, the GBParams and Timescale functions
    // are called here

    CalculateGBParams();                                                                // calculate giant branch parameters
    CalculateTimescales();                                                              // calculate timescales

    double dt = ChooseTimestep(m_Age);
    if (utils::Compare(TIMESTEP_REDUCTION_FACTOR, 1.0) != 0) {                          // timestep reduction factor == 1.0?
        if (!IsOneOf({ STELLAR_TYPE::MS_LTE_07, STELLAR_TYPE::MS_GT_07 })) {            // no - check stellar type
            dt /= TIMESTEP_REDUCTION_FACTOR;                                            // apply timestep reduction factor
        }
    }

    return LimitTimestep(dt);
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
        m_DtPrev          = m_Dt;
        m_MassPrev        = m_Mass;
        m_RadiusPrev      = m_Radius;
    }

    // the GBParams and Timescale calculations need to be done
    // before taking the timestep - since the binary code ultimately
    // calls this UpdateAttributesAndAgeOneTimestep, the GBParams and
    // Timescale functions are called here.
    //
    // JR: todo: we should revisit where and how often we recalulate
    // GBParams and Timescales.  The problem is that there are multiple
    // entry points into the calculate/take timestep code that it isn't
    // always abvious where we need to do this...  A project for another
    // time.

    CalculateGBParams();                                                                            // calculate giant branch parameters
    CalculateTimescales();                                                                          // calculate timescales
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
 *    - if required, the star's initial mass is changed by the amount passed as the p_DeltaMass0 parameter before
 *      other attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to
 *      the star's mass before the attributes of the star are calculated.  This should be used infrequently, and
 *      is really a kludge because the Mass0 attribute in Hurley et al. 2000 was overloaded after the introduction
 *      of mass loss (see section 7.1).  We should really separate the different uses of Mass0 in the code and
 *      use a different variable - initial mass shouldn't change (other than to initial mass upon entering a
 *      stellar phase - it doesn't make a lot of sense for initial mass to change during evolution through the
 *      phase).         JR: todo: action this paragraph.
 *
 *    - if required, the star is aged by the amount passed as the p_DeltaTime parameter, and the simulation time is
 *      advanced by the same amount, before other attributes are updated.  The p_deltaTime parameter may be zero,
 *      in which case no change is made to the star's age or the physical time attribute.
 *
 *
 * Checks whether the star:
 *    - is a massless remnant (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - has become a supernova (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - should skip this phase for this timestep (checked after applying p_DeltaMass, p_DeltaMass0 and p_DeltaTime)
 *
 * If none of the above are true the star evolves on phase for the specified timestep (which may be 0, in which case
 * the star's attributes other than age are re-calculated), then the need to evolve the star off phase is checked.
 *
 * If p_DeltaMass, p_DeltaMass0 and p_DeltaTime are all passed as zero the checks for massless remnant and supernova
 * are performed (and consequential changes made), but no other changes to the star's attributes are made - unless
 * the p_ForceRecalculate parameter is set true.        JR: todo: I'm not convinced p_ForceRecalculate is necessary - check it
 *
 * The functional return is the stellar type to which the star should evolve.  The returned stellar type is just the
 * stellar type of the star upon entry if it should remain on phase.  The star's stellar type is not changed here.
 *
 *
 * STELLAR_TYPE UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
 *                                                const double p_DeltaMass0,
 *                                                const double p_DeltaTime,
 *                                                const bool   p_ForceRecalculate)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_DeltaTime                 The timestep to evolve in Myr
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the three deltas are 0.0
 *                                              (optional, default = false)
 * @return                                      Stellar type to which star should evolve
 */
STELLAR_TYPE BaseStar::UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                         const double p_DeltaMass0,
                                                         const double p_DeltaTime,
                                                         const bool   p_ForceRecalculate) {
    STELLAR_TYPE stellarType = m_StellarType;                                               // default is no change


    UpdateAttributesAndAgeOneTimestepPreamble(p_DeltaMass, p_DeltaMass0, p_DeltaTime);      // apply mass changes and save current values if required

    if (ShouldBeMasslessRemnant()) {                                                        // ALEJANDRO - 02/12/2016 - Attempt to fix updating the star if it lost all of its mass
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;                                       // JR: should also pik up already massless remnant
    }
    else {
        stellarType = ResolveSupernova();                                                   // handle supernova     JR: moved this to start of timestep
        if (stellarType == m_StellarType) {                                                 // still on phase?

            if (p_ForceRecalculate                     ||                                   // force recalculate?
                utils::Compare(p_DeltaMass,  0.0) != 0 ||                                   // mass change? or...
                utils::Compare(p_DeltaMass0, 0.0) != 0 ||                                   // mass0 change? or...
                utils::Compare(p_DeltaTime,  0.0)  > 0) {                                   // age/time advance?
                                                                                            // yes - update attributes
                AgeOneTimestepPreamble(p_DeltaTime);                                        // advance dt, age, simulation time

                if (ShouldSkipPhase()) {                                                    // skip phase?
                    stellarType = ResolveSkippedPhase();                                    // phase skipped - do what's required
                }
                else {                                                                      // not skipped - execute phase
                    stellarType = EvolveOnPhase();                                          // evolve on phase

                    if (stellarType == m_StellarType) {                                     // still on phase?
                        stellarType = ResolveEndOfPhase();                                  // check for need to move off phase
                    }
                }
            }
        }
    }

    return stellarType;                                                                     // stellar type to which star should evolve
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

    EvolveOneTimestepPreamble();
}


/*
 * Evolve the star on it's current phase - take one timestep on the current phase
 *
 *
 * STELLAR_TYPE EvolveOnPhase()
 *
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::EvolveOnPhase() {

    STELLAR_TYPE stellarType = m_StellarType;

    if (ShouldEvolveOnPhase()) {                                                    // Evolve timestep on phase

        m_Tau         = CalculateTauOnPhase();

        m_COCoreMass  = CalculateCOCoreMassOnPhase();
        m_CoreMass    = CalculateCoreMassOnPhase();
        m_HeCoreMass  = CalculateHeCoreMassOnPhase();

        m_Luminosity  = CalculateLuminosityOnPhase();

        std::tie(m_Radius, stellarType) = CalculateRadiusAndStellarTypeOnPhase();   // Radius and possibly new stellar type

        ResolveEnvelopeMassOnPhase(m_Tau);

        m_Mu          = CalculatePerturbationMuOnPhase();

        PerturbLuminosityAndRadiusOnPhase();

        m_Temperature = CalculateTemperatureOnPhase();

        STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();                       // Resolve envelope loss if it occurs - possibly new stellar type
        if (thisStellarType != m_StellarType) {                                     // thisStellarType overrides stellarType (from CalculateRadiusAndStellarTypeOnPhase())
            stellarType = thisStellarType;
        }
    }

    return stellarType;
}


/*
 * Evolve the star onto the next phase if necessary - take one timestep at the end of the current phase
 *
 *
 * STELLAR_TYPE ResolveEndOfPhase()
 *
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::ResolveEndOfPhase() {

    STELLAR_TYPE stellarType = m_StellarType;

    if (IsEndOfPhase()) {                                                       // End of phase

        stellarType = ResolveEnvelopeLoss();                                    // Resolve envelope loss if it occurs
        if (stellarType == m_StellarType) {                                     // Staying on phase?

            m_Tau         = CalculateTauAtPhaseEnd();

            m_COCoreMass  = CalculateCOCoreMassAtPhaseEnd();
            m_CoreMass    = CalculateCoreMassAtPhaseEnd();
            m_HeCoreMass  = CalculateHeCoreMassAtPhaseEnd();

            m_Luminosity  = CalculateLuminosityAtPhaseEnd();
            
            m_Radius      = CalculateRadiusAtPhaseEnd();

            ResolveEnvelopeMassAtPhaseEnd(m_Tau);

            m_Mu          = CalculatePerturbationMuAtPhaseEnd();

            PerturbLuminosityAndRadiusAtPhaseEnd();

            m_Temperature = CalculateTemperatureAtPhaseEnd();

            stellarType   = EvolveToNextPhase();                                // determine the stellar type to which the star should evolve
        }
    }

    return stellarType;
}
