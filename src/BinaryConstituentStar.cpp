#include "BinaryConstituentStar.h"

#include "BaseBinaryStar.h"


class BaseBinaryStar;


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
 *    STAR_PROPERTY, STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY only.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.
 *
 *
 * COMPAS_VARIABLE StellarPropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BinaryConstituentStar::StellarPropertyValue(const T_ANY_PROPERTY p_Property) const {

    COMPAS_VARIABLE value;

    ANY_STAR_PROPERTY property;

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {

        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : { STAR_PROPERTY      prop = boost::get<STAR_PROPERTY>(p_Property);      property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : { STAR_1_PROPERTY    prop = boost::get<STAR_1_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : { STAR_2_PROPERTY    prop = boost::get<STAR_2_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: { SUPERNOVA_PROPERTY prop = boost::get<SUPERNOVA_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: { COMPANION_PROPERTY prop = boost::get<COMPANION_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;

        default:                                                                                                                // unexpected property type
            // the only ways this can happen are if someone added a stellar type property (into ANY_PROPERTY_TYPE)
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without a
            // warning - this is a code defect, so we flag it as an error and that will result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing property type or find and fix the code defect.

            THROW_ERROR(ERROR::UNEXPECTED_STELLAR_PROPERTY_TYPE);                                                               // throw error
    }

    switch (property) {
        case ANY_STAR_PROPERTY::BINDING_ENERGY_AT_COMMON_ENVELOPE:                  value = BindingEnergyAtCEE();                           break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE:                 value = BindingEnergyPreCEE();                          break;
        case ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMMON_ENVELOPE:                    value = COCoreMassAtCEE();                              break;
        case ANY_STAR_PROPERTY::CORE_MASS_AT_COMMON_ENVELOPE:                       value = CoreMassAtCEE();                                break;
        case ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE:           value = DynamicalTimescalePostCEE();                    break;
        case ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE:            value = DynamicalTimescalePreCEE();                     break;
        case ANY_STAR_PROPERTY::EXPERIENCED_RLOF:                                   value = ExperiencedRLOF();                              break;
        case ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMMON_ENVELOPE:                    value = HeCoreMassAtCEE();                              break;
        case ANY_STAR_PROPERTY::IS_RLOF:                                            value = IsRLOF();                                       break;
        case ANY_STAR_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE:                          value = LambdaAtCEE();                                  break;
        case ANY_STAR_PROPERTY::LUMINOSITY_POST_COMMON_ENVELOPE:                    value = LuminosityPostCEE();                            break;
        case ANY_STAR_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE:                     value = LuminosityPreCEE();                             break;
        case ANY_STAR_PROPERTY::MASS_LOSS_DIFF:                                     value = MassLossDiff();                                 break;
        case ANY_STAR_PROPERTY::MASS_TRANSFER_DIFF:                                 value = MassTransferDiff();                             break;
        case ANY_STAR_PROPERTY::ORBITAL_ENERGY_POST_SUPERNOVA:                      value = OrbitalEnergyPostSN();                          break;
        case ANY_STAR_PROPERTY::ORBITAL_ENERGY_PRE_SUPERNOVA:                       value = OrbitalEnergyPreSN();                           break;
        case ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE:    value = RadialExpansionTimescalePostCEE();              break;
        case ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE:     value = RadialExpansionTimescalePreCEE();               break;
        case ANY_STAR_PROPERTY::RECYCLED_NEUTRON_STAR:                              value = ExperiencedRecycledNS();                        break;
        case ANY_STAR_PROPERTY::TEMPERATURE_POST_COMMON_ENVELOPE:                   value = TemperaturePostCEE() * TSOL;                    break;
        case ANY_STAR_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE:                    value = TemperaturePreCEE() * TSOL;                     break;
        case ANY_STAR_PROPERTY::THERMAL_TIMESCALE_POST_COMMON_ENVELOPE:             value = ThermalTimescalePostCEE();                      break;
        case ANY_STAR_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE:              value = ThermalTimescalePreCEE();                       break;

        default:
            // not a constitient star property
            // try underlying star - any errors will be handled there
            value = Star::StellarPropertyValue(p_Property);
    }

    return value;
}


/*
 * Calculate (or set) pre common envelope values:
 *
 *    m_CEDetails.preCEE.bindingEnergy
 *    m_CEDetails.preCEE.dynamicalTimescale
 *    m_CEDetails.preCEE.eccentricity
 *    m_CEDetails.preCEE.luminosity
 *    m_CEDetails.preCEE.mass
 *    m_CEDetails.preCEE.radialExpansionTimescale
 *    m_CEDetails.preCEE.radius
 *    m_CEDetails.preCEE.semiMajorAxis
 *    m_CEDetails.preCEE.stellarType
 *    m_CEDetails.preCEE.temperature
 *    m_CEDetails.preCEE.thermalTimescale
 *
 *
 * void SetPreCEEValues()
 */
void BinaryConstituentStar::SetPreCEEValues() {

    m_CEDetails.preCEE.bindingEnergy            = m_CEDetails.bindingEnergy;
    m_CEDetails.preCEE.dynamicalTimescale       = CalculateDynamicalTimescale();
    m_CEDetails.preCEE.luminosity               = Luminosity();
    m_CEDetails.preCEE.mass                     = Mass();
    m_CEDetails.preCEE.radialExpansionTimescale = CalculateRadialExpansionTimescale();
    m_CEDetails.preCEE.radius                   = Radius();
    m_CEDetails.preCEE.stellarType              = StellarType();
    m_CEDetails.preCEE.temperature              = Temperature();
    m_CEDetails.preCEE.thermalTimescale         = CalculateThermalTimescale();
}


/*
 * Calculate (or set) post common envelope values:
 *
 *    m_CEDetails.postCEE.dynamicalTimescale
 *    m_CEDetails.postCEE.eccentricity
 *    m_CEDetails.postCEE.luminosity
 *    m_CEDetails.postCEE.mass
 *    m_CEDetails.postCEE.radialExpansionTimescale
 *    m_CEDetails.postCEE.radius
 *    m_CEDetails.postCEE.semiMajorAxis
 *    m_CEDetails.postCEE.stellarType
 *    m_CEDetails.postCEE.temperature
 *    m_CEDetails.postCEE.thermalTimescale
 *
 *
 * void SetPostCEEValues()
 */
void BinaryConstituentStar::SetPostCEEValues() {

    m_CEDetails.postCEE.dynamicalTimescale       = CalculateDynamicalTimescale();
    m_CEDetails.postCEE.luminosity               = Luminosity();
    m_CEDetails.postCEE.mass                     = Mass();
    m_CEDetails.postCEE.radialExpansionTimescale = CalculateRadialExpansionTimescale();
    m_CEDetails.postCEE.radius                   = Radius();
    m_CEDetails.postCEE.stellarType              = StellarType();
    m_CEDetails.postCEE.temperature              = Temperature();
    m_CEDetails.postCEE.thermalTimescale         = CalculateThermalTimescale();
}


/*
 * Calculate or set common envelope values:
 *
 *    m_CEDetails.HeCoreMass
 *    m_CEDetails.COCoreMass
 *    m_CEDetails.CoreMass
 *    m_CEDetails.bindingEnergy
 *    m_CEDetails.lambda
 *    m_CEDetails.convectiveEnvelopeMass
 *    m_CEDetails.radiativeIntershellMass
 *    m_CEDetails.convectiveEnvelopeBindingEnergy
 *
 * void CalculateCommonEnvelopeValues()
 */
void BinaryConstituentStar::CalculateCommonEnvelopeValues() {

    m_CEDetails.HeCoreMass = HeCoreMass();
    m_CEDetails.COCoreMass = COCoreMass();
    m_CEDetails.CoreMass   = CoreMass();

    m_CEDetails.lambda     = 0.0;                                                                               // default

    switch (OPTIONS->CommonEnvelopeLambdaPrescription()) {                                                      // which common envelope lambda prescription?

        case CE_LAMBDA_PRESCRIPTION::FIXED:
            m_CEDetails.lambda = OPTIONS->CommonEnvelopeLambda();
            break;

        case CE_LAMBDA_PRESCRIPTION::LOVERIDGE:
            m_CEDetails.lambda = CalculateLambdaLoveridge();
            break;

        case CE_LAMBDA_PRESCRIPTION::NANJING:
            m_CEDetails.lambda = CalculateLambdaNanjing();
            break;

        case CE_LAMBDA_PRESCRIPTION::KRUCKOW:
            m_CEDetails.lambda = CalculateLambdaKruckow();
            break;
            
        case CE_LAMBDA_PRESCRIPTION::DEWI:
            m_CEDetails.lambda = CalculateLambdaDewi();
            break;

        default:                                                                                                // unknown prescription
            // the only way this can happen is if someone added a CE_LAMBDA_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_CE_LAMBDA_PRESCRIPTION);                                                 // throw error            
    }

    if (utils::Compare(m_CEDetails.lambda, 0.0) <= 0) m_CEDetails.lambda = 0.0;                                 // force non-positive lambda to 0

    m_CEDetails.lambda *= OPTIONS->CommonEnvelopeLambdaMultiplier();                                            // multiply by constant (program option, default = 1.0)
                                                                        
    m_CEDetails.bindingEnergy = CalculateBindingEnergy(CoreMass(), Mass() - CoreMass(), Radius(), m_CEDetails.lambda);
    
    // properties relevant for the Hirai & Mandel (2022) formalism
    double maxConvectiveEnvelopeMass;
    std::tie(m_CEDetails.convectiveEnvelopeMass, maxConvectiveEnvelopeMass) = CalculateConvectiveEnvelopeMass();
    m_CEDetails.radiativeIntershellMass = Mass() - CoreMass() - m_CEDetails.convectiveEnvelopeMass;

    if (OPTIONS->CommonEnvelopeFormalism() == CE_FORMALISM::TWO_STAGE) {
        m_CEDetails.lambda = CalculateConvectiveEnvelopeLambdaPicker(std::tie(m_CEDetails.convectiveEnvelopeMass, maxConvectiveEnvelopeMass));
        m_CEDetails.bindingEnergy = CalculateConvectiveEnvelopeBindingEnergy(Mass(), m_CEDetails.convectiveEnvelopeMass, Radius(), m_CEDetails.lambda);
    }
}


/*
 * Calculate the circularisation timescale for the star
 *
 * Hurley et al. 2002, section 2.3
 *
 *
 * double CalculateCircularisationTimescale(const double p_SemiMajorAxis)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major Axis of the binary (Rsol)
 * @return                                      Circularisation timescale for the star (yr)
 */
double BinaryConstituentStar::CalculateCircularisationTimescale(const double p_SemiMajorAxis) {

    double q2     = m_Companion->Mass() / Mass();
    double rOverA = Radius() / p_SemiMajorAxis;

    double timescale;

	switch (DetermineEnvelopeType()) {

        case ENVELOPE::CONVECTIVE: {                                                                                                    // solve for stars with convective envelope, according to tides section (see Hurley et al. 2002, subsection 2.3.1)

	        double tauConv          = CalculateEddyTurnoverTimescale();
	        double fConv            = 1.0;                                                                                              // currently, as COMPAS doesn't have rotating stars tested, we set f_conv = 1 always.
            double fConvOverTauConv = fConv / tauConv;
            double rOverAPow8       = rOverA * rOverA * rOverA * rOverA * rOverA * rOverA * rOverA * rOverA;                            // use multiplication - pow() is slow

	        timescale               = 1.0 / (fConvOverTauConv * ((Mass() - CoreMass()) / Mass()) * q2 * (1.0 + q2) * rOverAPow8);
        } break;

        case ENVELOPE::RADIATIVE: {                                                                                                     // solve for stars with radiative envelope (see Hurley et al. 2002, subsection 2.3.2)

            double rInAU                  = Radius() * RSOL_TO_AU;
            double rInAUPow3              = rInAU * rInAU * rInAU;                                                                      // use multiplication - pow() is slow
            double rOverAPow10            = rOverA * rOverA * rOverA * rOverA * rOverA * rOverA * rOverA * rOverA * rOverA * rOverA;    // use multiplication - pow() is slow
            double rOverAPow21Over2       = rOverAPow10 * rOverA * std::sqrt(rOverA);                                                   // sqrt() is faster than pow()

		    double	secondOrderTidalCoeff = 1.592E-09 * PPOW(Mass(), 2.84);                                                             // aka E_2.    
		    double	freeFallFactor        = std::sqrt(G_AU_Msol_yr * Mass() / rInAUPow3);
		
		    timescale                     = 1.0 / ((21.0 / 2.0) * freeFallFactor * q2 * PPOW(1.0 + q2, 11.0/6.0) * secondOrderTidalCoeff * rOverAPow21Over2);
            } break;

        case ENVELOPE::REMNANT:                                                                                                         // remnants
            timescale = 0.0;
            break;

        default:                                                                                                                        // unknown envelope type
            // the only way this can happen is if someone added a new envelope type (to ENVELOPE)
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because DetermineEnvelopeType() returned an envelope type this code doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution of the
            // star or binary.
            // The correct fix for this is to add code for the missing envelope type or, if the missing envelope
            // type is incorrect/superfluous, remove it from ENVELOPE.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_TYPE);                                                                                  // throw error
        }

	return timescale;
}


/*
 * Calculate the synchronisation timescale for the star
 *
 * Hurley et al. 2002, section 2.3
 *
 *
 * double CalculateSynchronisationTimescale(const double p_SemiMajorAxis)
 *
 * @param   [IN]    p_SemiMajorAxis             Semi-major Axis of the binary (Rsol)
 * @return                                      Synchronisation timescale for the star (yr)
 */
double BinaryConstituentStar::CalculateSynchronisationTimescale(const double p_SemiMajorAxis) {

    double gyrationRadiusSquared   = CalculateMomentOfInertia() / Mass() / Radius() / Radius();
    double gyrationRadiusSquared_1 = 1.0 / gyrationRadiusSquared;
    double rOverA                  = Radius() / p_SemiMajorAxis;
    double rOverA_6                = rOverA * rOverA * rOverA * rOverA * rOverA * rOverA;
    double q2			           = m_Companion->Mass() / Mass();

	double timescale;

	switch (DetermineEnvelopeType()) {

            case ENVELOPE::CONVECTIVE: {                                                            // solve for stars with convective envelope, according to tides section (see Hurley et al. 2002, subsection 2.3.1)

                double tauConv = CalculateEddyTurnoverTimescale();
                double fConv   = 1.0;                                                               // currently, as COMPAS doesn't have rotating stars tested, we set f_conv = 1 always.
                double kOverTc = (2.0 / 21.0) * (fConv / tauConv) * ((Mass() - CoreMass()) / Mass());

                timescale      = 1.0 / (3.0 * kOverTc * q2 * gyrationRadiusSquared_1 * rOverA_6);
            } break;

            case ENVELOPE::RADIATIVE: {                                                             // solve for stars with radiative envelope (see Hurley et al. 2002, subsection 2.3.2)

                double coeff2         = 15.874010519681995;                                         // 5.0 * PPOW(2.0, 5.0 / 3.0) = 5.0 * 3.174802103936399
                double e2             = 1.592E-9 * PPOW(Mass(), 2.84);                              // second order tidal coefficient (a.k.a. E_2)
                double rAU            = Radius() * RSOL_TO_AU;
                double rAU_3          = rAU * rAU * rAU;
                double freeFallFactor = std::sqrt(G_AU_Msol_yr * Mass() / rAU_3);

		        timescale             = 1.0 / (coeff2 * freeFallFactor * gyrationRadiusSquared_1 * q2 * q2 * PPOW(1.0 + q2, 5.0 / 6.0) * e2 * PPOW(rOverA, 17.0 / 2.0));
            } break;

        case ENVELOPE::REMNANT:                                                                     // remnants
            timescale = 1.0 / ((1.0 / 1.3E7) * PPOW(Luminosity() / Mass(), 5.0 / 7.0) * rOverA_6);
            break;

        default:                                                                                    // unknown envelope type
            // the only way this can happen is if someone added a new envelope type (to ENVELOPE)
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because DetermineEnvelopeType() returned an envelope type this code doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution of the
            // star or binary.
            // The correct fix for this is to add code for the missing envelope type or, if the missing envelope
            // type is incorrect/superfluous, remove it from ENVELOPE.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_TYPE);                                              // throw error
	}

	return timescale;
}


/*
 * Set the Roche Lobe flags for a star based on its Roche Lobe radius
 *
 * Changes class member struct m_RLOFDetails
 *
 *
 * double SetRocheLobeFlags(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity)
 *
 * @param   [IN]    p_CommonEnvelope            Indicates whether a common envelope event is occurring
 * @param   [IN]    p_SemiMajorAxis             Semi major axis of the binary (in AU)
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary orbit
 */
void BinaryConstituentStar::SetRocheLobeFlags(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity) {

    m_RLOFDetails.isRLOF = false;                                                                                       // default - not overflowing Roche Lobe

    double starToRocheLobeRadiusRatio = StarToRocheLobeRadiusRatio(p_SemiMajorAxis, p_Eccentricity);

    if (utils::Compare(starToRocheLobeRadiusRatio, 1.0) >= 0) {                                                         // if star is equal to or larger than its Roche Lobe...
        m_RLOFDetails.isRLOF          = true;                                                                           // ... it is currently Roche Lobe overflowing
		m_RLOFDetails.experiencedRLOF = true;                                                                           // ... and for future checks, did Roche Lobe overflow
	}

	m_RLOFDetails.RLOFPostCEE = m_RLOFDetails.isRLOF && p_CommonEnvelope ? true : m_RLOFDetails.RLOFPostCEE;            // check for RLOF just after the CEE (if this flag was ever true for this system, it remains true)
}


/*
 * Ratio of star's radius to its Roche Lobe radius, calculated at the point of closest approach, periapsis
 *
 * double StarToRocheLobeRadiusRatio() const
 *
 * @param   [IN]    p_SemiMajorAxis             Semi major axis of the binary (in AU)
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary orbit
 * @return                                      Ratio of star's radius to its Roche lobe radius
 */
double BinaryConstituentStar::StarToRocheLobeRadiusRatio(const double p_SemiMajorAxis, const double p_Eccentricity) {
    // binary is unbound or semi-major axis is infinite (evolving single star as binary), so not in RLOF
    // JR: note, this will (probably) fail if option --fp-error-mode is not OFF (the calculation that resulted in p_SemiMajorAxis = inf will (probably) result in a trap)
    // Maybe use !std::isfinite(p_SemiMajorAxis) to ensure p_SemiMajorAxis is not NaN or inf
    if ((utils::Compare(p_SemiMajorAxis, 0.0) <= 0) || (utils::Compare(p_Eccentricity, 1.0) > 0) || isinf(p_SemiMajorAxis)) return 0.0;
    
    double rocheLobeRadius = BaseBinaryStar::CalculateRocheLobeRadius_Static(Mass(), m_Companion->Mass());
    return (Radius() * RSOL_TO_AU) / (rocheLobeRadius * p_SemiMajorAxis * (1.0 - p_Eccentricity));
}


/*
 * Initial calculations for mass transfer resolution
 *
 * Calculates the Roche Lobe radius based on mass of both stars
 * Sets Roche Lobe flags for the star
 * Set class member variable m_MassTransferDiff = 0.0
 *
 * void InitialiseMassTransfer(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity)
 *
 * @param   [IN]    p_CommonEnvelope            Indicates whether a common envelope event is occurring
 * @param   [IN]    p_SemiMajorAxis             Semi major axis of the binary (in AU)
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary orbit
 */
void BinaryConstituentStar::InitialiseMassTransfer(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity) {
    SetRocheLobeFlags(p_CommonEnvelope, p_SemiMajorAxis, p_Eccentricity);
    m_MassTransferDiff = 0.0;
}


/* 
 * Set class member variable m_MassTransferDiff based on the mass gained or lost during MT.
 *
 * Also modify changes to the H/He shell (relevant only for WDs)
 *
 * void SetMassTransferDiffAndResolveWDShellChange(const double p_MassTransferDiff) 
 *
 * @param   [IN]    p_MassTransferDiff          The amount of mass lost or gained by the star 
 */
void BinaryConstituentStar::SetMassTransferDiffAndResolveWDShellChange(const double p_MassTransferDiff) {
    m_MassTransferDiff = p_MassTransferDiff; 
    ResolveShellChange(p_MassTransferDiff);       // only applies to WDs
}

