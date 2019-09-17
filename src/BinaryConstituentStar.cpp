#include "BinaryConstituentStar.h"

#include "BaseBinaryStar.h"


class BaseBinaryStar;


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
COMPAS_VARIABLE BinaryConstituentStar::StellarPropertyValue(const T_ANY_PROPERTY p_Property) const {

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
            case ANY_STAR_PROPERTY::EXPERIENCED_RLOF:                                   value = ExperiencedRLOF();                                      break;
            case ANY_STAR_PROPERTY::IS_RLOF:                                            value = IsRLOF();                                               break;
            case ANY_STAR_PROPERTY::MASS_LOSS_DIFF:                                     value = MassLossDiff();                                         break;
            case ANY_STAR_PROPERTY::MASS_TRANSFER_CASE_INITIAL:                         value = static_cast<int>(MassTransferCaseInitial());            break;
            case ANY_STAR_PROPERTY::MASS_TRANSFER_DIFF:                                 value = MassTransferDiff();                                     break;
            case ANY_STAR_PROPERTY::ORBITAL_ENERGY_PRE_SUPERNOVA:                       value = PreSNeOrbitalEnergy();                                  break;
            case ANY_STAR_PROPERTY::ORBITAL_ENERGY_POST_SUPERNOVA:                      value = PostSNeOrbitalEnergy();                                 break;

            default:                                                                                                    // not a constitient star property - try underlying star
                std::tie(ok, value) = Star::StellarPropertyValue(p_Property);
        }
    }

    return std::make_tuple(ok, value);
}


/*
 * Calculate the mass accreted by a Neutron Star given mass and radius of companion
 *
 * JR: todo: flesh-out this documentation
 *
 *
 * double CalculateMassAccretedForNS(const double p_CompanionMass, const double p_CompanionRadius)
 *
 * @param   [IN]    p_CompanionMass             The mass of the companion star (Msol)
 * @param   [IN]    p_CompanionRadius           The radius of the companion star (Rsol)
 * @return                                      Mass accreted by the Neutron Star (Msol)
 */
double BinaryConstituentStar::CalculateMassAccretedForNS(const double p_CompanionMass, const double p_CompanionRadius) {

     double deltaMass;

     switch (OPTIONS->CommonEnvelopeMassAccretionPrescription()) {                                              // which prescription?

        case CE_ACCRETION_PRESCRIPTION::ZERO:                                                                   // ZERO
            deltaMass = 0.0;
            break;

        case CE_ACCRETION_PRESCRIPTION::CONSTANT:                                                               // CONSTANT
            deltaMass = OPTIONS->CommonEnvelopeMassAccretionConstant();                                         // use program option
            break;

        case CE_ACCRETION_PRESCRIPTION::UNIFORM:                                                                // UNIFROM
            deltaMass = RAND->Random(OPTIONS->CommonEnvelopeMassAccretionMin(), OPTIONS->CommonEnvelopeMassAccretionMax()); // uniform random distribution - Oslowski+ (2011)
            break;

        case CE_ACCRETION_PRESCRIPTION::MACLEOD: {                                                              // MACLEOD
                                                                                                                // linear regression estimated from Macleod+ (2015)
            double mm = -1.0714285714285712E-05;                                                                // gradient of the linear fit for gradient
            double cm =  0.00012057142857142856;                                                                // intercept of the linear fit for gradient
            double mc =  0.01588571428571428;                                                                   // gradient of the linear fit for intercept
            double cc = -0.15462857142857137;                                                                   // intercept of the linear fir for intercept
            double m  = mm * p_CompanionMass + cm ;                                                             // gradient of linear fit for mass
            double c  = mc * p_CompanionMass + cc ;                                                             // intercept of linear fit for mass

            // calculate mass to accrete and clamp to minimum and maximum from program options
            deltaMass = std::min(OPTIONS->CommonEnvelopeMassAccretionMax(), std::max(OPTIONS->CommonEnvelopeMassAccretionMin(), m * p_CompanionRadius + c));
            } break;

        default:                                                                                                // unknown common envelope accretion prescription - shouldn't happen
            deltaMass = 0.0;                                                                                    // default value
            SHOW_WARN(ERROR::UNKNOWN_CE_ACCRETION_PRESCRIPTION, "NS accreted mass = 0.0");                      // warn that an error occurred
    }

    return deltaMass;
}


/*
 * Resolve common envelope accretion
 *
 * For stellar types other than Neutron Star just set the star's mass to the parameter passed
 * For Neutron Stars calculate the mass accreted based on the companion's mass and radius
 *
 *
 * void ResolveCommonEnvelopeAccretion(const double p_FinalMass)
 *
 * @param   [IN]    p_FinalMass                 Mass of the star post mass transfer (Msol)
 */
void BinaryConstituentStar::ResolveCommonEnvelopeAccretion(const double p_FinalMass) {

    double deltaMass;

    if (IsOneOf({ STELLAR_TYPE::NEUTRON_STAR})) {           // only Neutron Star is different, so settled for this way rather than use class hierarchy...
        deltaMass = CalculateMassAccretedForNS(m_Companion->Mass(), m_Companion->Radius());
        m_MassTransferDiff = deltaMass;
    }
    else {
        deltaMass = p_FinalMass - Mass();
        // JR: todo: why isn't m_MassTransferDiff updated here (as it is for Neutron Stars)?
    }

    ResolveAccretion(deltaMass);
}


/*
 * Calculate the synchronisation timescale for the star
 *
 * Hurley et al. 2002, section 2.3
 *
 *
 * double CalculateSynchronisationTimescale(const double p_SemiMajorAxis)
 *
 * @return                                      synchronisation timescale for the star (yr)
 */
	// Function to calculate the synchronization timescale
	// For details see sec. 2.3 of Hurley+2002
	// [syncrhonizationTimescale] = yr
double BinaryConstituentStar::CalculateSynchronisationTimescale(const double p_SemiMajorAxis) {

    double gyrationRadiusSquared_1 = 1.0 / CalculateGyrationRadius();
    double rOverA                  = Radius() / p_SemiMajorAxis;
    double rOverA_6                = rOverA * rOverA * rOverA * rOverA * rOverA * rOverA;
    double q2			           = m_Companion->Mass() / Mass();

	double denominator;
	switch (DetermineEnvelopeTypeHurley2002()) {                                // JR: todo: this differs from envelopeType() in star.cpp and DetermineEnvelopeType() in new code - ok?

        case ENVELOPE::CONVECTIVE: {                                            // solve for stars with convective envelope, according to tides section (see Hurley et al. 2002, subsection 2.3.1)

                      double tauConv = CalculateEddyTurnoverTimescale();
                      double fConv   = 1.0;	                                    // currently, as COMPAS doesn't have rotating stars tested, we set f_conv = 1 always.
                      double kOverTc = (2.0 / 21.0) * (fConv / tauConv) * ((Mass() - CoreMass()) / Mass());

            denominator = 3.0 * kOverTc * q2 * gyrationRadiusSquared_1 * rOverA_6;
            } break;

        case ENVELOPE::RADIATIVE: {                                             // solve for stars with radiative envelope (see Hurley et al. 2002, subsection 2.3.2)

            constexpr double coeff2          = pow(52.0, 5.0 / 3.0);
                      double e2              = 1.592E-9 * pow(Mass(), 2.84);    // second order tidal coefficient (a.k.a. E_2)
                      double rAU             = Radius() * RSOL_TO_AU;
                      double rAU_3           = rAU * rAU * rAU;
                      double freeFallFactor  = sqrt(G1 * Mass() / rAU_3);

		    denominator = coeff2 * freeFallFactor * gyrationRadiusSquared_1 * q2 * q2 * pow(1.0 + q2, 5.0 / 6.0) * e2 * pow(rOverA, 17.0 / 2.0);
            } break;

        default:                                                                // all other envelope types (remnants?)
            denominator = (1.0 / 1.3E7) * pow(Luminosity() / Mass(), 5.0 / 7.0) * rOverA_6;
	}

	return 1.0 / denominator;
}


/*
 * Set the Roche Lobe flags for a star based on its Roche Lobe radius
 *
 * Changes class member struct m_Flags
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

    m_RocheLobeTracker = (Radius() * RSOL_TO_AU) / (m_RocheLobeRadius * p_SemiMajorAxis * (1.0 - p_Eccentricity));      // ratio of star's size to its Roche Lobe radius, calculated at the point of closest approach, periapsis

    if (utils::Compare(m_RocheLobeTracker, 1.0) >= 0) {                                                                 // if star is equal to or larger than its Roche Lobe...
		m_RLOFDetails.isRLOF          = true;                                                                           // ... it is currently Roche Lobe overflowing
		m_RLOFDetails.experiencedRLOF = true;                                                                           // ... and for future checks, did Roche Lobe overflow
	}

	m_RLOFDetails.RLOFPostCEE = m_RLOFDetails.isRLOF && p_CommonEnvelope ? true : m_RLOFDetails.RLOFPostCEE;            // check for RLOF just after the CEE     JR: todo: should the else part be false?
}


/*
 * Determine initial mass transfer case
 *
 * Three cases:
 *
 *    Case A: mass transfer while donor is on main sequence
 *    Case B: donor star is in (or evolving to) Red Giant phase
 *    Case C: SuperGiant phase
 *
 * Modifies class member variables m_MassTransferCaseInitial and m_FirstMassTransferEpisode
 *
 * void DetermineInitialMassTransferCase()
 */
void BinaryConstituentStar::DetermineInitialMassTransferCase() {
    if (!m_FirstMassTransferEpisode) m_MassTransferCaseInitial = DetermineMassTransferCase();       // JR: todo: are these actually used anywhere?  Only for printing perhaps...
    m_FirstMassTransferEpisode = true;                                                              // JR: todo: are these actually used anywhere?  Only for printing perhaps...
}


/*
 * Initial calculations for mass transfer resolution
 *
 * Calculates the Roche Lobe radius based on mass of both stars
 * Sets Roche Log flags for the star
 * Set class member variable m_MassTransferDiff = 0.0
 *
 * void InitialiseMassTransfer(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity)
 *
 * @param   [IN]    p_CommonEnvelope            Indicates whether a common envelope event is occurring
 * @param   [IN]    p_SemiMajorAxis             Semi major axis of the binary (in AU)
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary orbit
 */
void BinaryConstituentStar::InitialiseMassTransfer(const bool p_CommonEnvelope, const double p_SemiMajorAxis, const double p_Eccentricity) {
    m_RocheLobeRadius = BaseBinaryStar::CalculateRocheLobeRadius_Static(Mass(), m_Companion->Mass());
    SetRocheLobeFlags(p_CommonEnvelope, p_SemiMajorAxis, p_Eccentricity);
    m_MassTransferDiff = 0.0;
}
