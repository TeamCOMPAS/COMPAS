#include "Rand.h"
#include "NS.h"


/*
 * Calculate the luminosity of a Neutron Star
 *
 * Hurley et al. 2000, eq 93
 *
 * Called from CalculateCoreCollapseSNParams_Static(), so must be static.
 * 
 * 
 * double CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time since formation of the object in Myr
 * @return                                      Luminosity of the Neutron Star in Lsol
 */
double NS::CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time) {
    double t = std::max(p_Time, 0.1);
    return 0.02 * PPOW(p_Mass, 2.0 / 3.0) / (t * t);
}

/*
 * Choose timestep for Pulsar Evolution
 *
 * Pulsars evolve very fast when they are first born, and evolve slower as they age.
 * Hence, timestep is chosen to be small when pulsar is young, and is slowly increased
 * as the pulsar ages.
 *
 * 
 * double ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double NS::ChooseTimestep(const double p_Time) const {

    double result = 500.0;                                      // default value

         if (p_Time < 0.01 ) result = 0.001;
    else if (p_Time < 0.1  ) result = 0.01;
    else if (p_Time < 1.0  ) result = 0.1;
    else if (p_Time < 3000.0 ) result = 1.0;
    else if (p_Time > 3000.0 ) result = 100.0;
    // else if (p_Time < 500.0) {
    //     double slope      = 1.58859191006;                      // 1.58859191006 = log10(500.0) / (log10(500.0) - 1.0)
    //     double log10_step = slope * (log10(p_Time) - 1.0);
    //     result            = PPOW(10.0, log10_step);
    // }

    return result;
}


/*
 * Calculate Neutron Star radius according to selected equation of state (by commandline option)
 *
 * Indirectly called from CalculateCoreCollapseSNParams_Static(), so must be static.
 * 
 * 
 * double CalculateRadiusOnPhaseInKM(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius of Neutron Star in km
 */
double NS::CalculateRadiusOnPhaseInKM_Static(const double p_Mass) {

    double radius;

    switch (OPTIONS->NeutronStarEquationOfState()) {                                            // which equation-of-state?

        case NS_EOS::SSE:                                                                       // SSE
            radius = 10.0;
            break;

        case NS_EOS::ARP3: {                                                                    // ARP3

            // use table ARP3MassRadiusRelation defined in constants.h
            // don't extrapolate - masses outside table just set to extreme values

            std::map<double, double>::const_iterator iter = ARP3MassRadiusRelation.begin();     // first element;

            double ARP3MinimumMass         = iter->first;                                       // first element: mass
            double ARP3RadiusAtMinimumMass = iter->second;                                      // first element: radius

            iter = ARP3MassRadiusRelation.end();                                                // last + 1 element
            double ARP3MaximumMass = (--iter)->first;                                           // last element: mass
            double ARP3RadiusAtMaximumMass = iter->second;                                      // last element: radius

            if (utils::Compare(p_Mass, ARP3MinimumMass) < 0) {                                  // mass < minimum?
                radius = ARP3RadiusAtMinimumMass;                                               // yes, clamp to minimum
            }
            else if (utils::Compare(p_Mass, ARP3MaximumMass) > 0) {                             // not < minimum; > maximum?
                radius = ARP3RadiusAtMaximumMass;                                               // yes, clamp to maximum
            }
            else {
                radius = utils::SampleFromTabulatedCDF(p_Mass, ARP3MassRadiusRelation);         // no - mass in range - sample
            }
        } break;

        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added an NS_EOS
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNKNOWN_NS_EOS);                                          // throw error
	}

	return radius;
}


/*
 * Calculate core collapse Supernova parameters
 *
 * Called from GiantBranch, so must be static.
 * 
 * 
 * DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Tuple containing Luminosity, Radius and Temperature of Neutron Star
 */
DBL_DBL_DBL NS::CalculateCoreCollapseSNParams_Static(const double p_Mass) {
    double luminosity  = CalculateLuminosityOnPhase_Static(p_Mass, 0.0);                                        // luminosity of Neutron Star as it cools
    double radius      = CalculateRadiusOnPhase_Static(p_Mass);                                                 // radius of Neutron Star in Rsol
    double temperature = BaseStar::CalculateTemperatureOnPhase_Static(luminosity, radius);                      // temperature of NS

    return std::make_tuple(luminosity, radius, temperature);
}


/*
 * Calculate the spin period of a Pulsar at birth according to selected distribution (by commandline option)
 * Users should note that when choosing the NOSPIN option, 
 * pulsar spin frequency is set to 0 and spin period is infinity. 
 *
 * double CalculateBirthSpinPeriod()
 *
 * @return                                      Birth spin period of Pulsar in s
 */
double NS::CalculateBirthSpinPeriod() {

	double pSpin;

    if ((OPTIONS->MSPsFromAIC())&(ExperiencedAIC())){
        switch (OPTIONS->MSPBirthSpinPeriodDistribution()) {                                                     // which distribution?

            case MSP_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM: {                                                  // UNIFORM distribution between minimum and maximum value as in Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (default Pmin = and Pmax = )
                                                                                                                    // and also Kiel et al 2008 https://arxiv.org/abs/0805.0059 (default Pmin = 10 ms and Pmax 100 ms, section 3.4)
                double maximum = OPTIONS->MSPBirthSpinPeriodDistributionMax();
                double minimum = OPTIONS->MSPBirthSpinPeriodDistributionMin();

                pSpin = minimum + (RAND->Random() * (maximum - minimum));
                } break;

            case MSP_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL: {                                                   // NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

                double mean  = OPTIONS->MSPBirthSpinPeriodDistributionMean();
                double sigma = OPTIONS->MSPBirthSpinPeriodDistributionSigma();

                // this should terminate naturally, but just in case we add a guard
                std::size_t iterations = 0;
                do { pSpin = RAND->RandomGaussian(sigma) + mean;} while (iterations++ < PULSAR_SPIN_ITERATIONS && utils::Compare(pSpin, 0.0) <= 0);
                if (iterations >= PULSAR_SPIN_ITERATIONS) THROW_ERROR(ERROR::TOO_MANY_PULSAR_SPIN_ITERATIONS);

                } break;

            case MSP_BIRTH_SPIN_PERIOD_DISTRIBUTION::LOGNORMAL: {                                                   // NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

                double mean  = log10(OPTIONS->MSPBirthSpinPeriodDistributionMean());
                double sigma = log10(OPTIONS->MSPBirthSpinPeriodDistributionSigma());
                double log10pSpin;
                // this should terminate naturally, but just in case we add a guard
                std::size_t iterations = 0;
                do { log10pSpin = RAND->RandomGaussian(sigma) + mean;} while (iterations++ < PULSAR_SPIN_ITERATIONS && utils::Compare(log10pSpin, 0.0) <= 0);
                if (iterations >= PULSAR_SPIN_ITERATIONS) THROW_ERROR(ERROR::TOO_MANY_PULSAR_SPIN_ITERATIONS);
                pSpin = PPOW(log10pSpin, 10.0); 
                } break;

            default:                                                                                                // unknown prescription
                // the only way this can happen is if someone added a PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_MSP_BIRTH_SPIN_PERIOD_DISTRIBUTION);                                  // throw error
        }
    }
    else {
        switch (OPTIONS->PulsarBirthSpinPeriodDistribution()) {                                                     // which distribution?

            case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM: {                                                  // UNIFORM distribution between minimum and maximum value as in Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (default Pmin = and Pmax = )
                                                                                                                    // and also Kiel et al 2008 https://arxiv.org/abs/0805.0059 (default Pmin = 10 ms and Pmax 100 ms, section 3.4)
                double maximum = OPTIONS->PulsarBirthSpinPeriodDistributionMax();
                double minimum = OPTIONS->PulsarBirthSpinPeriodDistributionMin();

                pSpin = minimum + (RAND->Random() * (maximum - minimum));
                } break;

            case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL: {                                                   // NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

                double mean  = OPTIONS->PulsarBirthSpinPeriodDistributionMean();
                double sigma = OPTIONS->PulsarBirthSpinPeriodDistributionSigma();

                // this should terminate naturally, but just in case we add a guard
                std::size_t iterations = 0;
                do { pSpin = RAND->RandomGaussian(sigma) + mean;} while (iterations++ < PULSAR_SPIN_ITERATIONS && utils::Compare(pSpin, 0.0) <= 0);
                if (iterations >= PULSAR_SPIN_ITERATIONS) THROW_ERROR(ERROR::TOO_MANY_PULSAR_SPIN_ITERATIONS);

                } break;

            case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::LOGNORMAL: {                                                   // NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

                double mean  = log10(OPTIONS->PulsarBirthSpinPeriodDistributionMean());
                double sigma = log10(OPTIONS->PulsarBirthSpinPeriodDistributionSigma());
                double log10pSpin;
                // this should terminate naturally, but just in case we add a guard
                std::size_t iterations = 0;
                do { log10pSpin = RAND->RandomGaussian(sigma) + mean;} while (iterations++ < PULSAR_SPIN_ITERATIONS && utils::Compare(log10pSpin, 0.0) <= 0);
                if (iterations >= PULSAR_SPIN_ITERATIONS) THROW_ERROR(ERROR::TOO_MANY_PULSAR_SPIN_ITERATIONS);
                pSpin = PPOW(log10pSpin, 10.0); 
                } break;

            default:                                                                                                // unknown prescription
                // the only way this can happen is if someone added a PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION);                                  // throw error
        }
    }
    return pSpin * SECONDS_IN_MS;
}


/*
 * Calculate (log10 of) the magnetic field (in G) for a Pulsar at birth
 * according to selected distribution (by commandline option)
 *
 *
 * double CalculateBirthMagneticField()
 *
 * @return                                      log10 of the birth magnetic field in G
 */
double NS::CalculateBirthMagneticField() {

	double log10B;
    if ((OPTIONS->MSPsFromAIC())&(ExperiencedAIC())){

        switch (OPTIONS->MSPBirthMagneticFieldDistribution()) {                                                  // which distribution?

            case MSP_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG: {                                             // FLAT IN LOG distribution from Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (log10B0min = , log10B0max = )

                double maximum = OPTIONS->MSPBirthMagneticFieldDistributionMax();
                double minimum = OPTIONS->MSPBirthMagneticFieldDistributionMin();
                
                log10B = minimum + (RAND->Random() * (maximum - minimum));

                } break;

            case MSP_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM: {                                               // UNIFORM flat distribution used in Kiel et al 2008 https://arxiv.org/abs/0805.0059 (log10B0min = 11, log10B0max = 13.5 see section 3.4 and Table 1.)
                
                double maximum = PPOW(10.0, OPTIONS->MSPBirthMagneticFieldDistributionMax());
                double minimum = PPOW(10.0, OPTIONS->MSPBirthMagneticFieldDistributionMin());

                log10B = log10(minimum + (RAND->Random() * (maximum - minimum)));
                } break;

            case MSP_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL: {                                             // LOG NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

                double mean  = OPTIONS->MSPBirthMagneticFieldDistributionMean();
                double sigma = OPTIONS->MSPBirthMagneticFieldDistributionSigma();

                log10B = RAND->RandomGaussian(sigma) + mean;

                // add a guard to make sure magnetic field is always larger than the value set by --msp-minimum-magnetic-field
                std::size_t iterations = 0;
                do { log10B = RAND->RandomGaussian(sigma) + mean;} while (iterations++ < PULSAR_MAG_ITERATIONS && utils::Compare(log10B, log10(NS::NS_MAG_FIELD_LOWER_LIMIT)) <= 0);
                if (iterations >= PULSAR_MAG_ITERATIONS) THROW_ERROR(ERROR::TOO_MANY_PULSAR_MAG_ITERATIONS);
                } break;

            default:                                                                                                // unknown prescription
                // the only way this can happen is if someone added a MSP_BIRTH_MAGNETIC_FIELD_DISTRIBUTION
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_MSP_BIRTH_MAGNETIC_FIELD_DISTRIBUTION);                               // throw error
        }

    }
    else {

        switch (OPTIONS->PulsarBirthMagneticFieldDistribution()) {                                                  // which distribution?

            case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG: {                                             // FLAT IN LOG distribution from Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (log10B0min = , log10B0max = )

                double maximum = OPTIONS->PulsarBirthMagneticFieldDistributionMax();
                double minimum = OPTIONS->PulsarBirthMagneticFieldDistributionMin();
                
                log10B = minimum + (RAND->Random() * (maximum - minimum));

                } break;

            case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM: {                                               // UNIFORM flat distribution used in Kiel et al 2008 https://arxiv.org/abs/0805.0059 (log10B0min = 11, log10B0max = 13.5 see section 3.4 and Table 1.)
                
                double maximum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMax());
                double minimum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMin());

                log10B = log10(minimum + (RAND->Random() * (maximum - minimum)));
                } break;

            case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL: {                                             // LOG NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

                double mean  = OPTIONS->PulsarBirthMagneticFieldDistributionMean();
                double sigma = OPTIONS->PulsarBirthMagneticFieldDistributionSigma();

                log10B = RAND->RandomGaussian(sigma) + mean;

                // add a guard to make sure magnetic field is always larger than the value set by --pulsar-minimum-magnetic-field
                std::size_t iterations = 0;
                do { log10B = RAND->RandomGaussian(sigma) + mean;} while (iterations++ < PULSAR_MAG_ITERATIONS && utils::Compare(log10B, log10(NS::NS_MAG_FIELD_LOWER_LIMIT)) <= 0);
                if (iterations >= PULSAR_MAG_ITERATIONS) THROW_ERROR(ERROR::TOO_MANY_PULSAR_MAG_ITERATIONS);
                } break;

            default:                                                                                                // unknown prescription
                // the only way this can happen is if someone added a PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION);                               // throw error
        }
    }
    return log10B;
}

/*
 * Calculate the magnetic field decay due to accretion using user selected model
 * 
 * double CalculateMagneticFieldDecayAccretion_Static(const double p_initialMagField, const double p_MassGain)
 * 
 * @param   [IN]    p_initialMagField           Initial magnetic field of neutron star
 * @param   [IN]    p_MassGain                  Amount of mass accreted by neutron star
 * @return                                      Final magnetic field of neutron star
 * 
 */
double NS::CalculateMagneticFieldDecayAccretion_Static(const double p_initialMagField, const double p_MassGain){

    double magneticFieldNew = 0.0;

    switch (OPTIONS->PulsarMagneticFieldDecayAccretionModel()) {                                    // Which model to use for magnetic field decay due to accretion

        case PULSAR_MAGNETIC_FIELD_DECAY_ACCRETION_MODEL::NONE: {
            magneticFieldNew = p_initialMagField;
        } break;

        case PULSAR_MAGNETIC_FIELD_DECAY_ACCRETION_MODEL::EXPONENTIAL: {
            magneticFieldNew = CalculateMagneticFieldDecayAccretionExponential_Static(p_initialMagField, p_MassGain);
        } break;

        case PULSAR_MAGNETIC_FIELD_DECAY_ACCRETION_MODEL::SHIBAZAKI: {
            magneticFieldNew = CalculateMagneticFieldDecayAccretionShibazaki_Static(p_initialMagField, p_MassGain);
        } break;

        default:
            THROW_ERROR_STATIC(ERROR::UNKNOWN_PULSAR_MAGNETIC_FIELD_DECAY_ACCRETION_MODEL);               // Throw error
        break;
    }

    return magneticFieldNew;
}

/* 
 * Calculate the magnetic field decay due to accretion using the exponential model
 * 
 * See e.g., Equation 7 in Osłowski et al. 2011 (https://doi.org/10.1111/j.1365-2966.2010.18147.x)
 * 
 * double CalculateMagneticFieldDecayAccretionExponential_Static(const double p_initialMagField, const double p_MassGain)
 * 
 * @param   [IN]    p_initialMagField           Initial magnetic field of neutron star
 * @param   [IN]    p_MassGain                  Amount of mass accreted by neutron star
 * @return                                      Final magnetic field of neutron star
 */
double NS::CalculateMagneticFieldDecayAccretionExponential_Static(const double p_initialMagField, const double p_MassGain){

    double magneticFieldNew = (p_initialMagField - NS::NS_MAG_FIELD_LOWER_LIMIT) * std::exp(-p_MassGain / NS::NS_DECAY_MASS_SCALE) + NS::NS_MAG_FIELD_LOWER_LIMIT;

    return magneticFieldNew;
}

/* 
 * Calculate the magnetic field decay due to accretion using the model from Shibazaki et al. 1989
 * 
 * See Equation 1 in Shibazaki et al. 1989 (https://ui.adsabs.harvard.edu/abs/1989Natur.342..656S/abstract)
 * 
 * We also impose a minimum magnetic field strength here.
 * 
 * double CalculateMagneticFieldDecayAccretionShibazaki_Static(const double p_initialMagField, const double p_MassGain)
 * 
 * @param   [IN]    p_initialMagField           Initial magnetic field of neutron star
 * @param   [IN]    p_MassGain                  Amount of mass accreted by neutron star
 * @return                                      Final magnetic field of neutron star
 */
double NS::CalculateMagneticFieldDecayAccretionShibazaki_Static(const double p_initialMagField, const double p_MassGain){

    double magneticFieldNew = NS::NS_MAG_FIELD_LOWER_LIMIT + (p_initialMagField - NS::NS_MAG_FIELD_LOWER_LIMIT) / (1.0 + p_MassGain/NS::NS_DECAY_MASS_SCALE);

    return magneticFieldNew;
}

/*
 * Calculate the moment of inertia for a Neutron Star using a model independent relation between
 * the moment of inertia, mass and radius of a neutron star - return MoI in CGS.
 * 
 * Uses m_Mass and m_Radius to calculate moment of inertia.
 *
 * Raithel et al. 2016, eq 8 in  https://arxiv.org/abs/1603.06594
 *
 *
 * double CalculateMomentOfInertiaCGS_Static(const double p_Mass, const double p_Radius)
 *
 * @param   [IN]    p_Mass                      Mass of the Neutron Star (g)
 * @param   [IN]    p_Radius                    Radius of the Neutron Star in (cm)
 * @return                                      Moment of inertia (g cm^2)
 */
double NS::CalculateMomentOfInertiaCGS_Static(const double p_Mass, const double p_Radius) {
    double m_r = (p_Mass / MSOL_TO_G) / (p_Radius / KM_TO_CM);
    return 0.237 * p_Mass * p_Radius * p_Radius * (1.0 + (4.2 * m_r) + 90.0 * m_r * m_r * m_r * m_r);
}


/*
 * Calculate the spin down rate for isolated Neutron Stars in cgs
 *
 * See Equation 5 in https://arxiv.org/abs/2406.11428
 * Note that magnetic and rotational axes are orthogonal, leading to sin^2(alpha) = 1 in this equation. 
 * A model with evolving alpha will be implemented in a future version. 
 * 
 * Calculates spindown with P and Pdot, then converts to OmegaDot for recording in the log file.
 * Evolution of the inclination between pulsar magnetic and rotational axes will be considered in a future version. 
 *
 * double CalculateSpinDownRate(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius)
 *
 * @param   [IN]    p_Period                    Pulsar spin period (s). 
 * @param   [IN]    p_MomentOfInteria           Moment of Inertia of the Neutron Star (g cm^2)
 * @param   [IN]    p_MagField                  Magnetic field (Gauss )
 * @param   [IN]    p_Radius                    Radius of the Neutron Star (kilometres)
 * @return                                      Spin down rate (spin period derivative) of an isolated Neutron Star (s^-2)
 */
double NS::CalculateSpinDownRate(const double p_Period, const double p_MomentOfInteria, const double p_MagField, const double p_Radius) const {

   // pow() is slow - use multiplication

   double cgsRadius         = p_Radius * KM_TO_CM;                                                              // radius in cm
   double radius_6          = cgsRadius * cgsRadius * cgsRadius * cgsRadius * cgsRadius * cgsRadius;
   
   double magField_2        = p_MagField * p_MagField;
   constexpr double _8_PI_2 = 8.0 * PI_2;
   constexpr double _3_C_3  = 3.0E6 * C * C * C;                                                                // 3.0 * (C * 100.0) * (C * 100.0) * (C * 100.0)
   double pDotTop           = _8_PI_2 * radius_6 * magField_2;
   double pDotBottom        = _3_C_3 * p_MomentOfInteria * p_Period;
   double pDot              = pDotTop / pDotBottom;                                                             // period derivative 
   
   return pDot;                                                                           
}


/*
 * Calculate and set pulsar parameters at birth of pulsar.
 * Users should note that when choosing the NOSPIN option, 
 * pulsar spin frequency is set to 0 and spin period is infinity. 
 *
 * Modifies the following class member variables:
 *
 *    m_AngularMomentum_CGS
 *    m_MomentOfInertia_CGS
 *    m_PulsarDetails.birthPeriod
 *    m_PulsarDetails.birthSpinDownRate
 *    m_PulsarDetails.magneticField
 *    m_PulsarDetails.spinDownRate
 *    m_PulsarDetails.spinFrequency
 * 
 * 
 * void CalculateAndSetPulsarParameters()
 */
void NS::CalculateAndSetPulsarParameters() {

    m_PulsarDetails.magneticField     = PPOW(10.0, CalculateBirthMagneticField());                                                  // magnetic field in Gauss 
    m_PulsarDetails.spinPeriod        = CalculateBirthSpinPeriod();                                                                 // spin period in s
    m_MomentOfInertia_CGS             = CalculateMomentOfInertiaCGS();                                                              // MoI in CGS g cm^2
	
    m_PulsarDetails.spinFrequency     = _2_PI / m_PulsarDetails.spinPeriod;                              
    m_PulsarDetails.birthPeriod       = m_PulsarDetails.spinPeriod;                                         

    m_PulsarDetails.spinDownRate      = CalculateSpinDownRate(m_PulsarDetails.spinPeriod, m_MomentOfInertia_CGS, m_PulsarDetails.magneticField, m_Radius * RSOL_TO_KM);  
    m_PulsarDetails.birthSpinDownRate = m_PulsarDetails.spinDownRate; 
    m_AngularMomentum_CGS             = m_MomentOfInertia_CGS * m_PulsarDetails.spinFrequency;                              // in CGS g cm^2 s^-1
}

/*
 * Calculate the NS magnetic field decay timescale (in Myr)
 * 
 * taud = tauconst * (Bref/B)**alpha
 * where tauconst is taud(B=Bref). We assume Bref = 1E11 G.
 * 
 * From Dall'Osso et al. 2012 (https://academic.oup.com/mnras/article/422/4/2878/1048228) 
 * following Colpi et al. 2000 (https://iopscience.iop.org/article/10.1086/312448)
 * 
 * double CalculateMagneticFieldDecayTimescale
 * 
 * @return                                      Magnetic field decay timescale for an isolated neutron star in Myr
 *  
 * */
double NS::CalculateMagneticFieldDecayTimescale(){

    std::cout << "CalculateMagneticFieldDecayTimescale" << std::endl;

    double taud                 = 0.0;                                                          // Initialise variable to hold magnetic field decay timescale
    double Bref                 = 1E11;                                                         // Reference magnetic field (in G) at which OPTIONS->PulsarMagneticFieldDecayTimescale is defined
    double initialMagField_G    = m_PulsarDetails.magneticField * TESLA_TO_GAUSS;               // Convert T to G 

    std::cout << "Bref = " << Bref << std::endl;
    std::cout << "initialMagField_G = " << initialMagField_G << std::endl;

    if (OPTIONS->PulsarMagneticFieldDecayTimescalePower() == 0.0){                              // No scaling with magnetic field
        std::cout << "alpha = 0" << std::endl;
        taud = OPTIONS->PulsarMagneticFieldDecayTimescale();                                    // Decay timescale is just a constant
    }
    else{
        std::cout << "alpha != 0" << std::endl;
        std::cout << "B = " << m_PulsarDetails.magneticField << std::endl;
        taud = OPTIONS->PulsarMagneticFieldDecayTimescale() * PPOW(Bref/initialMagField_G, OPTIONS->PulsarMagneticFieldDecayTimescalePower());
    }
    
    std::cout << "tauconst = " << OPTIONS->PulsarMagneticFieldDecayTimescale() << std::endl;
    std::cout << "taud = " << taud << std::endl;

    return taud;
}

/*
 * Calculate the magnetic field strength as a function of time due to magnetic field decay
 * 
 * double CalculateMagneticFieldStrengthOnPhase(const double p_Time, const double p_initialMagField)
 * 
 * @param       [IN]    p_Time                  Time in seconds
 * @param       [IN]    p_initialMagField       Initial magnetic field strength in G
 * @return              Magnetic field strength timescale for an isolated neutron star in Myr
 * 
 */
double NS::CalculateMagneticFieldStrengthOnPhase(const double p_Time, const double p_initialMagField){

    std::cout << "CalculateMagneticFieldStrengthOnPhase" << std::endl;
    std::cout << "p_Time " << p_Time << std::endl;

    double magneticFieldStrength = 0.0;

    double magFieldLowerLimit    = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField()) * GAUSS_TO_TESLA; 
    double tau                   = CalculateMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR; 
    const double alpha           = OPTIONS->PulsarMagneticFieldDecayTimescalePower();                                

    if (alpha == 0.0){              // see Equation 6 in  arXiv:0903.3538v2    
        std::cout << "alpha == 0" << std::endl;
        magneticFieldStrength    = magFieldLowerLimit + (p_initialMagField - magFieldLowerLimit) * exp(-p_Time / tau);   // update pulsar magnetic field in SI. 
    }
    else{                   // Equation 8 in Dall'Osso et al. 2012 (https://ui.adsabs.harvard.edu/abs/2012MNRAS.422.2878D/abstract) but with a minimum magnetic field
        std::cout << "alpha != 0" << std::endl;
        magneticFieldStrength    = magFieldLowerLimit + (p_initialMagField - magFieldLowerLimit) * PPOW(1.0 + alpha*(p_Time/tau), -1.0/alpha);
        std::cout << "p_initialMagField, magneticFieldStrength = " << p_initialMagField << " " << magneticFieldStrength << std::endl;
    }

    return magneticFieldStrength;
}

/*
 * Update the magnetic field and spins of isolated pulsar
 * Users should note that when pulsar is not spinning, 
 * this function exits without changing anything. 
 * Note that we assume the rotational and magnetic axis 
 * are orthogonal and are not evolved in the current model.
 * A model with evolving alpha will be implemented in a future version. 
 *
 * Modifies the following class member variables:
 *
 *    m_AngularMomentum_CGS
 *    m_PulsarDetails.spinFrequency
 *    m_PulsarDetails.spinPeriod
 *    m_PulsarDetails.magneticField
 *    m_PulsarDetails.spinDownRate
 *
 *
 * void SpinDownIsolatedPulsar(const double p_Stepsize)
 *
 * @param   [IN]    p_Stepsize                  Timestep size for integration (in seconds)
 * @param   [IN]    p_RecycledNS                Boolean flag indicating whether this star is/was a recycled NS
 */
void NS::SpinDownIsolatedPulsar(const double p_Stepsize, const bool p_RecycledNS) {

    // Initialise variables, define constants
    double Psquared          = 0.0;

    double radius_IN_CM      = m_Radius * RSOL_TO_KM * KM_TO_CM;
    double radius_3          = radius_IN_CM * radius_IN_CM * radius_IN_CM;
    double radius_6          = radius_3 * radius_3;
    constexpr double _8_PI_2 = 8.0 * PI_2;
    constexpr double _3_C_3  = 3.0E6 * C * C * C;                                                                                           // 3.0 * (C * 100.0) * (C * 100.0) * (C * 100.0)
    double constant          = (_8_PI_2 * radius_6) / (_3_C_3 * m_MomentOfInertia_CGS);

    double initialMagField   = m_PulsarDetails.magneticField;     
    double initialSpinPeriod = m_PulsarDetails.spinPeriod;
 
    bool magFieldAboveMin = utils::Compare(initialMagField, NS::NS_MAG_FIELD_LOWER_LIMIT) > 0;

    std::cout << "initialMagField = " << initialMagField << std::endl;
    std::cout << "magFieldAboveMin = " << magFieldAboveMin << std::endl;

    // If the NS is not recycled, calculate magnetic field decay
    if (!p_RecycledNS && magFieldAboveMin){
        // calculate the decay of magnetic field for an isolated neutron star
        // see Equation 6 in  arXiv:0903.3538v2
        m_PulsarDetails.magneticField = NS::NS_MAG_FIELD_LOWER_LIMIT + (initialMagField - NS::NS_MAG_FIELD_LOWER_LIMIT) * std::exp(-p_Stepsize / NS::NS_DECAY_TIME_SCALE); // update pulsar magnetic field in cgs. 
        
        // calculate the spin period for an isolated neutron star 
        // with a decaying magnetic field
        // see Equation 3 in arxiv:2406.11428
        // Note that magnetic and rotational axes are orthogonal, leading to sin^2(alpha) = 1 in this equation. 
        // The rest of the calculations are carried out in cgs.   
        double term1           = NS::NS_MAG_FIELD_LOWER_LIMIT * NS::NS_MAG_FIELD_LOWER_LIMIT * p_Stepsize;
        double term2           = NS::NS_DECAY_TIME_SCALE * NS::NS_MAG_FIELD_LOWER_LIMIT * ( m_PulsarDetails.magneticField  - initialMagField);
        double term3           = (NS::NS_DECAY_TIME_SCALE / 2.0) * ((m_PulsarDetails.magneticField * m_PulsarDetails.magneticField) - (initialMagField * initialMagField));
        Psquared               = 2.0 * constant * (term1 - term2 - term3) + (initialSpinPeriod * initialSpinPeriod);
    
        std::cout << "term 1 = " << term1 << std::endl;
        std::cout << "term 2 = " << term2 << std::endl;
        std::cout << "term 3 = " << term3 << std::endl;
    
    }
    else{
        // Assume constant magnetic field strength
        m_PulsarDetails.magneticField = initialMagField;

        // calculate the final spin period for an isolated neutron star with a constant magnetic field
        // See Simon's note about pulsar spindown 
        Psquared = (initialSpinPeriod * initialSpinPeriod) + 2.0 * constant * initialMagField * initialMagField * p_Stepsize;
    }

    m_PulsarDetails.spinPeriod    = std::sqrt(Psquared);
    m_PulsarDetails.spinFrequency = _2_PI / m_PulsarDetails.spinPeriod;                                                                     // pulsar spin frequency

    std::cout << "initial spin period = " << initialSpinPeriod << std::endl;
    std::cout << "final spin period = " << m_PulsarDetails.spinPeriod << std::endl;

    m_PulsarDetails.spinDownRate  = CalculateSpinDownRate(m_PulsarDetails.spinPeriod, m_MomentOfInertia_CGS, m_PulsarDetails.magneticField, m_Radius * RSOL_TO_KM) ; 

    std::cout << "spinDownRate = " << m_PulsarDetails.spinDownRate << std::endl << std::endl;

    m_AngularMomentum_CGS         = m_PulsarDetails.spinFrequency * m_MomentOfInertia_CGS;                                                  // angular momentum of star in CGS
}


/*
 * Calculate the change in angular momentum wrt mass (dJ/dM) of a neutron star when it accretes mass through RLOF
 * Note: this function uses CGS units and requires parameters in CGS units where applicable
 * 
 * See sec. 2.2.1 in arxiv:1912.02415
 * 
 *
 * double DeltaJByAccretion_Static(const double p_Mass, const double p_Radius_6, const double p_MagField, const double p_SpinFrequency, const double p_mDot, const double p_Epsilon)
 * 
 * @param   [IN]    p_Mass                      Initial mass of the NS (g)
 * @param   [IN]    p_Radius_6                  (Radius of the NS (cm))^6 (for performance - so it isn't calculated at every integration step)
 * @param   [IN]    p_MagField                  NS magnetic field strength at the beginning of accretion (Gauss)
 * @param   [IN]    p_SpinFrequency             Angular frequency for the NS at the beginning of accretion (rad/s)
 * @param   [IN]    p_mDot                      Mass transfer rate (g s^-1)
 * @param   [IN]    p_Epsilon                   Efficiency factor allowing for uncertainties of coupling magnetic field and matter
 * @return                                      Change in angular momentum wrt mass (dJ/dM) of NS due to accretion
 */
double NS::DeltaJByAccretion_Static(const double p_Mass, const double p_Radius_6, const double p_MagField, const double p_SpinFrequency, const double p_mDot, const double p_Epsilon)  {

    // calculate the Alfven radius for an accreting neutron star
    // see eq 10 in arxiv:1912.02415 
    double spinPeriod     = _2_PI / p_SpinFrequency; 
    double p              = p_Radius_6 * p_Radius_6 / (p_mDot * p_mDot * p_Mass);
    double q              = PPOW(p, 1.0 / 7.0);
    double magneticRadius = ALFVEN_CONST * q * PPOW(p_MagField, 4.0 / 7.0) / 2.0;               // Alfven radius / 2.0 (cm)

    // calculate the difference in the keplerian angular velocity at the magnetic radius and surface angular velocity of the NS
    // see eq 2 in 1994MNRAS.269..455J / eq 9 in arxiv:1912.02415 
    double omegaK = std::sqrt(G_CGS * p_Mass / magneticRadius) / magneticRadius;                // rad/s
    double vDiff  = omegaK - p_SpinFrequency;                                                   // rad/s
    
    return p_Epsilon * vDiff * magneticRadius * magneticRadius;                                 // eq 12 in arXiv:0805.0059 / eq 8 in arxiv:1912.02415
}


/*
 * Update the magnetic field and spins of neutron stars in the following situations:
 * Note 1 : this function uses CGS units and requires parameters in CGS units where applicable
 * Note 2 : Users should note that when pulsar is not spinning (frequency or magnetic field == 0), 
 *          pulsar spin period is set to infinity and all other parameters to 0.
 * 
 * 1).  Isolated or post mass-transfer spin-down of neutron star
 * 2).  Neutron star in interacting binary system experiencing mass-transfer induced spin change for 
 *      2.1).  Roche Lobe overflow, and 
 *      2.2).  Common Envelope 
 * 
 * Modifies the following class member variables:
 *
 *    m_AngularMomentum_CGS
 *    m_MomentOfInertia_CGS            
 *    m_PulsarDetails.spinFrequency
 *    m_PulsarDetails.spinPeriod
 *    m_PulsarDetails.magneticField
 *    m_PulsarDetails.spinDownRate
 *
 *
 * void UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope, 
 *                                 const bool   p_RecycledNS, 
 *                                 const double p_Stepsize, 
 *                                 const double p_MassGain, 
 *                                 const double p_Epsilon)
 *
 * @param   [IN]    p_CommonEnvelope            Boolean flag indicating whether there there is a common envelope
 * @param   [IN]    p_RecycledNS                Boolean flag indicating whether this star is/was a recycled NS
 * @param   [IN]    p_Stepsize                  Timestep size for integration (seconds)
 * @param   [IN]    p_MassGain                  Required accretor mass gain for timestep (g)
 * @param   [IN]    p_Epsilon                   Uncertainty due to mass loss
 */
void NS::UpdateMagneticFieldAndSpin(const bool p_CommonEnvelope, const bool p_RecycledNS, double p_Stepsize, double p_MassGain, const double p_Epsilon) {

    if ((!p_RecycledNS && !p_CommonEnvelope) || (!p_RecycledNS && utils::Compare(p_MassGain, 0.0) == 0 )) {                                 // 'classical' isolated pulsars
        std::cout << "Isolated spindown 1" << std::endl;
        SpinDownIsolatedPulsar(p_Stepsize, p_RecycledNS);                                                                                                 // spin down
    }
    else if (p_CommonEnvelope && (OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::SURFACE)) {                                    // mass transfer through CE when accretion happens at the surface of the NS

        double massG           = m_Mass * MSOL_TO_G;                                                                                        // mass in g
        double radiusCM        = m_Radius * RSOL_TO_CM;                                                                                     // radius in cm
        m_MomentOfInertia_CGS  = CalculateMomentOfInertiaCGS(); 
        double jAcc            = std::sqrt(G_CGS * massG * radiusCM) * p_MassGain;    
        m_AngularMomentum_CGS += jAcc;                                                                                                      // angular momentum of the accreted material as it falls onto the surface of the NS
        if (utils::Compare(m_PulsarDetails.magneticField, NS::NS_MAG_FIELD_LOWER_LIMIT) < 0) {
            // if magnetic field is already lower than the lower limit, 
            // set it to the lower limit value.
            m_PulsarDetails.magneticField = NS::NS_MAG_FIELD_LOWER_LIMIT;
        }
        else {
            //m_PulsarDetails.magneticField = (m_PulsarDetails.magneticField - NS::NS_MAG_FIELD_LOWER_LIMIT) * std::exp(-p_MassGain / G_TO_KG / NS::NS_DECAY_MASS_SCALE) + NS::NS_MAG_FIELD_LOWER_LIMIT; // eq. 12 in arxiv:1912.02415 
            m_PulsarDetails.magneticField = CalculateMagneticFieldDecayAccretion_Static(m_PulsarDetails.magneticField, p_MassGain);
        }
        
        double previousSpinFrequency  = m_PulsarDetails.spinFrequency;
        m_PulsarDetails.spinFrequency = m_AngularMomentum_CGS / m_MomentOfInertia_CGS;
        m_PulsarDetails.spinPeriod    = _2_PI / m_PulsarDetails.spinFrequency;
        double fDot                   = (m_PulsarDetails.spinFrequency - previousSpinFrequency) / p_Stepsize;
        
        m_PulsarDetails.spinDownRate  = -fDot * m_PulsarDetails.spinPeriod * m_PulsarDetails.spinPeriod / _2_PI;
    } 
    else if (utils::Compare(p_MassGain, 0.0) > 0                                                                          &&
             utils::Compare(m_PulsarDetails.spinPeriod, 0.001) > 0                                                        && 
             (!p_CommonEnvelope || (p_CommonEnvelope && OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::DISK))) {                // pulsar recycling through accretion

        // recycling happens for pulsars with spin period larger than 1 ms and in a binary system with mass transfer
        // the pulsar being recycled is either in a common envelope, or should have started the recycling process in previous time steps
       
        // solve for the angular momentum of the NS after accretion
        // the accretor will gain p_MassGain g over p_Stepsize seconds

        double initialAngularMomentum_CGS = m_AngularMomentum_CGS;                                                                          // initial angular momentum
        double initialMagField            = m_PulsarDetails.magneticField;                                                                  // initial magnetic field
        double initialMass                = m_Mass * MSOL_TO_G;                                                                             // initial mass of NS in g

        double radius                     = m_Radius * RSOL_TO_CM;                                                                          // radius of NS in cm
        double radius_6                   = radius * radius * radius * radius * radius * radius;                                            // for performance - do it once
        double mDot                       = p_MassGain / p_Stepsize;                                                                        // required mass transfer rate (g s^-1)
        double massFinal                  = initialMass + p_MassGain;                                                                       // required final mass of NS (after accretion) in g

        // Calculate the co-rotation radius based on Eq. 10 in Li et al. 2021
        m_PulsarDetails.coRotationRadius  = 1.5E8 * PPOW(m_Mass * m_PulsarDetails.spinPeriod * m_PulsarDetails.spinPeriod , 1.0/3.0) ;

        // calculate the Magnetic/Alfven radius for an accreting neutron star
        // see eq 10 in arxiv:1912.02415 
        double p              = radius_6 * radius_6 / (mDot * mDot * initialMass);
        double q              = PPOW(p, 1.0 / 7.0);
        m_PulsarDetails.magneticRadius = ALFVEN_CONST * q * PPOW(initialMagField, 4.0 / 7.0) / 2.0;   
        // double vDiff3         = std::sqrt(G_CGS * initialMass / m_PulsarDetails.magneticRadius) / m_PulsarDetails.magneticRadius - m_PulsarDetails.spinFrequency;
        
        // std::cout <<m_PulsarDetails.coRotationRadius << " " << m_PulsarDetails.magneticRadius<< " " << vDiff3 << " " << m_Mass << " " << m_PulsarDetails.spinPeriod << std::endl;
        if (utils::Compare(m_PulsarDetails.magneticRadius, m_PulsarDetails.coRotationRadius) > 0 ) {m_PulsarDetails.propellerMode = true ;}
        else {m_PulsarDetails.propellerMode = false ;}
        // calculate initial mass slice size for integration
        double jAcc      = DeltaJByAccretion_Static(initialMass, radius_6, m_PulsarDetails.magneticField, m_PulsarDetails.spinFrequency, mDot, p_Epsilon);
        double massSlice = std::fabs(m_AngularMomentum_CGS / 1000.0 / jAcc);                                                                // abs(Jx10^-3 / dJdM)

        // use the boost ODE solver for speed and accuracy

        // define the ODE
        // initial state
        state_type x(1);
        x[0] = initialAngularMomentum_CGS;                                                                                                  // angular momentum

        // ODE
        struct ode {
            double p_Mass, p_Radius, p_Radius_6, p_MagField, p_Mdot, p_Epsilon, p_Propeller;
            ode(double mass, double radius, double radius_6, double magField, double mdot, double epsilon, bool propeller) :
                p_Mass(mass), p_Radius(radius), p_Radius_6(radius_6), p_MagField(magField), p_Mdot(mdot), p_Epsilon(epsilon), p_Propeller(propeller) { }

            // x is the current state of the ODE (x[0] = angular momentum J)
            // dxdm is the change of state wrt mass (dxdm[0] = dJdm)
            // p_MassDelta is the cumulative change in mass of the NS
            void operator () (const state_type& x, state_type& dxdm, double p_MassDelta ) const {
                double m = p_Mass + p_MassDelta;
                // double B = (p_MagField - NS::NS_MAG_FIELD_LOWER_LIMIT) * std::exp(-p_MassDelta / NS::NS_DECAY_MASS_SCALE) + NS::NS_MAG_FIELD_LOWER_LIMIT;
                double B = CalculateMagneticFieldDecayAccretion_Static(p_MagField, p_MassDelta);
                double f = x[0] / CalculateMomentOfInertiaCGS_Static(m, p_Radius);
                dxdm[0]  = DeltaJByAccretion_Static(m, p_Radius_6, B, f, p_Mdot, p_Epsilon);   
                // if (p_Propeller) {
                //     dxdm[0]  = -1.0 * abs(DeltaJByAccretion_Static(m, p_Radius_6, B, f, p_Mdot, p_Epsilon));         
                // }
                // else {
                //     dxdm[0]  = abs(DeltaJByAccretion_Static(m, p_Radius_6, B, f, p_Mdot, p_Epsilon));         
                // }                         
            }
        };

        // integrate
        controlled_stepper_type stepper;
        (void)integrate_adaptive(stepper, ode{ initialMass, radius, radius_6, initialMagField, mDot, p_Epsilon, m_PulsarDetails.propellerMode }, x, 0.0, p_MassGain, massSlice);
                
        // final values
        m_AngularMomentum_CGS = x[0];
        m_MomentOfInertia_CGS = CalculateMomentOfInertiaCGS_Static(massFinal, radius);
        if (utils::Compare(initialMagField, NS::NS_MAG_FIELD_LOWER_LIMIT) < 0) {
            // if magnetic field is already lower than the lower limit, 
            // set it to the lower limit value.
            m_PulsarDetails.magneticField = NS::NS_MAG_FIELD_LOWER_LIMIT;
        }
        else {
            //m_PulsarDetails.magneticField = (initialMagField - NS::NS_MAG_FIELD_LOWER_LIMIT) * std::exp(-p_MassGain / NS::NS_DECAY_MASS_SCALE) + NS::NS_MAG_FIELD_LOWER_LIMIT;
            m_PulsarDetails.magneticField = CalculateMagneticFieldDecayAccretion_Static(initialMagField, p_MassGain);
        }
        m_PulsarDetails.spinFrequency = m_AngularMomentum_CGS / m_MomentOfInertia_CGS;
        m_PulsarDetails.spinPeriod    = _2_PI / m_PulsarDetails.spinFrequency;
        double fDot                   = (m_AngularMomentum_CGS - initialAngularMomentum_CGS) / m_MomentOfInertia_CGS / p_Stepsize;          // eq. 11 in arxiv:1912.02415 
        m_PulsarDetails.spinDownRate  = -fDot * m_PulsarDetails.spinPeriod * m_PulsarDetails.spinPeriod / _2_PI;
    }      
    else {          
        std::cout << "Isolated spindown 2" << std::endl;                                                                                                                        // otherwise...    
        std::cout << "p_RecylcedNS = " << p_RecycledNS << std::endl;
        SpinDownIsolatedPulsar(p_Stepsize, p_RecycledNS);                                                                                                 // ...treat the pulsar as isolated - spin down
    }
}


/* 
 * Resolve common envelope accretion
 *
 * For stellar types other than Black hole or Neutron Star just set the star's mass to the parameter passed
 * For Black holes or Neutron Stars calculate the mass accreted during a CE
 *
 *
 * double ResolveCommonEnvelopeAccretion(const double p_FinalMass, 
 *                                       const double p_CompanionMass, 
 *                                       const double p_CompanionRadius, 
 *                                       const double p_CompanionEnvelope)
 *
 * @param   [IN]    p_FinalMass                 Mass of the accreting object post mass transfer (Msol) (not used here)
 * @param   [IN]    p_CompanionMass             Mass of the companion
 * @param   [IN]    p_CompanionRadius           Radius of the companion
 * @param   [IN]    p_CompanionEnvelope         Envelope of the companion pre-CE
 * @return                                      Mass delta                                      
 * 
 */
double NS::ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                          const double p_CompanionMass,
                                          const double p_CompanionRadius,
                                          const double p_CompanionEnvelope) {
    return CalculateMassAccretedForCO(Mass(), p_CompanionMass, p_CompanionRadius, p_CompanionEnvelope);
}

