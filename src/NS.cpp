#include "Rand.h"
#include "NS.h"
#include <random>


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
                pSpin = PPOW(10.0, log10pSpin); 
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
                pSpin = PPOW(10.0, log10pSpin); 
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
 * Calculate the maximum mass acceptance rate and accretion efficiency of a neutron star
 *
 * Mdot_acc = min(f * Mdot_donor, Mdot_Edd), where f is --neutron-star-accretion-efficiency-parameter
 *
 * The efficiency f describes inefficient accretion at sub-Eddington supply rates (e.g. propeller phases and
 * outflows in transient systems), so it applies to all stable mass transfer onto a NS, whether the donor
 * transfers mass on a nuclear or a thermal timescale.  At supply rates above Mdot_Edd / f the accretion rate
 * is Eddington-limited.  Implementing it here keeps the accretion efficiency consistent between the thermal
 * and nuclear timescale branches of mass transfer and the root solver used to find the donor's mass loss.
 *
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_AccretorMassRate          Thermal mass loss rate of the accretor (this star) - ignored here
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
DBL_DBL NS::CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate) {

    double efficiency       = OPTIONS->NeutronStarAccretionEfficiencyParameter();
    double acceptanceRate   = std::min(efficiency * p_DonorMassRate, CalculateEddingtonCriticalRate());
    double fractionAccreted = utils::Compare(p_DonorMassRate, 0.0) > 0 ? acceptanceRate / p_DonorMassRate : efficiency;

    return std::make_tuple(acceptanceRate, fractionAccreted);
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
double NS::CalculateMagneticFieldDecayAccretion_Static(const double p_MagField, const double p_MagFieldPreAccretion, const double p_AccretedMass, const double p_MassGain){

    double magneticFieldNew = 0.0;

    switch (OPTIONS->PulsarMagneticFieldDecayAccretionModel()) {                                    // Which model to use for magnetic field decay due to accretion

        case PULSAR_MAGNETIC_FIELD_DECAY_ACCRETION_MODEL::NONE: {
            magneticFieldNew = p_MagField;
        } break;

        case PULSAR_MAGNETIC_FIELD_DECAY_ACCRETION_MODEL::EXPONENTIAL: {                           // memoryless - only needs the current field and the mass gained this step
            magneticFieldNew = CalculateMagneticFieldDecayAccretionExponential_Static(p_MagField, p_MassGain);
        } break;

        case PULSAR_MAGNETIC_FIELD_DECAY_ACCRETION_MODEL::SHIBAZAKI: {                             // not memoryless - needs the pre-accretion field and the total accreted mass
            magneticFieldNew = CalculateMagneticFieldDecayAccretionShibazaki_Static(p_MagFieldPreAccretion, p_AccretedMass + p_MassGain);
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
 * Unlike the exponential model, this is not memoryless: B(dM1 + dM2) != B(B(dM1), dM2), so it must be
 * evaluated from the field at the onset of accretion and the total accreted mass, not applied per timestep
 * (applying it per timestep converges to exponential decay with mass scale NS_DECAY_MASS_SCALE).
 *
 * double CalculateMagneticFieldDecayAccretionShibazaki_Static(const double p_MagFieldPreAccretion, const double p_TotalAccretedMass)
 *
 * @param   [IN]    p_MagFieldPreAccretion      Magnetic field of neutron star at the onset of accretion
 * @param   [IN]    p_TotalAccretedMass         Total amount of mass accreted by neutron star since the onset of accretion
 * @return                                      Final magnetic field of neutron star
 */
double NS::CalculateMagneticFieldDecayAccretionShibazaki_Static(const double p_MagFieldPreAccretion, const double p_TotalAccretedMass){

    double x     = p_TotalAccretedMass / NS::NS_DECAY_MASS_SCALE;
    double index = OPTIONS->PulsarMagneticFieldDecayAccretionIndex();                                   // generalised: B = Bmin + (B0 - Bmin) / (1 + (dM / Mb)^index)
    if (utils::Compare(index, 1.0) != 0) x = PPOW(x, index);                                            // index 1 (default): Shibazaki et al. 1989

    double magneticFieldNew = NS::NS_MAG_FIELD_LOWER_LIMIT + (p_MagFieldPreAccretion - NS::NS_MAG_FIELD_LOWER_LIMIT) / (1.0 + x);

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
   double pDot              = pDotTop / pDotBottom;                                                             // period derivative (orthogonal vacuum dipole)
   
   return pDot * SpinDownTorqueFactor();                                                                         // scaled for the magnetic inclination (1 unless --pulsar-inclination-evolution is set)
}


/*
 * Spin-down torque relative to an orthogonal rotator
 *
 * --pulsar-inclination-evolution NONE: 1 (orthogonal vacuum dipole, the COMPAS default)
 * --pulsar-inclination-evolution MHD:  (k0 + k1 sin^2 alpha) / (k0 + k1), the force-free/MHD torque (Spitkovsky 2006; Philippov,
 *                                      Tchekhovskoy & Li 2014) normalised so that an orthogonal rotator spins down as in the default
 * --pulsar-inclination-evolution EXPONENTIAL: sin^2 alpha (--pulsar-spin-down-torque VACUUM) or the MHD factor above (MHD)
 *
 * double SpinDownTorqueFactor(const double p_Inclination)
 *
 * @param   [IN]    p_Inclination               Angle between the magnetic and rotation axes (radians)
 * @return                                      Torque relative to an orthogonal rotator with the same field and spin
 */
double NS::SpinDownTorqueFactor(const double p_Inclination) const {
    if (InclinationEvolution() == PULSAR_INCLINATION_EVOLUTION::NONE) return 1.0;
    double sinAlpha = std::sin(p_Inclination);
    if (InclinationEvolution() == PULSAR_INCLINATION_EVOLUTION::EXPONENTIAL &&
        OPTIONS->PulsarSpinDownTorque()       == PULSAR_SPIN_DOWN_TORQUE::VACUUM) return sinAlpha * sinAlpha;                  // vacuum dipole (e.g. Johnston & Karastergiou 2017)
    return (PULSAR_MHD_K0 + PULSAR_MHD_K1 * sinAlpha * sinAlpha) / (PULSAR_MHD_K0 + PULSAR_MHD_K1);
}


/*
 * Inclination evolution in force for this star: --pulsar-inclination-evolution, or after the onset of accretion
 * --pulsar-recycled-inclination-evolution (unless SAME)
 *
 * PULSAR_INCLINATION_EVOLUTION InclinationEvolution()
 */
PULSAR_INCLINATION_EVOLUTION NS::InclinationEvolution() const {
    if (utils::Compare(m_PulsarDetails.accretedMass, 0.0) > 0) {
        switch (OPTIONS->PulsarRecycledInclinationEvolution()) {
            case PULSAR_RECYCLED_INCLINATION_EVOLUTION::NONE: return PULSAR_INCLINATION_EVOLUTION::NONE;
            case PULSAR_RECYCLED_INCLINATION_EVOLUTION::MHD:  return PULSAR_INCLINATION_EVOLUTION::MHD;
            default: break;
        }
    }
    return OPTIONS->PulsarInclinationEvolution();
}


/*
 * Redraw the magnetic inclination (isotropic, cos alpha uniform) at the onset of accretion, if
 * --pulsar-reset-inclination-at-recycling is set and the inclination is evolved. Uses a generator seeded by the star
 * (not the COMPAS random number stream), so the rest of the evolution is unchanged.
 *
 * void ResetInclinationAtRecycling()
 */
void NS::ResetInclinationAtRecycling() {
    if (!OPTIONS->PulsarResetInclinationAtRecycling() || OPTIONS->PulsarInclinationEvolution() == PULSAR_INCLINATION_EVOLUTION::NONE) return;
    std::mt19937_64 generator(static_cast<std::uint64_t>(m_RandomSeed) ^ (std::hash<double>{}(m_MZAMS) + 0x51ED27E1C0FFEE01ULL));
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    m_PulsarDetails.inclination = std::acos(uniform(generator));
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
 *    m_PulsarDetails.birthMagneticField
 *    m_PulsarDetails.magneticField
 *    m_PulsarDetails.spinDownRate
 *    m_PulsarDetails.spinFrequency
 * 
 * 
 * void CalculateAndSetPulsarParameters()
 */
void NS::CalculateAndSetPulsarParameters() {

    m_PulsarDetails.magneticField     = PPOW(10.0, CalculateBirthMagneticField());                                                  // magnetic field in Gauss 
    m_PulsarDetails.birthMagneticField = m_PulsarDetails.magneticField;
    m_PulsarDetails.spinPeriod        = CalculateBirthSpinPeriod();                                                                 // spin period in s

    // birth inclination: isotropic (cos alpha uniform) for --pulsar-inclination-evolution MHD, otherwise orthogonal
    // drawn from a generator seeded by the star (not the COMPAS random number stream), so the rest of the evolution is unchanged
    if (OPTIONS->PulsarInclinationEvolution() == PULSAR_INCLINATION_EVOLUTION::NONE) {
        m_PulsarDetails.inclination = M_PI / 2.0;
    }
    else {
        std::mt19937_64 generator(static_cast<std::uint64_t>(m_RandomSeed) ^ (std::hash<double>{}(m_MZAMS) + 0x9E3779B97F4A7C15ULL));
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        m_PulsarDetails.inclination = std::acos(uniform(generator));
    }
    // field below which isolated decay stops: --pulsar-minimum-magnetic-field, or (GRABER2024) a lognormal per-star final field (Graber et al. 2024),
    // also drawn from a generator seeded by the star so that the rest of the evolution is unchanged
    m_PulsarDetails.finalMagneticField = NS::NS_MAG_FIELD_LOWER_LIMIT;
    if (OPTIONS->PulsarMagneticFieldDecayModel() == PULSAR_MAGNETIC_FIELD_DECAY_MODEL::GRABER2024) {
        std::mt19937_64 generator(static_cast<std::uint64_t>(m_RandomSeed) ^ (std::hash<double>{}(m_MZAMS) + 0xD1B54A32D192ED03ULL));
        std::normal_distribution<double> gaussian(OPTIONS->PulsarFinalMagneticFieldDistributionMean(), std::max(OPTIONS->PulsarFinalMagneticFieldDistributionSigma(), 1.0E-12));
        m_PulsarDetails.finalMagneticField = std::max(NS::NS_MAG_FIELD_LOWER_LIMIT, PPOW(10.0, gaussian(generator)));
    }
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

    // std::cout << "CalculateMagneticFieldDecayTimescale" << std::endl;

    double taud                 = 0.0;                                                          // Initialise variable to hold magnetic field decay timescale
    double Bref                 = 1E11;                                                         // Reference magnetic field (in G) at which OPTIONS->PulsarMagneticFieldDecayTimescale is defined
    double initialMagField_G    = m_PulsarDetails.magneticField;                                // m_PulsarDetails.magneticField is already in G

    // std::cout << "Bref = " << Bref << std::endl;
    // std::cout << "initialMagField_G = " << initialMagField_G << std::endl;

    if (OPTIONS->PulsarMagneticFieldDecayTimescalePower() == 0.0){                              // No scaling with magnetic field
        // std::cout << "alpha = 0" << std::endl;
        taud = OPTIONS->PulsarMagneticFieldDecayTimescale();                                    // Decay timescale is just a constant
    }
    else{
        // std::cout << "alpha != 0" << std::endl;
        // std::cout << "B = " << m_PulsarDetails.magneticField << std::endl;
        taud = OPTIONS->PulsarMagneticFieldDecayTimescale() * PPOW(Bref/initialMagField_G, OPTIONS->PulsarMagneticFieldDecayTimescalePower());
    }
    
    // std::cout << "tauconst = " << OPTIONS->PulsarMagneticFieldDecayTimescale() << std::endl;
    // std::cout << "taud = " << taud << std::endl;

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

    // std::cout << "CalculateMagneticFieldStrengthOnPhase" << std::endl;
    // std::cout << "p_Time " << p_Time << std::endl;

    double magneticFieldStrength = 0.0;

    double magFieldLowerLimit    = NS::NS_MAG_FIELD_LOWER_LIMIT;                                                       // G - same units as p_initialMagField
    double tau                   = CalculateMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR; 
    const double alpha           = OPTIONS->PulsarMagneticFieldDecayTimescalePower();                                

    if (alpha == 0.0){              // see Equation 6 in  arXiv:0903.3538v2    
        // std::cout << "alpha == 0" << std::endl;
        magneticFieldStrength    = magFieldLowerLimit + (p_initialMagField - magFieldLowerLimit) * exp(-p_Time / tau);   // update pulsar magnetic field in SI. 
    }
    else{                   // Equation 8 in Dall'Osso et al. 2012 (https://ui.adsabs.harvard.edu/abs/2012MNRAS.422.2878D/abstract) but with a minimum magnetic field
        // std::cout << "alpha != 0" << std::endl;
        magneticFieldStrength    = magFieldLowerLimit + (p_initialMagField - magFieldLowerLimit) * PPOW(1.0 + alpha*(p_Time/tau), -1.0/alpha);
        // std::cout << "p_initialMagField, magneticFieldStrength = " << p_initialMagField << " " << magneticFieldStrength << std::endl;
    }

    return magneticFieldStrength;
}


/*
 * Shape of the isolated pulsar field against age for --pulsar-magnetic-field-decay-model GRABER2024 (Graber et al. 2024, ApJ 968, 16, Appendix A):
 * broken power laws fitted to the Vigano et al. (2021) crustal magneto-thermal simulations, with timescales tau1 = A1 B0^b1 and tau2 = A2 B0^b2 set
 * by the birth field B0, followed by the late-time power law (1 + t/tau_late)^a_late:
 *
 *     B(t)/B0 = (1 + t/tau1)^a1 (1 + t/tau2)^(a2 - a1) (1 + t/tau_late)^(a_late - a2)    tau1 < tau2 < tau_late
 *             = (1 + t/tau1)^a1 (1 + t/tau_late)^(a_late - a1)                         tau1 < tau_late < tau2
 *             = (1 + t/tau_late)^a_late                                                tau_late < tau1
 *
 * The field is evolved multiplicatively, B(t + dt) = B(t) shape(t + dt) / shape(t), so that changes by accretion carry forward.
 *
 * double GraberFieldShape(const double p_Age)
 *
 * @param       [IN]    p_Age                   Age of the neutron star (s)
 * @return              B(t)/B0
 */
double NS::GraberFieldShape(const double p_Age) const {
    double tYr     = p_Age / SECONDS_IN_YEAR;
    double B0      = m_PulsarDetails.birthMagneticField;
    double tau1    = PULSAR_GRABER_A1 * PPOW(B0, PULSAR_GRABER_B1);                                                  // yr
    double tau2    = PULSAR_GRABER_A2 * PPOW(B0, PULSAR_GRABER_B2);                                                  // yr
    double tauLate = OPTIONS->PulsarLateMagneticFieldDecayTimescale() * MYR_TO_YEAR;                                  // yr
    double aLate   = OPTIONS->PulsarLateMagneticFieldDecayIndex();
    double f1      = std::log1p(tYr / tau1);
    double f2      = std::log1p(tYr / tau2);
    double fLate   = std::log1p(tYr / tauLate);
    double lnShape = (tau2 < tauLate) ? PULSAR_GRABER_ALPHA1 * f1 + (PULSAR_GRABER_ALPHA2 - PULSAR_GRABER_ALPHA1) * f2 + (aLate - PULSAR_GRABER_ALPHA2) * fLate
                   : (tau1 < tauLate) ? PULSAR_GRABER_ALPHA1 * f1 + (aLate - PULSAR_GRABER_ALPHA1) * fLate
                   : aLate * fLate;
    return std::exp(lnShape);
}


/*
 * Calculate the time integral of B^2 over a timestep, for B(t) given by CalculateMagneticFieldStrengthOnPhase()
 *
 * Used for the spin evolution of an isolated pulsar with a decaying field: P^2 = P0^2 + 2K * int_0^t B^2 dt'
 *
 * B(t) = Bmin + dB0 * f(t), with f(t) = exp(-t/tau) (alpha = 0) or f(t) = u^(-1/alpha), u = 1 + alpha t / tau (alpha != 0)
 *
 * double CalculateMagneticFieldSquaredIntegral(const double p_Time, const double p_initialMagField)
 *
 * @param       [IN]    p_Time                  Time in seconds
 * @param       [IN]    p_initialMagField       Initial magnetic field strength in G
 * @return              int_0^t B^2 dt' (G^2 s)
 */
double NS::CalculateMagneticFieldSquaredIntegral(const double p_Time, const double p_initialMagField){

    double Bmin  = NS::NS_MAG_FIELD_LOWER_LIMIT;
    double dB0   = p_initialMagField - Bmin;
    double tau   = CalculateMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR;
    double alpha = OPTIONS->PulsarMagneticFieldDecayTimescalePower();

    double result = Bmin * Bmin * p_Time;

    if (alpha == 0.0) {
        result += 2.0 * Bmin * dB0 * tau * (1.0 - std::exp(-p_Time / tau)) + 0.5 * dB0 * dB0 * tau * (1.0 - std::exp(-2.0 * p_Time / tau));
    }
    else {
        double U = 1.0 + alpha * p_Time / tau;
        auto F   = [U](const double p) { return utils::Compare(p, 1.0) == 0 ? std::log(U) : (PPOW(U, 1.0 - p) - 1.0) / (1.0 - p); };   // int_1^U u^-p du
        result  += (tau / alpha) * (2.0 * Bmin * dB0 * F(1.0 / alpha) + dB0 * dB0 * F(2.0 / alpha));
    }

    return result;
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
 
    bool   graber     = OPTIONS->PulsarMagneticFieldDecayModel() == PULSAR_MAGNETIC_FIELD_DECAY_MODEL::GRABER2024;
    double fieldFloor = graber ? m_PulsarDetails.finalMagneticField : NS::NS_MAG_FIELD_LOWER_LIMIT;                                       // G
    bool magFieldAboveMin = utils::Compare(initialMagField, fieldFloor) > 0;

    // std::cout << "initialMagField = " << initialMagField << std::endl;
    // std::cout << "magFieldAboveMin = " << magFieldAboveMin << std::endl;

    // GRABER2024: the field is a broken power law in age; tStart is the age at the start of the step (the star has already been aged)
    double tStart      = std::max(0.0, m_Age * MYR_TO_YEAR * SECONDS_IN_YEAR - p_Stepsize);                                                // s
    double shapeStart  = graber ? GraberFieldShape(tStart) : 1.0;
    auto   graberField = [&](const double p_T) { return std::max(fieldFloor, initialMagField * GraberFieldShape(tStart + p_T) / shapeStart); };  // field p_T (s) into the step
    // GRABER2024: composite Simpson over the step, uniform in ln(age + 1 yr), with PULSAR_GRABER_SUBSTEPS_PER_EFOLD substeps per e-folding of (age + 1 yr)
    // and at least PULSAR_ALIGNMENT_SUBSTEPS_PER_EFOLD per p_Timescale (s) in t (0: no extra requirement)
    auto   graberQuadrature = [&](auto &&p_G, const double p_Timescale) {
        double t0   = SECONDS_IN_YEAR;
        double u0   = std::log(tStart + t0);
        double u1   = std::log(tStart + p_Stepsize + t0);
        double nLog = PULSAR_GRABER_SUBSTEPS_PER_EFOLD * (u1 - u0);
        double nLin = p_Timescale > 0.0 ? PULSAR_ALIGNMENT_SUBSTEPS_PER_EFOLD * (u1 - u0) * (tStart + p_Stepsize + t0) / p_Timescale : 0.0;
        int    n    = 2 * std::max(1, std::min(100000, static_cast<int>(std::ceil(0.5 * std::max(nLog, nLin)))));
        double h    = (u1 - u0) / n;
        double sum  = 0.0;
        for (int i = 0; i <= n; i++) {
            double u = u0 + i * h;
            double t = std::min(std::max(std::exp(u) - t0 - tStart, 0.0), p_Stepsize);                                                       // time into the step
            sum     += ((i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0)) * p_G(t) * std::exp(u);                                         // dt = e^u du
        }
        return sum * h / 3.0;
    };

    // If the NS is not recycled, calculate magnetic field decay
    double magFieldSquaredIntegral;                                                                                                         // int B^2 dt over the step (G^2 s)
    bool   fieldDecays = !p_RecycledNS && magFieldAboveMin;
    if (fieldDecays && graber) {
        m_PulsarDetails.magneticField = graberField(p_Stepsize);
        magFieldSquaredIntegral       = graberQuadrature([&](const double p_T) { double b = graberField(p_T); return b * b; }, 0.0);
    }
    else if (fieldDecays){
        // calculate the decay of magnetic field for an isolated neutron star
        // see Equation 6 in  arXiv:0903.3538v2
        m_PulsarDetails.magneticField = CalculateMagneticFieldStrengthOnPhase(p_Stepsize, initialMagField); 
        magFieldSquaredIntegral       = CalculateMagneticFieldSquaredIntegral(p_Stepsize, initialMagField);
    }
    else{
        // Assume constant magnetic field strength
        m_PulsarDetails.magneticField = initialMagField;
        magFieldSquaredIntegral       = initialMagField * initialMagField * p_Stepsize;
    }

    if (InclinationEvolution() == PULSAR_INCLINATION_EVOLUTION::NONE) {
        // orthogonal rotator with a (possibly decaying) field: P^2 = P0^2 + 2K int B^2 dt
        // (generalises Equation 3 in arxiv:2406.11428, which assumes exponential decay with NS_DECAY_TIME_SCALE, to the power-law decay)
        Psquared = fieldDecays
                   ? 2.0 * constant * magFieldSquaredIntegral + (initialSpinPeriod * initialSpinPeriod)
                   : (initialSpinPeriod * initialSpinPeriod) + 2.0 * constant * initialMagField * initialMagField * p_Stepsize;      // (expressions as before, for identical results)
    }
    else if (InclinationEvolution() == PULSAR_INCLINATION_EVOLUTION::EXPONENTIAL) {
        // inclination decays on a fixed timescale, alpha(t) = alpha0 exp(-t/tau) (Tauris & Manchester 1998; Johnston & Karastergiou 2017),
        // independent of the spin, so P^2 is a quadrature:  P^2 = P0^2 + 2K int B(t)^2 f(alpha(t)) dt,  f = SpinDownTorqueFactor(alpha)
        // integrated with composite Simpson, with PULSAR_ALIGNMENT_SUBSTEPS_PER_EFOLD substeps per e-folding of alpha or of the field
        double tauAlign = OPTIONS->PulsarAlignmentTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR;                                              // s
        double alpha0   = m_PulsarDetails.inclination;
        if (fieldDecays && graber) {
            auto integrand = [&](const double p_T) { double b = graberField(p_T); return b * b * SpinDownTorqueFactor(alpha0 * std::exp(-p_T / tauAlign)); };
            Psquared       = (initialSpinPeriod * initialSpinPeriod) + 2.0 * constant * graberQuadrature(integrand, tauAlign);
        }
        else {
        double tauField = fieldDecays ? CalculateMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR : tauAlign;                  // s
        auto integrand  = [&](const double p_T) {
            double b = fieldDecays ? CalculateMagneticFieldStrengthOnPhase(p_T, initialMagField) : initialMagField;
            return b * b * SpinDownTorqueFactor(alpha0 * std::exp(-p_T / tauAlign));
        };
        int    nSub     = 2 * std::max(1, std::min(100000, static_cast<int>(std::ceil(0.5 * PULSAR_ALIGNMENT_SUBSTEPS_PER_EFOLD * p_Stepsize / std::min(tauAlign, tauField)))));
        double h        = p_Stepsize / nSub;
        double sum      = integrand(0.0) + integrand(p_Stepsize);
        for (int i = 1; i < nSub; i++) sum += (i % 2 == 1 ? 4.0 : 2.0) * integrand(i * h);
        Psquared                    = (initialSpinPeriod * initialSpinPeriod) + 2.0 * constant * sum * h / 3.0;
        }
        m_PulsarDetails.inclination = alpha0 * std::exp(-p_Stepsize / tauAlign);
    }
    else {
        // MHD torque with alignment; with y = P^2 and u = int B^2 dt:
        //     dy/du = 2K f(alpha),  dalpha/du = -K k2 sin(alpha) cos(alpha) / ((k0 + k1) y),  f(alpha) = (k0 + k1 sin^2 alpha) / (k0 + k1)
        // integrated with RK4, substeps chosen so that y changes by at most PULSAR_MHD_MAX_FRACTIONAL_STEP per substep
        double y     = initialSpinPeriod * initialSpinPeriod;
        double alpha = m_PulsarDetails.inclination;
        auto derivatives = [&](const double p_Y, const double p_Alpha, double &p_dY, double &p_dAlpha) {
            double sinAlpha = std::sin(p_Alpha);
            p_dY            = 2.0 * constant * (PULSAR_MHD_K0 + PULSAR_MHD_K1 * sinAlpha * sinAlpha) / (PULSAR_MHD_K0 + PULSAR_MHD_K1);
            p_dAlpha        = -constant * PULSAR_MHD_K2 * sinAlpha * std::cos(p_Alpha) / ((PULSAR_MHD_K0 + PULSAR_MHD_K1) * p_Y);
        };
        double growth = 2.0 * constant * magFieldSquaredIntegral / y;                                                                        // upper limit to the fractional change of y
        int    nSub   = std::max(1, std::min(1000000, static_cast<int>(std::ceil(growth / PULSAR_MHD_MAX_FRACTIONAL_STEP))));
        double du     = magFieldSquaredIntegral / nSub;
        for (int i = 0; i < nSub; i++) {
            double k1y, k1a, k2y, k2a, k3y, k3a, k4y, k4a;
            derivatives(y,                    alpha,                    k1y, k1a);
            derivatives(y + 0.5 * du * k1y,   alpha + 0.5 * du * k1a,   k2y, k2a);
            derivatives(y + 0.5 * du * k2y,   alpha + 0.5 * du * k2a,   k3y, k3a);
            derivatives(y + du * k3y,         alpha + du * k3a,         k4y, k4a);
            y     += du * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
            alpha += du * (k1a + 2.0 * k2a + 2.0 * k3a + k4a) / 6.0;
        }
        Psquared                    = y;
        m_PulsarDetails.inclination = std::max(0.0, alpha);
    }

    m_PulsarDetails.spinPeriod    = std::sqrt(Psquared);
    m_PulsarDetails.spinFrequency = _2_PI / m_PulsarDetails.spinPeriod;                                                                     // pulsar spin frequency

    // std::cout << "initial spin period = " << initialSpinPeriod << std::endl;
    // std::cout << "final spin period = " << m_PulsarDetails.spinPeriod << std::endl;

    m_PulsarDetails.spinDownRate  = CalculateSpinDownRate(m_PulsarDetails.spinPeriod, m_MomentOfInertia_CGS, m_PulsarDetails.magneticField, m_Radius * RSOL_TO_KM) ; 

    // std::cout << "spinDownRate = " << m_PulsarDetails.spinDownRate << std::endl << std::endl;

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
 * Note 3 : A neutron star formed during this timestep (age 0) is not updated: the timestep was taken by
 *          its progenitor, so the neutron star keeps its birth values until the next timestep.
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

    if (utils::Compare(m_Age, 0.0) <= 0) return;                                                                                                // formed during this timestep - keep birth values

    if ((!p_RecycledNS && !p_CommonEnvelope) || (!p_RecycledNS && utils::Compare(p_MassGain, 0.0) == 0 )) {                                 // 'classical' isolated pulsars
        // std::cout << "Isolated spindown 1" << std::endl;
        SpinDownIsolatedPulsar(p_Stepsize, p_RecycledNS);                                                                                                 // spin down
    }
    else if (p_CommonEnvelope && (OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::SURFACE)) {                                    // mass transfer through CE when accretion happens at the surface of the NS

        double massG           = m_Mass * MSOL_TO_G;                                                                                        // mass in g
        double radiusCM        = m_Radius * RSOL_TO_CM;                                                                                     // radius in cm
        m_MomentOfInertia_CGS  = CalculateMomentOfInertiaCGS(); 
        double jAcc            = std::sqrt(G_CGS * massG * radiusCM) * p_MassGain;    
        m_AngularMomentum_CGS += jAcc;                                                                                                      // angular momentum of the accreted material as it falls onto the surface of the NS
        if (utils::Compare(m_PulsarDetails.accretedMass, 0.0) <= 0) {                                                                       // onset of accretion
            m_PulsarDetails.magneticFieldPreAccretion = m_PulsarDetails.magneticField;                                                      // record field
            ResetInclinationAtRecycling();                                                                                                  // redraw inclination if requested
        }
        if (utils::Compare(m_PulsarDetails.magneticField, NS::NS_MAG_FIELD_LOWER_LIMIT) < 0) {
            // if magnetic field is already lower than the lower limit, 
            // set it to the lower limit value.
            m_PulsarDetails.magneticField = NS::NS_MAG_FIELD_LOWER_LIMIT;
        }
        else {
            m_PulsarDetails.magneticField = CalculateMagneticFieldDecayAccretion_Static(m_PulsarDetails.magneticField, m_PulsarDetails.magneticFieldPreAccretion, m_PulsarDetails.accretedMass, p_MassGain);
        }
        m_PulsarDetails.accretedMass += p_MassGain;
        
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

        if (utils::Compare(m_PulsarDetails.accretedMass, 0.0) <= 0) {                                                                       // onset of accretion
            m_PulsarDetails.magneticFieldPreAccretion = m_PulsarDetails.magneticField;                                                      // record field
            ResetInclinationAtRecycling();                                                                                                  // redraw inclination if requested
        }

        double initialAngularMomentum_CGS = m_AngularMomentum_CGS;                                                                          // initial angular momentum
        double initialMagField            = m_PulsarDetails.magneticField;                                                                  // initial magnetic field
        double preAccretionMagField       = m_PulsarDetails.magneticFieldPreAccretion;                                                      // magnetic field at onset of accretion
        double accretedMassBefore         = m_PulsarDetails.accretedMass;                                                                   // mass accreted before this step (g)
        double initialMass                = m_Mass * MSOL_TO_G - p_MassGain;                                                                // initial mass of NS in g (m_Mass has already been updated for this step's accretion)

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
            double p_Mass, p_Radius, p_Radius_6, p_MagField, p_MagFieldPreAccretion, p_AccretedMass, p_Mdot, p_Epsilon, p_Propeller;
            ode(double mass, double radius, double radius_6, double magField, double magFieldPreAccretion, double accretedMass, double mdot, double epsilon, bool propeller) :
                p_Mass(mass), p_Radius(radius), p_Radius_6(radius_6), p_MagField(magField), p_MagFieldPreAccretion(magFieldPreAccretion), p_AccretedMass(accretedMass), p_Mdot(mdot), p_Epsilon(epsilon), p_Propeller(propeller) { }

            // x is the current state of the ODE (x[0] = angular momentum J)
            // dxdm is the change of state wrt mass (dxdm[0] = dJdm)
            // p_MassDelta is the cumulative change in mass of the NS
            void operator () (const state_type& x, state_type& dxdm, double p_MassDelta ) const {
                double m = p_Mass + p_MassDelta;
                // double B = (p_MagField - NS::NS_MAG_FIELD_LOWER_LIMIT) * std::exp(-p_MassDelta / NS::NS_DECAY_MASS_SCALE) + NS::NS_MAG_FIELD_LOWER_LIMIT;
                double B = utils::Compare(p_MagField, NS::NS_MAG_FIELD_LOWER_LIMIT) < 0 ? NS::NS_MAG_FIELD_LOWER_LIMIT : CalculateMagneticFieldDecayAccretion_Static(p_MagField, p_MagFieldPreAccretion, p_AccretedMass, p_MassDelta);
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
        (void)integrate_adaptive(stepper, ode{ initialMass, radius, radius_6, initialMagField, preAccretionMagField, accretedMassBefore, mDot, p_Epsilon, m_PulsarDetails.propellerMode }, x, 0.0, p_MassGain, massSlice);
                
        // final values
        m_AngularMomentum_CGS = x[0];
        m_MomentOfInertia_CGS = CalculateMomentOfInertiaCGS_Static(massFinal, radius);
        if (utils::Compare(initialMagField, NS::NS_MAG_FIELD_LOWER_LIMIT) < 0) {
            // if magnetic field is already lower than the lower limit, 
            // set it to the lower limit value.
            m_PulsarDetails.magneticField = NS::NS_MAG_FIELD_LOWER_LIMIT;
        }
        else {
            m_PulsarDetails.magneticField = CalculateMagneticFieldDecayAccretion_Static(initialMagField, preAccretionMagField, accretedMassBefore, p_MassGain);
        }
        m_PulsarDetails.accretedMass += p_MassGain;
        m_PulsarDetails.spinFrequency = m_AngularMomentum_CGS / m_MomentOfInertia_CGS;
        m_PulsarDetails.spinPeriod    = _2_PI / m_PulsarDetails.spinFrequency;
        double fDot                   = (m_AngularMomentum_CGS - initialAngularMomentum_CGS) / m_MomentOfInertia_CGS / p_Stepsize;          // eq. 11 in arxiv:1912.02415 
        m_PulsarDetails.spinDownRate  = -fDot * m_PulsarDetails.spinPeriod * m_PulsarDetails.spinPeriod / _2_PI;
    }      
    else {          
        // std::cout << "Isolated spindown 2" << std::endl;                                                                                                                        // otherwise...    
        // std::cout << "p_RecylcedNS = " << p_RecycledNS << std::endl;
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

