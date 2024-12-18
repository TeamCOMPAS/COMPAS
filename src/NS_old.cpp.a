#include "Rand.h"
#include "NS.h"

#include "gsl/gsl_poly.h"

/*
 * Calculate the luminosity of a Neutron Star
 *
 * Hurley et al. 2000, eq 93
 *
 * Called (indirectly) from GiantBranch, so must be static.
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
 * double ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double NS::ChooseTimestep(const double p_Time) const {
    double result = 500.0;                                                                                                      // default value

    if (p_Time < 0.01) {
        result = 0.001;
    }
    else if (p_Time < 0.1) {
        result = 0.01;
    }
    else if (p_Time < 1.0) {
        result = 0.1;
    }
    else if (p_Time < 10.0) {
        result = 1.0;
    }
    else if (p_Time < 500.0) {
        double slope      = 1.58859191006;                                                                                      // 1.58859191006 = log10(500.0) / (log10(500.0) - 1.0)
        double log10_step = slope * (log10(p_Time) - 1.0);
        result = PPOW(10.0, log10_step);
    }

    return result;
}


/*
 * Calculate Neutron Star radius according to selected equation of state (by commandline option)
 *
 *
 * double CalculateRadiusOnPhaseInKM_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius of Neutron Star in km
 */
double NS::CalculateRadiusOnPhaseInKM_Static(const double p_Mass) {

    double radius;

    switch (OPTIONS->NeutronStarEquationOfState()) {                                                                            // which equation-of-state?

        case NS_EOS::SSE:                                                                                                       // SSE
            radius = 10.0;
            break;

        case NS_EOS::ARP3: {                                                                                                    // ARP3

            // We don't extrapolate so masses outside table just set to extreme values

            std::map<double, double>::const_iterator iter;

            iter = ARP3MassRadiusRelation.begin();
            double ARP3MinimumMass = iter->first;
            double ARP3RadiusAtMinimumMass = iter->second;

            iter = ARP3MassRadiusRelation.end();
            double ARP3MaximumMass = (--iter)->first;
            double ARP3RadiusAtMaximumMass = iter->second;

            if (utils::Compare(p_Mass, ARP3MinimumMass) < 0) {
                radius = ARP3RadiusAtMinimumMass;
            }
            else if (utils::Compare(p_Mass, ARP3MaximumMass) > 0) {
                radius = ARP3RadiusAtMaximumMass;
            }
            else {
                radius = utils::SampleFromTabulatedCDF(p_Mass, ARP3MassRadiusRelation);
            }
        } break;

        default:                                                                                                                // unknown equation-of-state
            SHOW_WARN_STATIC(ERROR::UNKNOWN_NS_EOS,                                                                             // show warning
                             "Using default NS radius = 10.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            radius = 10.0;
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
    double luminosity  = CalculateLuminosityOnPhase_Static(p_Mass, 0.0);                                                        // Luminosity of Neutron Star as it cools
    double radius      = CalculateRadiusOnPhase_Static(p_Mass);                                                                 // Radius of Neutron Star in Rsol
    double temperature = BaseStar::CalculateTemperatureOnPhase_Static(luminosity, radius);                                      // Temperature of NS

    return std::make_tuple(luminosity, radius, temperature);
}


/*
 * Calculate the spin period of a Pulsar at birth according to selected distribution (by commandline option)
 *
 *
 * double CalculatePulsarBirthSpinPeriod()
 *
 * @return                                      Birth spin period of Pulsar in ms
 */
double NS::CalculatePulsarBirthSpinPeriod() {

	double pSpin;

    switch (OPTIONS->PulsarBirthSpinPeriodDistribution()) {                                                                     // which distribution?

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO:                                                                       // ZERO
            pSpin = 0.0;
            break;

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::FIXED:                                                                      // FIXED  constant value as used in default model in Oslowski et al 2011 https://arxiv.org/abs/0903.3538
            SHOW_WARN_STATIC(ERROR::UNSUPPORTED_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                                          // show warning
                             "Using spin = 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            pSpin = 0.0;
            break;

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM: {                                                                  // UNIFORM distribution between minimum and maximum value as in Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (default Pmin = and Pmax = )
                                                                                                                                // and also Kiel et al 2008 https://arxiv.org/abs/0805.0059 (default Pmin = 10 ms and Pmax 100 ms, section 3.4)

            double maximum = OPTIONS->PulsarBirthSpinPeriodDistributionMax();
            double minimum = OPTIONS->PulsarBirthSpinPeriodDistributionMin();

            pSpin = minimum + (RAND->Random() * (maximum - minimum));
            } break;

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL: {                                                                   // NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

            // Values hard-coded for now, can make them options if necessary
            // pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Mean = 300.0;
            // pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Std = 150.0;

            double mean  = 300.0;
            double sigma = 150.0;

            do { pSpin = RAND->RandomGaussian(sigma) + mean;} while (utils::Compare(pSpin, 0.0) < 0);

            } break;

        default:                                                                                                                // unknown distribution
            SHOW_WARN_STATIC(ERROR::UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                                              // show warning
                             "Using spin = 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            pSpin = 0.0;
    }

    return pSpin;
}


/*
 * Calculate (log10 of) the magnetic field (in G) for a Pulsar at birth
 * according to selected distribution (by commandline option)
 *
 *
 * double CalculatePulsarBirthMagneticField()
 *
 * @return                                      log10 of the birth magnetic field in G
 */
double NS::CalculatePulsarBirthMagneticField() {

	double log10B;

    switch (OPTIONS->PulsarBirthMagneticFieldDistribution()) {                                                                  // which distribution?

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO:                                                                    // ZERO
            log10B = 0.0;
            break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FIXED:                                                                   // FIXED - set to a fixed constant value
            SHOW_WARN_STATIC(ERROR::UNSUPPORTED_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,                                       // show warning
                             "Using 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            log10B = 0.0;
            break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG: {                                                             // FLAT IN LOG distribution from Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (log10B0min = , log10B0max = )

            double maximum = OPTIONS->PulsarBirthMagneticFieldDistributionMax();
            double minimum = OPTIONS->PulsarBirthMagneticFieldDistributionMin();

            log10B = minimum + (RAND->Random() * (maximum - minimum));

            } break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM: {                                                               // UNIFORM flat distribution used in Kiel et al 2008 https://arxiv.org/abs/0805.0059 (log10B0min = 11, log10B0max = 13.5 see section 3.4 and Table 1.)
            
      
            double maximum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMax());
            double minimum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMin());

            log10B = log10(minimum + (RAND->Random() * (maximum - minimum)));
            } break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL: {                                                             // LOG NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

            // Values hard-coded for now, can make them options if necessary
            // pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Mean = 12.65
            // pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Std = 0.55

            double mean  = 12.65;
            double sigma = 0.55;

            log10B = RAND->RandomGaussian(sigma) + mean;
            } break;

        default:                                                                                                                // unknown distribution
            SHOW_WARN_STATIC(ERROR::UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,                                           // show warning
                             "Using 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            log10B = 0.0;
    }

    return log10B;
}


/*
 * Calculate the moment of inertia for a Neutron Star using a model independent relation between
 * the moment of inertia, mass and radius of a neutron star - return MoI in CGS.
 * 
 * Uses m_Mass and m_Radius to calculate moment of inertia.
 *
 * Raithel et al. 2016, eq 8 in  https://arxiv.org/abs/1603.06594
 * https://tap.arizona.edu/sites/tap.arizona.edu/files/Raithel_2015_manuscript.pdf
 *
 *
 * double CalculateMomentOfInertiaCGS()
 *
 * @return                                      Moment of inertia in g cm^2
 */
double NS::CalculateMomentOfInertiaCGS() const {

    // pow() is slow - use multiplication

    double r_km   = m_Radius * RSOL_TO_KM;
    double m_r    = m_Mass / r_km;
    double m_r_4  = m_r * m_r * m_r * m_r;

    double r_cm   = m_Radius * RSOL_TO_CM;
    double r_cm_2 = r_cm * r_cm;

    return 0.237 * m_Mass * MSOL_TO_G * r_cm_2 * (1.0 + (4.2 * m_r) + 90.0 * m_r_4);
}


/*
 * Calculate the spin down rate for isolated Neutron Stars in cgs
 *
 * See Equation 2 in https://arxiv.org/pdf/1912.02415.pdf
 *
 * This is changed to the form of calculating spindown with P and Pdot, then convert to OmegaDot and to be recorded in the output file.
 * Evolution of the inclination between pulsar magnetic and rotational axes will be considered in a future version. 
 *
 * double CalculateSpinDownRate(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius)
 *
 * @param   [IN]    p_Omega                     Pulsar spin frequency. 
 * @param   [IN]    p_MomentOfInteria           Moment of Interia of the Neutron Star in kg m^2
 * @param   [IN]    p_MagField                  Magnetic field in Tesla
 * @param   [IN]    p_Radius                    Radius of the Neutron Star in kilometres
 * @return                                      Spin down rate (spin frequency derivative) of an isolated Neutron Star in s^(-2)
 */
double NS::CalculateSpinDownRate(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius) const {

   // pow() is slow - use multiplication

   double period            = _2_PI / p_Omega;                                                                                  // convert frequency to period
   double cgsRadius         = p_Radius * KM_TO_CM;                                                                              // radius in cm
   double radius_6          = cgsRadius * cgsRadius * cgsRadius * cgsRadius * cgsRadius * cgsRadius;
   double cgsMagField       = p_MagField * TESLA_TO_GAUSS;                                                                      // B field in G
   double magField_2        = cgsMagField * cgsMagField;
   constexpr double _8_PI_2 = 8.0 * PI_2;
   constexpr double _3_C_3  = 3.0E6 * C * C * C;                                                                                // 3.0 * (C * 100.0) * (C * 100.0) * (C * 100.0)
   double pDotTop           = _8_PI_2 * radius_6 * magField_2;
   double pDotBottom        = _3_C_3 * p_MomentOfInteria * period;
   double pDot              = pDotTop / pDotBottom;                                                                             // period derivative 
   
   return(-pDot * p_Omega / period);                                                                                            // convert period derivative to frequency derivative, which is what is recorded in the output
}


/*
 * Calculates and sets pulsar parameters at birth of pulsar
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

    m_PulsarDetails.magneticField     = PPOW(10.0, CalculatePulsarBirthMagneticField()) * GAUSS_TO_TESLA;                       // magnetic field in Gauss -> convert to Tesla
    m_PulsarDetails.spinPeriod        = CalculatePulsarBirthSpinPeriod();                                                       // spin period in ms
    m_PulsarDetails.spinFrequency     = _2_PI / (m_PulsarDetails.spinPeriod * SECONDS_IN_MS);
    m_PulsarDetails.birthPeriod       = m_PulsarDetails.spinPeriod * SECONDS_IN_MS;                                             // convert from ms to s 
    m_MomentOfInertia_CGS             = CalculateMomentOfInertiaCGS();                                                          // in CGS g cm^2
	
    // Note we convert neutronStarMomentOfInertia from CGS to SI here
    m_PulsarDetails.spinDownRate      = CalculateSpinDownRate(m_PulsarDetails.spinFrequency, m_MomentOfInertia_CGS, m_PulsarDetails.magneticField, m_Radius * RSOL_TO_KM);  
    m_PulsarDetails.birthSpinDownRate = m_PulsarDetails.spinDownRate; 
    m_AngularMomentum_CGS             = m_MomentOfInertia_CGS * m_PulsarDetails.spinFrequency;                                  // in CGS g cm^2 s^-1
}


/*
 * Update the magnetic field and spins of neutron stars when it's deemed as an isolated pulsar. 
 *
 * This function is called in multiple situations in the NS::UpdateMagneticFieldAndSpin() function
 *
 * Modifies the following class member variables:
 *
 *    m_AngularMomentum_CGS
 *    m_PulsarDetails.spinFrequency
 *    m_PulsarDetails.magneticField
 *    m_PulsarDetails.spinDownRate
 *
 *
 * void SpinDownIsolatedPulsar(const double p_Stepsize)
 *
 * @param   [IN]    p_Stepsize                  Timestep size for integration (in seconds)
 */
void NS::SpinDownIsolatedPulsar(const double p_Stepsize) {
    
    double NSradius_IN_CM         = m_Radius * RSOL_TO_KM * KM_TO_CM;
    double NSradius_3             = NSradius_IN_CM * NSradius_IN_CM * NSradius_IN_CM;
    double NSradius_6             = NSradius_3 * NSradius_3;
    constexpr double _8_PI_2      = 8.0 * PI_2;
    constexpr double _3_C_3       = 3.0 * C * C * C * 1000000.0;
    
    double initialMagField        = m_PulsarDetails.magneticField;                                                              // (in T)
    double initialMagField_G      = initialMagField * TESLA_TO_GAUSS;
    double initialSpinPeriod      = _2_PI / m_PulsarDetails.spinFrequency;
    double magFieldLowerLimit     = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField()) * GAUSS_TO_TESLA;    
    double magFieldLowerLimit_G   = magFieldLowerLimit * TESLA_TO_GAUSS;                                   
    double tau                    = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR;                                 

    // calculate isolated decay of the magnetic field for a neutron star
    // see Equation 6 in  arXiv:0903.3538v2       
    m_PulsarDetails.magneticField = magFieldLowerLimit + (initialMagField - magFieldLowerLimit) * exp(-p_Stepsize / tau);       // update pulsar magnetic field in SI. 
    
    // calculate the spin frequency for isolated neutron stars
    // see Equation 6 in arxiv:1912.02415
    // The rest of the calculations are carried out in cgs.   
    double constant2              = (_8_PI_2 * NSradius_6) / (_3_C_3 * m_MomentOfInertia_CGS);
    double term1                  = magFieldLowerLimit_G * magFieldLowerLimit_G * p_Stepsize;
    double term2                  = tau * magFieldLowerLimit_G * ( m_PulsarDetails.magneticField * TESLA_TO_GAUSS - initialMagField_G);
    double term3                  = (tau / 2.0) * (TESLA_TO_GAUSS * TESLA_TO_GAUSS * (m_PulsarDetails.magneticField * m_PulsarDetails.magneticField) - (initialMagField_G * initialMagField_G));
    double Psquared               = constant2 * (term1 - term2 - term3) + (initialSpinPeriod * initialSpinPeriod);
    
    double P_f                    = std::sqrt(Psquared);
    m_PulsarDetails.spinFrequency = _2_PI / P_f;                                                                                // pulsar spin frequency

    // calculate the spin down rate for isolated neutron stars
    // see Equation 4 in arXiv:0903.3538v2 (Our version is in cgs)      
    m_PulsarDetails.spinDownRate  = CalculateSpinDownRate(m_PulsarDetails.spinFrequency,m_MomentOfInertia_CGS  , m_PulsarDetails.magneticField, NSradius_IN_CM / KM_TO_CM);
    m_AngularMomentum_CGS         = m_PulsarDetails.spinFrequency * m_MomentOfInertia_CGS;                                      // angular momentum of star in CGS
}


/*
 * Calculate the magnetic field, spin, spin-down rate and angular momentum of a neutron star,
 * when it's accreting mass from companion through stable mass transfer. 
 * The calculations used in this function follow closely Sec. 2.2.1 in arxiv:1912.02415 
 * We carry out the calculations in this function using cgs units. 
 * This function is called in the NS::UpdateMagneticFieldAndSpin() function
 *
 * Returns, in a tuple, these four quantities of a neutron star:
 *
 *    magnetic field strength
 *    spin frequency
 *    spin-down rate
 *    angular momentum
 *
 * DBL_DBL_DBL_DBL PulsarAccretion(const double p_MagField, const double p_SpinFrequency, const double p_AngularMomentum, const double p_Stepsize, const double p_MassGainPerTimeStep, const double kappa, const double p_Epsilon)
 * @param   [IN]    p_MagField                  NS magnetic field strength at the beginning of accretion (in Tesla)
 * @param   [IN]    p_SpinFrequency             Spin frequency for the NS at the beginning of accretion (in Hz)
 * @param   [IN]    p_AngularMomentum           Angular momentum of the NS at the beginning of accretion (g cm-2 / s)
 * @param   [IN]    p_Stepsize                  Timestep size for integration (in seconds)
 * @param   [IN]    p_MassGainPerTimeStep       Mass transferred from the secondary for each iteration (in g)
 * @param   [IN]    p_Kappa                     Magnetic field mass decay scale (in g)
 * @param   [IN]    p_Epsilon                   Efficiency factor allowing for uncertainties of coupling magnetic field and matter.
 * @return                                      Tuple containing the updated magnetic field strength, spin frequency, spin-down rate and angular momentum of neutron star
 */
DBL_DBL_DBL_DBL NS::PulsarAccretion(const double p_MagField, const double p_SpinFrequency, const double p_AngularMomentum, const double p_Stepsize, const double p_MassGainPerTimeStep, const double p_Kappa, const double p_Epsilon) {
    double initialMagField_G      = p_MagField * TESLA_TO_GAUSS;         
    double magFieldLowerLimit_G   = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField()) ;                               // convert to Gauss                                
    double mass_g                 = m_Mass * MSOL_TO_G;                                                                    // in g
    double r_cm                   = m_Radius * RSOL_TO_KM * KM_TO_CM;    
    double angularMomentum        = p_AngularMomentum;

    //Magnetic field decay due to mass transfer. 
    //Follows Eq. 12 in arxiv:1912.02415 
    double newPulsarMagneticField = (initialMagField_G - magFieldLowerLimit_G) * exp(-p_MassGainPerTimeStep / p_Kappa) + magFieldLowerLimit_G;
    
    // calculate the Alfven radius for an accreting neutron star
    // see Equation 10 in arxiv:1912.02415 
    double mDot           =  p_MassGainPerTimeStep / p_Stepsize;
    double R_CM_6         =  r_cm * r_cm * r_cm * r_cm * r_cm * r_cm;
    double p              =  R_CM_6 * R_CM_6 / (mass_g * mass_g * mDot);
    double q              =  PPOW(p, 1.0/7.0);
    double alfvenConst    =  PPOW(2 * PI_2 / G_CGS, 1.0/7.0) ; 
    double alfvenRadius   =  alfvenConst * q * PPOW(initialMagField_G, 4.0/7.0); 
    double magneticRadius =  alfvenRadius / 2.0;
    // calculate the difference in the keplerian angular velocity at the magnetic radius
    // and surface angular velocity of the neutron star
    // magnetic radius is half of alfven radius 
    // see Equation 2 in 1994MNRAS.269..455J / Equation 9 in arxiv:1912.02415 
    double keplerianVelocityAtMagneticRadius        = std::sqrt((G_CGS) * mass_g / magneticRadius); 
    double keplerianAngularVelocityAtMagneticRadius = keplerianVelocityAtMagneticRadius / magneticRadius;
    double omegaDifference                          = keplerianAngularVelocityAtMagneticRadius - p_SpinFrequency;

    // calculate the change in angular momentum due to accretion
    // see Equation 12 in arXiv:0805.0059/ Equation 8 in arxiv:1912.02415 
    double Jdot        =  p_Epsilon * omegaDifference * magneticRadius * magneticRadius * mDot;
    angularMomentum    = angularMomentum + Jdot * p_Stepsize;
    return std::make_tuple(newPulsarMagneticField * GAUSS_TO_TESLA, angularMomentum / m_MomentOfInertia_CGS, Jdot / m_MomentOfInertia_CGS, angularMomentum);
}



/*
 * Calculate the magnetic field, spin, spin-down rate and angular momentum of a neutron star,
 * when it's accreting mass from companion through stable mass transfer. 
 * The calculations used in this function follow closely Sec. 2.2.1 in arxiv:1912.02415 
 * We carry out the calculations in this function using cgs units. 
 * This function is called in the NS::UpdateMagneticFieldAndSpin() function
 *
 * Returns the updated: 
 *    dJ/dm: change in angular momentum per mass transferred. 
 *
 * DBL_DBL_DBL_DBL PulsarAccretion(const double p_MagField, const double p_SpinFrequency, const double p_AngularMomentum, const double p_Stepsize, const double p_MassGainPerTimeStep, const double kappa, const double p_Epsilon)
 * @param   [IN]    p_MagField                  NS magnetic field strength at the beginning of accretion (in Gauss)
 * @param   [IN]    p_SpinFrequency             Spin frequency for the NS at the beginning of accretion (in Hz)
 * @param   [IN]    p_AngularMomentum           Angular momentum of the NS at the beginning of accretion (g cm-2 / s)
 * @param   [IN]    p_Stepsize                  Timestep size for integration (in seconds)
 * @param   [IN]    p_MassGainPerTimeStep       Mass transferred from the secondary for each iteration (in g)
 * @param   [IN]    p_Kappa                     Magnetic field mass decay scale (in g)
 * @param   [IN]    p_Epsilon                   Efficiency factor allowing for uncertainties of coupling magnetic field and matter.
 * @param   [IN]    p_MoI                       Moment of Inertia
 * @return                                      Tuple containing the updated magnetic field strength, spin frequency, spin-down rate and angular momentum of neutron star
 */
DBL_DBL_DBL_DBL NS::deltaAngularMomentumByPulsarAccretion(const double p_MassGainPerTimeStep, const double p_Mass, const double p_Radius, const double p_MagField, const double p_SpinFrequency, const double p_AngularMomentum, const double p_Stepsize, const double p_Kappa, const double p_Epsilon, const double p_MoI)  {
    if (utils::Compare(p_MassGainPerTimeStep, 0.0) <= 0) {
        return std::make_tuple(0.0, 0.0, 0.0, 0.0);
    }
    double initialMagField_G      = p_MagField;         
    double magFieldLowerLimit_G   = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField()) ;                               // convert to Gauss                                
    double mass_g                 = p_Mass * MSOL_TO_G;                                                                    // in g
    double r_cm                   = p_Radius * RSOL_TO_KM * KM_TO_CM;    
    double angularMomentum        = p_AngularMomentum;
    //std::cout << "in Func: " << p_Kappa <<" " << p_SpinFrequency << " " << mass_g << " " << r_cm  << " " << angularMomentum << std::endl;
    //Magnetic field decay due to mass transfer. 
    //Follows Eq. 12 in arxiv:1912.02415 
    double newPulsarMagneticField = (initialMagField_G - magFieldLowerLimit_G) * exp(-p_MassGainPerTimeStep / p_Kappa) + magFieldLowerLimit_G;
    
    // calculate the Alfven radius for an accreting neutron star
    // see Equation 10 in arxiv:1912.02415 
    double mDot           =  p_MassGainPerTimeStep / p_Stepsize;
    double R_CM_6         =  r_cm * r_cm * r_cm * r_cm * r_cm * r_cm;
    double p              =  R_CM_6 * R_CM_6 / (mass_g * mass_g * mDot);
    double q              =  PPOW(p, 1.0/7.0);
    double alfvenConst    =  PPOW(2 * PI_2 / G_CGS, 1.0/7.0) ; 
    double alfvenRadius   =  alfvenConst * q * PPOW(initialMagField_G, 4.0/7.0); 
    double magneticRadius =  alfvenRadius / 2.0;
    // std::cout<< "MagRad: " << mDot << " "  << r_cm << " " << p << " " << mass_g << " " << q << " " << alfvenConst << std::endl;
    // calculate the difference in the keplerian angular velocity at the magnetic radius
    // and surface angular velocity of the neutron star
    // magnetic radius is half of alfven radius 
    // see Equation 2 in 1994MNRAS.269..455J / Equation 9 in arxiv:1912.02415 
    double keplerianVelocityAtMagneticRadius        = std::sqrt((G_CGS) * mass_g / magneticRadius); 
    double keplerianAngularVelocityAtMagneticRadius = keplerianVelocityAtMagneticRadius / magneticRadius;
    double omegaDifference                          = keplerianAngularVelocityAtMagneticRadius - p_SpinFrequency;
    // std::cout<< "Denomenators: " << p_Kappa << " " << mass_g << " " << p_Stepsize <<  " " << magneticRadius << std::endl;
    // calculate the change in angular momentum due to accretion
    // see Equation 12 in arXiv:0805.0059/ Equation 8 in arxiv:1912.02415 
    double deltaJ        =  p_Epsilon * omegaDifference * magneticRadius * magneticRadius ;
    double Jdot          =  p_Epsilon * omegaDifference * magneticRadius * magneticRadius * mDot; 
    return std::make_tuple(newPulsarMagneticField, 
                           (angularMomentum + Jdot * p_Stepsize) / p_MoI, 
                           deltaJ * mDot / p_MoI, 
                           deltaJ);
}


/*
 * Update the magnetic field and spins of neutron stars. 
 * Modifies the following class member variables:
 *
 *    m_AngularMomentum_CGS
 *    m_PulsarDetails.spinFrequency
 *    m_PulsarDetails.magneticField
 *    m_PulsarDetails.spinDownRate
 *
 *
 * void UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope, 
 *                                 const bool   p_RecycledNS, 
 *                                 const double p_Stepsize, 
 *                                 const double p_MassGainPerTimeStep, 
 *                                 const double p_Epsilon)
 *
 * @param   [IN]    p_CommonEnvelope            Indicates whether there there is a common envelope - true or false
 * @param   [IN]    p_RecycledNS                Indicates whether this star is/was a recycled neutron star - true or false
 * @param   [IN]    p_Stepsize                  Timestep size for integration (in seconds)
 * @param   [IN]    p_MassGainPerTimeStep       Mass transferred from the secondary for each iteration (in kg)
 * @param   [IN]    p_Epsilon                   Efficiency factor allowing for uncertainties of coupling magnetic field and matter.
 */
void NS::UpdateMagneticFieldAndSpin(const bool p_CommonEnvelope, const bool p_RecycledNS, const double p_Stepsize, const double p_MassGainPerTimeStep, const double p_Epsilon) {
    double kappa                  = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;     
    double magFieldLowerLimit_G   = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField()) ; 
    double massNS                 = m_Mass;
    double radiusNS               = m_Radius;  
    if ((!p_RecycledNS && !p_CommonEnvelope) || (p_RecycledNS && utils::Compare(p_MassGainPerTimeStep, 0.0) == 0 )) {
        // These are the ''classical'' isolated pulsars. They experience spin-down.  
        SpinDownIsolatedPulsar(p_Stepsize);        
    }

    else if ( utils::Compare(p_MassGainPerTimeStep, 0.0) > 0) {
        // Next, we update pulsar that goes through mass transfer. 
        // Note that we do not need to set hard lower limit on the spin period,
        // as propeller effect should be included in the calculation, 
        // which means pulsars will start spinning down when the AM is no longer transferred along mass transfer. 
        if ((!p_CommonEnvelope) || (p_CommonEnvelope && OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::DISK)){
            // Pulsar evolution during stable mass transfer (RLOF), or if using the 'DISK' option 
            // in the CE accretion option

            //double newAM = -1.0; //updated angular momentum
            double thisB = m_PulsarDetails.magneticField * TESLA_TO_GAUSS ;
            double thisF = m_PulsarDetails.spinFrequency  ;
            double thisFdot = m_PulsarDetails.spinDownRate  ;
            double thisM = m_Mass ;
            double divideTimestepBy = 10.0;
            
            //while (utils::Compare(newAM, 0.0) < 0) {   
            double thisMassGain = p_MassGainPerTimeStep / G_TO_KG / divideTimestepBy; 
            double thisTimestepSize = p_Stepsize / divideTimestepBy;
            controlled_stepper_type controlled_stepper;
            //boost::numeric::odeint::runge_kutta4<state_type> stepper;

            state_type x(5);
            x[0] = m_AngularMomentum_CGS;
            x[1] = thisB;
            x[2] = thisF;
            x[3] = thisFdot;
            x[4] = thisM;
            std::cout<< "the original AM is " << x[0] << std::endl;
            std::cout<< massNS << ", "<< radiusNS << ", "<< thisB<< ", " << thisF<< ", " << x[0]<< ", " << thisTimestepSize<< ", "<< kappa<< ", " << p_Epsilon << std::endl;
            
            // Solve for the angular momentum of the NS after accretion. 
            // Use boost adaptive ODE solver for speed and accuracy and ensure the result is always positive. 
            struct DynamicODE
            {
                double mass, radius, magField, spinFrequency, angularMomentum, stepsize, kkappa, epsilon, momentI;
                DynamicODE( double mass0, double radius0, double magField0, double spinFreq0, double AM0, 
                double stepsize0, double kkappa0, double epsilon0, double momentI0) : mass(mass0), radius(radius0), magField(magField0), 
                spinFrequency(spinFreq0), angularMomentum(AM0), stepsize(stepsize0), kkappa(kkappa0), epsilon(epsilon0), momentI(momentI0)  {}

                void operator()( state_type& x , state_type& dxdt , double p_MassChange ) const {
                    DBL_DBL_DBL_DBL results = NS::deltaAngularMomentumByPulsarAccretion(
                        p_MassChange, mass, radius, magField, spinFrequency, angularMomentum, stepsize, kkappa, epsilon, momentI);
                    dxdt[0] = std::get<3>(results) ;
                    x[1] = std::get<0>(results);
                    x[2] = std::get<1>(results);
                    x[3] = std::get<2>(results);
                    x[4] = mass + p_MassChange/MSOL_TO_G;
            }
            };
            std::cout << thisMassGain << " " << p_MassGainPerTimeStep << std::endl;
            //DynamicODE ode(massNS, radiusNS, thisB, thisF, x[0], thisTimestepSize, kappa, p_Epsilon);
            integrate_adaptive(controlled_stepper, DynamicODE{x[4], radiusNS, x[1], x[2], x[0], thisTimestepSize, 
            kappa, p_Epsilon, m_MomentOfInertia_CGS}, x, 0.0, p_MassGainPerTimeStep / G_TO_KG, thisMassGain);
            // for (double dm = 0.0; dm < p_MassGainPerTimeStep / G_TO_KG; dm += thisMassGain){
            //     stepper.do_step( ode,  x, dm, thisMassGain);
            
            //     thisB = (thisB - magFieldLowerLimit_G) * exp(-thisMassGain / kappa) + magFieldLowerLimit_G;;
            //     thisB = thisB ; 
            //     thisF = x[0] / m_MomentOfInertia_CGS; 
            //     thisFdot = (x[0]-newAM)/ thisTimestepSize / m_MomentOfInertia_CGS;
            //     newAM = x[0];
            //     std::cout << thisB << ", " << thisF << ", " << newAM <<std::endl;
            // }
            std::cout<< "the updated AM is " << x[0] << " and others " << x[1] << " " << x[2]
            << " " << x[3] << " " << x[4] << std::endl;
            // newAM = x[0];
            // if (utils::Compare(newAM,0) < 0) {
            //     std::cout << "before stopping " << thisMassGain  << " " << p_MassGainPerTimeStep / G_TO_KG << 
            //     " " << newAM << " " <<  x[1] << " " << x[2] << " " << x[3] << "  " << x[4] << std::endl;
            //     divideTimestepBy *= 2.0 ; 
            //     thisB = m_PulsarDetails.magneticField * TESLA_TO_GAUSS ;
            //     thisF = m_PulsarDetails.spinFrequency  ;
            //     thisFdot = m_PulsarDetails.spinDownRate  ;
            //     std::exit(1);
            // } else {
            //     std::cout << "DB "  << divideTimestepBy << " SZ " << p_Stepsize <<  std::endl;
            // }
            m_AngularMomentum_CGS         = x[0] ;
            m_PulsarDetails.magneticField = x[1] * GAUSS_TO_TESLA;
            m_PulsarDetails.spinFrequency = x[2] ;
            m_PulsarDetails.spinDownRate  = x[3] ;

            //}
            // m_AngularMomentum_CGS = newAM;

            // m_PulsarDetails.magneticField = thisB * GAUSS_TO_TESLA;
            // m_PulsarDetails.spinFrequency = thisF; 
            // m_PulsarDetails.spinDownRate  = thisFdot;   
        }
        
        
        else if (p_CommonEnvelope && OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::SURFACE) {
            // Mass transfer through CE when accretion happens at the surface of the NS
            double initialMagField        = m_PulsarDetails.magneticField;
            double initialMagField_G      = initialMagField * TESLA_TO_GAUSS;        
            
            double mass_g                 = m_Mass * MSOL_TO_G;                                                                    // in g
            double r_cm                   = m_Radius * RSOL_TO_KM * KM_TO_CM;                                                        // in cm
            m_MomentOfInertia_CGS         = CalculateMomentOfInertiaCGS() ; 
            double newPulsarMagneticField = (initialMagField_G - magFieldLowerLimit_G) * exp(-p_MassGainPerTimeStep / G_TO_KG / kappa) + magFieldLowerLimit_G;
            double r_cm_3                 = r_cm * r_cm * r_cm ;
            // angular momentum of the accreted materials as they fall onto the surface of the NS. 
            double Jacc                   = m_MomentOfInertia_CGS * PPOW(G_CGS * mass_g /r_cm_3, 0.5) * p_MassGainPerTimeStep / G_TO_KG / mass_g ;
            
            m_AngularMomentum_CGS         = m_AngularMomentum_CGS + Jacc ;
            m_PulsarDetails.magneticField = newPulsarMagneticField * GAUSS_TO_TESLA;
            m_PulsarDetails.spinFrequency = m_AngularMomentum_CGS / m_MomentOfInertia_CGS ;
            m_PulsarDetails.spinDownRate  = Jacc / p_Stepsize  / m_MomentOfInertia_CGS;
        }
   
        else if (p_CommonEnvelope && OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::ZERO) {
            // If mass transfer is happening through CE and chosen 'ZERO' in CE accretion mode,
            // treat the pulsar as isolated and should be spun down.  
            SpinDownIsolatedPulsar(p_Stepsize);
        }
    }
}
