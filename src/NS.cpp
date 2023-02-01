#include "Rand.h"

#include "NS.h"


/*
 * Calculate the luminosity of a Neutron Star
 *
 * Hurley et al. 2000, eq 93
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
    return 0.02 * p_Mass / (t * t);
}

/*
 * Choose timestep for pulsarevolution
 *
 * Can choose whatever way you want.
 * Given in the form discussed by Simon & Robert  
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double NS::ChooseTimestep(const double p_Time) const {
    double result;
    result = 1;
    return result;
    // if (p_Time < 50.0) {
    //     result = 0.1;
    // }
    // else if (p_Time <= 500) {
    //     double slope      = (log10(500) - log10(0.1)) / (log10(500) - log10(50));
    //     double intersect  = log10(0.1) - slope * log10(50);
    //     double log10_step = intersect + slope * log10(p_Time);
    //     result = pow(10, log10_step);
    // }
    // else {
    //     result = 500.0;
    // }

    // //below was the Hurley+2000 discussion. 
    // // double dtk = std::min(std::max(0.1, 10.0 * std::max(0.1, 10.0 * p_Time)), 5.0E2);
    // // double dte = dtk;

    // // return std::max(dte, NUCLEAR_MINIMUM_TIMESTEP);
    // std::cout << "NS Time step = " << result << " at " << p_Time << "\n";
    // std::cout.flush();

    // return result;
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

    switch (OPTIONS->NeutronStarEquationOfState()) {                                                    // which equation-of-state?

        case NS_EOS::SSE:                                                                               // SSE
            radius = 10.0;
            break;

        case NS_EOS::ARP3: {                                                                            // ARP3

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
            else{
                radius = utils::SampleFromTabulatedCDF(p_Mass, ARP3MassRadiusRelation);
            }
        } break;

        default:                                                                                        // unknown equation-of-state
            SHOW_WARN_STATIC(ERROR::UNKNOWN_NS_EOS,                                                     // show warning
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
 *
 * DBL_DBL_DBL CalculateCoreCollapseSNParams_Static(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Tuple containing Luminosity, Radius and Temperature of Neutron Star
 */
DBL_DBL_DBL NS::CalculateCoreCollapseSNParams_Static(const double p_Mass) {
    double luminosity  = CalculateLuminosityOnPhase_Static(p_Mass, 0.0);                    // Luminosity of Neutron Star as it cools
    double radius      = CalculateRadiusOnPhase_Static(p_Mass);                             // Radius of Neutron Star in Rsol
    double temperature = BaseStar::CalculateTemperatureOnPhase_Static(luminosity, radius);  // Temperature of NS

    return std::make_tuple(luminosity, radius, temperature);
}


/*
 * Calculate the spin period of a Pulsar at birth according to selected distribution (by commandline option)
 *
 *
 * double CalculatePulsarBirthSpinPeriod_Static()
 *
 * @return                                      Birth spin period of Pulsar in ms
 */
double NS::CalculatePulsarBirthSpinPeriod_Static() {

	double pSpin;

    switch (OPTIONS->PulsarBirthSpinPeriodDistribution()) {                                                                 // which distribution?

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::ZERO:                                                                   // ZERO
            pSpin = 0.0;
            break;

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::FIXED:                                                                  // FIXED  constant value as used in default model in Oslowski et al 2011 https://arxiv.org/abs/0903.3538
            SHOW_WARN_STATIC(ERROR::UNSUPPORTED_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                                      // show warning
                             "Using spin = 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            pSpin = 0.0;
            break;

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM: {                                                              // UNIFORM distribution between minimum and maximum value as in Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (default Pmin = and Pmax = )
                                                                                                                            // and also Kiel et al 2008 https://arxiv.org/abs/0805.0059 (default Pmin = 10 ms and Pmax 100 ms, section 3.4)

            double maximum = OPTIONS->PulsarBirthSpinPeriodDistributionMax();
            double minimum = OPTIONS->PulsarBirthSpinPeriodDistributionMin();

            pSpin = minimum + (RAND->Random() * (maximum - minimum));
            } break;

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL: {                                                               // NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

            // Values hard-coded for now, can make them options if necessary
            // pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Mean = 300.0;
            // pulsarBirthSpinPeriodDistributionFaucherGiguereKaspi2006Std = 150.0;

            double mean  = 300.0;
            double sigma = 150.0;

            do { pSpin = RAND->RandomGaussian(sigma) + mean;} while (utils::Compare(pSpin, 0.0) < 0);

            } break;

        default:                                                                                                            // unknown distribution
            SHOW_WARN_STATIC(ERROR::UNKNOWN_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION,                                          // show warning
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
 * double CalculatePulsarBirthMagneticField_Static()
 *
 * @return                                      log10 of the birth magnetic field in G
 */
double NS::CalculatePulsarBirthMagneticField_Static() {

	double log10B;

    switch (OPTIONS->PulsarBirthMagneticFieldDistribution()) {                                                          // which distribution?

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::ZERO:                                                            // ZERO
            log10B = 0.0;
            break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FIXED:                                                           // FIXED - set to a fixed constant value
            SHOW_WARN_STATIC(ERROR::UNSUPPORTED_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,                               // show warning
                             "Using 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            log10B = 0.0;
            break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG: {                                                     // FLAT IN LOG distribution from Oslowski et al 2011 https://arxiv.org/abs/0903.3538 (log10B0min = , log10B0max = )

            double maximum = OPTIONS->PulsarBirthMagneticFieldDistributionMax();
            double minimum = OPTIONS->PulsarBirthMagneticFieldDistributionMin();

            log10B = minimum + (RAND->Random() * (maximum - minimum));

            } break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM: {                                                       // UNIFORM flat distribution used in Kiel et al 2008 https://arxiv.org/abs/0805.0059 (log10B0min = 11, log10B0max = 13.5 see section 3.4 and Table 1.)
            
      
            double maximum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMax());
            double minimum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMin());

            log10B = log10(minimum + (RAND->Random() * (maximum - minimum)));
            } break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL: {                                                     // LOG NORMAL distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585

            // Values hard-coded for now, can make them options if necessary
            // pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Mean = 12.65
            // pulsarBirthMagneticFieldDistributionFaucherGiguereKaspi2006Std = 0.55

            double mean  = 12.65;
            double sigma = 0.55;

            log10B = RAND->RandomGaussian(sigma) + mean;
            } break;

        default:                                                                                                        // unknown distribution
            SHOW_WARN_STATIC(ERROR::UNKNOWN_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION,                                   // show warning
                             "Using 0.0",
                             OBJECT_TYPE::BASE_STAR,
                             STELLAR_TYPE::NEUTRON_STAR);
            log10B = 0.0;
    }

    return log10B;
}


/*
 * Calculate the moment of inertia for a Neutron Star using a model independent relation between
 * the moment of inertia, mass and radius of a neutron star
 *
 * Raithel et al. 2016, eq 8 in  https://arxiv.org/abs/1603.06594
 * https://tap.arizona.edu/sites/tap.arizona.edu/files/Raithel_2015_manuscript.pdf
 *
 *
 * double CalculateMomentOfInertia_Static(const double p_Mass, const double p_Radius)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Radius                    Radius in km
 * @return                                      Moment of inertia in g cm^2
 */
double NS::CalculateMomentOfInertia_Static(const double p_Mass, const double p_Radius) {

    // pow() is slow - use multiplication
    double m_r    = p_Mass / p_Radius;
    double m_r_4  = m_r * m_r * m_r * m_r;
    double r_km   = p_Radius * KM_TO_CM;
    double r_km_2 = r_km * r_km;

    return 0.237 * p_Mass * MSOL_TO_G * r_km_2 * (1.0 + (4.2 * m_r) + 90 * m_r_4);
}


/*
 * Calculate the spin down rate for isolated Neutron Stars in SI
 *
 * See Equation X in arxiv.xxxx  JR: Find reference
 *
 *
 * double CalculateSpinDownRate_Static(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius, double const p_Alpha)
 *
 * @param   [IN]    p_Omega                     Mass in Msol
 * @param   [IN]    p_MomentOfInteria           Moment of Interia of the Neutron Star in kg m^2
 * @param   [IN]    p_MagField                  Magnetic field in Tesla
 * @param   [IN]    p_Radius                    Radius of the Neutron Star in metres
 * @param   [IN]    p_Alpha                     Mass in Msol
 * @return                                      Spin down rate of an isolated Neutron Star in SI
 */
double NS::CalculateSpinDownRate_Static(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius, double const p_Alpha) {

   // pow() is slow - use multiplication
//    double omega_3              = p_Omega * p_Omega * p_Omega;
//    double radius_6             = p_Radius * p_Radius * p_Radius * p_Radius * p_Radius * p_Radius;
//    double magField_2           = p_MagField * p_MagField;
//    constexpr double PI_8       = 8.0 * M_PI;
//    constexpr double MU_0_3_C_3 = 3.0 * MU_0 * C * C * C;
//    double omegaDotTop          = PI_8 * omega_3 * radius_6 * magField_2;
//    double omegaDotBottom       = MU_0_3_C_3 * p_MomentOfInteria;
//    double omega_3              = p_Omega * p_Omega * p_Omega;
//    double cgs_Radius           = p_Radius * 1e3;
//    double radius_6             = cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius;
//    double cgs_MagField         = p_MagField * 1e4; 
//    double magField_2           = cgs_MagField * cgs_MagField;
//    constexpr double PI_8       = 8.0 * M_PI;
//    constexpr double _3_C_3     = 3.0 * C * C * C;
//    double omegaDotTop          = PI_8 * omega_3 * radius_6 * magField_2;
//    double omegaDotBottom       = _3_C_3 * p_MomentOfInteria / 1e3 / 1e2 / 1e2;
//   double omega_3              = p_Omega * p_Omega * p_Omega;
// CalculateSpinDownRate_Static(m_PulsarDetails.spinFrequency, m_MomentOfInertia * factor, m_PulsarDetails.magneticField,  NEUTRON_STAR_RADIUS, 1.0)
   double Period               = 2 * M_PI / p_Omega ;
   double cgs_Radius           = p_Radius * KM_TO_CM;
   double radius_6             = cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius;
   double cgs_MagField         = p_MagField * TESLA_TO_GAUSS; 
   double magField_2           = cgs_MagField * cgs_MagField;
   constexpr double _8_PI_2    = 8.0 * M_PI * M_PI;
   constexpr double _3_C_3     = 3.0 * (C * 100) * (C * 100) * (C * 100) ;
   double pDotTop              = _8_PI_2 * radius_6 * magField_2;
   double pDotBottom           = _3_C_3 * p_MomentOfInteria * Period ;
   double pDot                 = pDotTop / pDotBottom ;
   return(  -1 * pDot * p_Omega / Period);
//    string name1 = "omega_3: ";
//    std::cout <<  name1 ;
//    std::cout << omega_3 << "\n";
//    string name2 = " radius_6: ";
//    std::cout << name2;
//    std::cout << radius_6 << "\n";
//    string name3 = " magField_2: ";
//    std::cout << name3 ;
//    std::cout << magField_2;
//    string name4 = " PI_8: ";
//    std::cout << name4 ;
//    std::cout << PI_8 << "\n";
//    string name5 = " MU_0_3_C_3: ";
//    std::cout << name5 ;
//    std::cout << MU_0_3_C_3 << "\n";
//    string name6 = " omegaDotTop: ";
//    std::cout << name6;
//    std::cout << omegaDotTop << "\n";
//    string name7 = " omegaDotBottom: ";
//    std::cout << name7;
//    std::cout << omegaDotBottom << "\n";
//    std::cout.flush();
//   return (-omegaDotTop / omegaDotBottom);
}


/*
 * Calculates and sets pulsar parameters
 *
 * JR: todo: description?
 *
 *
 *
 * void CalculateAndSetPulsarParameters()
 */
void NS::CalculateAndSetPulsarParameters() {

    m_PulsarDetails.magneticField = PPOW(10.0, CalculatePulsarBirthMagneticField_Static()) * GAUSS_TO_TESLA ;                                                                       // magnetic field in Gauss -> convert to Tesla
    m_PulsarDetails.spinPeriod    = CalculatePulsarBirthSpinPeriod_Static();                                                                                                        // spin period in ms
    m_PulsarDetails.spinFrequency = _2_PI / (m_PulsarDetails.spinPeriod * SECONDS_IN_MS);
    m_PulsarDetails.birthPeriod   = m_PulsarDetails.spinPeriod; 
    
    m_MomentOfInertia             = CalculateMomentOfInertia_Static(m_Mass, NEUTRON_STAR_RADIUS * RSOL_TO_KM) ;      //m_Radius * RSOL_TO_KM) ;                                                                                // in CGS g cm^2

	// Note we convert neutronStarMomentOfInertia from CGS to SI here
	constexpr double factor       = G_TO_KG * CM_TO_M * CM_TO_M;
    m_PulsarDetails.spinDownRate  = CalculateSpinDownRate_Static(m_PulsarDetails.spinFrequency, m_MomentOfInertia, m_PulsarDetails.magneticField,  NEUTRON_STAR_RADIUS * RSOL_TO_KM, 1.0);   // was using m_radius, alpha = 1.0
    m_PulsarDetails.birthSpindown = m_PulsarDetails.spinDownRate; 
    m_AngularMomentum             = _2_PI * m_MomentOfInertia / (m_PulsarDetails.spinPeriod * SECONDS_IN_MS) * factor;                                                              // in kg m^2 sec^-1
    //std::cout << "Mass: " << m_Mass << "(MSol) Radius: " <<  NEUTRON_STAR_RADIUS << "(RSol) I: " << m_MomentOfInertia << "g cm2" << " Spin: " << m_PulsarDetails.spinPeriod * SECONDS_IN_MS << "\n";
    //std::cout.flush();
}


/*
 * Update the magnetic field and spins of neutron stars
 *
 * Modifies the following class member variables:
 *
 *    m_AngularMomentum
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
 * @param   [IN]    p_RecycledNS                Indicates whether this star is/was a recyled neutron star - true or false
 * @param   [IN]    p_Stepsize                  Timestep size for integration (in seconds)
 * @param   [IN]    p_MassGainPerTimeStep       Mass loss from the secondary for each iteration (in kg)
 * @param   [IN]    p_Epsilon                   Uncertainty due to mass loss
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
void NS::UpdateMagneticFieldAndSpin(const bool p_CommonEnvelope, const bool p_RecycledNS, const double p_Stepsize, const double p_MassGainPerTimeStep, const double p_Epsilon) {

    //constexpr double NSRadius_IN_M = NEUTRON_STAR_RADIUS * RSOL_TO_KM * KM_TO_M ;
    constexpr double NSRadius_IN_CM = NEUTRON_STAR_RADIUS * RSOL_TO_KM * KM_TO_CM ;
    constexpr double NSRadius_3    = NSRadius_IN_CM * NSRadius_IN_CM * NSRadius_IN_CM;
    constexpr double NSRadius_6    = NSRadius_3 * NSRadius_3;
    constexpr double _8_PI_2          = 8.0 * M_PI * M_PI;
    //constexpr double MU_0_3_C_3    = 3.0 * MU_0 * C * C * C;
    constexpr double _3_C_3        = 3.0 * C * C * C * 1000000;
    //constexpr double unitsMoI      = G_TO_KG * CM_TO_M * CM_TO_M;
    
    //double mass               = m_Mass * MSOL_TO_KG;
    double mass               = m_Mass * MSOL_TO_G;
    double radius             = NEUTRON_STAR_RADIUS * RSOL_TO_CM;// m_Radius * RSOL * 100; // in meters
    double initialMagField    = m_PulsarDetails.magneticField; // (in T)
    double initialMagField_G  = initialMagField * TESLA_TO_GAUSS;
    double initialSpinPeriod  = 2 * M_PI / m_PulsarDetails.spinFrequency;
    double magFieldLowerLimit = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField()) * GAUSS_TO_TESLA;    
    double magFieldLowerLimit_G = magFieldLowerLimit * TESLA_TO_GAUSS;                                   
    double momentOfInertia    = m_MomentOfInertia; //* unitsMoI;
    double tau                = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR;                                 
    double kappa              = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;      
    //std::cout << "tau:" << tau <<  " Initial Period:  " << initialSpinPeriod << " Initial Pdot: " << m_PulsarDetails.spinDownRate << "\n";// (-1 * m_PulsarDetails.spinDownRate * initialSpinPeriod / m_PulsarDetails.spinFrequency) << "\n";                                                    

    if ((!p_RecycledNS && !p_CommonEnvelope) || (!p_RecycledNS && utils::Compare(p_MassGainPerTimeStep, 0.0) == 0 )) {

        // calculate isolated decay of the magnetic field for a neutron star see Equation 6 in  arXiv:0903.3538v2       
        //std::cout << magFieldLowerLimit << ' ' << p_Stepsize << ' ' << tau << ' ' << initialMagField << "\n";
        //m_PulsarDetails.magneticField = magFieldLowerLimit + exp(-p_Stepsize / tau) * (initialMagField - magFieldLowerLimit);                   // pulsar magnetic field
        m_PulsarDetails.magneticField = magFieldLowerLimit + (initialMagField - magFieldLowerLimit) * exp(-p_Stepsize / tau);
        // calculate the spin down rate for isolated neutron stars, see Equation 6 in arxiv:1912.02415   
        //std::cout << "for K: _8_PI_2 = " << _8_PI_2 << " NSRadius_6 = " << NSRadius_6 << " _3_C_3 = " << _3_C_3 << " MOI " << momentOfInertia << "\n";   
        double constant_2             = (_8_PI_2 * NSRadius_6) / (_3_C_3 * momentOfInertia);
        double term1                  = magFieldLowerLimit_G * magFieldLowerLimit_G * p_Stepsize;
        double term2                  = tau * magFieldLowerLimit_G * ( m_PulsarDetails.magneticField * TESLA_TO_GAUSS - initialMagField_G);
        double term3                  = (tau / 2.0) * (TESLA_TO_GAUSS * TESLA_TO_GAUSS * (m_PulsarDetails.magneticField * m_PulsarDetails.magneticField) - (initialMagField_G * initialMagField_G));
        //double oneOverOmegaSquared    = constant_2 * (term1 - term2 - term3) + (1.0 / (m_PulsarDetails.spinFrequency * m_PulsarDetails.spinFrequency));
        double PSquared    = 2 * constant_2 * (term1 - term2 - term3) + (initialSpinPeriod * initialSpinPeriod);
        //double oneOverOmegaSquared    = constant_2 * (term1 - term2 - term3) + (1.0 / (initialSpinFreq * initialSpinFreq));
        // string name00 = "MOI: ";
        // std::cout << name00;
        // std::cout << momentOfInertia << "\n";
        // string name0 = "B field:";
        // std::cout << name0;
        // std::cout << m_PulsarDetails.magneticField << "\n";
        //string name1 = "constant_2: ";
        // std::cout << name1;
        // std::cout << constant_2 << "\n";
        // string name2 = "term1: ";
        // std::cout << name2;
        // std::cout << term1 << "\n";
        // string name3 = "term2: ";
        // std::cout << name3;
        // std::cout << term2 << "\n";
        // string name4 = "term3: ";
        // std::cout << name4;
        // std::cout << term3 << "\n";
        // string name5 = "oneOverOmegaSquared : ";
        // std::cout << name5;
        // std::cout << oneOverOmegaSquared  << "\n";
        // string name6 = "Spin Frequency: ";
        // std::cout << name6;
        // std::cout << 1.0 / std::sqrt(oneOverOmegaSquared) << "\n";
        // std::cout.flush();
        double   P_f                  = std::sqrt(PSquared);
        m_PulsarDetails.spinFrequency = 2 * M_PI / P_f;                                                                        // pulsar spin frequency

        // calculate the spin down rate for isolated neutron stars, see Equation 4 in arXiv:0903.3538v2 (Our version is in SI)      
        //double omegaDotTop           = PI_8 * m_PulsarDetails.spinFrequency * m_PulsarDetails.spinFrequency * m_PulsarDetails.spinFrequency * NSRadius_6 * m_PulsarDetails.magneticField * m_PulsarDetails.magneticField;
        //double omegaDotBottom        = MU_0_3_C_3 * momentOfInertia;
        double pDotTop               = constant_2 * TESLA_TO_GAUSS * TESLA_TO_GAUSS * m_PulsarDetails.magneticField * m_PulsarDetails.magneticField;
        //double pDotBot               = _3_C_3 * momentOfInertia * P_f;
        double pDot                  = pDotTop / P_f;// / pDotBot;
        m_PulsarDetails.spinDownRate = -2 * M_PI * pDot / (P_f * P_f);  
        //std::cout << "Stepsize: " << p_Stepsize << " Radius: " << m_Radius << " Period: " << P_f << " Pdot: " << pDot << " B: " << m_PulsarDetails.magneticField << "\n";
        //std::cout.flush();
        //m_PulsarDetails.spinDownRate = -omegaDotTop / omegaDotBottom;                                                                           // pulsar spin down rate
   
        m_AngularMomentum            = m_PulsarDetails.spinFrequency * momentOfInertia;                                                         // angular momentum of star
    }
    else if ((p_RecycledNS || p_CommonEnvelope) && utils::Compare(p_MassGainPerTimeStep, 0.0) > 0) {

        // calculate the Alfven radius for an accreting neutron star, see Equation 8 in  arXiv:0903.3538v2       
        double mDot         = p_MassGainPerTimeStep / p_Stepsize ;
        double p            = ((radius * radius * radius * radius * radius * radius) / (std::sqrt(mass) * mDot));
        double q            = PPOW(p, 2.0 / 7.0);
        double constant     = PPOW((2.0 * M_PI * M_PI) / (G * MU_0 * MU_0), 1.0 / 7.0);                                                         
        double alfvenRadius = constant * q * PPOW(m_PulsarDetails.magneticField, 4.0 / 7.0);
  
        // calculate the difference in the keplerian angular velocity and surface angular velocity of the neutron star in m - see Equation 2 in 1994MNRAS.269..455J       
        double keplerianVelocityAtAlfvenRadius        = std::sqrt(G * mass) / std::sqrt(alfvenRadius / 2.0);
        double keplerianAngularVelocityAtAlfvenRadius = 2.0 * (keplerianVelocityAtAlfvenRadius / alfvenRadius);
        double velocityDifference                     = keplerianAngularVelocityAtAlfvenRadius - m_PulsarDetails.spinFrequency;

        // calculate accretion induced magnetic field decay for an accreting neutron star, see Equation 7 in arXiv:0903.3538v2       
        m_PulsarDetails.magneticField = magFieldLowerLimit + exp(-(1.0 / kappa) * p_MassGainPerTimeStep) * (m_PulsarDetails.magneticField - magFieldLowerLimit);  // pulsar magnetic field

        // calculate the change in angular momentum due to accretion, see Equation 12 in arXiv:0805.0059/ Equation 8 in arxiv:1912.02415 
        double deltaAngularMomentum   = 0.25 * p_Epsilon * alfvenRadius * alfvenRadius * p_MassGainPerTimeStep * velocityDifference;

        m_AngularMomentum             = m_AngularMomentum + deltaAngularMomentum;                                                               // angular momentum of star
        m_PulsarDetails.spinFrequency = m_AngularMomentum / momentOfInertia;                                                                    // pulsar spin frequency
        m_PulsarDetails.spinDownRate  = (deltaAngularMomentum / p_Stepsize) / momentOfInertia;                                                  // pulsar spin down rate
    }
}
