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
 * Choose timestep for Pulsar Evolution
 *
 * Can choose whatever way you want.
 * Given in the form discussed by Simon & Robert  
 *
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double NS::ChooseTimestep(const double p_Time) const {
    double result;
    if (p_Time < 0.01) {
        result = 0.001;
    }
    else if (p_Time < 0.1) {
        result = 0.01;
    }
    else if (p_Time < 1) {
        result = 0.1;
    }
    else if (p_Time < 10) {
        result = 1;
    }
    else if (p_Time <= 500) {
        double slope      = (log10(500) - log10(1)) / (log10(500) - log10(10));
        double intersect  = log10(1) - slope * log10(10);
        double log10_step = intersect + slope * log10(p_Time);
        result = pow(10, log10_step);
    }
    else {
        result = 500.0;
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
 * Calculate the spin down rate for isolated Neutron Stars in cgs
 *
 * See Equation 2 in https://arxiv.org/pdf/1912.02415.pdf
 *
 * This is changed to the form 
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

   double Period               = 2 * M_PI / p_Omega ;           //Convert frequency to period
   double cgs_Radius           = p_Radius * KM_TO_CM;           // radius in cm
   double radius_6             = cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius * cgs_Radius;
   double cgs_MagField         = p_MagField * TESLA_TO_GAUSS;   //B field in G
   double magField_2           = cgs_MagField * cgs_MagField;
   constexpr double _8_PI_2    = 8.0 * M_PI * M_PI;
   constexpr double _3_C_3     = 3.0 * (C * 100) * (C * 100) * (C * 100) ;
   double pDotTop              = _8_PI_2 * radius_6 * magField_2;
   double pDotBottom           = _3_C_3 * p_MomentOfInteria * Period ;
   double pDot                 = pDotTop / pDotBottom ;         //Period derivative 
   return(  -1 * pDot * p_Omega / Period);                      //convert period derivative to frequency derivative

}


/*
 * Calculates and sets pulsar parameters
 *
 * JR: todo: description?
 * Robert: I think this function calculates the pulsar parameters at birth. 
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

}



/*
 * Update the magnetic field and spins of neutron stars when it's deemed as an isolated pulsar. 
 *
 * This function is called in multiple situations in the NS::UpdateMagneticFieldAndSpin() function
 *
 * Modifies the following class member variables:
 *
 *    m_AngularMomentum
 *    m_PulsarDetails.spinFrequency
 *    m_PulsarDetails.magneticField
 *    m_PulsarDetails.spinDownRate
 *
 *
 * void spinDownIsolatePulsar(const double p_Stepsize)
 *
 * @param   [IN]    p_Stepsize                  Timestep size for integration (in seconds)
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
void NS::spinDownIsolatePulsar(const double p_Stepsize) {

    constexpr double NSRadius_IN_CM = NEUTRON_STAR_RADIUS * RSOL_TO_KM * KM_TO_CM ;
    constexpr double NSRadius_3    = NSRadius_IN_CM * NSRadius_IN_CM * NSRadius_IN_CM;
    constexpr double NSRadius_6    = NSRadius_3 * NSRadius_3;
    constexpr double _8_PI_2       = 8.0 * M_PI * M_PI;
    constexpr double _3_C_3        = 3.0 * C * C * C * 1000000;
    constexpr double unitsMoI      = G_TO_KG * CM_TO_M * CM_TO_M;
    
    double mass               = m_Mass * MSOL_TO_G;
    double radius             = NEUTRON_STAR_RADIUS * RSOL_TO_CM;//  in cm
    double initialMagField    = m_PulsarDetails.magneticField; // (in T)
    double initialMagField_G  = initialMagField * TESLA_TO_GAUSS;
    double initialSpinPeriod  = 2 * M_PI / m_PulsarDetails.spinFrequency;
    double magFieldLowerLimit = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField()) * GAUSS_TO_TESLA;    
    double magFieldLowerLimit_G = magFieldLowerLimit * TESLA_TO_GAUSS;                                   
    double momentOfInertia    = m_MomentOfInertia; //* unitsMoI;
    double tau                = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR;                                 

    
    // calculate isolated decay of the magnetic field for a neutron star see Equation 6 in  arXiv:0903.3538v2       
    m_PulsarDetails.magneticField = magFieldLowerLimit + (initialMagField - magFieldLowerLimit) * exp(-p_Stepsize / tau);
    // calculate the spin down rate for isolated neutron stars, see Equation 6 in arxiv:1912.02415   
    double constant_2             = (_8_PI_2 * NSRadius_6) / (_3_C_3 * momentOfInertia);
    double term1                  = magFieldLowerLimit_G * magFieldLowerLimit_G * p_Stepsize;
    double term2                  = tau * magFieldLowerLimit_G * ( m_PulsarDetails.magneticField * TESLA_TO_GAUSS - initialMagField_G);
    double term3                  = (tau / 2.0) * (TESLA_TO_GAUSS * TESLA_TO_GAUSS * (m_PulsarDetails.magneticField * m_PulsarDetails.magneticField) - (initialMagField_G * initialMagField_G));
    double PSquared    = 2 * constant_2 * (term1 - term2 - term3) + (initialSpinPeriod * initialSpinPeriod);
    
    double   P_f                  = std::sqrt(PSquared);
    m_PulsarDetails.spinFrequency = 2 * M_PI / P_f;                                                                        // pulsar spin frequency

    // calculate the spin down rate for isolated neutron stars, see Equation 4 in arXiv:0903.3538v2 (Our version is in SI)      
    double pDotTop               = constant_2 * TESLA_TO_GAUSS * TESLA_TO_GAUSS * m_PulsarDetails.magneticField * m_PulsarDetails.magneticField;
    double pDot                  = pDotTop / P_f;// / pDotBot;
    m_PulsarDetails.spinDownRate = -2 * M_PI * pDot / (P_f * P_f);  

    m_AngularMomentum            = m_PulsarDetails.spinFrequency * momentOfInertia;                                                         // angular momentum of star
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
    constexpr double unitsMoI      = G_TO_KG * CM_TO_M * CM_TO_M;
    
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
    double kappa              = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_KG;      
    double unitsG             = 100 * 100 * 100 / 1000 ; //units of G in cgs 
    //std::cout<<  "Updating pulsar parameters " << std::endl;
    if ((!p_RecycledNS && !p_CommonEnvelope) || (!p_RecycledNS && utils::Compare(p_MassGainPerTimeStep, 0.0) == 0 )) {
        //This is the ''classical'' isolated pulsars
        spinDownIsolatePulsar(p_Stepsize);
        
    }
    else if ((m_PulsarDetails.spinFrequency < 2 * M_PI * 1000.0) && (p_RecycledNS || p_CommonEnvelope) && utils::Compare(p_MassGainPerTimeStep, 0.0) > 0) {
        //This is partially recycling the pulsars in a CE system 
        double mass_kg              = m_Mass * MSOL_TO_KG; //in kg
        double r_m                  = NEUTRON_STAR_RADIUS * RSOL_TO_KM * 1e3; //in meters
        
        double MOI_SI               = m_MomentOfInertia * unitsMoI;
        double angularMomentum_SI   = m_AngularMomentum * unitsMoI;
       
        double newPulsarMagneticField = (initialMagField - magFieldLowerLimit) * exp(-1 * p_MassGainPerTimeStep / 1000 / kappa) + magFieldLowerLimit ;
        

        //Calculate the Alfven radius for an accreting neutron star, see Equation 8 in  arXiv:0903.3538v2
        double mDot          =  p_MassGainPerTimeStep / 1000 / p_Stepsize ;
        double R_M_6         =  r_m * r_m * r_m * r_m * r_m * r_m;
        double B_4           =  newPulsarMagneticField * newPulsarMagneticField * newPulsarMagneticField * newPulsarMagneticField;
        double R_a_top       = 8 * R_M_6 * R_M_6 * B_4;
        double R_a_bot       = mass_kg * mDot * mDot * G;
        double alfvenRadius  = PPOW(R_a_top / R_a_bot, 1.0/7.0);
        
        // calculate the difference in the keplerian angular velocity and surface angular velocity of the neutron star in m - see Equation 2 in 1994MNRAS.269..455J       
        double keplerianVelocityAtAlfvenRadius          = std::sqrt(2.0 * G * mass_kg / alfvenRadius); 
        double keplerianAngularVelocityAtAlfvenRadius   = 4.0 * M_PI * keplerianVelocityAtAlfvenRadius / alfvenRadius;
        double velocityDifference                       = keplerianAngularVelocityAtAlfvenRadius - m_PulsarDetails.spinFrequency;

        // calculate the change in angular momentum due to accretion, see Equation 12 in arXiv:0805.0059/ Equation 8 in arxiv:1912.02415 
        double Jdot                     =  p_Epsilon * velocityDifference * alfvenRadius * alfvenRadius * mDot ;
        
        angularMomentum_SI              = angularMomentum_SI + Jdot * p_Stepsize  ;
        
        if (angularMomentum_SI / MOI_SI > 0) {
            m_PulsarDetails.magneticField   = newPulsarMagneticField  ;
            m_PulsarDetails.spinFrequency   = angularMomentum_SI / MOI_SI;
            m_PulsarDetails.spinDownRate    = Jdot / MOI_SI;
            m_AngularMomentum               = angularMomentum_SI / unitsMoI;
        } 
        else {
            spinDownIsolatePulsar(p_Stepsize);
            //std::cout << m_PulsarDetails.magneticField << std::endl;
        }
        

    }
    else if  (p_RecycledNS && utils::Compare(p_MassGainPerTimeStep, 0.0) == 0 ) {
        //After recycling is completed, pulsar keeps evolving as an isolated pulsar. 
        //std::cout << "Post-recycling evolution " << std::endl;
        spinDownIsolatePulsar(p_Stepsize);

    }
}
