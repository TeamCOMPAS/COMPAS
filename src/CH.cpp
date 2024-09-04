#include "CH.h"

/*
 * Calculate the helium abundance in the core of the star
 * 
 * Currently just a simple linear model from the initial helium abundance to 
 * the maximum helium abundance (assuming that all hydrogen is converted to
 * helium). 
 * 
 * Should be updated to match detailed models.
 *
 * double CalculateHeliumAbundanceCoreOnPhase(p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 *
 * @return                                      Helium abundance in the core (Y_c)
 */
double CH::CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const {

    double heliumAbundanceCoreMax = 1.0 - m_Metallicity;

    double heliumAbundanceCore = ((heliumAbundanceCoreMax - m_initialHeliumAbundance) * p_Tau) + m_initialHeliumAbundance;
    
    return heliumAbundanceCore;
}

/*
 * Calculate the helium abundance at the surface of the star
 * 
 * Since the star is homogeneous, the core and the surface have the same abundances
 * 
 * double CalculateHeliumAbundanceSurfaceOnPhase(p_Tau)
 * 
 * @param   [IN]    p_Tau                      Fraction of main sequence lifetime
 *
 * @return                                     Helium abundance at the surface (Y_s)
 */
double CH::CalculateHeliumAbundanceSurfaceOnPhase(const double p_Tau) const {

    double heliumAbundanceSurface = CalculateHeliumAbundanceCoreOnPhase(p_Tau);

    return heliumAbundanceSurface;
}

/*
 * Calculate the hydrogen abundance in the core of the star
 * 
 * Currently just a simple linear model. Assumes that hydrogen in the core of 
 * the star is burned to helium at a constant rate throughout the lifetime. 
 * 
 * Should be updated to match detailed models.
 *
 * double CalculateHydrogenAbundanceCoreOnPhase(p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 *
 * @return                                      Hydrogen abundance in the core (X_c)
 */
double CH::CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const {

    double hydrogenAbundanceCore = m_initialHydrogenAbundance * (1.0 - p_Tau);

    return hydrogenAbundanceCore;
}

/*
 * Calculate the hydrogen abundance at the surface of the star
 * 
 * Since the star is homogeneous, the core and the surface have the same abundances
 *
 * double CalculateHydrogenAbundanceSurfaceOnPhase(p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 * 
 * @return                                      Hydrogen abundance at the surface (X_s)
 */
double CH::CalculateHydrogenAbundanceSurfaceOnPhase(const double p_Tau) const {

    double hydrogenAbundanceSurface = CalculateHydrogenAbundanceCoreOnPhase(p_Tau);

    return hydrogenAbundanceSurface;
}

/*
 * Calculate the log of the ratio of the lifetimes of a CH star and a normal star of the same mass (log(t_CHE/t_MS))
 * 
 * Polynomial fit derived using IZw18 and IZw18CHE BOoST models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203)
 *
 * double CalculateLogLifetimeRatio(p_Mass)
 * 
 * @param   [IN]    p_Mass                      Mass in Msol
 *
 * @return                                      Log of the ratio of the lifetimes (log(t_CHE/t_MS)) of a CHE star to a normal star of the same mass
 */
double CH::CalculateLogLifetimeRatio(const double p_Mass) const {

    // Define variables
    double x = log10(p_Mass);

    // Coefficients obtained from np.polynomial.Polynomial.fit
    double term1 = -0.15929168474199387;
    double term2 = 1.050069750483549 * x;
    double term3 = -0.8233601359988406 * x * x;
    double term4 = 0.17772610259473764 * x * x * x;

    double logLifetimeRatio = term1 + term2 + term3 + term4;

    return logLifetimeRatio;
}

/*
 * Calculate the ratio of the MS lifetime of a CH star and a normal star of the same mass/metallicity etc
 *
 * Polynomial fit derived using IZw18 and IZw18CHE BOoST models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203)
 *
 * double CalculateLifetimeRatio(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 *
 * @return                                      Ratio of MS lifetimes of a CHE star to a normal star of the same mass
 */
double CH::CalculateLifetimeRatio(const double p_Mass) const {

    double lifetimeRatio = 1.0;

    // If user wants to increase CH MS lifetime, calculate the ratio of CH to MS lifetimes
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {
        double logLifetimeRatio = CalculateLogLifetimeRatio(p_Mass);
        lifetimeRatio = PPOW(10.0, logLifetimeRatio);
    }

    return lifetimeRatio;
}

/*
 * Calculate the ratio of the log luminosities (logL_CH / logL_MS) of a CH star and a normal star of the same mass
 *
 * Polynomial fit derived using IZw18 and IZw18CHE BOoST models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203)
 * 
 * double CalculateLogLuminosityRatio(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       Dimensionless main-sequence age
 *
 * @return                                      Ratio of log luminosities of a CHE star to a normal star of the same mass and age
 */
double CH::CalculateLogLuminosityRatio(const double p_Mass, const double p_Tau) const {
    
    // Define some variables
    double logLuminosityRatio = 1.0;   // log(L_CH) / log(L_MS)

    // If user wants to increase CH MS luminosity, calculate the ratio of CH to MS luminosity
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {
        
        // Polynomial fits to the ratio of logL in models from Szecsi et al., derived using np.polynomial.Polynomial
        double x = log10(p_Mass);

        double term1 =  1.8261540986808193;
        double term2 = -1.1636822407341694 * x;
        double term3 =  0.5876329884434304 * x * x;
        double term4 = -0.10236336828026288 * x * x * x;

        logLuminosityRatio = term1 + term2 + term3 + term4;

        // Make sure we don't reduce the luminosity
        if (logLuminosityRatio < 1.0) {
            logLuminosityRatio = 1.0;
        }

    }

    // Make enhancement grow from 1 to logLuminosityRatio over main-sequence
    logLuminosityRatio = 1.0 + (logLuminosityRatio - 1.0) * p_Tau * p_Tau;
    
    return logLuminosityRatio;
}

/*
 * Calculate the luminosity of a CH star on the main sequence
 * 
 * Start with the luminosity of a 'normal' (non-CH) star of the same mass/metallicity etc.
 * Then multiply that luminosity by the ratio of luminosities 
 * 
 * double CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS)
 *
 * @param   [IN]    p_Time                      Time (after ZAMS) in Myr
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_LZAMS                     Zero Age Main Sequence (ZAMS) Luminosity
 * @return                                      Luminosity of a CH star as a function of time
 */
double CH::CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const {

    // Firstly, calculate the standard main sequence luminosity 
    double MS_Luminosity = MainSequence::CalculateLuminosityOnPhase(p_Time, p_Mass, p_LZAMS);

    // Then, calculate the ratio of the standard luminosity to the CHE luminosity [log(L_CH) / log(L_MS)]
    double logLuminosityRatio = CalculateLogLuminosityRatio(p_Mass, m_Tau);

    // logLuminosityRatio is log(L_CH) / log(L_MS). Hence, log(L_CH) = log(L_MS) * logLuminosityRatio;
    double CH_logLuminosity = log10(MS_Luminosity) * logLuminosityRatio;

    // Exponentiate the logLuminosity to get the luminosity L_CH = 10**log(L_CH)
    double CH_Luminosity = PPOW(10.0, CH_logLuminosity);

    return CH_Luminosity;

}

/*
 * Calculate luminosity at the end of the CH main sequence
 *
 * Start with the luminosity of a 'normal' (non-CH) star of the same mass/metallicity etc.
 * Then multiply that luminosity by the ratio of luminosities of a CH star to a normal star
 *
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the Main Sequence in Lsol
 */
double CH::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {
    // Firstly, calculate the luminosity at the end of the main sequence for a normal star
    double TAMS_Luminosity = MainSequence::CalculateLuminosityAtPhaseEnd(p_Mass);

    // Then, calculate the ratio of the standard luminosity to the CHE luminosity [log(L_CH) / log(L_MS)]
    double logLuminosityRatio = CalculateLogLuminosityRatio(p_Mass, 1.0); // Tau = 1.0 at phase end by construction

    // logLuminosityRatio is log(L_CH) / log(L_MS). Hence, log(L_CH) = log(L_MS) * logLuminosityRatio;
    double CH_logLuminosity = log10(TAMS_Luminosity) * logLuminosityRatio;

    // Exponentiate the logLuminosity to get the luminosity L_CH = 10**log(L_CH)
    double CH_Luminosity = PPOW(10.0, CH_logLuminosity);

    // Finally, multiply the luminosity by the ratio above to get the CH star luminosity
    return CH_Luminosity;
}

/*
 * Calculate timescales in units of Myr
 *
 * Timescales depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * void CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_Timescales            Timescales
 */
void CH::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    // Calculate the standard lifetimes in the normal way
    timescales(tBGB)   = CalculateLifetimeToBGB(p_Mass);
    timescales(tMS)    = CalculateLifetimeOnPhase(p_Mass, timescales(tBGB));

    // Calculate the ratio of the lifetime of a CH star to a normal MS star
    double lifetimeRatio = CalculateLifetimeRatio(p_Mass);

    // Multiply the main sequence lifetime by this ratio
    timescales(tBGB) *= lifetimeRatio;
    timescales(tMS)  *= lifetimeRatio;

#undef timescales
}

/*
 * Recalculates the star's age after mass loss
 *
 * Hurley et al. 2000, section 7.1
 *
 * Modifies attribute m_Age
 *
 *
 * UpdateAgeAfterMassLoss()
 *
 */
void CH::UpdateAgeAfterMassLoss() {

    double tMS       = m_Timescales[static_cast<int>(TIMESCALE::tMS)];

    // Recalculate tBGB and tMS
    double tBGBprime = CalculateLifetimeToBGB(m_Mass);
    double tMSprime  = MainSequence::CalculateLifetimeOnPhase(m_Mass, tBGBprime);

    // Calculate the ratio of the lifetime of a CH star to a normal MS star
    double lifetimeRatio = CalculateLifetimeRatio(m_Mass);

    tBGBprime *= lifetimeRatio;
    tMSprime  *= lifetimeRatio;

    m_Age *= tMSprime / tMS;
}

/*
 * Calculate the weight for the OB and WR mass loss prescriptions
 *
 * According to the prescription described in Yoon et al. 2006 (and Szecsi et al. 2015)
 * Use OB mass loss rates until Y1 = 0.55, WR mass loss rates above Y2 = 0.7, and 
 * linearly interpolate for Y1 < Ys < Y2 where Ys is the surface helium fraction.
 * 
 * Since we don't have Ys by default in our models, and detailed models show Ys
 * rises ~linearly from 0.2 to 1.0 across the main-sequence, here we simply 
 * use tau = t / tMS as a proxy for Ys.
 *
 * CalculateMassLossRateWeightOB()
 *
 * @param   [IN]        p_HeliumAbundanceSurface           Surface helium abundance Ys
 *
 * @return                                   weight for OB mass loss rate
 */
double CH::CalculateMassLossRateWeightOB(const double p_HeliumAbundanceSurface){

    // Define variables
    double weight = 0.0;
    const double y1 = 0.55;
    const double y2 = 0.70;
    const double weight_bot = y2 - y1;

    double Ys = p_HeliumAbundanceSurface; // Get surface helium abundance

    // Calculate the weight as a function of Ys according to the prescription from Yoon et al. 2006
    if (Ys < y1){
        weight = 1.0;
    }
    else if (Ys > y2){
        weight = 0.0;
    }
    else{
        weight = (y2 - Ys) / weight_bot;
    }

    return weight;
}

/*
 * Calculate the dominant mass loss mechanism and associated rate for the star 
 * at the current evolutionary phase
 * 
 * According to Belczynski et al. 2010 prescription - based on implementation in StarTrack
 *
 * Modifications for CH stars
 * 
 * double CalculateMassLossRateBelczynski2010()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double CH::CalculateMassLossRateBelczynski2010() {

    // Define variables
    double Mdot    = 0.0;
    double Mdot_OB = 0.0;
    double Mdot_WR = 0.0;
    double weight  = 1.0;           // Initialised to 1.0 to allow us to use the OB mass loss rate by default 

    // Calculate OB mass loss rate according to the Vink et al. formalism
    Mdot_OB = BaseStar::CalculateMassLossRateBelczynski2010();

    // If user wants to transition between OB and WR mass loss rates
    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()){

        // Calculate WR mass loss rate
        Mdot_WR = BaseStar::CalculateMassLossRateWolfRayetZDependent(0.0);

        // Calculate weight for combining these into total mass-loss rate
        weight = CalculateMassLossRateWeightOB(m_HeliumAbundanceSurface);

    }

    // Set dominant mass loss rate
    m_DominantMassLossRate = weight == 0 ? MASS_LOSS_TYPE::WR : MASS_LOSS_TYPE::OB;   

    // Finally, combine each of these prescriptions according to the weight
    Mdot = (weight * Mdot_OB) + ((1.0 - weight) * Mdot_WR);

    // Enhance mass loss rate due to rotation
    Mdot *= CalculateMassLossRateEnhancementRotation();

    return Mdot;
}

/*
 * Calculate the dominant mass loss mechanism and associated rate for the star 
 * at the current evolutionary phase
 * 
 * According to Merritt et al. 2024 prescription
 *
 * Modifications for CH stars
 * 
 * double CalculateMassLossRateMerritt2024()
 * 
 * @return                                      Mass loss rate in Msol per year
 */
double CH::CalculateMassLossRateMerritt2024() {
    
    // Define variables
    double Mdot    = 0.0;
    double Mdot_OB = 0.0;
    double Mdot_WR = 0.0;
    double weight  = 1.0;           // Initialised to 1.0 to allow us to use the OB mass loss rate by default
    
    // Calculate OB mass loss rate according to the chosen prescription
    Mdot_OB = BaseStar::CalculateMassLossRateOB(OPTIONS->OBMassLossPrescription());  

    // If user wants to transition between OB and WR mass loss rates
    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()){

        // Here we are going to pretend that this CH star is an HeMS star by
        // cloning it, so that we can ask it what its mass loss rate would be if it were
        // a HeMS star
        //HeMS *clone = HeMS::Clone(static_cast<BaseStar&>(*this), OBJECT_PERSISTENCE::EPHEMERAL, false);        // Do not initialise so that we can use same mass, luminosity, radius etc
        HeMS *clone = HeMS::Clone((HeMS&)static_cast<const CH&>(*this), OBJECT_PERSISTENCE::EPHEMERAL, false);   // Do not initialise so that we can use same mass, luminosity, radius etc
        
        Mdot_WR = clone->CalculateMassLossRateMerritt2024();       // Calculate WR mass loss rate              
        delete clone; clone = nullptr;    

        // Calculate weight for combining these into total mass-loss rate
        weight = CalculateMassLossRateWeightOB(m_HeliumAbundanceSurface);

    }

    // Set dominant mass loss rate
    m_DominantMassLossRate = weight == 0 ? MASS_LOSS_TYPE::WR : MASS_LOSS_TYPE::OB;   

    // Finally, combine each of these prescriptions according to the weight
    Mdot = (weight * Mdot_OB) + ((1.0 - weight) * Mdot_WR);

    // Enhance mass loss rate due to rotation
    Mdot *= CalculateMassLossRateEnhancementRotation();

    return Mdot;
}

STELLAR_TYPE CH::EvolveToNextPhase() {

    STELLAR_TYPE stellarType = STELLAR_TYPE::MS_GT_07;

    if (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]) {           // evolving off because of age?
        stellarType = STELLAR_TYPE::MS_GT_07;                               // no - must have spun down - evolve as MS star now
        m_CHE       = false;                                                // evolved CH->MS
    }
    else {                                                                  // yes
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                   // evolve as HeMS star now
        m_Age       = 0.0;                                                  // can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        m_Tau       = 0.0;
        m_CHE       = true;                                                 // stayed on MS as CH
    }

    return stellarType;
}

