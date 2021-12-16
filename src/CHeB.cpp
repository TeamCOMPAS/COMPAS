#include "CHeB.h"
#include "EAGB.h"
#include "HeMS.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate timescales in units of Myr
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster.  This function isn't
 * called too often, but the pattern is the same for others that are called many, many times.
 *
 *
 * void CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_Timescales            Timescales - calculated here
 */
void CHeB::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]                      // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]                  // for convenience and readability - undefined at end of function

    GiantBranch::CalculateTimescales(p_Mass, p_Timescales);                             // calculate common values

    timescales(tHe)     = CalculateLifetimeOnPhase(p_Mass);
	timescales(tau_BL)  = CalculateLifetimeOnBluePhase(p_Mass);

    // JR: The blue loop on CHeB can be 0-length/duration (see Hurley et al., 2000, section 5.3, 
    // particularly eq 58 and beyond).  COMPAS does not allow for a 0-length blue loop - some of 
    // the equations used (e.g. to calculate Radius) result in nan or inf.  As a temporary workaround 
    // until we work out how to skip the blue loop (when it is 0-length) we will set the length of a 
    // 0-length blue loop to the absolute minimum timestep (currently 100 seconds).
    //
    // Note that this works around a long-standing problem, which was worked around in legacy COMPAS
    // by the following code in calculateBluePhaseFBL() in star.cpp:
	//
    //    if(brackets ==0){brackets = 1e-12;}  //If zero gives R=NaN Coen Neijssel 10-01-2017
    //
    // The workaround implemented here is closer to the source of the problem (the blue loop does not
    // actually exist for some stars), and maybe a bit more meaningful (we're just using a very short 
    // duration blue loop instead of a no duration (non-existent) one)

    if (timescales(tau_BL) <= 0.0) timescales(tau_BL) = ABSOLUTE_MINIMUM_TIMESTEP;      // don't use utils::Compare() here

    // Calculate the relative age at the start of the blue phase of Core Helium Burning
    // Hurley et al. 2000, just before eq 59
    // Naturally clamped to [0, 1]
	timescales(tauX_BL) = (utils::Compare(p_Mass, massCutoffs(MHeF)) >= 0 && utils::Compare(p_Mass, massCutoffs(MFGB)) < 0)
                                    ? 1.0 - timescales(tau_BL)                          // intermediate mass stars
                                    : 0.0;                                              // low and high mass stars

    // Calculate the relative age at the end of the blue phase of Core Helium Burning
    // Hurley et al. 2000, just before eq 64
    // Naturally clamped to [0, 1]
	timescales(tauY_BL) = utils::Compare(p_Mass, massCutoffs(MFGB)) >= 0
                                    ? timescales(tau_BL)                                // high mass stars
                                    : 1.0;                                              // intermediate and low mass stars

#undef massCutoffs
#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                LAMBDA CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the common envelope lambda parameter using the "Dewi" prescription.
 *
 * Fit from Appendix A in Claeys+2014, based on the method described in Dewi & Tauris 2000
 * arXiv:1401.2895 for Claeys+2014
 * arXiv:0007034 for Dewi and Tauris 2000
 *
 * ALEJANDRO - 17/05/2017 - Not fully tested, nor fully coded, nor fully trusted. Any other \lambda prescription is personally preferred.
 * Missing (A.6),(A.7),(A.8),(A.9),(A.10),(A.11) and (A.12), which have to do with ionization energy
 *
 *
 * double CalculateLambdaDewi()
 *
 * @return                                      Dewi lambda for use in common envelope
 */
double CHeB::CalculateLambdaDewi() const {

    double lambda3 = std::min(-0.9, 0.58 + (0.75 * log10(m_Mass))) - (0.08 * log10(m_Luminosity));                          // (A.4) Claeys+2014
    double lambda1 = std::min(lambda3, std::min(0.8, 1.25 - (0.15 * log10(m_Luminosity))));                                 // (A.5) Top, Claeys+2014
	double lambda2 = 0.42 * PPOW(m_RZAMS / m_Radius, 0.4);                                                                  // (A.2) Claeys+2014          // RTW - Consider replacing this with a 2/5 root function (somehow) to avoid NaNs if the base is negative
	double envMass = utils::Compare(m_CoreMass, 0.0) > 0 && utils::Compare(m_Mass, m_CoreMass) > 0 ? m_Mass - m_CoreMass : 0.0;

    double lambdaCE;

         if (utils::Compare(envMass, 1.0) >= 0) lambdaCE = 2.0 * lambda1;                                                   // (A.1) Bottom, Claeys+2014
	else if (utils::Compare(envMass, 0.0) >  0) lambdaCE = 2.0 * (lambda2 + (std::sqrt(envMass) * (lambda1 - lambda2)));         // (A.1) Mid, Claeys+2014
	else                                        lambdaCE = 2.0 * lambda2;	                                                // (A.1) Top, Claeys+2014

	return	lambdaCE;
}


/*
 * Calculate the common envelope lambda parameter using the "Nanjing" prescription
 * from X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010) as implemented in STARTRACK
 *
 * This implementation adapted from the STARTRACK implementation (STARTRACK courtesy Chris Belczynski)
 *
 *
 * JR: todo: the coefficients and factors here are hard-coded until I figure out an
 * efficient way of putting them in constants.h.  Because they are indexed by a few
 * things: stellar type, metallicity and ZAMS mass the easiet thing would be to  put
 * them in a map - but because they can be re-calculated at every timestep the hashing
 * overhead becomes a performance concern.  Vectors are a good alternative, but I need
 * to figure out how best to structure them for reasonable (and inuitive) access.
 *
 * This function good for CHeB stars.
 *
 *
 * double CalculateLambdaNanjing()
 *
 * @return                                      Nanjing lambda for use in common envelope
 */
double CHeB::CalculateLambdaNanjing() const {

	DBL_VECTOR maxBG    = {};                                                       // [0] = maxB, [1] = maxG
	DBL_VECTOR lambdaBG = {};                                                       // [0] = lambdaB, [1] = lambdaG
	DBL_VECTOR a        = {};                                                       // 0..5 a_coefficients
	DBL_VECTOR b        = {};                                                       // 0..5 b_coefficients

    if (utils::Compare(m_Metallicity, LAMBDA_NANJING_ZLIMIT) > 0) {                 // Z>0.5 Zsun: popI
        if (utils::Compare(m_MZAMS, 1.5) < 0) {
            maxBG = { 2.5, 1.5 };
            if (utils::Compare(m_Radius, 200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else {
                a = { 46.00978, -298.64993, 727.40936, -607.66797, 0.0, 0.0 };
                b = { 63.61259, -399.89494, 959.62055, -795.20699, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
                 if (utils::Compare(m_Radius, 340.0) > 0)                                     lambdaBG = { 3.589970, 0.514132 };
            else if (utils::Compare(m_Radius, 8.5) > 0 && utils::Compare(m_Radius, 60.0) < 0) lambdaBG = { 3.0, 1.2 };
            else {
                a = { 34.41826, -6.65259, 0.43823, -0.00953, 0.0, 0.0 };
                b = { 13.66058, -2.48031, 0.15275, -0.00303, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 3.5) < 0) {
            maxBG = { 500.0, 10.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 116.935557, 0.848808 };
            else {
                maxBG = { 2.5, 1.5 };
                a     = { -42.98513, 7.90134, -0.54646, 0.01863,  3.13101E-04, 2.07468E-06 };
                b     = { -6.73842 , 1.06656, -0.05344, 0.00116, -9.34446E-06, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 4.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 410.0) > 0) lambdaBG = { 52.980056, 1.109736 };
            else {
                maxBG = { 2.5, 1.5 };
                a     = { -7.3098 , 0.56647, -0.01176, 7.90112E-05, 0.0, 0.0 };
                b     = { -3.80455, 0.29308, -0.00603, 4.00471E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 5.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 430.0) > 0) lambdaBG = { 109.593522, 1.324248 };
            else {
                a = { -9.93647, 0.42831, -0.00544, 2.25848E-05, 0.0, 0.0 };
                b = { -5.33279, 0.22728, -0.00285, 1.16408E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 6.5) < 0) {
            maxBG = { 25.5, 5.0 };
            if (utils::Compare(m_Radius, 440.0) > 0) lambdaBG = { 16.279603, 1.352166 };
            else {
                a = { 13.91465, -0.55579, 0.00809, -4.94872E-05, 1.08899E-07, 0.0 };
                b = {  7.68768, -0.30723, 0.00445, -2.70449E-05, 5.89712E-08, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 7.5) < 0) {
            maxBG = { 9.0, 3.0 };
            if (utils::Compare(m_Radius, 420.0) > 0) lambdaBG = { 5.133959, 1.004036 };
            else {
                a = { 4.12387, -0.12979, 0.00153    , -7.43227E-06, 1.29418E-08, 0.0 };
                b = { 2.18952, -0.06892, 8.00936E-04, -3.78092E-06, 6.3482E-09 , 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 8.5) < 0) {
            maxBG = { 7.0, 3.0 };
            if (utils::Compare(m_Radius, 490.0) > 0) lambdaBG = { 4.342985, 0.934659 };
            else {
                maxBG = { 1.0, 0.5 };
                a     = { -3.89189, 0.19378, -0.0032 , 2.39504E-05, -8.28959E-08, 1.07843E-10 };
                b     = { -2.24354, 0.10918, -0.00179, 1.33244E-05, -4.57829E-08, 5.90313E-11 };
            }
        }
        else if (utils::Compare(m_MZAMS, 9.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 530.0) > 0) lambdaBG = { 2.441672, 0.702310 };
            else {
                a = { 0.86369, -0.00995,  4.80837E-05, -6.10454E-08, -2.79504E-12, 0.0 };
                b = { -0.7299,  0.0391 , -5.78132E-04,  3.7072E-06 , -1.07036E-08, 1.14833E-11 };
            }
        }
        else if (utils::Compare(m_MZAMS, 11.0) < 0) {
            maxBG = { 3.0, 1.5 };
            if (utils::Compare(m_Radius, 600.0) > 0) lambdaBG = { 1.842314, 0.593854 };
            else {
                a = { 0.74233, -0.00623, 2.04197E-05, -1.30388E-08, 0.0, 0.0 };
                b = { 0.36742, -0.00344, 1.27838E-05, -1.0722E-08 , 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 13.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 850.0) > 0) lambdaBG = { 0.392470, 0.176660 };
            else if (utils::Compare(m_Radius, 0.0) > 0 && utils::Compare(m_Radius, 350.0) <= 0) {
                a = { 1.28593, -0.02209, 1.79764E-04, -6.21556E-07, 7.59444E-10, 0.0 };
                b = { 0.68544, -0.01394, 1.20845E-04, -4.29071E-07, 5.29169E-10, 0.0 };
            }
            else if (utils::Compare(m_Radius, 350.0) > 0 && utils::Compare(m_Radius, 600.0) <= 0) {
                a = { -11.99537,  0.0992, -2.8981E-04,  3.62751E-07, -1.65585E-10, 0.0 };
                b = {   0.46156, -0.0066,  3.9625E-05, -9.98667E-08, -8.84134E-11, 0.0 };
            }
            else {
                a = { -58.03732, 0.23633, -3.20535E-04, 1.45129E-07, 0.0, 0.0 };
                b = { -15.11672, 0.06331, -8.81542E-05, 4.0982E-08 , 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 15.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 1000.0) > 0)                                      lambdaBG = { 0.414200, 0.189008 };
            else if (utils::Compare(m_Radius, 69.0) > 0 && utils::Compare(m_Radius, 126.0) < 0) lambdaBG = { 0.5 - (m_Radius * 8.77E-04), 0.18 };
            else {
                a = { 1.12889, -0.00901, 3.04077E-05, -4.31964E-08, 2.14545E-11, 0.0 };
                b = { 0.568  , -0.0047 , 1.57818E-05, -2.21207E-08, 1.08472E-11, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1050.0) > 0) lambdaBG = { 0.2, 0.1 };
            else {
                a = { 0.84143, -0.00576, 1.68854E-05, -2.0827E-08 , 8.97813E-12, 0.0 };
                b = { 0.36014, -0.00254, 7.49639E-06, -9.20103E-09, 3.93828E-12, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else {
                a = { 0.48724, -0.00177   , 2.60254E-06, -1.25824E-09, 0.0, 0.0 };
                b = { 0.22693, -8.7678E-04, 1.28852E-06, -6.12912E-10, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 75.0) < 0) {
            maxBG = { 1.0, 0.5 };
            a     = { 0.31321, -7.50384E-04, 5.38545E-07, -1.16946E-10, 0.0, 0.0 };
            b     = { 0.159  , -3.94451E-04, 2.88452E-07, -6.35132E-11, 0.0, 0.0 };
        }
        else {
            maxBG = { 1.0, 0.5 };
            a     = { 0.376 , -0.0018 , 2.81083E-06, -1.67386E-09, 3.35056E-13, 0.0 };
            b     = { 0.2466, -0.00121, 1.89029E-06, -1.12066E-09, 2.2258E-13 , 0.0 };
        }
    }
    else {                                                                  // Z<=0.5 Zsun: popI and popII
        if (utils::Compare(m_MZAMS, 1.5) < 0) {
            maxBG = { 2.0, 1.5 };
            if (utils::Compare(m_Radius, 160.0) > 0) lambdaBG = { 0.05, 0.05 };
            else {
                a = { 0.37294, -0.05825, 0.00375, -7.59191E-05, 0.0, 0.0 };
                b = { 0.24816, -0.04102, 0.0028 , -6.20419E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
                 if (utils::Compare(m_Radius, 350.0) > 0)                                     lambdaBG = { 2.868539, 0.389991 };
            else if (utils::Compare(m_Radius, 6.0) > 0 && utils::Compare(m_Radius, 50.0) < 0) lambdaBG = { 0.8, 0.35 };
            else {
                a = { -103.92538, 25.37325, -2.03273, 0.0543 , 0.0, 0.0 };
                b = {  -56.03478, 13.6749 , -1.09533, 0.02925, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 3.5) < 0) {
            maxBG = { 600.0, 2.0 };
                 if (utils::Compare(m_Radius, 400.0) > 0)                               lambdaBG = { 398.126442, 0.648560 };
            else if (utils::Compare(m_Radius, 36.0) > 0 && utils::Compare(m_Radius, 53.0) < 0) lambdaBG = { 1.0, 1.0 };
            else {
                a = { -12.40832, 1.59021, -0.06494, 8.69587E-04, 0.0, 0.0 };
                b = { -6.47476 , 0.8328 , -0.03412, 4.58399E-04, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 4.5) < 0) {
            maxBG = { 600.0, 2.0 };
                 if (utils::Compare(m_Radius, 410.0) > 0)                                      lambdaBG = { 91.579093, 1.032432 };
            else if (utils::Compare(m_Radius, 19.0) > 0 && utils::Compare(m_Radius, 85.0) < 0) lambdaBG = { 0.255, 0.115 };
            else {
                a = { -5.89253, 0.54296, -0.01527, 1.38354E-04, 0.0, 0.0 };
                b = { -3.21299, 0.29583, -0.00833, 7.55646E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 5.5) < 0) {
            maxBG = { 10.0, 3.0 };
                 if (utils::Compare(m_Radius, 320.0) > 0)                                       lambdaBG = { 7.618019, 1.257919 };
            else if (utils::Compare(m_Radius, 85.0) > 0 && utils::Compare(m_Radius, 120.0) < 0) lambdaBG = { 0.4, 0.1 };
            else {
                a = { -0.67176, 0.07708, -0.00175   , 1.1991E-05 , 0.0, 0.0 };
                b = { -0.38561, 0.0427 , -9.6948E-04, 6.64455E-06, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 6.5) < 0) {
            maxBG = { 4.0, 1.5 };
                 if (utils::Compare(m_Radius, 330.0) > 0)                                        lambdaBG = { 2.390575, 0.772091 };
            else if (utils::Compare(m_Radius, 115.0) > 0 && utils::Compare(m_Radius, 165.0) < 0) lambdaBG = { 0.2, 0.1 };
            else {
                a = { 0.30941, 0.00965, -2.31975E-04, 1.26273E-06, 0.0, 0.0 };
                b = { 0.14576, 0.00562, -1.30273E-04, 7.06459E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 7.5) < 0) {
            maxBG = { 2.5, 1.0 };
                 if (utils::Compare(m_Radius, 360.0) > 0)                                        lambdaBG = { 1.878174, 0.646353 };
            else if (utils::Compare(m_Radius, 150.0) > 0 && utils::Compare(m_Radius, 210.0) < 0) lambdaBG = { 0.2, 0.1 };
            else {
                a = { 0.44862, 0.00234, -9.23152E-05, 4.67797E-07, 0.0, 0.0 };
                b = { 0.21873, 0.00154, -5.18806E-05, 2.60283E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 8.5) < 0) {
            maxBG = { 2.0, 1.0 };
                 if (utils::Compare(m_Radius, 400.0) > 0)                                        lambdaBG = { 1.517662, 0.553169 };
            else if (utils::Compare(m_Radius, 190.0) > 0 && utils::Compare(m_Radius, 260.0) < 0) lambdaBG = { 0.2, 0.1 };
            else {
                a = { 0.50221, -3.19021E-04, -3.81717E-05, 1.80726E-07, 0.0, 0.0 };
                b = { 0.24748, -9.9338E-05 , -1.99272E-05, 9.47504E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 9.5) < 0) {
            maxBG = { 1.6, 1.0 };
                 if (utils::Compare(m_Radius, 440.0) > 0)                                        lambdaBG = { 1.136394, 0.478963 };
            else if (utils::Compare(m_Radius, 180.0) > 0 && utils::Compare(m_Radius, 300.0) < 0) lambdaBG = { 0.15, 0.08 };
            else {
                a = { 0.39342, 0.00259    , -4.97778E-05, 1.69533E-07, 0.0, 0.0 };
                b = { 0.20796, 6.62921E-04, -1.84663E-05, 6.58983E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 11.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 500.0) > 0) lambdaBG = { 1.068300, 0.424706 };
            else {
                a = { 0.75746, -0.00852, 3.51646E-05, -4.57725E-08, 0.0, 0.0 };
                b = { 0.35355, -0.00388, 1.56573E-05, -1.98173E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 13.0) < 0) {
            maxBG = { 1.6, 1.0 };
                 if (utils::Compare(m_Radius, 600.0) > 0)                                        lambdaBG = { 0.537155, 0.211105 };
            else if (utils::Compare(m_Radius, 200.0) > 0 && utils::Compare(m_Radius, 410.0) < 0) lambdaBG = { 0.08, 0.05 };
            else {
                a = { 0.85249, -0.00861, 2.99246E-05, -3.21416E-08, 0.0, 0.0 };
                b = { 0.37188, -0.00365, 1.24944E-05, -1.32388E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 15.0) < 0) {
            maxBG = { 1.6, 1.0 };
                 if (utils::Compare(m_Radius, 650.0) > 0)                                        lambdaBG = { 0.3, 0.160696 };
            else if (utils::Compare(m_Radius, 250.0) > 0 && utils::Compare(m_Radius, 490.0) < 0) lambdaBG = { 0.06, 0.05 };
            else {
                a = { 0.85271, -0.00793, 2.5174E-05 , -2.4456E-08 , 0.0, 0.0 };
                b = { 0.36163, -0.00328, 1.03119E-05, -9.92712E-09, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 750.0) > 0)                                        lambdaBG = { 0.5, 0.204092 };
            else if (utils::Compare(m_Radius, 200.0) > 0 && utils::Compare(m_Radius, 570.0) < 0) lambdaBG = { 0.1, 0.05 };
            else {
                a = { 0.83254, -0.00696, 1.9597E-05 , -1.67985E-08, 0.0, 0.0 };
                b = { 0.34196, -0.0028 , 7.82865E-06, -6.66684E-09, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 900.0) > 0)                                        lambdaBG = { 0.2, 0.107914 };
            else if (utils::Compare(m_Radius, 230.0) > 0 && utils::Compare(m_Radius, 755.0) < 0) lambdaBG = { 0.1, 0.05 };
            else {
                a = { 0.69746, -0.0043 , 8.97312E-06, -5.83402E-09, 0.0, 0.0 };
                b = { 0.26691, -0.00161, 3.3378E-06 , -2.1555E-09 , 0.0, 0.0 };
            }
        }
        else if (utils::Compare(m_MZAMS, 75.0) < 0) {
            maxBG = { 20.0, 3.0 };
            a     = { 0.821  , -0.00669, 1.57665E-05, -1.3427E-08 , 3.74204E-12, 0.0 };
            b     = { 0.49287, -0.00439, 1.06766E-05, -9.22015E-09, 2.58926E-12, 0.0 };
        }
        else {
            maxBG = { 4.0, 2.0 };
            a     = { 1.25332, -0.02065, 1.3107E-04 , -3.67006E-07, 4.58792E-10, -2.09069E-13 };
            b     = { 0.81716, -0.01436, 9.31143E-05, -2.6539E-07 , 3.30773E-10, -1.51207E-13 };
        }
    }

    if (lambdaBG.empty()) {                                                 // calculate lambda B & G - not approximated by hand
        if (utils::Compare(m_Metallicity, LAMBDA_NANJING_ZLIMIT) > 0 && utils::Compare(m_MZAMS, 1.5) < 0) {
            double x  = (m_Mass - m_CoreMass) / m_Mass;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x3 * x2;

            double y1 = a[0] + (a[1] * x) + (a[2] * x2) + (a[3] * x3) + (a[4] * x4) + (a[5] * x5);
            double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

            lambdaBG = { 1.0 / y1, 1.0 / y2 };
        }
        else {
            double x  = m_Radius;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x3 * x2;

            double y1 = a[0] + (a[1] * x) + (a[2] * x2) + (a[3] * x3) + (a[4] * x4) + (a[5] * x5);
            double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

            lambdaBG = { y1, y2 };
        }
    }

    // Limit lambda to some 'reasonable' range
    lambdaBG[0] = std::min(std::max(0.05, lambdaBG[0]), maxBG[0]);                      // clamp lambda B to [0.05, maxB]
    lambdaBG[1] = std::min(std::max(0.05, lambdaBG[1]), maxBG[1]);                      // clamp lambda G to [0.05, maxG]

    // Calculate lambda as some combination of lambda_b and lambda_g by
    // lambda = alpha_th • lambda_b    +  (1-alpha_th) • lambda_g
    // Note that this is different from STARTRACK
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaBG[0]) + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaBG[1]);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate minimum luminosity during Core Helium Burning
 * for intermediate mass (IM) stars (JR: is there a check for IM?)
 *
 * Hurley et al. 2000, eq 51
 *
 *
 * double CalculateMinimumLuminosityOnPhase(const double      p_Mass,
 *                                          const double      p_Alpha1,
 *                                          const double      p_MHeF,
 *                                          const double      p_MFGB,
 *                                          const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurly et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_MFGB                      Maximum initial mass for which helium ignites on the First Giant Branch
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Minimum luminosity during Core Helium Burning in Lsol
 *
 * p_MHef, p_MFGB and p_BnCoefficients passed as parameters so function can be declared static
 */
double CHeB::CalculateMinimumLuminosityOnPhase(const double      p_Mass,
                                               const double      p_Alpha1,
                                               const double      p_MHeF,
                                               const double      p_MFGB,
                                               const DBL_VECTOR &p_BnCoefficients) const {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double LHeI = GiantBranch::CalculateLuminosityAtHeIgnition_Static(p_Mass, p_Alpha1, p_MHeF, p_BnCoefficients);
    double c    = (b[17] / PPOW(p_MFGB, 0.1)) + (((b[16] * b[17]) - b[14]) / (PPOW(p_MFGB, (b[15] + 0.1))));

    return  LHeI * ((b[14] + (c * PPOW(p_Mass, (b[15] + 0.1)))) / (b[16] + PPOW(p_Mass, b[15])));

#undef b
}


/*
 * Calculate luminosity at the start of the blue phase of Core Helium Burning
 *
 * Hurley et al. 2000, eq 59
 *
 *
 * double CalculateLuminosityAtBluePhaseStart(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the start of the blue phase of core helium burning in Lsol
 */
double CHeB::CalculateLuminosityAtBluePhaseStart(const double p_Mass) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double Lx;
    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        Lx = GiantBranch::CalculateLuminosityOnZAHB_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
    }
    else if (utils::Compare(p_Mass, massCutoffs(MFGB)) < 0) {
        Lx = CalculateMinimumLuminosityOnPhase(p_Mass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_BnCoefficients);
    }
    else {
        Lx = GiantBranch::CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);
    }

    return Lx;

#undef massCutoffs
}


/*
 * Calculate luminosity at the end of the blue phase of Core Helium Burning
 *
 * Hurley et al. 2000, eq 61 (see discussion just before eq 64)
 *
 *
 * double CalculateLuminosityAtBluePhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the blue phase of Core Helium Burning in Lsol
 */
double CHeB::CalculateLuminosityAtBluePhaseEnd(const double p_Mass) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double Ly;

    double tx = timescales(tauX_BL);
    double ty = timescales(tauY_BL);

    double Lx = CalculateLuminosityAtBluePhaseStart(p_Mass);

    if (utils::Compare(ty, tx) >= 0) {
        double Rx      = CalculateRadiusAtBluePhaseStart(p_Mass);
        double RMinHe  = CalculateMinimumRadiusOnPhase_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
        double epsilon = std::min(2.5, std::max(0.4, RMinHe / Rx));
        double lambda  = (utils::Compare(ty, tx) == 0) ? 0.0 : PPOW(((ty - tx) / (1.0 - tx)), epsilon);     // JR: tx can be 1.0 here - if so, lambda = 0.0
        Ly             = Lx * PPOW(CalculateLuminosityAtBAGB(p_Mass) / Lx, lambda);
    }
    else {
        // pow() is slow - use multiplication
        double tmp         = (tx - ty) / tx;                                                        // JR: tx cannot be 0.0 here - so safe (tx > ty, ty = [0, 1])
        double lambdaPrime = tmp * tmp * tmp;
        Ly                 = Lx * PPOW(GiantBranch::CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients) / Lx, lambdaPrime);
    }

    return Ly;

#undef massCutoffs
#undef timescales
}


/*
 * Calculate luminosity during Core Helium Burning
 *
 * Hurley et al. 2000, eqs 61, 62 & 63
 *
 *
 * double CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       CHeB relative age
 * @return                                      Luminosity during Core Helium Burning in Lsol
 */
double CHeB::CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double lCHeB;

    double tx = timescales(tauX_BL);                                                                                                        // 0 for LM and HM stars, non-zero for IM stars
    double Lx = CalculateLuminosityAtBluePhaseStart(p_Mass);

    if (utils::Compare(p_Tau, tx) >= 0) {
        double Rx      = CalculateRadiusAtBluePhaseStart(p_Mass);
        double RmHe    = CalculateMinimumRadiusOnPhase_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
        double LBAGB   = CalculateLuminosityAtBAGB(p_Mass);
        double epsilon = std::min(2.5, std::max(0.4, (RmHe / Rx)));
        double lambda  = (utils::Compare(p_Tau, tx) == 0) ? 0.0 : PPOW((p_Tau - tx) / (1.0 - tx), epsilon);                                  // JR: tx can be 1.0 here - if so, lambda = 0.0
        lCHeB          = Lx * PPOW(LBAGB / Lx, lambda);
    }
    else {
        double LHeI        = GiantBranch::CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);    // pow() is slow - use multiplication
        double tmp         = (tx - p_Tau) / tx;                                                                                             // JR: tx cannot be 0.0 here, so safe (tx > tau, tau = [0, 1])
        double lambdaPrime = tmp * tmp * tmp;
        lCHeB              = Lx * PPOW((LHeI / Lx), lambdaPrime);
    }

    return lCHeB;

#undef masCutoffs
#undef timescales
}


/*
 * Calculate luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity()
 *
 * @return                                      Luminosity of remnant core in Lsol
 */
double CHeB::CalculateRemnantLuminosity() const {
    return HeMS::CalculateLuminosityOnPhase_Static(m_CoreMass, CalculateTauOnPhase());
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


double CHeB::CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function
    return EAGB::CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, massCutoffs(MHeF), m_BnCoefficients);
#undef massCutoffs
}


/*
 * Calculate minimum radius during Core Helium Burning (on the blue loop)
 *
 *  Hurley et al. 2000, eq 55
 *
 *
 * double CalculateMinimumRadiusOnPhase_Static(const double      p_Mass,
 *                                             const double      p_CoreMass,
 *                                             const double      p_Alpha1,
 *                                             const double      p_MHeF,
 *                                             const double      p_MFGB,
 *                                             const double      p_MinimumLuminosityOnPhase,
 *                                             const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_CoreMass                  Core Mass in Msol
 * @param   [IN]    p_Alpha1                    alpha1 in Hurly et al. 2000 (just after eq 49)
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_MFGB                      Maximum initial mass for which helium ignites on the First Giant Branch
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Minimum radius during Core Helium Burning in Rsol
 *
 * p_Alpha1, p_MHeF, p_MFGB and p_BnCoefficients passed as parameters so function can be declared static
 */
double CHeB::CalculateMinimumRadiusOnPhase_Static(const double      p_Mass,
                                                  const double      p_CoreMass,
                                                  const double      p_Alpha1,
                                                  const double      p_MHeF,
                                                  const double      p_MFGB,
                                                  const double      p_MinimumLuminosityOnPhase,
                                                  const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    double minR = 0.0;  // Minimum radius

    if (utils::Compare(p_MHeF, p_Mass) < 0) {
        double m_b28 = PPOW(p_Mass, b[28]);     // pow() is slow - do it once only
        minR = ((b[24] * p_Mass) + (PPOW((b[25] * p_Mass), b[26]) * m_b28)) / (b[27] + m_b28);
    }
    else {
        double LZAHB_MHeF = GiantBranch::CalculateLuminosityOnZAHB_Static(p_MHeF, p_CoreMass, p_Alpha1, p_MHeF, p_MFGB, p_MinimumLuminosityOnPhase, p_BnCoefficients);
        double LZAHB      = GiantBranch::CalculateLuminosityOnZAHB_Static(p_Mass, p_CoreMass, p_Alpha1, p_MHeF, p_MFGB, p_MinimumLuminosityOnPhase, p_BnCoefficients);
        double mu         = p_Mass / p_MHeF;
        double MHeF_b28   = PPOW(p_MHeF, b[28]);     // pow() is slow - do it once only
        double top        = ((b[24] * p_MHeF) + (PPOW((b[25] * p_MHeF), b[26]) * MHeF_b28)) / (b[27] + MHeF_b28);
        double bottom     = GiantBranch::CalculateRadiusOnPhase_Static(p_MHeF, LZAHB_MHeF, p_BnCoefficients);

        minR = GiantBranch::CalculateRadiusOnPhase_Static(p_Mass, LZAHB, p_BnCoefficients) * PPOW(top / bottom, mu);
    }

    return minR;

#undef b
}


/*
 * Calculate the radius at the start of the blue phase of Core Helium Burning
 *
 * Hurley et al. 2000, eq 60
 *
 *
 * double CalculateRadiusAtBluePhaseStart(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at the start of the blue phase of Core Helium Burning in Rsol
 */
double CHeB::CalculateRadiusAtBluePhaseStart(const double p_Mass) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double Rx;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        Rx = GiantBranch::CalculateRadiusOnZAHB_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
    }
    else if (utils::Compare(p_Mass, massCutoffs(MFGB)) < 0) {
        double luminosity = CalculateMinimumLuminosityOnPhase(p_Mass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_BnCoefficients);
        Rx = GiantBranch::CalculateRadiusOnPhase(p_Mass, luminosity);
    }
    else {
        Rx = CalculateRadiusAtHeIgnition(p_Mass);
    }

    return Rx;

#undef massCutoffs
}


/*
 * Calculate the radius at the end of the blue phase of Core Helium Burning
 *
 * Hurley et al. 2000, just before eq 64
 *
 *
 * double CalculateRadiusAtBluePhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at the end of the blue phase of Core Helium Burning in Rsol
 */
double CHeB::CalculateRadiusAtBluePhaseEnd(const double p_Mass) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    // JR: in all cases in the original code (function radiusCHeBY):
    //      - tau = ty
    //      - tau (ty) can never be < tx (so tau (ty) must be >= tx)
    //      - Ry = RAGB(Lx) (which is correct according to Hurley et al. 2000)
    return EAGB::CalculateRadiusOnPhase_Static(p_Mass, CalculateLuminosityAtBluePhaseEnd(m_Mass0), massCutoffs(MHeF), m_BnCoefficients);

#undef massCutoffs
}


/*
 * Calculate the parameter Rho for radius during Core Helium Burning
 *
 * Hurley et al. 2000, eq 65
 *
 *
 * double CalculateRadiusRho(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       CHeB relative age
 * @return                                      Rho
 */
double CHeB::CalculateRadiusRho(const double p_Mass, const double p_Tau) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double tx = timescales(tauX_BL);
    double ty = timescales(tauY_BL);

    double Rx   = CalculateRadiusAtBluePhaseStart(p_Mass);
    double Ry   = CalculateRadiusAtBluePhaseEnd(p_Mass);
    double RmHe = CalculateMinimumRadiusOnPhase_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
    double Rmin = std::min(RmHe, Rx);

    double ty_tx = ty - tx;

    double one   = std::cbrt(log((Ry / Rmin)));
    double two   = (p_Tau - tx) / ty_tx;
    double three = std::cbrt(log((Rx / Rmin)));
    double four  = (ty - p_Tau) / ty_tx;

    return (one * two) - (three * four);

#undef massCutoffs
#undef timescales
}


/*
 * Calculate the radius during Core Helium Burning
 *
 * Hurley et al. 2000, eq 64
 *
 *
 * double CalculateRadiusOnPhase(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       CHeB relative age
 * @return                                      Radius during Core Helium Burning in Rsol
 */
double CHeB::CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity, const double p_Tau) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double RCHeB;

    double tx  = timescales(tauX_BL);
    double ty  = timescales(tauY_BL);
    double rho = 0.0;

    if (utils::Compare(p_Tau, tx) < 0) {
        RCHeB = GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity);
    }
    else if (utils::Compare(p_Tau, ty) > 0) {
        RCHeB = EAGB::CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, massCutoffs(MHeF), m_BnCoefficients);
    }
    else  {
        double RmHe = CalculateMinimumRadiusOnPhase_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
        double Rx   = CalculateRadiusAtBluePhaseStart(p_Mass);
        double Rmin = std::min(RmHe, Rx);

        rho         = std::abs(CalculateRadiusRho(p_Mass, p_Tau));
        RCHeB       = Rmin * exp(rho * rho * rho);
    }

    return RCHeB;

#undef massCutoffs
#undef timescales
}


/*
 * Calculate radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantRadius()
 *
 * @return                                      Radius of remnant core in Rsol
 */
double CHeB::CalculateRemnantRadius() const {
    return HeMS::CalculateRadiusOnPhase_Static(m_CoreMass, CalculateTauOnPhase());
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the core mass between helium ignition and the Base of the Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 67
 *
 *
 * double CalculateCoreMassOnPhase(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       Relative age on phase
 * @return                                      Core mass on the First Giant Branch in Msol
 */
double CHeB::CalculateCoreMassOnPhase(const double p_Mass, const double p_Tau) const {
    return std::min(((1.0 - p_Tau) * CalculateCoreMassAtHeIgnition(p_Mass)) + (p_Tau * CalculateCoreMassAtBAGB(p_Mass)), m_Mass);               //He mass capped at total mass (should become HeMS star)
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate relative age during Core Helium Burning
 *
 * Hurley et al. 2000, just before eq 59
 * Naturally bounded by [0, 1], but clamp here anyway
 *
 *
 * double CalculateTauOnPhase()
 *
 * @return                                      CHeB relative age, clamped to [0, 1]
 */

double CHeB::CalculateTauOnPhase() const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
    return std::max(0.0, std::min(1.0, (m_Age - timescales(tHeI)) / timescales(tHe)));
#undef timescales
}


/*
 * Calculate the lifetime of Core Helium Burning
 *
 * Hurley et al. 2000, eq 57
 *
 *
 * double CalculateLifetimeOnPhase(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Lifetime of Core Helium Burning in MYRs (tHe)
 *
 * JR: changed this to use m_Timescales[TS::tBGB] instead of parameter
 */
double CHeB::CalculateLifetimeOnPhase(const double p_Mass) {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]      // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double tHe;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        double tHeMS = HeMS::CalculateLifetimeOnPhase_Static(m_CoreMass);   // can't use Timescales here - calculated using Mass not CoreMass
        double mu    = p_Mass / massCutoffs(MHeF);

        tHe = (b[39] + ((tHeMS - b[39]) * PPOW((1.0 - mu), b[40]))) * (1.0 + (m_Alpha4 * exp(15.0 * (p_Mass - massCutoffs(MHeF)))));
    }
    else {
        double m_5 = p_Mass * p_Mass * p_Mass * p_Mass * p_Mass;            // pow() is slow - use multiplication (sqrt() is much faster than pow())

        tHe = timescales(tBGB) * (((b[41] * PPOW(p_Mass, b[42])) + (b[43] * m_5)) / (b[44] + m_5));
    }

    return tHe;

#undef massCutoffs
#undef timescales
#undef b
}


/*
 * Calculate the parameter f_bl
 *
 * Hurley et al. 2000, just after eq 58
 *
 *
 * double CalculateBluePhaseFBL(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Blue phase fbl - see Hurley et al. 2000, eq 58 and just after
 */
double CHeB::CalculateBluePhaseFBL(const double p_Mass) {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    // Calculate RmHe for M > MFGB > MHeF
    double m_b28 = PPOW(p_Mass, b[28]);  // pow() is slow - do it once only
    double top   = ((b[24] * p_Mass) + (PPOW((b[25] * p_Mass), b[26]) * m_b28)) / (b[27] + m_b28);

    // Might be that we are supposed to use min(RmHe, Rx=RHeI)
    double RHeI = CalculateRadiusAtHeIgnition(p_Mass);
    double LHeI = GiantBranch::CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);

    top = std::min(top, RHeI);

    // Calculate RAGB(LHeI(M)) for M > MFGB > MHeF
    double bottom   = EAGB::CalculateRadiusOnPhase_Static(p_Mass, LHeI, massCutoffs(MHeF), m_BnCoefficients);
    double brackets = 1.0 - (top / bottom);

    return PPOW(p_Mass, b[48]) * PPOW(brackets, b[49]);

#undef massCutoffs
#undef b
}


/*
 * Calculate tHe relative lifetime of the blue phase of Core Helium Burning
 *
 * Hurley et al. 2000, eq 58
 *
 *
 * double CalculateLifetimeOnBluePhase(double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Relative lifetime of blue phase of Core Helium Burning, clamped to [0, 1]
 */
double CHeB::CalculateLifetimeOnBluePhase(const double p_Mass) {
#define b m_BnCoefficients                                              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double tbl;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        tbl = 1.0;
    }
    else if (utils::Compare(p_Mass, massCutoffs(MFGB)) <= 0) {
        double m_MFGB    = p_Mass / massCutoffs(MFGB);
        double firstTerm = (b[45] * PPOW(m_MFGB, 0.414));
        tbl              = firstTerm + ((1.0 - firstTerm) * PPOW((log10(m_MFGB) / log10(massCutoffs(MHeF) / massCutoffs(MFGB))), b[46]));
    }
    else {
        double fblM    = CalculateBluePhaseFBL(p_Mass);
        double fblMFGB = CalculateBluePhaseFBL(massCutoffs(MFGB));

        tbl = (1.0 - b[47]) * (fblM / fblMFGB);
    }

    return std::min(1.0, std::max(0.0, tbl));

#undef massCutoffs
#undef b
}

/*
 * Determine whether star should continue to evolve on phase
 *
 *
 * bool            ShouldEvolveOnPhase()
 *
 * @return         true if evolution should continue on phase, false otherwise
 */

bool CHeB::ShouldEvolveOnPhase() const {
    bool afterHeIgnition = (m_Age >= m_Timescales[static_cast<int>(TIMESCALE::tHeI)]);
    bool beforeEndOfHeBurning = (m_Age < (m_Timescales[static_cast<int>(TIMESCALE::tHeI)] + m_Timescales[static_cast<int>(TIMESCALE::tHe)]));
    bool coreIsNotTooMassive = (m_HeCoreMass < m_Mass);
    // Evolve on CHeB phase if age after He Ign and while He Burning and He core mass does not exceed total mass (could happen due to mass loss)
    return (afterHeIgnition && beforeEndOfHeBurning && coreIsNotTooMassive);
}

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Determine the star's envelope type.
 *
 *
 *
 * ENVELOPE DetermineEnvelopeType()
 *
 * @return                                      ENVELOPE::{ RADIATIVE, CONVECTIVE, REMNANT }
 */
ENVELOPE CHeB::DetermineEnvelopeType() const {
    
    ENVELOPE envelope = ENVELOPE::CONVECTIVE;                                                       // default envelope type  is CONVECTIVE
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                                 // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:
            envelope = ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:
            envelope =  utils::Compare(Temperature() * TSOL, CONVECTIVE_BOUNDARY_TEMPERATURE) ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;  // Envelope is radiative if temperature exceeds fixed threshold, otherwise convective
            break;
            
        default:                                                                                    // unknown prescription - use default envelope type
            SHOW_WARN(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION, "Using Envelope = CONVECTIVE");   // show warning
    }
    
    return envelope;
}


/*
 * Choose timestep for evolution
 *
 * Can obviously do this your own way
 * Given in the discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double CHeB::ChooseTimestep(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double dtk = 2.0E-3 * timescales(tHe);
    double dte = timescales(tHeI) + timescales(tHe) - p_Time;

// JR: todo: I took this next IF out - was it left in by mistake after some debugging?
//    if (dtk < dte) {
//        dtk = 2E-3 * m_Timescales[TS::tHe];
//    }

    double timestep = std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);

    return timestep;

#undef timescales
}



/*
 * Modify the star after it loses its envelope
 *
 * Hurley et al. 2000, section 6 just before eq 76 and after Eq. 105
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBParams
 *     - m_Luminosity
 *     - m_Radius
 *     - m_Mass
 *     - m_Mass0
 *     - m_CoreMass
 *     - m_HeCoreMass
 *     - m_COCoreMass
 *     - m_Age
 *
 *
 * STELLAR_TYPE ResolveEnvelopeLoss()
 *
 * @return                                      Stellar Type to which star shoule evolve after losing envelope
 */
STELLAR_TYPE CHeB::ResolveEnvelopeLoss(bool p_NoCheck) {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_NoCheck || utils::Compare(m_CoreMass, m_Mass) >= 0) {                     // Envelope loss

        m_CoreMass   = m_Mass;
        m_Mass0      = m_Mass;
        m_HeCoreMass = m_CoreMass;
        m_COCoreMass = 0.0;

 		// set evolved time for naked helium star since already has some core mass.
		double tHeIPrime = timescales(tHeI);
		double tHePrime  = timescales(tHe);

        m_Tau = ((m_Age - tHeIPrime) / tHePrime);                                   // Hurley et al. 2000, just after eq 81: tau = t/tHeMS
        m_Age = m_Tau * HeMS::CalculateLifetimeOnPhase_Static(m_Mass0);             // see Hurley et al. 2000, eq 76 and following discussion

        CalculateTimescales(m_Mass0, m_Timescales);
        CalculateGBParams(m_Mass0, m_GBParams);

        m_Luminosity = HeMS::CalculateLuminosityOnPhase_Static(m_Mass, m_Tau);
        m_Radius     = HeMS::CalculateRadiusOnPhase_Static(m_Mass, m_Tau);

        stellarType  = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                          // will evolve to an evolved helium star
    }

    return stellarType;

#undef timescales
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE CHeB::EvolveToNextPhase() {
    return STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH;
}
