#include "HG.h"
#include "HeMS.h"
#include "HeWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                     COEFFICIENT AND CONSTANT CALCULATIONS ETC.                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the parameter rho for the Hertzsprung Gap
 *
 * Hurley et al. 2000, eq 29
 *
 *
 * double CalculateRho(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Rho - Constant such that MCTMS = Rho * MCEHG
 */
double HG::CalculateRho(const double p_Mass) const {

    double m_5_25 = p_Mass * p_Mass * p_Mass * p_Mass * p_Mass * std::sqrt(std::sqrt(p_Mass));    // pow() is slow - use multiplication (sqrt() is much faster than pow())

    return (1.586 + m_5_25) / (2.434 + (1.02 * m_5_25));
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
double HG::CalculateLambdaDewi() const {

	double lambda1 = std::min(0.80, (3.0 / (2.4 + PPOW(m_Mass,-3.0 / 2.0))) - (0.15 * log10(m_Luminosity)));                // (A.3) Claeys+2014
	double lambda2 = 0.42 * PPOW(m_RZAMS / m_Radius, 0.4);                                                                  // (A.2) Claeys+2014
	double envMass = utils::Compare(m_CoreMass, 0.0) > 0 && utils::Compare(m_Mass, m_CoreMass) > 0 ? m_Mass - m_CoreMass : 0.0;

    double lambdaCE;

         if (utils::Compare(envMass, 1.0) >= 0) lambdaCE = 2.0 * lambda1;                                                   // (A.1) Bottom, Claeys+2014
	else if (utils::Compare(envMass, 0.0) >  0) lambdaCE = 2.0 * (lambda2 + (std::sqrt(envMass) * (lambda1 - lambda2)));         // (A.1) Mid, Claeys+2014
	else                                        lambdaCE = 2.0 * lambda2;			                                        // (A.1) Top, Claeys+2014

	return	lambdaCE;
}


/*
 * Calculate the common envelope lambda parameter using the enhanced "Nanjing" prescription
 * from X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010)
 *
 * This function good for HG and FGB stars.
 *
 *
 * double CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind)
 *
 * @return                                      Nanjing lambda for use in common envelope
 */
double HG::CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const {

	DBL_VECTOR maxBG    = {};                                                           // [0] = maxB, [1] = maxG
	DBL_VECTOR lambdaBG = {};                                                           // [0] = lambdaB, [1] = lambdaG
	DBL_VECTOR a        = {};                                                           // 0..5 a_coefficients
	DBL_VECTOR b        = {};                                                           // 0..5 b_coefficients
    double     Rmax     = std::numeric_limits<double>::max();                           // Upper R limit to applicability of Nanjing polynomials.

    switch(p_Zind) {
        // Pop. I metallicity
        case 1:
            switch(p_MassInd) {
                case 0: {
                    maxBG = { 2.5, 1.5 };
                    Rmax = 200.0;
                    double R_in = std::min(Rmax, m_Radius);
                    if (utils::Compare(R_in, 2.7) > 0) {
                        lambdaBG = { 2.33 - (R_in * 9.18E-03), 1.12 - (R_in * 4.59E-03) };
                    }
                    else {
                        a = {  8.35897, -18.89048, 10.47651, 0.99352, 0.0, 0.0 };
                        b = { 17.58328, -34.84355, 10.70536, 8.49042, 0.0, 0.0 };
                    }
                    break;
                }
                case 1:
                    maxBG = { 4.0, 2.0 };
                    Rmax = 340.0;
                    a = { 2.05363, -0.00685, -3.42739E-04, 3.93987E-06, -1.18237E-08, 0.0 };
                    b = { 1.07658, -0.01041, -4.90553E-05, 1.13528E-06, -3.91609E-09, 0.0 };
                    break;
                case 2:
                    Rmax = 400.0;
                    maxBG = { 2.5, 1.5 };
                    a     = { 2.40831, -0.42459, 0.03431, -9.26879E-04, 8.24522E-06, 0.0 };
                    b     = { 1.30705, -0.22924, 0.01847, -5.06216E-04, 4.57098E-06, 0.0 };
                    break;
                case 3:
                    Rmax = 410.0;
                    maxBG = { 2.5, 1.5 };
                    a     = { 1.8186 , -0.17464, 0.00828, -1.31727E-04, 7.08329E-07, 0.0 };
                    b     = { 1.02183, -0.1024 , 0.00493, -8.16343E-05, 4.55426E-07, 0.0 };
                    break;
                case 4:
                    maxBG = { 1000.0, 8.0 };
                    Rmax = 430.0;
                    a = { 1.52581, -0.08125, 0.00219, -2.0527E-05 , 6.79169E-08, 0.0 };
                    b = { 0.85723, -0.04922, 0.00137, -1.36163E-05, 4.68683E-08, 0.0 };
                    break;
                case 5:
                    maxBG = { 25.5, 5.0 };
                    Rmax = 440.0;
                    a = { 1.41601, -0.04965, 8.51527E-04, -5.54384E-06, 1.32336E-08, 0.0 };
                    b = { 0.78428, -0.02959, 5.2013E-04 , -3.45172E-06, 8.17248E-09, 0.0 };
                    break;
                case 6:
                    maxBG = { 9.0, 3.0 };
                    Rmax = 420.0;
                    a = { 1.38344, -0.04093, 5.78952E-04, -3.19227E-06, 6.40902E-09, 0.0 };
                    b = { 0.76009, -0.02412, 3.47104E-04, -1.92347E-06, 3.79609E-09, 0.0 };
                    break;
                case 7:
                    maxBG = { 7.0, 3.0 };
                    Rmax = 490.0;
                    a = { 1.35516, -0.03414, 4.02065E-04, -1.85931E-06, 3.08832E-09, 0.0 };
                    b = { 0.73826, -0.01995, 2.37842E-04, -1.09803E-06, 1.79044E-09, 0.0 };
                    break;
                case 8:
                    maxBG = { 4.0, 2.0 };
                    Rmax = 530.0;
                    a = { 1.32549, -0.02845, 2.79097E-04, -1.07254E-06, 1.46801E-09, 0.0 };
                    b = { 0.71571, -0.01657, 1.64607E-04, -6.31935E-07, 8.52082E-10, 0.0 };
                    break;
                case 9:
                    Rmax = 600.0;
                    maxBG = { 1.0, 0.6 };
                    a     = { 1.29312, -0.02371, 1.93764E-04, -6.19576E-07, 7.04227E-10, 0.0 };
                    b     = { 0.69245, -0.01398, 1.17256E-04, -3.81487E-07, 4.35818E-10, 0.0 };
                    break;
                case 10:
                    maxBG = { 1.5, 1.0 };
                    Rmax = 850.0;
                    if (utils::Compare(m_Radius, 0.0) > 0 && utils::Compare(m_Radius, 350.0) <= 0) {
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
                    break;
                case 11:
                    maxBG = { 1.5, 1.0 };
                    Rmax = 1000.0;
                    if (utils::Compare(m_Radius, 190.0) > 0 && utils::Compare(m_Radius, 600.0) < 0) lambdaBG = { 0.15, 0.15 };
                    else {
                        a = { 1.39332, -0.0318 , 3.95917E-04, -2.23132E-06, 4.50831E-09, 0.0 };
                        b = { 0.78215, -0.02326, 3.25984E-04, -1.94991E-06, 4.08044E-09, 0.0 };
                    }
                    break;
                case 12:
                    maxBG = { 1.5, 1.0 };
                    Rmax = 1050.0;
                    if (utils::Compare(m_Radius, 120.0) > 0 && utils::Compare(m_Radius, 170.0) < 0) lambdaBG = { 0.2, 0.2 };
                    else {
                        a = { 1.43177, -0.03533, 5.11128E-04, -3.57633E-06, 9.36778E-09, 0.0 };
                        b = { 0.85384, -0.03086, 5.50878E-04, -4.37671E-06, 1.25075E-08, 0.0 };
                    }
                    break;
                case 13: {
                    maxBG = { 1.5, 1.0 };
                    Rmax = 1200.0;
                    double R_in = std::min(Rmax, m_Radius);
                    lambdaBG = { 1.2 * exp(-R_in / 90.0), 0.55 * exp(-R_in / 160.0) };
                    break;
                }
                case 14:
                    maxBG = { 1.0, 0.5 };
                    a     = { 0.31321, -7.50384E-04, 5.38545E-07, -1.16946E-10, 0.0, 0.0 };
                    b     = { 0.159  , -3.94451E-04, 2.88452E-07, -6.35132E-11, 0.0, 0.0 };
                    break;
                case 15:
                    maxBG = { 1.0, 0.5 };
                    a     = { 0.376 , -0.0018 , 2.81083E-06, -1.67386E-09, 3.35056E-13, 0.0 };
                    b     = { 0.2466, -0.00121, 1.89029E-06, -1.12066E-09, 2.2258E-13 , 0.0 };
                    break;
                }
                break;

        // Pop. II metallicity
        case 0:
            switch(p_MassInd) {
                case 0: {
                    maxBG = { 2.0, 1.5 };
                    Rmax = 160.0;
                    double R_in = std::min(Rmax, m_Radius);
                    if (utils::Compare(R_in, 12.0)  > 0) {
                        lambdaBG = { 1.8 * exp(-R_in / 80.0), exp(-R_in / 45.0) };
                    }
                    else {
                        a = { 0.24012, -0.01907, 6.09529E-04, -8.17819E-06, 4.83789E-08, -1.04568E-10 };
                        b = { 0.15504, -0.01238, 3.96633E-04, -5.3329E-06 , 3.16052E-08, -6.84288E-11 };
                    }
                    break;
                }
                case 1:
                    maxBG = { 4.0, 2.0 };
                    Rmax = 350.0;
                    if (utils::Compare(m_Radius, 22.0) > 0 && utils::Compare(m_Radius, 87.0) < 0) lambdaBG = { 1.95, 0.85 };
                    else {
                        a = { 2.56108, -0.75562, 0.1027 , -0.00495, 8.05436E-05, 0.0 };
                        b = { 1.41896, -0.4266 , 0.05792, -0.00281, 4.61E-05   , 0.0 };
                    }
                    break;
                case 2:
                    maxBG = { 600.0, 2.0 };
                    Rmax = 400.0;
                    a = { 1.7814 , -0.17138, 0.00754, -9.02652E-05, 0.0, 0.0 };
                    b = { 0.99218, -0.10082, 0.00451, -5.53632E-05, 0.0, 0.0 };
                    break;
                case 3:
                    maxBG = { 600.0, 2.0 };
                    Rmax = 410.0;
                    a = { 1.65914, -0.10398, 0.0029 , -2.24862E-05, 0.0, 0.0 };
                    b = { 0.92172, -0.06187, 0.00177, -1.42677E-05, 0.0, 0.0 };
                    break;
                case 4:
                    maxBG = { 10.0, 3.0 };
                    Rmax = 320.0;
                    a = { 1.58701, -0.06897, 0.00129    , -6.99399E-06, 0.0, 0.0 };
                    b = { 0.87647, -0.04103, 7.91444E-04, -4.41644E-06, 0.0, 0.0 };
                    break;
                case 5:
                    maxBG = { 4.0, 1.5 };
                    Rmax = 330.0;
                    a = { 1.527  , -0.04738, 6.1373E-04 , -2.36835E-06, 0.0, 0.0 };
                    b = { 0.83636, -0.02806, 3.73346E-04, -1.47016E-06, 0.0, 0.0 };
                    break;
                case 6:
                    maxBG = { 2.5, 1.0 };
                    Rmax = 360.0;
                    a = { 1.49995, -0.03921, 4.2327E-04, -1.37646E-06, 0.0, 0.0 };
                    b = { 0.81688, -0.02324, 2.5804E-04, -8.54696E-07, 0.0, 0.0 };
                    break;
                case 7:
                    maxBG = { 2.0, 1.0 };
                    Rmax = 400.0;
                    a = { 1.46826, -0.03184, 2.85622E-04, -7.91228E-07, 0.0, 0.0 };
                    b = { 0.79396, -0.01903, 1.77574E-04, -5.04262E-07, 0.0, 0.0 };
                    break;
                case 8:
                    maxBG = { 1.6, 1.0 };
                    Rmax = 440.0;
                    a = { 1.49196, -0.03247, 3.08066E-04, -9.53247E-07, 0.0, 0.0 };
                    b = { 0.805  , -0.02   , 2.01872E-04, -6.4295E-07 , 0.0, 0.0 };
                    break;
                case 9: {
                    maxBG = { 1.6, 1.0 };
                    Rmax = 500.0;
                    double R_in = std::min(Rmax, m_Radius);
                    lambdaBG = { 1.75 * exp(-R_in / 35.0), 0.9 * exp(-R_in /35.0) };
                    break;
                }
                case 10:
                    maxBG = { 1.6, 1.0 };
                    Rmax = 600.0;
                    a = { 1.63634, -0.04646, 7.49351E-04, -5.23622E-06, 0.0, 0.0 };
                    b = { 1.17934, -0.08481, 0.00329    , -4.69096E-05, 0.0, 0.0 };
                    break;
                case 11:
                    maxBG = { 1.6, 1.0 };
                    Rmax = 650.0;
                    a = { 1.45573, -0.00937, -0.00131,  3.07004E-05, 0.0, 0.0 };
                    b = { 1.19526, -0.08503,  0.00324, -4.58919E-05, 0.0, 0.0 };
                    break;
                case 12:
                    maxBG = { 1.5, 1.0 };
                    Rmax = 750.0;
                    a = { 1.33378,  0.01274, -0.00234,  4.6036E-05 , 0.0, 0.0 };
                    b = { 1.17731, -0.07834,  0.00275, -3.58108E-05, 0.0, 0.0 };
                    break;
                case 13:
                    maxBG = { 1.5, 1.0 };
                    Rmax = 900.0;
                    a = { 1.27138,  0.00538, -0.0012 ,  1.80776E-05, 0.0, 0.0 };
                    b = { 1.07496, -0.05737,  0.00153, -1.49005E-05, 0.0, 0.0 };
                    break;
                case 14:
                    maxBG = { 20.0, 3.0 };
                    a     = { 0.821  , -0.00669, 1.57665E-05, -1.3427E-08 , 3.74204E-12, 0.0 };
                    b     = { 0.49287, -0.00439, 1.06766E-05, -9.22015E-09, 2.58926E-12, 0.0 };
                    break;
                case 15:
                    maxBG = { 4.0, 2.0 };
                    a     = { 1.25332, -0.02065, 1.3107E-04 , -3.67006E-07, 4.58792E-10, -2.09069E-13 };
                    b     = { 0.81716, -0.01436, 9.31143E-05, -2.6539E-07 , 3.30773E-10, -1.51207E-13 };
                    break;
            }
            break;
        }

    if (lambdaBG.empty()) {
        if ( (p_Zind == 1) && (p_MassInd == 0) ) {                        // Pop. I metallicity and M = 1 Msun
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
            double x  = std::min(m_Radius, Rmax);
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
    lambdaBG[1] = std::min( std::max(0.05, lambdaBG[1]), std::min(1.0, maxBG[1]) );         // clamp lambda G to [0.05, min(1,maxG)]
    lambdaBG[0] = std::max( std::min(lambdaBG[0],maxBG[0]), std::max(0.05, lambdaBG[1]) );  // clamp lambda B to [ max(0.05,lambdaG), maxB]

    // Calculate lambda as some combination of lambda_b and lambda_g by
    // lambda = alpha_th • lambda_b    +  (1-alpha_th) • lambda_g
    // STARTRACK uses alpha_th = 1/2
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaBG[0]) + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaBG[1]);
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
 * This function good for HG and FGB stars.
 *
 *
 * double CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity)
 *
 * @return                                      Nanjing lambda for use in common envelope
 */
double HG::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const {

	DBL_VECTOR maxBG    = {};                                                           // [0] = maxB, [1] = maxG
	DBL_VECTOR lambdaBG = {};                                                           // [0] = lambdaB, [1] = lambdaG
	DBL_VECTOR a        = {};                                                           // 0..5 a_coefficients
	DBL_VECTOR b        = {};                                                           // 0..5 b_coefficients

    if (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT) > 0) {                     // Z>0.5 Zsun: popI
        if (utils::Compare(p_Mass, 1.5) < 0) {
            maxBG = { 2.5, 1.5 };
                 if (utils::Compare(m_Radius, 200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else if (utils::Compare(m_Radius, 2.7  ) > 0) lambdaBG = { 2.33 - (m_Radius * 9.18E-03), 1.12 - (m_Radius * 4.59E-03) };
            else {
                a = {  8.35897, -18.89048, 10.47651, 0.99352, 0.0, 0.0 };
                b = { 17.58328, -34.84355, 10.70536, 8.49042, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 340.0) > 0) lambdaBG = { 3.589970, 0.514132 };
            else {
                a = { 2.05363, -0.00685, -3.42739E-04, 3.93987E-06, -1.18237E-08, 0.0 };
                b = { 1.07658, -0.01041, -4.90553E-05, 1.13528E-06, -3.91609E-09, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 3.5) < 0) {
            maxBG = { 500.0, 10.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 116.935557, 0.848808 };
            else {
                maxBG = { 2.5, 1.5 };
                a     = { 2.40831, -0.42459, 0.03431, -9.26879E-04, 8.24522E-06, 0.0 };
                b     = { 1.30705, -0.22924, 0.01847, -5.06216E-04, 4.57098E-06, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 4.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 410.0) > 0) lambdaBG = { 52.980056, 1.109736 };
            else {
                maxBG = { 2.5, 1.5 };
                a     = { 1.8186 , -0.17464, 0.00828, -1.31727E-04, 7.08329E-07, 0.0 };
                b     = { 1.02183, -0.1024 , 0.00493, -8.16343E-05, 4.55426E-07, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 5.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 430.0) > 0) lambdaBG = { 109.593522, 1.324248 };
            else {
                a = { 1.52581, -0.08125, 0.00219, -2.0527E-05 , 6.79169E-08, 0.0 };
                b = { 0.85723, -0.04922, 0.00137, -1.36163E-05, 4.68683E-08, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 6.5) < 0) {
            maxBG = { 25.5, 5.0 };
            if (utils::Compare(m_Radius, 440.0) > 0) lambdaBG = { 16.279603, 1.352166 };
            else {
                a = { 1.41601, -0.04965, 8.51527E-04, -5.54384E-06, 1.32336E-08, 0.0 };
                b = { 0.78428, -0.02959, 5.2013E-04 , -3.45172E-06, 8.17248E-09, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 7.5) < 0) {
            maxBG = { 9.0, 3.0 };
            if (utils::Compare(m_Radius, 420.0) > 0) lambdaBG = { 5.133959, 1.004036 };
            else {
                a = { 1.38344, -0.04093, 5.78952E-04, -3.19227E-06, 6.40902E-09, 0.0 };
                b = { 0.76009, -0.02412, 3.47104E-04, -1.92347E-06, 3.79609E-09, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 8.5) < 0) {
            maxBG = { 7.0, 3.0 };
            if (utils::Compare(m_Radius, 490.0) > 0) lambdaBG = { 4.342985, 0.934659 };
            else {
                a = { 1.35516, -0.03414, 4.02065E-04, -1.85931E-06, 3.08832E-09, 0.0 };
                b = { 0.73826, -0.01995, 2.37842E-04, -1.09803E-06, 1.79044E-09, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 9.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 530.0) > 0) lambdaBG = { 2.441672, 0.702310 };
            else {
                a = { 1.32549, -0.02845, 2.79097E-04, -1.07254E-06, 1.46801E-09, 0.0 };
                b = { 0.71571, -0.01657, 1.64607E-04, -6.31935E-07, 8.52082E-10, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 11.0) < 0) {
            maxBG = { 3.0, 1.5 };
            if (utils::Compare(m_Radius, 600.0) > 0) lambdaBG = { 1.842314, 0.593854 };
            else {
                maxBG = { 1.0, 0.6 };
                a     = { 1.29312, -0.02371, 1.93764E-04, -6.19576E-07, 7.04227E-10, 0.0 };
                b     = { 0.69245, -0.01398, 1.17256E-04, -3.81487E-07, 4.35818E-10, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 13.0) < 0) {
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
        else if (utils::Compare(p_Mass, 15.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 1000.0) > 0)                                       lambdaBG = { 0.414200, 0.189008 };
            else if (utils::Compare(m_Radius, 190.0) > 0 && utils::Compare(m_Radius, 600.0) < 0) lambdaBG = { 0.15, 0.15 };
            else {
                a = { 1.39332, -0.0318 , 3.95917E-04, -2.23132E-06, 4.50831E-09, 0.0 };
                b = { 0.78215, -0.02326, 3.25984E-04, -1.94991E-06, 4.08044E-09, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 1050.0) > 0)                                       lambdaBG = { 0.2, 0.1 };
            else if (utils::Compare(m_Radius, 120.0) > 0 && utils::Compare(m_Radius, 170.0) < 0) lambdaBG = { 0.2, 0.2 };
            else {
                a = { 1.43177, -0.03533, 5.11128E-04, -3.57633E-06, 9.36778E-09, 0.0 };
                b = { 0.85384, -0.03086, 5.50878E-04, -4.37671E-06, 1.25075E-08, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else                                      lambdaBG = { 1.2 * exp(-m_Radius / 90.0), 0.55 * exp(-m_Radius / 160.0) };
        }
        else if (utils::Compare(p_Mass, 75.0) < 0) {
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
        if (utils::Compare(p_Mass, 1.5) < 0) {
            maxBG = { 2.0, 1.5 };
                 if (utils::Compare(m_Radius, 160.0) > 0) lambdaBG = { 0.05, 0.05 };
            else if (utils::Compare(m_Radius, 12.0)  > 0) lambdaBG = { 1.8 * exp(-m_Radius / 80.0), exp(-m_Radius / 45.0) };
            else {
                a = { 0.24012, -0.01907, 6.09529E-04, -8.17819E-06, 4.83789E-08, -1.04568E-10 };
                b = { 0.15504, -0.01238, 3.96633E-04, -5.3329E-06 , 3.16052E-08, -6.84288E-11 };
            }
        }
        else if (utils::Compare(p_Mass, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
                 if (utils::Compare(m_Radius, 350.0) > 0)                                      lambdaBG = { 2.868539, 0.389991 };
            else if (utils::Compare(m_Radius, 22.0) > 0 && utils::Compare(m_Radius, 87.0) < 0) lambdaBG = { 1.95, 0.85 };
            else {
                a = { 2.56108, -0.75562, 0.1027 , -0.00495, 8.05436E-05, 0.0 };
                b = { 1.41896, -0.4266 , 0.05792, -0.00281, 4.61E-05   , 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 3.5) < 0) {
            maxBG = { 600.0, 2.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 398.126442, 0.648560 };
            else {
                a = { 1.7814 , -0.17138, 0.00754, -9.02652E-05, 0.0, 0.0 };
                b = { 0.99218, -0.10082, 0.00451, -5.53632E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 4.5) < 0) {
            maxBG = { 600.0, 2.0 };
            if (utils::Compare(m_Radius, 410.0) > 0) lambdaBG = { 91.579093, 1.032432 };
            else {
                a = { 1.65914, -0.10398, 0.0029 , -2.24862E-05, 0.0, 0.0 };
                b = { 0.92172, -0.06187, 0.00177, -1.42677E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 5.5) < 0) {
            maxBG = { 10.0, 3.0 };
            if (utils::Compare(m_Radius, 320.0) > 0) lambdaBG = { 7.618019, 1.257919 };
            else {
                a = { 1.58701, -0.06897, 0.00129    , -6.99399E-06, 0.0, 0.0 };
                b = { 0.87647, -0.04103, 7.91444E-04, -4.41644E-06, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 6.5) < 0) {
            maxBG = { 4.0, 1.5 };
            if (utils::Compare(m_Radius, 330.0) > 0) lambdaBG = { 2.390575, 0.772091 };
            else {
                a = { 1.527  , -0.04738, 6.1373E-04 , -2.36835E-06, 0.0, 0.0 };
                b = { 0.83636, -0.02806, 3.73346E-04, -1.47016E-06, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 7.5) < 0) {
            maxBG = { 2.5, 1.0 };
            if (utils::Compare(m_Radius, 360.0) > 0) lambdaBG = { 1.878174, 0.646353 };
            else {
                a = { 1.49995, -0.03921, 4.2327E-04, -1.37646E-06, 0.0, 0.0 };
                b = { 0.81688, -0.02324, 2.5804E-04, -8.54696E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 8.5) < 0) {
            maxBG = { 2.0, 1.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 1.517662, 0.553169 };
            else {
                a = { 1.46826, -0.03184, 2.85622E-04, -7.91228E-07, 0.0, 0.0 };
                b = { 0.79396, -0.01903, 1.77574E-04, -5.04262E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 9.5) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 440.0) > 0) lambdaBG = { 1.136394, 0.478963 };
            else {
                a = { 1.49196, -0.03247, 3.08066E-04, -9.53247E-07, 0.0, 0.0 };
                b = { 0.805  , -0.02   , 2.01872E-04, -6.4295E-07 , 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 11.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 500.0) > 0) lambdaBG = { 1.068300, 0.424706 };
            else                                     lambdaBG = { 1.75 * exp(-m_Radius / 35.0), 0.9 * exp(-m_Radius /35.0) };
        }
        else if (utils::Compare(p_Mass, 13.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 600.0) > 0) lambdaBG = { 0.537155, 0.211105 };
            else {
                a = { 1.63634, -0.04646, 7.49351E-04, -5.23622E-06, 0.0, 0.0 };
                b = { 1.17934, -0.08481, 0.00329    , -4.69096E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 15.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 650.0) > 0) lambdaBG = { 0.3, 0.160696 };
            else {
                a = { 1.45573, -0.00937, -0.00131,  3.07004E-05, 0.0, 0.0 };
                b = { 1.19526, -0.08503,  0.00324, -4.58919E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 750.0) > 0) lambdaBG = { 0.5, 0.204092 };
            else {
                a = { 1.33378,  0.01274, -0.00234,  4.6036E-05 , 0.0, 0.0 };
                b = { 1.17731, -0.07834,  0.00275, -3.58108E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 900.0) > 0) lambdaBG = { 0.2, 0.107914 };
            else {
                a = { 1.27138,  0.00538, -0.0012 ,  1.80776E-05, 0.0, 0.0 };
                b = { 1.07496, -0.05737,  0.00153, -1.49005E-05, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 75.0) < 0) {
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
        if (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT) > 0 && utils::Compare(p_Mass, 1.5) < 0) {
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
    lambdaBG[0] = std::min(std::max(0.05, lambdaBG[0]), maxBG[0]);          // clamp lambda B to [0.05, maxB]
    lambdaBG[1] = std::min(std::max(0.05, lambdaBG[1]), maxBG[1]);          // clamp lambda G to [0.05, maxG]

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
 * Calculate luminosity at the end of the Hertzsprung Gap
 *
 * Hurley et al. 2000, just before eq 8
 *
 *
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the Hertzsprung Gap in Lsol
 */
double HG::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    return (utils::Compare(p_Mass, massCutoffs(MFGB)) < 0)
            ? GiantBranch::CalculateLuminosityAtPhaseBase_Static(p_Mass, m_AnCoefficients)
            : GiantBranch::CalculateLuminosityAtHeIgnition_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);

#undef massCutoffs
}


/*
 * Calculate luminosity on the Hertzsprung Gap
 *
 * Hurley et al. 2000, eq 26
 *
 *
 * double CalculateLuminosityOnPhase(const double p_Age, const double p_Mass)
 *
 * @param   [IN]    p_Age                       Effective age in Myr
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity on the HG in Lsol
 */
double HG::CalculateLuminosityOnPhase(const double p_Age, const double p_Mass) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double LTMS = MainSequence::CalculateLuminosityAtPhaseEnd(p_Mass);
    double LEHG = CalculateLuminosityAtPhaseEnd(p_Mass);
    double tMS  = timescales(tMS);
    double tBGB = timescales(tBGB);
    double tau  = (p_Age - tMS) / (tBGB - tMS);

    return LTMS * PPOW((LEHG / LTMS), tau);

#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate radius at the end of the Hertzsprung Gap
 *
 * Hurley et al. 2000, eqs 7 & 8
 *
 *
 * double CalculateRadiusAtPhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at the end of the Hertzsprung Gap in Rsol
 */
double HG::CalculateRadiusAtPhaseEnd(const double p_Mass) const {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    return (utils::Compare(p_Mass, massCutoffs(MFGB)) < 0)
            ? GiantBranch::CalculateRadiusOnPhase(p_Mass, GiantBranch::CalculateLuminosityAtPhaseBase_Static(p_Mass, m_AnCoefficients))
            : GiantBranch::CalculateRadiusAtHeIgnition(p_Mass);

#undef massCutoffs
}


/*
 * Calculate radius on the Hertzsprung Gap
 *
 * Hurley et al. 2000, eq 27
 *
 *
 * double CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       Relative age
 * @param   [IN]    p_RZAMS                     Zero Age Main Sequence (ZAMS) Radius
 * @return                                      Radius on the Hertzsprung Gap in Rsol
 */
double HG::CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const {

    double RTMS = MainSequence::CalculateRadiusAtPhaseEnd(p_Mass, p_RZAMS);
    double REHG = CalculateRadiusAtPhaseEnd(p_Mass);

    return RTMS * PPOW((REHG / RTMS), p_Tau);
}


/*
 * Calculate the radial extent of the star's convective envelope (if it has one)
 *
 * Hurley et al. 2000, sec. 2.3, particularly subsec. 2.3.1, eqs 36-40
 *
 * (Technically not a radius calculation I suppose, but "radial extent" is close enough to put it with the radius calculations...)
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the star's convective envelope in Rsol
 */
double HG::CalculateRadialExtentConvectiveEnvelope() const {

	BaseStar clone = *this;                         // clone this star so can manipulate without changes persisiting
	clone.ResolveEnvelopeLoss(true);                // update clone's attributes after envelope is lost

    return std::sqrt(m_Tau) * (m_Radius - clone.Radius());
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



/*
 * Calculate core mass at the end of the Hertzsprung Gap
 *
 * Hurley et al. 2000, eq 28
 *
 *
 * double CalculateCoreMassAtPhaseEnd(const double p_Mass) {
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Core mass at the end of the Hertzsprung Gap (Base of the Giant Branch) in Msol
 */
double HG::CalculateCoreMassAtPhaseEnd(const double p_Mass) const {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]                // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    double coreMass;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        double LBGB = GiantBranch::CalculateLuminosityAtPhaseBase_Static(p_Mass, m_AnCoefficients);
        coreMass    = BaseStar::CalculateCoreMassGivenLuminosity_Static(LBGB, m_GBParams);
    }
    else if (utils::Compare(p_Mass, massCutoffs(MFGB)) < 0) {
        coreMass = gbParams(McBGB);
    }
    else {
        coreMass = CalculateCoreMassAtHeIgnition(p_Mass);
    }

    return coreMass;

#undef massCutoffs
#undef gbParams
}


/*
 * Calculate core mass on the Hertzsprung Gap
 *
 * Hurley et al. 2000, eq 30
 *
 *
 * double CalculateCoreMassOnPhase(const double p_Mass, const double p_Time)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time after ZAMS in Myr (tBGB <= time <= tHeI)
 * @return                                      Core mass on the Hertzsprung Gap in Msol
 */
double HG::CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double McEHG = CalculateCoreMassAtPhaseEnd(p_Mass);
    double rhoHG = CalculateRho(p_Mass);
    double tau   = (p_Time - timescales(tMS)) / (timescales(tBGB) - timescales(tMS));

    // If the star is losing mass, choose core mass as the maximum of the core mass
    // at the previous time-step and the value given by Hurley et al. 2000, eq 30
    return std::max((((1.0 - tau) * rhoHG) + tau) * McEHG, m_CoreMass);

#undef timescales
}


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * Always returns 1.0 for HG
 * The rest of the code is here so sanity checks can be made and an error displayed if a bad prescription was specified in the program options
 *
 *
 * double CalculateMassTransferRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double HG::CalculateMassTransferRejuvenationFactor() const {

    double fRej = 1.0;

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {          // which rejuvenation prescription?

        case MT_REJUVENATION_PRESCRIPTION::NONE:                        // use default Hurley et al. 2000 prescription = 1.0
        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                   // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf

            if (utils::Compare(m_Mass, m_MassPrev) <= 0) {              // Rejuvenation factor is unity for mass losing stars
                fRej = 1.0;
            }
            break;

        default:                                                        // unknow prescription - use default Hurley et al. 2000 prescription = 1.0
            SHOW_WARN(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION);     // show warning
    }

    return fRej;
}


/*
 * Determines if mass transfer is unstable according to the critical mass ratio.
 *
 * See e.g de Mink et al. 2013, Claeys et al. 2014, and Ge et al. 2010, 2015, 2020 for discussions.
 *
 * Assumes this star is the donor; relevant accretor details are passed as parameters.
 * Critical mass ratio is defined as qCrit = mAccretor/mDonor.
 *
 * double HG::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate) 
 *
 * @param   [IN]    p_AccretorIsDegenerate      Boolean indicating if accretor in degenerate (true = degenerate)
 * @return                                      Critical mass ratio for unstable MT 
 */
double HG::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate) const {

    double qCrit;

    qCrit = p_AccretorIsDegenerate
                ? OPTIONS->MassTransferCriticalMassRatioHGDegenerateAccretor()              // degenerate accretor
                : OPTIONS->MassTransferCriticalMassRatioHGNonDegenerateAccretor();          // non-degenerate accretor
                                                                                                                        
    return qCrit;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate relative age on the Hertzsprung Gap
 *
 * Hurley et al. 2000, eq 25
 * Naturally bounded by [0, 1], but clamp here anyway
 *
 * double CalculateTauOnPhase()
 *
 * @return                                      HG relative age, clamped to [0, 1]
 */
double HG::CalculateTauOnPhase() const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
    return std::max(0.0, std::min(1.0, (m_Age - timescales(tMS)) / (timescales(tBGB) - timescales(tMS))));
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
void HG::UpdateAgeAfterMassLoss() {

    double tBGB      = m_Timescales[static_cast<int>(TIMESCALE::tBGB)];
    double tBGBprime = CalculateLifetimeToBGB(m_Mass);

    double tMS       = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
    double tMSprime  = MainSequence::CalculateLifetimeOnPhase(m_Mass, tBGBprime);

    m_Age = tMSprime + (((tBGBprime - tMSprime) / (tBGB - tMS)) * (m_Age - tMS));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



/*
 * Determine the star's envelope type.
 *
 * Some calculations on this can be found in sec. 2.3.4 of Belczynski et al. 2008.  For now, we will only do the calculation using stellarType.
 *
 *
 * ENVELOPE DetermineEnvelopeType()
 *
 * @return                                      ENVELOPE::{ RADIATIVE, CONVECTIVE, REMNANT }
 */
ENVELOPE HG::DetermineEnvelopeType() const {
 
    ENVELOPE envelope = ENVELOPE::RADIATIVE;                                                         // default envelope type is RADIATIVE
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                                  // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:
            envelope = ENVELOPE::RADIATIVE;                                                          // default treatment
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:                                                   // Eq. (39,40) of Hurley+ (2002) and end of section 7.2 of Hurley+ (2000) describe gradual growth of convective envelope over HG; we approximate it as convective here
            envelope = ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:
            envelope =  utils::Compare(Temperature() * TSOL, OPTIONS->ConvectiveEnvelopeTemperatureThreshold()) > 0 ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;  // Envelope is radiative if temperature exceeds fixed threshold, otherwise convective
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
double HG::ChooseTimestep(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double dtk      = 0.05 * (timescales(tBGB) - timescales(tMS));
    double dte      = timescales(tBGB) - p_Time;
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
 * STELLAR_TYPE ResolveEnvelopeLoss()
 *
 * @return                                      Stellar Type to which star shoule evolve after losing envelope
 */
STELLAR_TYPE HG::ResolveEnvelopeLoss(bool p_NoCheck) {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]          // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_NoCheck || utils::Compare(m_CoreMass, m_Mass) > 0) {                  // envelope loss

        m_Mass       = std::min(m_CoreMass, m_Mass);
        
        if (utils::Compare(m_Mass0, massCutoffs(MHeF)) < 0) {                   // star evolves to Helium White Dwarf

            stellarType  = STELLAR_TYPE::HELIUM_WHITE_DWARF;

            m_Radius     = HeWD::CalculateRadiusOnPhase_Static(m_Mass);
            m_Age        = 0.0;                                                 // see Hurley et al. 2000, discussion after eq 76
        }
        else {                                                                  // star evolves to Zero age Naked Helium Main Star

            stellarType  = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;

            m_Mass0      = m_Mass;
            m_Radius     = HeMS::CalculateRadiusAtZAMS_Static(m_Mass);          
            m_Luminosity = HeMS::CalculateLuminosityAtZAMS_Static(m_Mass);
            m_Age        = 0.0;                                                 // can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        }

    }

    return stellarType;

#undef massCutoffs
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
STELLAR_TYPE HG::EvolveToNextPhase() {
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]  // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType;

    if (utils::Compare(m_Mass0, massCutoffs(MFGB)) < 0) {
        stellarType = STELLAR_TYPE::FIRST_GIANT_BRANCH;
    }
    else {
        stellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;
    }    
    return stellarType;

#undef massCutoffs
}

/*
 * Update effective initial mass
 *
 * Per Hurley et al. 2000, section 7.1, the effective initial mass on the HG tracks the stellar mass -- unless it would yield an unphysical decrease in the core mass
 *
 *
 * void UpdateInitialMass()
 *
 */
void HG::UpdateInitialMass() {
    if (utils::Compare(m_CoreMass,HG::CalculateCoreMassOnPhase(m_Mass, m_Age)) < 0 ) {        //The current mass would yield a core mass larger than the current core mass -- i.e., no unphysical core mass decrease would ensue
        m_Mass0 = m_Mass;
    }
}
