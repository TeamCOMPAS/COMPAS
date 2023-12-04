#include "EAGB.h"
#include "HeMS.h"
#include "HeGB.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
void EAGB::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
#define timescales(x) p_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

    double p1   = gbParams(p) - 1.0;
    double q1   = gbParams(q) - 1.0;
    double p1_p = p1 / gbParams(p);
    double q1_q = q1 / gbParams(q);

    double tBAGB = CalculateLifetimeToBAGB(timescales(tHeI), timescales(tHe));
    double LBAGB = CalculateLuminosityAtBAGB(p_Mass);

    timescales(tinf1_FAGB) = tBAGB + ((1.0 / (p1 * gbParams(AHe) * gbParams(D))) * PPOW((gbParams(D) / LBAGB), p1_p));
    timescales(tMx_FAGB)   = timescales(tinf1_FAGB) - ((timescales(tinf1_FAGB) - tBAGB) * PPOW((LBAGB / gbParams(Lx)), p1_p));
    timescales(tinf2_FAGB) = timescales(tMx_FAGB) + ((1.0 / (q1 * gbParams(AHe) * gbParams(B))) * PPOW((gbParams(B) / gbParams(Lx)), q1_q));

#undef gbParams
#undef timescales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                LAMBDA CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the common envelope lambda parameter using the enhanced "Nanjing" prescription
 * from X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010)
 *
 * This function good for EAGB stars.
 *
 *
 * double CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind)
 *
 * @return                                      Nanjing lambda for use in common envelope
 */
double EAGB::CalculateLambdaNanjingEnhanced(const int p_MassInd, const int p_Zind) const {

	DBL_VECTOR maxBG    = {};                                                       // [0] = maxB, [1] = maxG
	DBL_VECTOR lambdaBG = {};                                                       // [0] = lambdaB, [1] = lambdaG
	DBL_VECTOR a        = {};                                                       // 0..5 a_coefficients
	DBL_VECTOR b        = {};                                                       // 0..5 b_coefficients
    double     Rmax     = std::numeric_limits<double>::max();                       // Upper R limit to applicability of Nanjing polynomials.

    switch(p_Zind) {
        // Pop. I metallicity
        case 1:
            switch(p_MassInd) {
                case 0: {
                    maxBG = { 2.5, 1.5 };
                    Rmax = 200.0;
                    double R_in = std::min(Rmax, m_Radius);
                    double tmp = 0.1 - ( R_in * 3.57E-04);
                    lambdaBG   = { tmp, tmp };
                    break;
                }
                case 1:
                    maxBG = { 4.0, 2.0 };
                    Rmax = 340.0;
                    a = { 0.88954, 0.0098 , -3.1411E-05 , 7.66979E-08,  0.0       , 0.0 };
                    b = { 0.48271, 0.00584, -6.22051E-05, 2.41531E-07, -3.1872E-10, 0.0 };
                    break;
                case 2:
                    maxBG = { 500.0, 10.0 };
                    Rmax = 400.0;
                    a = { -0.04669, 0.00764, -4.32726E-05, 9.31942E-08, 0.0        ,  0.0 };
                    b = {  0.44889, 0.01102, -6.46629E-05, 5.66857E-09, 7.21818E-10, -1.2201E-12 };
                    break;
                case 3:
                    maxBG = { 1000.0, 8.0 };
                    Rmax = 410.0;
                    a = { -0.37322, 0.00943, -3.26033E-05, 5.37823E-08, 0.0, 0.0 };
                    b = {  0.13153, 0.00984, -2.89832E-05, 2.63519E-08, 0.0, 0.0 };
                    break;
                case 4:
                    maxBG = { 1000.0, 8.0 };
                    Rmax = 430.0;
                    a = { -0.80011, 0.00992, -3.03247E-05,  5.26235E-08, 0.0, 0.0 };
                    b = { -0.00456, 0.00426,  4.71117E-06, -1.72858E-08, 0.0, 0.0 };
                    break;
                case 5:
                    maxBG = { 25.5, 5.0 };
                    Rmax = 440.0;
                    a = { -2.7714 ,  0.06467, -4.01537E-04,  7.98466E-07, 0.0, 0.0 };
                    b = {  0.23083, -0.00266,  2.21788E-05, -2.35696E-08, 0.0, 0.0 };
                    break;
                case 6:
                    maxBG = { 9.0, 3.0 };
                    Rmax = 420.0;
                    a = { -0.63266,  0.02054, -1.3646E-04 ,  2.8661E-07 , 0.0, 0.0 };
                    b = {  0.26294, -0.00253,  1.32272E-05, -7.12205E-09, 0.0, 0.0 };
                    break;
                case 7:
                    maxBG = { 7.0, 3.0 };
                    Rmax = 490.0;
                    a = { -0.1288 ,  0.0099 , -6.71455E-05,  1.33568E-07, 0.0, 0.0 };
                    b = {  0.26956, -0.00219,  7.97743E-06, -1.53296E-09, 0.0, 0.0 };
                    break;
                case 8:
                    maxBG = { 4.0, 2.0 };
                    Rmax = 530.0;
                    a = { 1.19804, -0.01961, 1.28222E-04, -3.41278E-07, 3.35614E-10, 0.0 };
                    b = { 0.40587, -0.0051 , 2.73866E-05, -5.74476E-08, 4.90218E-11, 0.0 };
                    break;
                case 9:
                    maxBG = { 3.0, 1.5 };
                    Rmax = 600.0;
                    a = { 0.3707 ,  2.67221E-04, -9.86464E-06, 2.26185E-08, 0.0, 0.0 };
                    b = { 0.25549, -0.00152    ,  3.35239E-06, 2.24224E-10, 0.0, 0.0 };
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
                    a = { -106.90553, 0.36469, -4.1472E-04 , 1.57349E-07, 0.0, 0.0 };
                    b = {  -39.93089, 0.13667, -1.55958E-04, 5.94076E-08, 0.0, 0.0 };
                    break;
                case 12:
                    maxBG = { 1.5, 1.0 };
                    Rmax = 1050.0;
                    a = { -154.70559, 0.46718, -4.70169E-04, 1.57773E-07, 0.0, 0.0 };
                    b = {  -65.39602, 0.19763, -1.99078E-04, 6.68766E-08, 0.0, 0.0 };
                    break;
                case 13:
                    maxBG = { 1.5, 1.0 };
                    Rmax = 1200.0;
                    a = { -260484.85724, 4.26759E+06, -2.33016E+07, 4.24102E+07, 0.0, 0.0 };
                    b = { -480055.67991, 7.87484E+06, -4.30546E+07, 7.84699E+07, 0.0, 0.0 };
                    break;
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
                case 0:
                    maxBG = { 2.0, 1.5 };
                    Rmax = 160.0;
                    a = { 0.24012, -0.01907, 6.09529E-04, -8.17819E-06, 4.83789E-08, -1.04568e-10 };
                    b = { 0.15504, -0.01238, 3.96633E-04, -5.3329E-06 , 3.16052E-08, -6.84288e-11 };
                    break;
                case 1:
                    maxBG = { 4.0, 2.0 };
                    Rmax = 350.0;
                    a = { 0.5452 ,  0.00212    , 6.42941E-05, -1.46783E-07, 0.0       ,  0.0 };
                    b = { 0.30594, -9.58858E-04, 1.12174E-04, -1.04079E-06, 3.4564E-09, -3.91536e-12 };
                    break;
                case 2:
                    maxBG = { 600.0, 2.0 };
                    Rmax = 400.0;
                    if (utils::Compare(m_Radius, 36.0) > 0 && utils::Compare(m_Radius, 53.0) < 0) lambdaBG = { 1.0, 1.0 };
                    else {
                        a = { -0.475  , -0.00328, 1.31101E-04, -6.03669E-07, 8.49549E-10, 0.0 };
                        b = {  0.05434,  0.0039 , 9.44609E-06, -3.87278E-08, 0.0        , 0.0 };
                    }
                    break;
                case 3:
                    maxBG = { 600.0, 2.0 };
                    Rmax = 410.0;
                    a = { -0.2106 , -0.01574, 2.01107E-04, -6.90334E-07, 7.92713E-10, 0.0 };
                    b = {  0.36779, -0.00991, 1.19411E-04, -3.59574E-07, 3.33957E-10, 0.0 };
                    break;
                case 4:
                    maxBG = { 10.0, 3.0 };
                    Rmax = 320.0;
                    a = { -0.12027,  0.01981, -2.27908E-04,  7.55556E-07, 0.0, 0.0 };
                    b = {  0.31252, -0.00527,  3.60348E-05, -3.22445E-08, 0.0, 0.0 };
                    break;
                case 5:
                    maxBG = { 4.0, 1.5 };
                    Rmax = 330.0;
                    a = { 0.26578,  0.00494, -7.02203E-05, 2.25289E-07, 0.0, 0.0 };
                    b = { 0.26802, -0.00248,  6.45229E-06, 1.69609E-08, 0.0, 0.0 };
                    break;
                case 6:
                    maxBG = { 2.5, 1.0 };
                    Rmax = 360.0;
                    a = { 0.8158 , -0.01633, 1.46552E-04, -5.75308E-07, 8.77711E-10, 0.0 };
                    b = { 0.26883, -0.00219, 4.12941E-06,  1.33138E-08, 0.0        , 0.0 };
                    break;
                case 7:
                    maxBG = { 2.0, 1.0 };
                    Rmax = 400.0;
                    a = { 0.74924, -0.01233, 9.55715E-05, -3.37117E-07, 4.67367E-10, 0.0 };
                    b = { 0.25249, -0.00161, 8.35478E-07,  1.25999E-08, 0.0        , 0.0 };
                    break;
                case 8:
                    maxBG = { 1.6, 1.0 };
                    Rmax = 440.0;
                    a = { 0.73147, -0.01076, 7.54308E-05, -2.4114E-07 , 2.95543E-10, 0.0 };
                    b = { 0.31951, -0.00392, 2.31815E-05, -6.59418E-08, 7.99575E-11, 0.0 };
                    break;
                case 9:
                    maxBG = { 1.6, 1.0 };
                    Rmax = 500.0;
                    a = { -9.26519,  0.08064, -2.30952E-04, 2.21986E-07, 0.0, 0.0 };
                    b = {  0.81491, -0.00161, -8.13352E-06, 1.95775E-08, 0.0, 0.0 };
                    break;
                case 10:
                    maxBG = { 1.6, 1.0 };
                    Rmax = 600.0;
                    if (utils::Compare(m_Radius, 390.0) > 0 && utils::Compare(m_Radius, 460.0) < 0) lambdaBG = { 0.08, 0.05 };
                    else {
                        a = { -51.15252, 0.30238, -5.95397E-04, 3.91798E-07, 0.0, 0.0 };
                        b = { -13.44   , 0.08141, -1.641E-04  , 1.106E-07  , 0.0, 0.0 };
                    }
                    break;
                case 11: {
                    maxBG = { 1.6, 1.0 };
                    Rmax = 650.0;
                    double R_in = std::min(Rmax, m_Radius);
                    if (utils::Compare(R_in, 480.0) >  0 && utils::Compare(R_in, 540.0) <  0) lambdaBG = { 0.06, 0.05 };
                    else if (utils::Compare(R_in, 540.0) >= 0 && utils::Compare(R_in, 650.0) <= 0) {
                        lambdaBG = { (R_in * 1.8E-03) - 0.88, (R_in* 9.1E-04) - 0.43 };   // JR: todo: inserted "=" in both cases here
                    }
                    else {
                        a = { -140.0   , 0.7126 , -0.00121    , 6.846E-07  , 0.0, 0.0 };
                        b = {  -44.1964, 0.22592, -3.85124E-04, 2.19324E-07, 0.0, 0.0 };
                    }
                    break;
                }
                case 12: {
                    maxBG = { 1.5, 1.0 };
                    Rmax = 750.0;
                    double R_in = std::min(Rmax, m_Radius);
                    if (utils::Compare(R_in, 560.0) >  0 && utils::Compare(R_in, 650.0) <  0) lambdaBG = { 0.1, 0.05 };
                    else if (utils::Compare(R_in, 650.0) >= 0 && utils::Compare(R_in, 750.0) <= 0) {
                        lambdaBG = { (R_in * 4.0E-03) - 2.5, (R_in * 1.5E-03) - 0.93 };    // JR: todo: inserted "=" in both cases here
                    }
                    else {
                        a = { -358.4    , 1.599  , -0.00238   , 1.178E-06  , 0.0, 0.0 };
                        b = { -118.13757, 0.52737, -7.8479E-04, 3.89585E-07, 0.0, 0.0 };
                    }
                    break;
                }
                case 13: {
                    maxBG = { 1.5, 1.0 };
                    Rmax = 900.0;
                    double R_in = std::min(Rmax, m_Radius);
                    if (utils::Compare(R_in, 725.0) >  0 && utils::Compare(R_in, 850.0) <  0) lambdaBG = { 0.1, 0.05 };
                    else if (utils::Compare(R_in, 850.0) >= 0 && utils::Compare(R_in, 900.0) <= 0) {
                        lambdaBG = { (R_in * 2.0E-03) - 1.6, (R_in * 1.0E-03) - 0.8 };    // JR: todo: inserted "=" in both cases here
                    }
                    else {
                        a = { -436.00777, 1.41375, -0.00153    , 5.47573E-07, 0.0, 0.0 };
                        b = { -144.53456, 0.46579, -4.99197E-04, 1.78027E-07, 0.0, 0.0 };
                    }
                    break;
                }
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
            double x  = std::min(m_Radius, Rmax);                           // Evaluate lambda with maximum allowed radius to prevent exceeding domain of the polynomial fits
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
 * This function good for EAGB stars.
 *
 *
 * double CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity)
 *
 * @return                                      Nanjing lambda for use in common envelope
 */
double EAGB::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const {

	DBL_VECTOR maxBG    = {};                                                       // [0] = maxB, [1] = maxG
	DBL_VECTOR lambdaBG = {};                                                       // [0] = lambdaB, [1] = lambdaG
	DBL_VECTOR a        = {};                                                       // 0..5 a_coefficients
	DBL_VECTOR b        = {};                                                       // 0..5 b_coefficients

    if (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT) > 0) {                 // Z>0.5 Zsun: popI
        if (utils::Compare(p_Mass, 1.5) < 0) {                                     // Should probably use effective mass m_Mass0 instead for Lambda calculations
            maxBG = { 2.5, 1.5 };
            if (utils::Compare(m_Radius, 200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else  {
                double tmp = 0.1 - (m_Radius * 3.57E-04);
                lambdaBG   = { tmp, tmp };
            }
        }
        else if (utils::Compare(p_Mass, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 340.0) > 0) lambdaBG = { 3.589970, 0.514132 };
            else {
                a = { 0.88954, 0.0098 , -3.1411E-05 , 7.66979E-08,  0.0       , 0.0 };
                b = { 0.48271, 0.00584, -6.22051E-05, 2.41531E-07, -3.1872E-10, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 3.5) < 0) {
            maxBG = { 500.0, 10.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 116.935557, 0.848808 };
            else {
                a = { -0.04669, 0.00764, -4.32726E-05, 9.31942E-08, 0.0        ,  0.0 };
                b = {  0.44889, 0.01102, -6.46629E-05, 5.66857E-09, 7.21818E-10, -1.2201E-12 };
            }
        }
        else if (utils::Compare(p_Mass, 4.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 410.0) > 0) lambdaBG = { 52.980056, 1.109736 };
            else {
                a = { -0.37322, 0.00943, -3.26033E-05, 5.37823E-08, 0.0, 0.0 };
                b = {  0.13153, 0.00984, -2.89832E-05, 2.63519E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 5.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 430.0) > 0) lambdaBG = { 109.593522, 1.324248 };
            else {
                a = { -0.80011, 0.00992, -3.03247E-05,  5.26235E-08, 0.0, 0.0 };
                b = { -0.00456, 0.00426,  4.71117E-06, -1.72858E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 6.5) < 0) {
            maxBG = { 25.5, 5.0 };
            if (utils::Compare(m_Radius, 440.0) > 0) lambdaBG = { 16.279603, 1.352166 };
            else {
                a = { -2.7714 ,  0.06467, -4.01537E-04,  7.98466E-07, 0.0, 0.0 };
                b = {  0.23083, -0.00266,  2.21788E-05, -2.35696E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 7.5) < 0) {
            maxBG = { 9.0, 3.0 };
            if (utils::Compare(m_Radius, 420.0) > 0) lambdaBG = { 5.133959, 1.004036 };
            else {
                a = { -0.63266,  0.02054, -1.3646E-04 ,  2.8661E-07 , 0.0, 0.0 };
                b = {  0.26294, -0.00253,  1.32272E-05, -7.12205E-09, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 8.5) < 0) {
            maxBG = { 7.0, 3.0 };
            if (utils::Compare(m_Radius, 490.0) > 0) lambdaBG = { 4.342985, 0.934659 };
            else {
                a = { -0.1288 ,  0.0099 , -6.71455E-05,  1.33568E-07, 0.0, 0.0 };
                b = {  0.26956, -0.00219,  7.97743E-06, -1.53296E-09, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 9.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 530.0) > 0) lambdaBG = { 2.441672, 0.702310 };
            else {
                a = { 1.19804, -0.01961, 1.28222E-04, -3.41278E-07, 3.35614E-10, 0.0 };
                b = { 0.40587, -0.0051 , 2.73866E-05, -5.74476E-08, 4.90218E-11, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 11.0) < 0) {
            maxBG = { 3.0, 1.5 };
            if (utils::Compare(m_Radius, 600.0) > 0) lambdaBG = { 1.842314, 0.593854 };
            else {
                a = { 0.3707 ,  2.67221E-04, -9.86464E-06, 2.26185E-08, 0.0, 0.0 };
                b = { 0.25549, -0.00152    ,  3.35239E-06, 2.24224E-10, 0.0, 0.0 };
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
            if (utils::Compare(m_Radius, 1000.0) > 0) lambdaBG = { 0.414200, 0.189008 };
            else {
                a = { -106.90553, 0.36469, -4.1472E-04 , 1.57349E-07, 0.0, 0.0 };
                b = {  -39.93089, 0.13667, -1.55958E-04, 5.94076E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1050.0) > 0) lambdaBG = { 0.2, 0.1 };
            else {
                a = { -154.70559, 0.46718, -4.70169E-04, 1.57773E-07, 0.0, 0.0 };
                b = {  -65.39602, 0.19763, -1.99078E-04, 6.68766E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else {
                a = { -260484.85724, 4.26759E+06, -2.33016E+07, 4.24102E+07, 0.0, 0.0 };
                b = { -480055.67991, 7.87484E+06, -4.30546E+07, 7.84699E+07, 0.0, 0.0 };
            }
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
    else {                                                                          // Z<=0.5 Zsun: popI and popII
        if (utils::Compare(p_Mass, 1.5) < 0) {
            maxBG = { 2.0, 1.5 };
            if (utils::Compare(m_Radius, 160.0) > 0) lambdaBG = { 0.05, 0.05 };
            else {
                a = { 0.24012, -0.01907, 6.09529E-04, -8.17819E-06, 4.83789E-08, -1.04568e-10 };
                b = { 0.15504, -0.01238, 3.96633E-04, -5.3329E-06 , 3.16052E-08, -6.84288e-11 };
            }
        }
        else if (utils::Compare(p_Mass, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 350.0) > 0) lambdaBG = { 2.868539, 0.389991 };
            else {
                a = { 0.5452 ,  0.00212    , 6.42941E-05, -1.46783E-07, 0.0       ,  0.0 };
                b = { 0.30594, -9.58858E-04, 1.12174E-04, -1.04079E-06, 3.4564E-09, -3.91536e-12 };
            }
        }
        else if (utils::Compare(p_Mass, 3.5) < 0) {
            maxBG = { 600.0, 2.0 };
                 if (utils::Compare(m_Radius, 400.0) > 0)                                      lambdaBG = { 398.126442, 0.648560 };
            else if (utils::Compare(m_Radius, 36.0) > 0 && utils::Compare(m_Radius, 53.0) < 0) lambdaBG = { 1.0, 1.0 };
            else {
                a = { -0.475  , -0.00328, 1.31101E-04, -6.03669E-07, 8.49549E-10, 0.0 };
                b = {  0.05434,  0.0039 , 9.44609E-06, -3.87278E-08, 0.0        , 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 4.5) < 0) {
            maxBG = { 600.0, 2.0 };
            if (utils::Compare(m_Radius, 410.0) > 0) lambdaBG = { 91.579093, 1.032432 };
            else {
                a = { -0.2106 , -0.01574, 2.01107E-04, -6.90334E-07, 7.92713E-10, 0.0 };
                b = {  0.36779, -0.00991, 1.19411E-04, -3.59574E-07, 3.33957E-10, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 5.5) < 0) {
            maxBG = { 10.0, 3.0 };
            if (utils::Compare(m_Radius, 320.0) > 0) lambdaBG = { 7.618019, 1.257919 };
            else {
                a = { -0.12027,  0.01981, -2.27908E-04,  7.55556E-07, 0.0, 0.0 };
                b = {  0.31252, -0.00527,  3.60348E-05, -3.22445E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 6.5) < 0) {
            maxBG = { 4.0, 1.5 };
            if (utils::Compare(m_Radius, 330.0) > 0) lambdaBG = { 2.390575, 0.772091 };
            else {
                a = { 0.26578,  0.00494, -7.02203E-05, 2.25289E-07, 0.0, 0.0 };
                b = { 0.26802, -0.00248,  6.45229E-06, 1.69609E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 7.5) < 0) {
            maxBG = { 2.5, 1.0 };
            if (utils::Compare(m_Radius, 360.0) > 0) lambdaBG = { 1.878174, 0.646353 };
            else {
                a = { 0.8158 , -0.01633, 1.46552E-04, -5.75308E-07, 8.77711E-10, 0.0 };
                b = { 0.26883, -0.00219, 4.12941E-06,  1.33138E-08, 0.0        , 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 8.5) < 0) {
            maxBG = { 2.0, 1.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 1.517662, 0.553169 };
            else {
                a = { 0.74924, -0.01233, 9.55715E-05, -3.37117E-07, 4.67367E-10, 0.0 };
                b = { 0.25249, -0.00161, 8.35478E-07,  1.25999E-08, 0.0        , 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 9.5) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 440.0) > 0) lambdaBG = { 1.136394, 0.478963 };
            else {
                a = { 0.73147, -0.01076, 7.54308E-05, -2.4114E-07 , 2.95543E-10, 0.0 };
                b = { 0.31951, -0.00392, 2.31815E-05, -6.59418E-08, 7.99575E-11, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 11.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 500.0) > 0) lambdaBG = { 1.068300, 0.424706 };
            else {
                a = { -9.26519,  0.08064, -2.30952E-04, 2.21986E-07, 0.0, 0.0 };
                b = {  0.81491, -0.00161, -8.13352E-06, 1.95775E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 13.0) < 0) {
            maxBG = { 1.6, 1.0 };
                 if (utils::Compare(m_Radius, 600.0) > 0)                                        lambdaBG = { 0.537155, 0.211105 };
            else if (utils::Compare(m_Radius, 390.0) > 0 && utils::Compare(m_Radius, 460.0) < 0) lambdaBG = { 0.08, 0.05 };
            else {
                a = { -51.15252, 0.30238, -5.95397E-04, 3.91798E-07, 0.0, 0.0 };
                b = { -13.44   , 0.08141, -1.641E-04  , 1.106E-07  , 0.0, 0.0 };
            }
        }
        else if (p_Mass < 15.0) {
            maxBG = { 1.6, 1.0 };
                 if (utils::Compare(m_Radius, 650.0) >  0)                                         lambdaBG = { 0.3, 0.160696 };
            else if (utils::Compare(m_Radius, 480.0) >  0 && utils::Compare(m_Radius, 540.0) <  0) lambdaBG = { 0.06, 0.05 };
            else if (utils::Compare(m_Radius, 540.0) >= 0 && utils::Compare(m_Radius, 650.0) <= 0) lambdaBG = { (m_Radius * 1.8E-03) - 0.88, (m_Radius * 9.1E-04) - 0.43 };   // JR: todo: inserted "=" in both cases here
            else {
                a = { -140.0   , 0.7126 , -0.00121    , 6.846E-07  , 0.0, 0.0 };
                b = {  -44.1964, 0.22592, -3.85124E-04, 2.19324E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 750.0) >  0)                                         lambdaBG = { 0.5, 0.204092 };
            else if (utils::Compare(m_Radius, 560.0) >  0 && utils::Compare(m_Radius, 650.0) <  0) lambdaBG = { 0.1, 0.05 };
            else if (utils::Compare(m_Radius, 650.0) >= 0 && utils::Compare(m_Radius, 750.0) <= 0) lambdaBG = { (m_Radius * 4.0E-03) - 2.5, (m_Radius * 1.5E-03) - 0.93 };    // JR: todo: inserted "=" in both cases here
            else {
                a = { -358.4    , 1.599  , -0.00238   , 1.178E-06  , 0.0, 0.0 };
                b = { -118.13757, 0.52737, -7.8479E-04, 3.89585E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 900.0) >  0)                                         lambdaBG = { 0.2, 0.107914 };
            else if (utils::Compare(m_Radius, 725.0) >  0 && utils::Compare(m_Radius, 850.0) <  0) lambdaBG = { 0.1, 0.05 };
            else if (utils::Compare(m_Radius, 850.0) >= 0 && utils::Compare(m_Radius, 900.0) <= 0) lambdaBG = { (m_Radius * 2.0E-03) - 1.6, (m_Radius * 1.0E-03) - 0.8 };    // JR: todo: inserted "=" in both cases here
            else {
                a = { -436.00777, 1.41375, -0.00153    , 5.47573E-07, 0.0, 0.0 };
                b = { -144.53456, 0.46579, -4.99197E-04, 1.78027E-07, 0.0, 0.0 };
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
        if (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT) > 0 &&
            (utils::Compare(p_Mass, 1.5) < 0 || (utils::Compare(p_Mass, 25.0) < 0 && utils::Compare(p_Mass, 18.0) >= 0)) ) {
            double x  = (m_Mass - m_CoreMass) / m_Mass;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x3 * x2;

            double y1 = a[0] + (a[1] * x) + (a[2] * x2) + (a[3] * x3) + (a[4] * x4) + (a[5] * x5);
            double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

            lambdaBG = { 1.0 / y1, 1.0 / y2 };
        }
        else if ( (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT) > 0  && utils::Compare(p_Mass, 2.5) >= 0 && utils::Compare(p_Mass, 5.5) < 0) ||
                  (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT) <= 0 && utils::Compare(p_Mass, 2.5) >= 0 && utils::Compare(p_Mass, 4.5) < 0)) {
            double x  = m_Radius;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x3 * x2;

            double y1 = a[0] + (a[1] * x) + (a[2] * x2) + (a[3] * x3) + (a[4] * x4) + (a[5] * x5);
            double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

            lambdaBG = { pow(10.0, y1), y2 };
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
 * Calculate luminosity on the Early Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 37
 *
 *
 * double CalculateLuminosityOnPhase(const double p_CoreMass)
 *
 * @param   [IN]    p_CoreMass                  Core mass in Msol
 * @return                                      Luminosity on the Early Asymptotic Giant Branch in Lsol
 */
double EAGB::CalculateLuminosityOnPhase(const double p_CoreMass) const {
    return CalculateLuminosityGivenCoreMass(p_CoreMass);
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
double EAGB::CalculateRemnantLuminosity() const {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function
    return HeGB::CalculateLuminosityOnPhase_Static(m_COCoreMass, gbParams(B), gbParams(D));
#undef gbParams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate radius on the Early Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 74
 *
 *
 * double CalculateRadiusOnPhase_Static(const double      p_Mass,
 *                                      const double      p_Luminosity,
 *                                      const double      p_MHeF,
 *                                      const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Radius on the Early Asymptotic Giant Branch in Rsol
 *
 * p_MHeF and p_BnCoefficients passed as parameters so function can be declared static
 */
double EAGB::CalculateRadiusOnPhase_Static(const double      p_Mass,
                                           const double      p_Luminosity,
                                           const double      p_MHeF,
                                           const DBL_VECTOR &p_BnCoefficients) {
#define b p_BnCoefficients  // for convenience and readability - undefined at end of function

    // calculate radius constant A (Hurley et al. 2000, eq 74)
    // and coefficient b(50)
    double A;
    double b50;

    if (utils::Compare(p_Mass, p_MHeF) >= 0) {
        A = std::min((b[51] * PPOW(p_Mass, -b[52])), (b[53] * PPOW(p_Mass, -b[54])));
        b50 = b[55] * b[3];
    }
    else if (utils::Compare(p_Mass, (p_MHeF - 0.2)) <= 0) {
        A = b[56] + (b[57] * p_Mass);
        b50 = b[3];
    }
    else {  // Linear interpolation between end points
        double x1        = p_MHeF - 0.2;
        double x2        = p_MHeF;
        double x2_x1     = x2 - x1;

        double y1        = b[56] + (b[57] * x1);
        double y2        = std::min((b[51] * PPOW(x2, -b[52])), (b[53] * PPOW(x2, -b[54])));
        double gradient  = (y2 - y1) / x2_x1;
        double intercept = y2 - (gradient * x2);
               A         = (gradient * p_Mass) + intercept;

               y1        = b[3];
               y2        = b[55] * b[3];
               gradient  = (y2 - y1) / x2_x1;
               intercept = y2 - (gradient * x2);
               b50       = (gradient * p_Mass) + intercept;
    }

    // now calculate the radius
    return A * (PPOW(p_Luminosity, b[1]) + (b[2] * PPOW(p_Luminosity, b50)));

#undef b
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
double EAGB::CalculateRemnantRadius() const {
    double R1, R2;
    std::tie(R1, R2) = HeGB::CalculateRadiusOnPhase_Static(m_HeCoreMass, CalculateRemnantLuminosity());

    return std::min(R1, R2);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate CO core mass on the Early Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 39, modified as described in Section 5.4
 *
 *
 * double CalculateCOCoreMassOnPhase(const double p_Time)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Time                      Time in Myr
 * @return                                      Core mass on the Early Asymptotic Giant Branch in Msol
 */
double EAGB::CalculateCOCoreMassOnPhase(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)] // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    return utils::Compare(p_Time, timescales(tMx_FAGB)) <= 0
            ? PPOW((gbParams(p) - 1.0) * gbParams(AHe) * gbParams(D) * (timescales(tinf1_FAGB) - p_Time), 1.0 / (1.0 - gbParams(p)))
            : PPOW((gbParams(q) - 1.0) * gbParams(AHe) * gbParams(B) * (timescales(tinf2_FAGB) - p_Time), 1.0 / (1.0 - gbParams(q)));

#undef gbParams
#undef timescales
}


/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * According to Hurley et al. 2000
 *
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double EAGB::CalculateMassLossRateHurley() {
    double rateNJ = CalculateMassLossRateNieuwenhuijzenDeJager();
    double rateKR = CalculateMassLossRateKudritzkiReimers();
    double rateVW = CalculateMassLossRateVassiliadisWood();
    double rateWR = CalculateMassLossRateWolfRayet(m_Mu);
    double dominantRate;
    m_DominantMassLossRate = MASS_LOSS_TYPE::GB;

    if (utils::Compare(rateNJ, rateKR) > 0) {
        dominantRate = rateNJ;
    } else {
        dominantRate = rateKR;
    }

    if (utils::Compare(rateVW, dominantRate) > 0) {
        dominantRate = rateVW;
    }

    if (utils::Compare(rateWR, dominantRate) > 0) {
        m_DominantMassLossRate = MASS_LOSS_TYPE::WR;
        dominantRate = rateWR;
    }
    return dominantRate;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the lifetime to second dredge up (nb there is no explicit first)
 * This is the time of transition between the EAGB and the TPAGB
 *
 * Hurley et al. 2000, eqs 70, 71 & 72
 *
 *
 * double CalculateLifetimeTo2ndDredgeUp(const double p_Tinf1_FAGB, const double p_Tinf2_FAGB)
 *
 * @param   [IN]    p_Tinf1_FAGB                tinf1_FAGB (Timescales[TS::tinf1_FAGB]) - Early Asymptotic Giant Branch tinf1
 * @param   [IN]    p_Tinf2_FAGB                tinf2_FAGB (Timescales[TS::tinf2_FAGB]) - Early Asymptotic Giant Branch tinf2
 * @return                                      Lifetime to second dredge up (tDU)
 */
double EAGB::CalculateLifetimeTo2ndDredgeUp(const double p_Tinf1_FAGB, const double p_Tinf2_FAGB) const {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    double p1  = gbParams(p) - 1.0;
    double q1  = gbParams(q) - 1.0;

    double LDU = CalculateLuminosityGivenCoreMass(gbParams(McDU));

    return utils::Compare(LDU, gbParams(Lx)) <= 0
            ? p_Tinf1_FAGB - (1.0 / (p1 * gbParams(AHe) * gbParams(D))) * PPOW((gbParams(D) / LDU), (p1 / gbParams(p)))
            : p_Tinf2_FAGB - (1.0 / (q1 * gbParams(AHe) * gbParams(B))) * PPOW((gbParams(B) / LDU), (q1 / gbParams(q)));

#undef gbParams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
double EAGB::ChooseTimestep(const double p_Time) const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    double dtk = utils::Compare(p_Time, timescales(tMx_FAGB)) <= 0
                    ? 0.02 * (timescales(tinf1_FAGB) - p_Time)
                    : 0.02 * (timescales(tinf2_FAGB) - p_Time);

    double dte      = dtk;                                          // FINISH ME        JR: todo: ? Finish how?
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
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE EAGB::ResolveEnvelopeLoss(bool p_NoCheck) {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType = m_StellarType;

    if(ShouldEnvelopeBeExpelledByPulsations())          { m_EnvelopeJustExpelledByPulsations = true; }

    if (p_NoCheck || utils::Compare(m_HeCoreMass, m_Mass) >= 0 || m_EnvelopeJustExpelledByPulsations) {              // Envelope lost, form an evolved naked helium giant

        m_Mass     = std::min(m_HeCoreMass, m_Mass);
        m_HeCoreMass  = m_Mass;
        m_Mass0    = m_Mass;
        m_CoreMass = m_COCoreMass;

        double p1   = gbParams(p) - 1.0;
        double q1   = gbParams(q) - 1.0;
        double p1_p = p1 / gbParams(p);
        double q1_q = q1 / gbParams(q);

        timescales(tHeMS) = HeMS::CalculateLifetimeOnPhase_Static(m_Mass);  // calculate common values

        double LTHe = HeMS::CalculateLuminosityAtPhaseEnd_Static(m_Mass);

        timescales(tinf1_HeGB) = timescales(tHeMS) + (1.0 / ((p1 * gbParams(AHe) * gbParams(D))) * PPOW((gbParams(D) / LTHe), p1_p));
        timescales(tx_HeGB) = timescales(tinf1_HeGB) - (timescales(tinf1_HeGB) - timescales(tHeMS)) * PPOW((LTHe / gbParams(Lx)), p1_p);
        timescales(tinf2_HeGB) = timescales(tx_HeGB) + ((1.0 / (q1 * gbParams(AHe) * gbParams(B))) * PPOW((gbParams(B) / gbParams(Lx)), q1_q));

        m_Age      = HeGB::CalculateAgeOnPhase_Static(m_Mass, m_COCoreMass, timescales(tHeMS), m_GBParams);
        HeHG::CalculateGBParams_Static(m_Mass0, m_Mass, LogMetallicityXi(), m_MassCutoffs, m_AnCoefficients, m_BnCoefficients, m_GBParams);  // IM: order of type change and parameter updates to be revisited (e.g., why not just CalculateGBParams(m_Mass0, m_GBParams)?)  JR: static function has no access to class variables
        m_Luminosity = HeGB::CalculateLuminosityOnPhase_Static(m_COCoreMass, gbParams(B), gbParams(D));

        double R1, R2;
        std::tie(R1, R2) = HeGB::CalculateRadiusOnPhase_Static(m_Mass, m_Luminosity);
        if (utils::Compare(R1, R2) < 0) {
            m_Radius    = R1;
            stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;
        }
        else {
            m_Radius    = R2;
            stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH;     // Has a deep convective envelope
        }
    }

    return stellarType;

#undef gbParams
#undef timescales
}


/*
 * Determine if this phase should be skipped entirely
 *
 *
 * bool ShouldSkipPhase()
 *
 * @return                                      Boolean flag: true if this phase should be skipped, false if not
 */
bool EAGB::ShouldSkipPhase() const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

    double McCOBAGB = CalculateCOCoreMassOnPhase(timescales(tHeI) + timescales(tHe));
    double McSN     = std::max(gbParams(McSN), 1.05 * McCOBAGB);                                // hack from Hurley fortran code, doesn't seem to be in the paper   JR: do we know why?

    return (utils::Compare(McSN, m_COCoreMass) < 0);                                            // skip phase if core is heavy enough to go supernova

#undef gbParams
#undef timescales
}


/*
 * Determine if evolution should continue on this phase, or whether evolution
 * on this phase should end (and so evolve to next phase)
 *
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                                      Boolean flag: true if evolution on this phase should continue, false if not
 */
bool EAGB::ShouldEvolveOnPhase() const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

	double tDU = CalculateLifetimeTo2ndDredgeUp(timescales(tinf1_FAGB), timescales(tinf2_FAGB));
    return ((utils::Compare(m_Age, tDU) < 0 || utils::Compare(gbParams(McBAGB), MCBUR2) >= 0) && !ShouldEnvelopeBeExpelledByPulsations());

#undef gbParams
#undef timescales
}


/*
 * Determine if star should continue evolution as a Supernova
 *
 *
 * bool IsSupernova()
 *
 * @return                                      Boolean flag: true if star has gone Supernova, false if not
 */
bool EAGB::IsSupernova() const {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

    double McCOBAGB = CalculateCOCoreMassOnPhase(timescales(tHeI) + timescales(tHe));
    double McSN     = std::max(gbParams(McSN), 1.05 * McCOBAGB);                                // hack from Hurley fortran code, doesn't seem to be in the paper   JR: do we know why?

    return (utils::Compare(McSN, m_COCoreMass) < 0);                                            // core is heavy enough to go Supernova

#undef gbParams
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
STELLAR_TYPE EAGB::EvolveToNextPhase() {
    return STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH;
}
