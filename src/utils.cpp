#include <iostream>
#include <stdarg.h>
#include <fstream>
#include <algorithm>
#include <cstring>
#include "utils.h"
#include "Rand.h"
#include "changelog.h"

/*
 * utility functions that don't belong in any class
 *
 */
namespace utils {

    /*
     * Announce COMPAS
     * 
     * Constructs and returns a splash string.  Prints string to stdout if required.
     *
     *
     * void SplashScreen(const bool p_Print)
     * 
     * @param   [IN]    p_Print             Boolean indicating whether splash string should be printed.  Default is TRUE
     * @return                              Splash string
     */
    std::string SplashScreen(const bool p_Print) {

        // Construct the splash string
        std::string splashString = "\nCOMPAS v" + 
                                   VERSION_STRING + 
                                   "\nCompact Object Mergers: Population Astrophysics and Statistics"
                                   "\nby Team COMPAS (http://compas.science/index.html)"
                                   "\nA binary star simulator\n";

        if (p_Print) std::cout << splashString << std::endl;    // print the splash string if required

        return splashString;                                    // return the splash string
    }


    /*
     * Compare floating-point numbers with tolerance
     *
     * Absolute and relative tolerance can be different - see constants.h
     * Set relative tolerance = 0.0 to always use absolute
     * Set absolute tolerance = 0.0 to always use relative
     * Set both to zero for no tolerance - or #undef COMPARE_WITH_TOLERANCE for performance
     *
     *
     * int Compare(const double p_X, const double p_Y)
     *
     * @param   [IN]    p_X                 Floating-point value to be compared
     * @param   [IN]    p_Y                 Floating-point value to be compared
     * @return                              Integer indicating result of comparison:
     *                                         -1 indicates p_X is less than p_Y
     *                                          0 indicates equality
     *                                          1 indicates p_X is greater than p_Y
     */
    int Compare(const double p_X, const double p_Y) {
    #ifdef COMPARE_WITH_TOLERANCE
        return (std::abs(p_X - p_Y) <= std::max(FLOAT_TOLERANCE_ABSOLUTE, FLOAT_TOLERANCE_RELATIVE * std::max(std::abs(p_X), fabs(p_Y)))) ? 0 : (p_X < p_Y ? -1 : 1);
    #else
        return (p_X == p_Y) ? 0 : (p_X < p_Y ? -1 : 1);
    #endif
    }


    /*
     * Determine if filename exists - input parameter is character array
     *
     *
     * bool FileExists(const char *p_Filename)
     *
     * @param   [IN]    p_Filename                  Filename to be checked - should be fully qualified
     * @return                                      Boolean indicating whether file exists
     */
    bool FileExists(const char *p_Filename) {
        std::ifstream ifile(p_Filename);
        return (bool)ifile;
    }


    /*
     * Determine if filename exists - input parameter is std::string
     *
     *
     * bool FileExists(const std::string& p_Filename)
     *
     * @param   [IN]    p_Filename                  Filename to be checked - should be fully qualified
     * @return                                      Boolean indicating whether file exists
     */
    bool FileExists(const std::string& p_Filename) {
        return FileExists(p_Filename.c_str());
    }


    /*
     * Upshift lowercase characters in string to uppercase
     *
     * This only works with ASCII data, but I think that's all we need
     *
     *
     * bool std::string ToUpper(std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to be upshifted
     * @return                                      Upshifted string
     */
    std::string ToUpper(std::string p_Str) {
        std::transform(p_Str.begin(), p_Str.end(), p_Str.begin(), ::toupper);
        return p_Str;
    }


    /*
     * Downshift uppercase characters in string to lowercase
     *
     * This only works with ASCII data, but I think that's all we need
     *
     *
     * bool std::string ToLower(std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to be downshifted
     * @return                                      Downshifted string
     */
    std::string ToLower(std::string p_Str) {
        std::transform(p_Str.begin(), p_Str.end(), p_Str.begin(), ::tolower);
        return p_Str;
    }


    /*
     * Case-insensitive comparison of strings
     *
     * This only works with ASCII data, but I think that's all we need
     * Note that std::string has an == operator to test for equality (actually calls std::strcmp)
     *
     *
     * bool Equals(std::string p_Str1, std::string p_Str2)
     *
     * @param   [IN]    p_Str1                      String to be compared
     * @param   [IN]    p_Str2                      String to be compared
     * @return                                      Boolean indicating equality (true = equal)
     */
    bool Equals(std::string p_Str1, std::string p_Str2) {
        std::transform(p_Str1.begin(), p_Str1.end(), p_Str1.begin(), ::tolower);
        std::transform(p_Str2.begin(), p_Str2.end(), p_Str2.begin(), ::tolower);
        return (std::strcmp(p_Str1.c_str(), p_Str2.c_str()) == 0);
    }


    /*
     * Pads string to specified length by prepending the string with "0"
     *
     * This only works with ASCII data, but I think that's all we need
     * Note that std::string has an == operator to test for equality (actually calls std::strcmp)
     *
     *
     * std::string PadLeadingZeros(const std::string p_Str, const std::size_t p_MaxLength)
     *
     * @param   [IN]    p_Str                       String to be padded with leading "0"s
     * @param   [IN]    p_MaxLength                 The required length of the resultant string compared
     * @return                                      String padded with leading "0"s - will be unchanged from input string if length alread >= required length
     */
    std::string PadLeadingZeros(const std::string p_Str, const std::size_t p_MaxLength) {
        return (p_Str.length() < p_MaxLength) ? std::string(p_MaxLength - p_Str.length(), '0') + p_Str : p_Str;
    }


    /*
     * Centre-justfies string to specified width by prepending and appending spaces.
     * Extra space will be at the end of the string if necessary.
     *
     *
     * std::string CentreJustify(const std::string p_Str, const std::size_t p_Width)
     *
     * @param   [IN]    p_Str                       String to be centre-justified
     * @param   [IN]    p_Width                     The required width of the resultant string
     * @return                                      String padded with leading and trailing spaces so as to (as close as possible) centre-justify p_Str
     *                                              The string returned will always be p_Width characters in length
     */
     std::string CentreJustify(const std::string p_Str, const std::size_t p_Width) {

        std::string result = p_Str;                                                                 // default is no change

        if (p_Str.length() < p_Width) {                                                             // p_Str < field width
            int numLeadingSpaces = (p_Width - p_Str.length()) / 2;                                  // number of spaces to add at start - half the deficit
            int numTralingSpaces = p_Width - p_Str.length() - numLeadingSpaces;                     // number of spaces to adda t the end - whatever is left
            std::string leadingSpaces(numLeadingSpaces, ' ');                                       // blank string to prepend
            std::string trailingSpaces(numTralingSpaces, ' ');                                      // blank string to append
            result = leadingSpaces + p_Str + trailingSpaces;                                        // add leading and trailing spaces to p_Str
        }

        return result;
    }


    /*
     * Formats value per sprintf() and returns string
     *
     * From https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf/49812018#49812018
     *
     *
     * string vFormat(const char* const p_zcFormat, ...)
     *
     * @param   [IN]    p_zcFormat                  Format string
     * @param   [IN]    ...                         Parameters to be formatted
     * @return                                      Formatted string
     */
    const std::string vFormat(const char* const p_zcFormat, ...) {

        // initialize use of the variable argument array
        va_list vaArgs;
        va_start(vaArgs, p_zcFormat);

        // reliably acquire the size from a copy of the variable argument array
        // and a functionally reliable call to mock the formatting
        va_list vaArgsCopy;
        va_copy(vaArgsCopy, vaArgs);
        const int iLen = std::vsnprintf(NULL, 0, p_zcFormat, vaArgsCopy);
        va_end(vaArgsCopy);

        // return a formatted string without risking memory mismanagement
        // and without assuming any compiler or platform specific behavior
        std::vector<char> zc(iLen + 1);
        std::vsnprintf(zc.data(), zc.size(), p_zcFormat, vaArgs);
        va_end(vaArgs);

        return std::string(zc.data(), iLen);
    }


    /*
     * Solve quadratic Ax^2 + Bx + C
     *
     * Returns either root, depending on discriminant will return:
     *
     *    0.0               if 0 roots
     *    root              if 1 root
     *    max(root1, root2) if 2 roots
     *
     *
     * std::tuple<ERROR, double> SolveQuadratic(const double p_A, const double p_B, double p_C)
     *
     * @param   [IN]    p_A                       Coefficient of x^2
     * @param   [IN]    p_B                       Coefficient of x^1
     * @param   [IN]    p_C                       Coefficient of x^0 (Constant)
     * @return                                    Root found (see above)
     * @return                                    Tuple containing (in order): error value, root found (see above)
     *                                            The error value returned will be:
     *                                                ERROR::NONE if no error occurred
     *                                                ERROR::NO_REAL_ROOTS if the equation has no real roots
     *                                              If the error returned is not ERROR:NONE, use root returned at your own risk
     */
    std::tuple<ERROR, double> SolveQuadratic(const double p_A, const double p_B, double p_C) {

        ERROR error = ERROR::NONE;                                  // error

        double discriminant = (p_B * p_B) - (4.0 * p_A * p_C);      // d = B^2 - 4AC

        double root = 0.0;                                          // root found

        // JR: check < 0 first so don't have to check = 0.0 (will almost never happen after calculation - need epsilon)
        if (discriminant < 0.0) {                                   // no real roots? (leave this as an absolute compare)
            error = ERROR::NO_REAL_ROOTS;                           // no real roots - set error
        }
        else if (discriminant > 0.0) {                              // 2 real roots (leae this as an absolute compare)
            double sqrtD = sqrt(discriminant);
            double A2    = p_A + p_A;

            double root1 = (-p_B + sqrtD) / A2;                     // (-B + SQRT(B^2 - 4AC)) / 2A
            double root2 = (-p_B - sqrtD) / A2;                     // (-B - SQRT(B^2 - 4AC)) / 2A

            root         = std::max(root1, root2);
        }
        else {                                                      // 1 real root
            root = -p_B / (p_A + p_A);                              // -B / 2A,discriminant = 0.0
        }

        return std::make_tuple(error, root);
    }


    /*
     * Sample from tabulated CDF
     *
     * Finds Y given X and an ordered map<double, double> of (X, Y) pairs
     * Uses simple linear interpolation
     *
     *
     * double SampleFromTabulatedCDF(const double p_X, const std::map<double, double> p_Table)
     *
     * @param   [IN]    p_X                 The X value for which Y is to be calculated
     * @param   [IN]    p_Table             The table to interpolate on (Ordered map<double, double> of (X, Y) pairs)
     * @return                              Interpolated Y value
     */
    double SampleFromTabulatedCDF(const double p_X, const std::map<double, double> pTable) {

        double yInterp = 0.0;

        std::map<double, double>::const_iterator iter;

        iter        = pTable.begin();
        double xMin = iter->first;

        iter        = pTable.end();
        double xMax = (--iter)->first;

        if (pTable.size() > 0 && (p_X >= xMin && p_X <= xMax)) {                                // sanity check - map must not be empty, and p_X must be in [xMin, xMax]  (Leave these as absolute compares)

            for (iter = pTable.begin(); (iter != pTable.end() && iter->first < p_X); ++iter);   // leave this as an absolute compare

            if (iter == pTable.begin()) {
                yInterp = iter->second;
            }
            else if (iter == pTable.end()) {
                yInterp = (--iter)->second;
            }
            else {

                double xAbove = iter->first;
            double yAbove = iter->second;
                    double xBelow = (--iter)->first;
                double yBelow = iter->second;

                double gradient = (yAbove - yBelow) / (xAbove - xBelow);

                yInterp = yBelow + ((p_X - xBelow) * gradient);
            }
        }

        return yInterp;
    }


    /*
     * Inverse sample from tabulated CDF
     *
     * Finds X given Y and an ordered map<double, double> of (X, Y) pairs
     * Uses simple linear interpolation
     *
     *
     * double InverseSampleFromTabulatedCDF(const double p_Y, const std::map<double, double> p_Table)
     *
     * @param   [IN]    p_Y                 The Y value for which X is to be calculated
     * @param   [IN]    p_Table             The table to interpolate on (Ordered map<double, double> of (X, Y) pairs)
     * @return                              Interpolated X value
     */
    double InverseSampleFromTabulatedCDF(const double p_Y, const std::map<double, double> p_Table) {

        double xInterp = 0.0;

        if (p_Table.size() > 0 && (p_Y >= 0.0 && p_Y < 1.0)) {                                  // sanity check - map must not be empty, and p_Y must in [0.0, 1.0)      (Levae these as absolute compares)

            std::map<double, double>::const_iterator iter;

            iter        = p_Table.begin();
            double yMin = iter->second;

            iter        = p_Table.end();
            double yMax = (--iter)->second;

            double y = yMin + (p_Y * (std::max(yMin, yMax) - std::min(yMax, yMin)));            // normalise y - clamp to [yMin, yMax)

            for (iter = p_Table.begin(); (iter != p_Table.end() && iter->second < y); ++iter);  // leave this as an absolute compare

            double xAbove = iter->first;
            double yAbove = iter->second;

            if (iter == p_Table.begin()) {
                xInterp = xAbove;
            }
            else {
                iter--;
                double xBelow = iter->first;
                double yBelow = iter->second;

                double gradient = (yAbove - yBelow) / (xAbove - xBelow);

                xInterp = xBelow + ((y - yBelow) / gradient);
            }
        }

        return xInterp;
    }


    /*
     * Draw sample from a power law distribution p(x) ~ x^(n) between p_Xmin and p_Xmax
     *
     * double InverseSampleFromPowerLaw(const double p_Power, const double p_Xmax, const double p_Xmin)
     *
     * @param   [IN]    p_Power             The power for the power law
     * @param   [IN]    p_Xmax              Maximum of the X-interval from which to sample
     * @param   [IN]    p_Xmin              Minimum of the X-interval from which to sample
     * @return                              Drawn sample
     */
    double InverseSampleFromPowerLaw(const double p_Power, const double p_Xmax, const double p_Xmin) {

        double result;

        double rand = RAND->Random();                   // Draw a random number between 0 and 1

        if (p_Power == -1.0) {                          // JR: todo: find a better way of doing this
            result = exp(rand * log(p_Xmax / p_Xmin)) * p_Xmin;
        }
        else {
            double powerPlus1     = p_Power + 1.0;
            double min_powerPlus1 = utils::POW(p_Xmin, powerPlus1);

            result = utils::POW((rand * (utils::POW(p_Xmax, powerPlus1) - min_powerPlus1) + min_powerPlus1), 1.0 / powerPlus1);
        }

        return result;
    }


    /*
     * Convert a period in days to a semi-major axis in AU
     *
     *
     * double ConvertPeriodInDaysToSemiMajorAxisInAU(const double p_Mass1, const double p_Mass2, const double p_Period)
     *
     * @param   [IN]    p_Mass1             The mass of star 1 (Msol)
     * @param   [IN]    p_Mass2             The mass of star 2 (Msol)
     * @param   [IN]    p_Period            Orbital period (days)
     * @return                              Semi-major axis in AU
     */
    double ConvertPeriodInDaysToSemiMajorAxisInAU(const double p_Mass1, const double p_Mass2, const double p_Period) {

        double a_cubed_SI_top    = G * ((p_Mass1 * MSOL_TO_KG) + (p_Mass2 * MSOL_TO_KG)) * p_Period * p_Period * SECONDS_IN_DAY * SECONDS_IN_DAY;
        double a_cubed_SI_bottom = 4.0 * M_PI * M_PI;
        double a_cubed_SI        = a_cubed_SI_top / a_cubed_SI_bottom;
        double a_SI              = utils::POW(a_cubed_SI, 1.0 / 3.0);

        return a_SI / AU;
    }


    /*
     * Calculate x^y where x is double and y is an integer
     *
     * faster than utils::POW() for integer exponent
     *
     *
     * double intPow(const double p_Base, const int p_Exponent)
     *
     * @param   [IN]    p_Base              Base - number to be raised to integer power
     * @param   [IN]    p_Exponent          Exponent - integer to wich base should be raised
     * @return                              Base ^ Exponent
     */
    double intPow(const double p_Base, const int p_Exponent) {

        double result = 1.0;                                            // for exponent = 0

        if (p_Exponent != 0) {                                          // exponent not zero
            int times = abs(p_Exponent);                                // number of times to multiply
            for (int i = 0; i < times; i++) result *= p_Base;           // multiply
        }

        return p_Exponent < 0 ? 1.0 / result : result;                  // invert if negative exponent
    }


    /*
     * Determines if the stellar type passed is one of a list of stellar types passed
     *
     *
     * bool IsOneOf(const STELLAR_TYPE p_StellarType, const std::initializer_list<ST> p_List)
     *
     * @param   [IN]    p_StellarType               Stellar type to check
     * @param   [IN]    p_List                      List of stellar types
     * @return                                      Boolean - true if p_StellarType is in list, false if not
     */
    bool IsOneOf(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE_LIST p_List) {
        for (auto elem: p_List)
            if (p_StellarType == elem) return true;
        return false;
    }


    /*
     * Returns a single SN type based on the SN_EVENT parameter passed
     * 
     * Returns (in priority order):
     *
     *    SN_EVENT::CCSN  iff CCSN  bit is set and USSN bit is not set
     *    SN_EVENT::ECSN  iff ECSN  bit is set
     *    SN_EVENT::PISN  iff PISN  bit is set
     *    SN_EVENT::PPISN iff PPISN bit is set
     *    SN_EVENT::USSN  iff USSN  bit is set
     *    SN_EVENT::NONE  otherwise
     * 
     *
     * @param   [IN]    p_SNEvent                   SN_EVENT mask to check for SN event type
     * @return                                      SN_EVENT
     */
    SN_EVENT SNEventType(const SN_EVENT p_SNEvent) {

        if ((p_SNEvent & (SN_EVENT::CCSN | SN_EVENT::USSN)) == SN_EVENT::CCSN ) return SN_EVENT::CCSN;
        if ((p_SNEvent & SN_EVENT::ECSN )                   == SN_EVENT::ECSN ) return SN_EVENT::ECSN;
        if ((p_SNEvent & SN_EVENT::PISN )                   == SN_EVENT::PISN ) return SN_EVENT::PISN;
        if ((p_SNEvent & SN_EVENT::PPISN)                   == SN_EVENT::PPISN) return SN_EVENT::PPISN;
        if ((p_SNEvent & SN_EVENT::USSN )                   == SN_EVENT::USSN ) return SN_EVENT::USSN;
        
        return SN_EVENT::NONE;
    }


    /*
     * Solve Kepler's Equation using root finding techniques. Here we use Newton-Raphson.
     *
     * For a definition of all the anomalies see here:
     *
     *    https://en.wikipedia.org/wiki/Mean_anomaly
     *    https://en.wikipedia.org/wiki/True_anomaly
     *    https://en.wikipedia.org/wiki/Eccentric_anomaly
     *
     *
     * std::tuple<ERROR, double, double> SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity)
     *
     * @param   [IN]    p_MeanAnomaly               The mean anomaly
     * @param   [IN]    p_Eccentricity              Eccentricity of the binary
     * @return                                      Tuple containing (in order): error value, eccentric anomaly, true anomaly
     *                                              The error value returned will be:
     *                                                  ERROR::NONE if no error occurred
     *                                                  ERROR::NO_CONVERGENCE if the Newton-Raphson iteration did not converge
     *                                                  ERROR::OUT_OF_BOUNDS if the eccentric anomaly returned is < 0 or > 2pi
     *                                              If the error returned is not ERROR:NONE, use the eccentric anomaly and 
     *                                              true anomaly returned at your own risk
     */
    std::tuple<ERROR, double, double> SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity) {

        ERROR  error = ERROR::NONE;                                                                                                     // error

        double e = p_Eccentricity;
        double M = p_MeanAnomaly;
        double E = p_MeanAnomaly;                                                                                                       // initial guess at E is M - correct for e = 0

        double kepler = E - (e * sin(E)) - M;                                                                                           // let f(E) = 0.  Equation (92) in my "A simple toy model" document

        int iteration = 0;
        while (std::abs(kepler) >= NEWTON_RAPHSON_EPSILON && iteration++ < MAX_KEPLER_ITERATIONS) {                                     // repeat the approximation until E is within the specified error of the true value, or max iterations exceeded
            double keplerPrime = 1.0 - (e * cos(E));                                                                                    // derivative of f(E), f'(E).  Equation (94) in my "A simple toy model" document
            E = E - kepler / keplerPrime;
            kepler = E - (e * sin(E)) - M;                                                                                              // let f(E) = 0.  Equation (92) in my "A simple toy model" document
        }

        if (iteration >= MAX_KEPLER_ITERATIONS) error = ERROR::NO_CONVERGENCE;                                                          // no convergence - set error

        double nu = 2.0 * atan((sqrt((1.0 + e) / (1.0 - e))) * tan(0.5*E));                                                             // convert eccentric anomaly into true anomaly.  Equation (96) in my "A simple toy model" document

             if (utils::Compare(E, M_PI) >= 0 && utils::Compare(E, _2_PI) <= 0) nu += _2_PI;                                            // add 2PI if necessary
        else if (utils::Compare(E, 0.0)  <  0 || utils::Compare(E, _2_PI) >  0) error = ERROR::OUT_OF_BOUNDS;                           // out of bounds - set error

        return std::make_tuple(error, E, nu);
    }

    /*
     * Custom Caller/Callee reporting for repeated calls to "utils::POW(a,b)" with the same arguments
     */
    bool isprofiling=false;
    double _uamin=__DBL_MAX__;
    double _uamax=-__DBL_MAX__;
    double _ubmin=__DBL_MAX__;
    double _ubmax=-__DBL_MAX__;
    int ucount=0;
    std::unordered_map<std::string, int> umap;
    
    void setProfiling(bool yesno) {
        isprofiling = yesno;
    }
    
    bool checkDuplicatesLessThan10(std::string &thisone,int count) {
        bool b=false;
        if (count>=5) {
            for (int i=0;i<count;i++) {
                b = umap[thisone]==count;
            }
        }
        return b;
    }
    
    
    bool checkDuplicates10_100(std::string &thisone,int count) {
        bool b=false;
        for (int i=0;i<10;i++) {
            b = umap[thisone]==count;
        }
        return b;
    }
    
    
    bool checkDuplicates100_500(std::string &thisone,int count) {
        bool b=false;
        for (int i=0;i<20;i++) {
            b = umap[thisone]==count;
        }
        return b;
    }
    
    
    bool checkDuplicates500_1000(std::string &thisone,int count) {
        bool b=false;
        for (int i=0;i<100;i++) {
            b = umap[thisone]==count;
        }
        return b;
    }
    
    
    bool checkDuplicates1000_5000(std::string &thisone,int count) {
        bool b=false;
        for (int i=0;i<200;i++) {
            b = umap[thisone]==count;
        }
        return b;
    }
    
    
    bool checkDuplicates5000_10000(std::string &thisone,int count) {
        bool b=false;
        for (int i=0;i<1000;i++) {
            b = umap[thisone]==count;
        }
        return b;
    }
    
    
    bool checkDuplicatesMoreThan10000(std::string &thisone,int count) {
        bool b=false;
        for (int i=0;i<2000;i++) {
            b = umap[thisone]==count;
        }
        return b;
    }
    
    
    bool checkDuplicateCount(std::string &thisone,int count) {        
        if (count>1 && count<10) return checkDuplicatesLessThan10(thisone, count);
        else if (count <100) return checkDuplicates10_100(thisone, count);
        else if (count <500) return checkDuplicates100_500(thisone, count);
        else if (count <1000) return checkDuplicates500_1000(thisone, count);
        else if (count <5000) return checkDuplicates1000_5000(thisone, count);
        else if (count <10000) return checkDuplicates5000_10000(thisone, count);
        else return checkDuplicatesMoreThan10000(thisone, count);
        return false;
    }
    
    
    double POW(double a, double b) {
	    if (isprofiling) {
            _uamin = std::min(a,_uamin);
            _uamax = std::max(a,_uamax);
            _ubmin = std::min(b,_ubmin);
            _ubmax = std::max(b,_ubmax);
            ucount++;
        
            std::string thisone= "POW"+std::to_string(a)+","+std::to_string(b);
            auto index=umap.find(thisone);
            if (index==umap.end()) {
                umap[thisone] = 1;
            }
            else {
                int count = umap[thisone]++;
                checkDuplicateCount(thisone, count);
            }
	    }
        return pow(a,b);
    }

    
    void finalisePOW() {
        std::cout << "===" <<ucount<<std::endl;
        std::cout << "===amin " <<_uamin<<std::endl;
        std::cout << "===amax " <<_uamax<<std::endl;
        std::cout << "===bmin " <<_ubmin<<std::endl;
        std::cout << "===bmax " <<_ubmax<<std::endl;
        std::unordered_map<std::string, int>::iterator it = umap.begin();
        while (it != umap.end())
        {
            std::string word = it->first;                    // accessing KEY from element pointed by it.
            int count = it->second;                          // accessing VALUE from element pointed by it.
            std::cout << word << "," << count << std::endl;
            it++;                                            // increment the Iterator to point to next entry
        }
    }
}
