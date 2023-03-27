#include <iostream>
#include <stdarg.h>
#include <algorithm>
#include <cstring>
#include "profiling.h"
#include "utils.h"
#include "Rand.h"
#include "changelog.h"

/*
 * utility functions that don't belong in any class
 *
 */
namespace utils {

    // Alphabetical - so I can find them...

     /*
     * Iterative binary search
     *
     * For a given number x and a sorted array arr, return the lower and upper bin edges of x in arr.
     *
     *
     * std::vector<int> binarySearch(const std::vector<double> p_Arr, const double p_x)
     *
     * @param   [IN]    p_Array             Sorted array to search over
     * @param   [IN]    p_x                 Value to search for

     * @return                              Vector containing indices of the lower and upper
                                            bin edges containing x. If x < min(Arr), return
                                            {-1, 0}. If x > max(Arr), return {0, -1}. If x
                                            is equal to an array element, return index of that
                                            element.
     */
    std::vector<int> binarySearch(const std::vector<double> p_Arr, const double p_x) {
        int low = 0;
        int up = p_Arr.size() - 1;
        int mid = 0;

        // If x is not within array limits...
        if      (p_x < p_Arr[low]) { return {-1, 0}; }
        else if (p_x > p_Arr[up])  { return {0, -1}; }

        while(1) {
            mid = roundl( 0.5*(up + low) );
            if (std::abs(low - up) == 1) { return {low, low+1}; }    // arr(low) < x < arr(up), so return low
            else if (p_x == p_Arr[low])  { return {low, low}; }      // arr(low) = x. In this case, return low = up
            else if (p_x == p_Arr[up])   { return {up, up}; }        // arr(up) = x. In this case, return low = up
            else if (p_x == p_Arr[mid])  { return {mid, mid}; }      // arr(mid) = x. In this case, return low = up = mid
            else if (p_x < p_Arr[mid])   { up = mid; }               // Bring down upper bound
            else                         { low = mid; }              // Bring up lower bound
        }
    }


    /*
     * Calculate the value of the CDF of the Kroupa (2001) IMF at p_Mass
     *
     * If p_Mass is outside the bounds of the IMF (< p_Min or >= p_Max), the returned CDF value will be 0.0
     * 
     * 
     * double CalculateCDFKroupa(const double p_Mass, const double p_Max, const double p_Min)
     *
     * @param   [IN]    p_Mass                      Mass value (in Msol) at which to calculate the CDF
     * @param   [IN]    p_Max                       IMF maximum
     * @param   [IN]    p_Min                       IMF minimum
     * @return                                      CDF value
     */
    double CalculateCDFKroupa(const double p_Mass, const double p_Max, const double p_Min) {

        if ((p_Mass < p_Min) || (p_Mass >= p_Max)) return 0.0;      // return 0.0 if mass is out of bounds of function
    
        double CDF = 0.0;

        if (p_Min <= KROUPA_BREAK_1 &&
            p_Max >  KROUPA_BREAK_1 &&
            p_Max <= KROUPA_BREAK_2) {

            double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
            double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (PPOW(p_Max, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);

            double C1 = 1.0 / (term1 + term2);
            double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;

            if (p_Mass >= p_Min && p_Mass < KROUPA_BREAK_1) {
                CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_1) - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
            }
            else if (p_Mass >= KROUPA_BREAK_1 && p_Mass < KROUPA_BREAK_2) {
                CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                      ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
            }
        }
        else if (p_Min <= KROUPA_BREAK_1 &&
                 p_Max >  KROUPA_BREAK_2) {

            double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
            double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
            double term3 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

            double C1 = 1.0 / (term1 + term2 + term3);
            double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;
            double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

            if (p_Mass >= p_Min && p_Mass < KROUPA_BREAK_1) {
                CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_1) - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
            }
            else if (p_Mass >= KROUPA_BREAK_1 && p_Mass < KROUPA_BREAK_2) {
                CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                      ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
            }
            else if (p_Mass >= KROUPA_BREAK_2 && p_Mass < p_Max) {
                CDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                      ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2) +
                      ONE_OVER_KROUPA_POWER_3_PLUS1 * C3 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
            }
        }
        else if (p_Min >  KROUPA_BREAK_1 &&
                 p_Min <= KROUPA_BREAK_2 &&
                 p_Max >  KROUPA_BREAK_2) {

            double term1 = ONE_OVER_KROUPA_POWER_2_PLUS1 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
            double term2 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

            double C2 = 1.0 / (term1 + term2);
            double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

            if (p_Mass >= p_Min && p_Mass < KROUPA_BREAK_2) {
                CDF = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
            }
            else if (p_Mass >= KROUPA_BREAK_2 && p_Mass < p_Max) {
                CDF = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2)) +
                      ONE_OVER_KROUPA_POWER_3_PLUS1 * C3 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
            }
        }

        return CDF;
    }


    /*
     * Centre-justifies string to specified width by prepending and appending spaces.
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
        double a_SI              = std::cbrt(a_cubed_SI); 

        return a_SI / AU;
    }


    /*
     * Draw the angular components of the supernova kick theta and phi.
     *
     * 
     * DBL_DBL DrawKickDirection(const KICK_DIRECTION_DISTRIBUTION p_KickDirectionDistribution, const double p_KickDirectionPower)
     * 
     * @param   [IN]    p_KickDirectionDistribution The kick direction distribution to use - program option
     * @param   [IN]    p_KickDirectionPower        Exponent for power law - program option
     * @return                                      Tuple containing theta and phi 
     */
    DBL_DBL DrawKickDirection(const KICK_DIRECTION_DISTRIBUTION p_KickDirectionDistribution, const double p_KickDirectionPower) {

        double theta = 0.0;                                                                                         // theta, angle out of the plane
        double phi   = 0.0;                                                                                         // phi, angle in the plane

        double delta = 1.0 * DEGREE;                                                                                // small angle () in radians - could be set by user in options

        double rand = RAND->Random();                                                                               // do this here to be consistent with legacy code - allows comparison tests (won't work for long - soon there will be too many changes to the code...)

        switch (p_KickDirectionDistribution) {                                                                      // which kick direction distribution?

            case KICK_DIRECTION_DISTRIBUTION::ISOTROPIC:                                                            // ISOTROPIC: Draw theta and phi isotropically
                theta = acos(1.0 - (2.0 * RAND->Random())) - M_PI_2;
                phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
                break;

            case KICK_DIRECTION_DISTRIBUTION::POWERLAW: {                                                           // POWERLAW: Draw phi uniform in [0,2pi], theta according to a powerlaw
                                                                                                                    // (power law power = 0 = isotropic, +infinity = kick along pole, -infinity = kick in plane)
                // Choose magnitude of power law distribution -- if using a negative power law that blows up at 0,
                // need a lower cutoff (currently set at 1E-6), check it doesn't affect things too much
                // JR: todo: should these be in constants.h?
                double magnitude_of_cos_theta = utils::InverseSampleFromPowerLaw(p_KickDirectionPower, 1.0, 1E-6);
                if (p_KickDirectionPower < 0.0) magnitude_of_cos_theta = 1.0 - magnitude_of_cos_theta;              // don't use utils::Compare() here

                double actual_cos_theta = magnitude_of_cos_theta;

                if (RAND->Random() < 0.5) actual_cos_theta = -magnitude_of_cos_theta;                               // By taking the magnitude of cos theta we lost the information about whether it was up or down, so we put that back in randomly here

                actual_cos_theta = std::min(1.0, std::max(-1.0, actual_cos_theta));                                 // clamp to [-1.0, 1.0]

                theta = acos(actual_cos_theta);                                                                     // set the kick angle out of the plane theta
                phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
                } break;

            case KICK_DIRECTION_DISTRIBUTION::INPLANE:                                                              // INPLANE: Force the kick to be in the plane theta = 0
                theta = 0.0;                                                                                        // force the kick to be in the plane - theta = 0
                phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
                break;

            case KICK_DIRECTION_DISTRIBUTION::PERPENDICULAR:                                                        // PERPENDICULAR: Force kick to be along spin axis
                theta = M_PI_2;                                                                                     // pi/2 UP
                if (rand >= 0.5) theta = -theta;                                                                    // switch to DOWN
                phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
                break;

            case KICK_DIRECTION_DISTRIBUTION::POLES:                                                                // POLES: Direct the kick in a small cone around the poles
                if (rand < 0.5) theta = M_PI_2 - fabs(RAND->RandomGaussian(delta));                                 // UP - slightly less than or equal to pi/2
                else            theta = fabs(RAND->RandomGaussian(delta)) - M_PI_2;                                 // DOWN - slightly more than or equal to -pi/2

                phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
                break;

            case KICK_DIRECTION_DISTRIBUTION::WEDGE:                                                                // WEDGE: Direct kick into a wedge around the equator (theta = 0)
                theta = RAND->RandomGaussian(delta);                                                                // Gaussian around 0 with a deviation delta
                phi   = RAND->Random() * _2_PI;                                                                     // allow to randomly take an angle 0 - 2pi in the plane
                break;

            default:                                                                                                // unknown kick direction distribution - use ISOTROPIC
                // the kick direction distribution is set by a commandline option that
                // has already been checked by the time we get to this function - the
                // default case should never be taken.
                theta = acos(1.0 - (2.0 * RAND->Random())) - M_PI_2;
                phi   = RAND->Random() * _2_PI;
        }

        return std::make_tuple(theta, phi);
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
     * Calculate x^y where x is double and y is an integer
     *
     * faster than pow() for integer exponent
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
            double min_powerPlus1 = PPOW(p_Xmin, powerPlus1);

            result = PPOW((rand * (PPOW(p_Xmax, powerPlus1) - min_powerPlus1) + min_powerPlus1), 1.0 / powerPlus1);
        }

        return result;
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

        if (p_Table.size() > 0 && (p_Y >= 0.0 && p_Y < 1.0)) {                                  // sanity check - map must not be empty, and p_Y must in [0.0, 1.0)      (Leave these as absolute compares)

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
     * Determines if the string passed as p_Str is a valid BOOL
     * (as defined by Boost)
     *
     * In this context (the Boost context), a valid boolean is one of:
     * 
     *     - 0|1        ("0" or "1")
     *     - true|false ("true" or "false" - case insensitive)
     *     - yes|no     ("yes" or "no" - case insensitive)
     *     - on|off     ("on" or "off" - case insensitive)
     *
     * The function will retiurn one of {0, 1, 2, 3, 4, -1, -2, -3, -4} to indicate the result:
     * 
     *     0 = not a valid boolean
     *     1 = valid: 0|1
     *     2 = valid: true|false
     *     3 = valid: yes|no
     *     4 = valid: on|off
     * 
     *     A positive return value indicates the boolean value is TRUE; a negative, FALSE
     * 
     * 
     * int IsBOOL(const std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to check
     * @return                                      Result - as described above
     */
    int IsBOOL(const std::string p_Str) {

        if (p_Str.empty()) return 0;                    // not valid: empty string

        if (utils::Equals(p_Str, "0")    ) return  1;   // valid: 0|1       : TRUE
        if (utils::Equals(p_Str, "1")    ) return -1;   // valid: 0|1       : FALSE
        if (utils::Equals(p_Str, "true") ) return  2;   // valid: true|false: TRUE
        if (utils::Equals(p_Str, "false")) return -2;   // valid: true|false: FALSE
        if (utils::Equals(p_Str, "yes")  ) return  3;   // valid: yes|no    : TRUE
        if (utils::Equals(p_Str, "no")   ) return -3;   // valid: yes|no    : FALSE
        if (utils::Equals(p_Str, "on")   ) return  4;   // valid: on|off    : TRUE
        if (utils::Equals(p_Str, "off")  ) return -4;   // valid: on|off    : FALSE

        return 0;                                       // not valid
    }


    /*
     * Determines if the string passed as p_Str is a valid DOUBLE
     *
     * In this context, to be a valid DOUBLE the string must convert to a
     * double succussfully via the std::stod() function
     * 
     * 
     * int IsDOUBLE(const std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to check
     * @return                                      Result - TRUE if string is a valid DOUBLE, else FALSE
     */
    bool IsDOUBLE(const std::string p_Str) {

        bool result = false;                        // default result

        try {
            size_t lastChar;
            (void)std::stod(p_Str, &lastChar);      // try conversion

            result = lastChar == p_Str.size();      // valid DOUBLE if p_Str completely consumed
        }
        catch (const std::out_of_range& e) {        // conversion failed
            result = false;                         // not a valid DOUBLE
        }
        catch (const std::invalid_argument& e) {    // conversion failed
            result = false;                         // not a valid DOUBLE
        }
        return result;
    }


    /*
     * Determines if the string passed as p_Str is a valid FLOAT
     *
     * In this context, to be a valid FLOAT the string must convert to a
     * double succussfully via the std::stof() function
     * 
     * 
     * int IsFLOAT(const std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to check
     * @return                                      Result - TRUE if string is a valid FLOAT, else FALSE
     */
    bool IsFLOAT(const std::string p_Str) {

        bool result = false;                        // default result

        try {
            size_t lastChar;
            (void)std::stof(p_Str, &lastChar);      // try conversion

            result = lastChar == p_Str.size();      // valid FLOAT if p_Str completely consumed
        }
        catch (const std::out_of_range& e) {        // conversion failed
            result = false;                         // not a valid FLOAT
        }
        catch (const std::invalid_argument& e) {    // conversion failed
            result = false;                         // not a valid FLOAT
        }
        return result;
    }


    /*
     * Determines if the string passed as p_Str is a valid INT
     *
     * In this context, to be a valid INT the string must convert to an
     * integer succussfully via the std::stoi() function
     * 
     * 
     * int IsINT(const std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to check
     * @return                                      Result - TRUE if string is a valid INT, else FALSE
     */
    bool IsINT(const std::string p_Str) {

        bool result = false;                        // default result

        try {
            size_t lastChar;
            (void)std::stoi(p_Str, &lastChar);      // try conversion

            result = lastChar == p_Str.size();      // valid INT if p_Str completely consumed
        }
        catch (const std::out_of_range& e) {        // conversion failed
            result = false;                         // not a valid INT
        }
        catch (const std::invalid_argument& e) {    // conversion failed
            result = false;                         // not a valid INT
        }
        return result;
    }


    /*
     * Determines if the string passed as p_Str is a valid LONG DOUBLE
     *
     * In this context, to be a valid LONG DOUBLE the string must convert to a
     * long double succussfully via the std::stold() function
     * 
     * 
     * int IsLONGDOUBLE(const std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to check
     * @return                                      Result - TRUE if string is a valid LONG DOUBLE, else FALSE
     */
    bool IsLONGDOUBLE(const std::string p_Str) {

        bool result = false;                        // default result

        try {
            size_t lastChar;
            (void)std::stold(p_Str, &lastChar);     // try conversion

            result = lastChar == p_Str.size();      // valid LONG DOUBLE if p_Str completely consumed
        }
        catch (const std::out_of_range& e) {        // conversion failed
            result = false;                         // not a valid LONG DOUBLE
        }
        catch (const std::invalid_argument& e) {    // conversion failed
            result = false;                         // not a valid LONG DOUBLE
        }
        return result;
    }


    /*
     * Determines if the string passed as p_Str is a valid LONG INT
     *
     * In this context, to be a valid LONG INT the string must convert to a
     * long integer succussfully via the std::stol() function
     * 
     * 
     * int IsLONGINT(const std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to check
     * @return                                      Result - TRUE if string is a valid LONGINT, else FALSE
     */
    bool IsLONGINT(const std::string p_Str) {

        bool result = false;                        // default result

        try {
            size_t lastChar;
            (void)std::stol(p_Str, &lastChar);      // try conversion

            result = lastChar == p_Str.size();      // valid LONG INT if p_Str completely consumed
        }
        catch (const std::out_of_range& e) {        // conversion failed
            result = false;                         // not a valid LONG INT
        }
        catch (const std::invalid_argument& e) {    // conversion failed
            result = false;                         // not a valid LONG INT
        }
        return result;
    }


    /*
     * Determines if the string passed as p_Str is a valid UNSIGNED LONG INT
     *
     * In this context, to be a valid UNSIGNED LONG INT the string must convert to a
     * unsigned long integer succussfully via the std::stoul() function
     * 
     * 
     * int IsULONGINT(const std::string p_Str)
     *
     * @param   [IN]    p_Str                       String to check
     * @return                                      Result - TRUE if string is a valid UNSIGNED LONGINT, else FALSE
     */
    bool IsULONGINT(const std::string p_Str) {

        bool result = false;                        // default result

        try {
            size_t lastChar;
            (void)std::stoul(p_Str, &lastChar);     // try conversion

            result = lastChar == p_Str.size();      // valid UNSIGNED LONG INT if p_Str completely consumed
        }
        catch (const std::out_of_range& e) {        // conversion failed
            result = false;                         // not a valid UNSIGNED LONG INT
        }
        catch (const std::invalid_argument& e) {    // conversion failed
            result = false;                         // not a valid UNSIGNED LONG INT
        }
        return result;
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
     * Pads string to specified length by prepending the string with "0"
     *
     * This only works with ASCII data, but I think that's all we need
     *
     *
     * std::string PadLeadingZeros(const std::string p_Str, const std::size_t p_MaxLength)
     *
     * @param   [IN]    p_Str                       String to be padded with leading "0"s
     * @param   [IN]    p_MaxLength                 The required length of the resultant string
     * @return                                      String padded with leading "0"s - will be unchanged from input string if length alread >= required length
     */
    std::string PadLeadingZeros(const std::string p_Str, const std::size_t p_MaxLength) {
        return (p_Str.length() < p_MaxLength) ? std::string(p_MaxLength - p_Str.length(), '0') + p_Str : p_Str;
    }


    /*
     * Pads string to specified length by appending the string with " "
     *
     * This only works with ASCII data, but I think that's all we need
     *
     *
     * std::string PadTrailingSpaces(const std::string p_Str, const std::size_t p_MaxLength)
     *
     * @param   [IN]    p_Str                       String to be padded with trailing " "s
     * @param   [IN]    p_MaxLength                 The required length of the resultant string
     * @return                                      String padded with leading "0"s - will be unchanged from input string if length alread >= required length
     */
    std::string PadTrailingSpaces(const std::string p_Str, const std::size_t p_MaxLength) {
        return (p_Str.length() < p_MaxLength) ? p_Str + std::string(p_MaxLength - p_Str.length(), ' ') : p_Str;
    }


    /*
     * Trim leading whitespace characters from a string.
     *
     *
     * std::string& ltrim(std::string& p_Str)
     *
     * @param   [IN]    p_Str                       String to be trimmed of whitespace
     * @return                                      Trimmed string
     */
    std::string& ltrim(std::string& p_Str) {
        p_Str.erase(0, p_Str.find_first_not_of("\t\n\v\f\r "));
        return p_Str;
    }


    /*
     * Trim trailing whitespace characters from a string.
     *
     *
     * std::string& rtrim(std::string& p_Str)
     *
     * @param   [IN]    p_Str                       String to be trimmed of whitespace
     * @return                                      Trimmed string
     */
    std::string& rtrim(std::string& p_Str) {
        p_Str.erase(p_Str.find_last_not_of("\t\n\v\f\r ") + 1);
        return p_Str;
    }


    /*
     * Trim both leading and trailing whitespace characters from a string.
     *
     *
     * std::string& trim(std::string& p_Str)
     *
     * @param   [IN]    p_Str                       String to be trimmed of whitespace
     * @return                                      Trimmed string
     */
    std::string& trim(std::string& p_Str) {
        return utils::ltrim(utils::rtrim(p_Str));
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
     * Draw eccentricity from the distribution specified by the user
     *
     *
     * double SampleEccentricity()
     *
     * @param   [IN]    p_Edist                     The eccentricity distribution to use to draw the eccentricity
     * @param   [IN]    p_Max                       Distribution maximum
     * @param   [IN]    p_Min                       Distribution minimum
     * @return                                      Eccentricity
     */
    double SampleEccentricity(const ECCENTRICITY_DISTRIBUTION p_Edist, const double p_Max, const double p_Min) {

        double eccentricity;

        switch (p_Edist) {                                                                              // which distribution?

            case ECCENTRICITY_DISTRIBUTION::ZERO:                                                       // ZERO - all systems are initially circular i.e. have zero eccentricity
                eccentricity = 0.0;
                break;

            case ECCENTRICITY_DISTRIBUTION::FLAT:                                                       // FLAT
                eccentricity = utils::InverseSampleFromPowerLaw(0.0, p_Max, p_Min);
                break;

            case ECCENTRICITY_DISTRIBUTION::THERMAL:                                                    // THERMAL eccentricity distribution p(e) = 2e
                eccentricity = utils::InverseSampleFromPowerLaw(1.0, p_Max, p_Min);
                break;

            case ECCENTRICITY_DISTRIBUTION::GELLER_2013:                                                // M35 eccentricity distribution from Geller, Hurley and Mathieu 2013
                // Gaussian with mean 0.38 and sigma 0.23
                // http://iopscience.iop.org/article/10.1088/0004-6256/145/1/8/pdf
                // Sampling function taken from binpop.f in NBODY6

                do {
                    eccentricity = 0.23 * std::sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 0.38;
                } while (eccentricity < p_Min || eccentricity > p_Max);                                 // JR: don't use utils::Compare() here
                break;

            case ECCENTRICITY_DISTRIBUTION::DUQUENNOYMAYOR1991:                                        // eccentricity distribution from Duquennoy & Mayor (1991)
                // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
                // Sampling function taken from binpop.f in NBODY6

                do {
                    eccentricity = 0.15 * std::sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 0.3;
                } while (eccentricity < p_Min or eccentricity > p_Max);                                 // JR: don't use utils::Compare() here
                break;

            case ECCENTRICITY_DISTRIBUTION::SANA2012:                                                   // Sana et al 2012
                // (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
                // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
                // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

                eccentricity = utils::InverseSampleFromPowerLaw(-0.42, p_Max, p_Min);
                break;

            default:                                                                                    // unknown distribution - should not be possible (options code should prevent this)
                eccentricity = 0.0;
        }

        return eccentricity;
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
     * Draw mass from the distribution specified by the user
     *
     *
     * double SampleInitialMass(const INITIAL_MASS_FUNCTION p_IMF, const double p_Max, const double p_Min, const double p_Power)
     *
     * @param   [IN]    p_IMF                       The IMF to use to draw the mass
     * @param   [IN]    p_Max                       IMF maximum
     * @param   [IN]    p_Min                       IMF minimum
     * @param   [IN]    p_Power                     IMF power (for IMF::POWERLAW)
     * @return                                      Mass
     */
    double SampleInitialMass(const INITIAL_MASS_FUNCTION p_IMF, const double p_Max, const double p_Min, const double p_Power) {

        double thisMass = 0.0;

        switch (p_IMF) {                                                                                            // which IMF?

            case INITIAL_MASS_FUNCTION::SALPETER:                                                                   // SALPETER

                thisMass = utils::InverseSampleFromPowerLaw(SALPETER_POWER, p_Max, p_Min);
                break;

            case INITIAL_MASS_FUNCTION::POWERLAW:                                                                   // POWER LAW

                thisMass = utils::InverseSampleFromPowerLaw(p_Power, p_Max, p_Min);
                break;

            case INITIAL_MASS_FUNCTION::UNIFORM:                                                                    // UNIFORM - convienience function for POWERLAW with slope of 0

                thisMass = RAND->Random(p_Min, p_Max);
                break;

            case INITIAL_MASS_FUNCTION::KROUPA:                                                                     // KROUPA

                // find out where the user specificed their minimum and maximum masses to generate
                if (utils::Compare(p_Min, KROUPA_BREAK_1) <= 0 && utils::Compare(p_Max, KROUPA_BREAK_1) <= 0) {
                    thisMass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_1, p_Max, p_Min);                      // draw mass using inverse sampling
                }
                else if (utils::Compare(p_Min, KROUPA_BREAK_1) > 0 && utils::Compare(p_Min, KROUPA_BREAK_2) <= 0 &&
                         utils::Compare(p_Max, KROUPA_BREAK_1) > 0 && utils::Compare(p_Max, KROUPA_BREAK_2) <= 0) {

                    thisMass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_2, p_Max, p_Min);                      // draw mass using inverse sampling
                }
                else if (utils::Compare(p_Min, KROUPA_BREAK_2) > 0 && utils::Compare(p_Max, KROUPA_BREAK_2) > 0) {

                    thisMass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_3, p_Max, p_Min);                      // draw mass using inverse sampling
                }
                else if (utils::Compare(p_Min, KROUPA_BREAK_1) <= 0 &&
                         utils::Compare(p_Max, KROUPA_BREAK_1)  > 0 && utils::Compare(p_Max, KROUPA_BREAK_2) <= 0) {

                    double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
                    double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (PPOW(p_Max, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);

                    double C1    = 1.0 / (term1 + term2);
                    double C2    = C1 * KROUPA_BREAK_1_POWER_1_2;
                    double A     = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));

                    double rand  = RAND->Random();                                                                  // draw a random number between 0 and 1
                    thisMass = utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_1, p_Max, p_Min)) < 0
                                ? PPOW(rand * (KROUPA_POWER_PLUS1_1 / C1) + PPOW(p_Min, KROUPA_POWER_PLUS1_1), ONE_OVER_KROUPA_POWER_1_PLUS1)
                                : PPOW((rand - A) * (KROUPA_POWER_PLUS1_2 / C2) + KROUPA_BREAK_1_PLUS1_2, ONE_OVER_KROUPA_POWER_2_PLUS1);
                }
                else if (utils::Compare(p_Min, KROUPA_BREAK_1) <= 0 && utils::Compare(p_Max, KROUPA_BREAK_2_POWER_2_3) > 0) {

                    double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
                    double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
                    double term3 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
                    
                    double C1    = 1.0 / (term1 + term2 + term3);
                    double C2    = C1 * KROUPA_BREAK_1_POWER_1_2;
                    double C3    = C2 * KROUPA_BREAK_2_POWER_2_3;

                    double A     = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
                    double B     = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);

                    double rand  = RAND->Random();                                                                  // draw a random number between 0 and 1

                    if (utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_1, p_Max, p_Min)) < 0)
                        thisMass = PPOW(rand * (KROUPA_POWER_PLUS1_1 / C1) + PPOW(p_Min, KROUPA_POWER_PLUS1_1), ONE_OVER_KROUPA_POWER_1_PLUS1);
                    else if (utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_2, p_Max, p_Min)) < 0)
                        thisMass = PPOW((rand - A) * (KROUPA_POWER_PLUS1_2 / C2) + KROUPA_BREAK_1_PLUS1_2, ONE_OVER_KROUPA_POWER_2_PLUS1);
                    else
                        thisMass = PPOW((rand - A - B) * (KROUPA_POWER_PLUS1_3 / C3) + KROUPA_BREAK_2_PLUS1_3, ONE_OVER_KROUPA_POWER_3_PLUS1);
                }
                else if (utils::Compare(p_Min, KROUPA_BREAK_1)  > 0 &&
                         utils::Compare(p_Min, KROUPA_BREAK_2) <= 0 && utils::Compare(p_Max, KROUPA_BREAK_2) > 0) {

                    double term1 = ONE_OVER_KROUPA_POWER_2_PLUS1 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
                    double term2 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

                    double C2    = 1.0 / (term1 + term2);
                    double C3    = C2 * KROUPA_BREAK_2_POWER_2_3;
                    double B     = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2));

                    double rand  = RAND->Random();                                                                  // draw a random number between 0 and 1

                    thisMass = utils::Compare(rand, CalculateCDFKroupa(KROUPA_BREAK_2, p_Max, p_Min)) < 0
                                ? PPOW(rand * (KROUPA_POWER_PLUS1_2 / C2) + PPOW(p_Min, KROUPA_POWER_PLUS1_2), ONE_OVER_KROUPA_POWER_2_PLUS1)
                                : PPOW((rand - B) * (KROUPA_POWER_PLUS1_3 / C3) + KROUPA_BREAK_2_PLUS1_3, ONE_OVER_KROUPA_POWER_3_PLUS1);
                }
                // JR: no other case possible - as long as p_Min < p_Max (currently enforced in Options.cpp)
                break;

            default:                                                                                                // unknown IMF
                thisMass = utils::InverseSampleFromPowerLaw(SALPETER_POWER, SALPETER_MAXIMUM, SALPETER_MINIMUM);    // calculate mass using power law with default values
        }

        return thisMass;
    }


    /*
     * Draw mass ratio q from the distribution specified by the user
     *
     *
     * double SampleMassRatio(const MASS_RATIO_DISTRIBUTION p_Qdist, const double p_Max, const double p_Min)
     *
     * @param   [IN]    p_IMF                       The distribution to use to draw the ratio
     * @param   [IN]    p_Max                       Distribution maximum
     * @param   [IN]    p_Min                       Distribution minimum
     * @return                                      Mass ratio q
     */
    double SampleMassRatio(const MASS_RATIO_DISTRIBUTION p_Qdist, const double p_Max, const double p_Min) {

        double q;

        switch (p_Qdist) {

            case MASS_RATIO_DISTRIBUTION::FLAT:                                                                 // FLAT mass ratio distriution
                q = utils::InverseSampleFromPowerLaw(0.0, p_Max, p_Min);
                break;

            case MASS_RATIO_DISTRIBUTION::DUQUENNOYMAYOR1991:                                                   // mass ratio distribution from Duquennoy & Mayor (1991) (http://adsabs.harvard.edu/abs/1991A%26A...248..485D)

                do {                                                                                            // JR: todo: catch non-convergence?
                    q = 0.42 * std::sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 0.23;
                } while (q < p_Min || q > p_Max);                                                               // JR: don't use utils::Compare() here
                break;

            case MASS_RATIO_DISTRIBUTION::SANA2012:                                                             // Sana et al 2012 (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
                // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
                // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

                q = utils::InverseSampleFromPowerLaw(-0.1, p_Max, p_Min);                                       // de Mink and Belczynski use min = 0.1, max = 1.0
                break;

            default:                                                                                            // unknown q-distribution
                q = utils::InverseSampleFromPowerLaw(0.0, 1.0, 0.0);                                            // calculate q using power law with default values
        }

        return std::min(std::max(p_Min, q), p_Max);                                                             // clamp to [min, max]
    }


    /*
     * Draw eccentricity from the distribution specified by the user
     *
     *
     * double SampleMetallicity()
     *
     * @param   [IN]    p_Zdist                     The metallicity distribution to use to draw the metallicity
     * @param   [IN]    p_Max                       Distribution maximum
     * @param   [IN]    p_Min                       Distribution minimum
     * @return                                      Metallicity
     */
    double SampleMetallicity(const METALLICITY_DISTRIBUTION p_Zdist, const double p_Max, const double p_Min) {

        double metallicity;

        switch (p_Zdist) {                                                                              // which distribution?

            case METALLICITY_DISTRIBUTION::ZSOLAR:                                                      // ZSOLAR - all systems have Z = ZSOLAR
                metallicity = ZSOL_ASPLUND;
                break;

            case METALLICITY_DISTRIBUTION::LOGUNIFORM: {                                                // LOGUNIFORM - sample Z uniformly in the log
                double logMin = log10(p_Min);
                double logMax = log10(p_Max);
                metallicity = PPOW(10, (logMin + ((logMax - logMin) * RAND->Random())));
            } break;

            default:                                                                                    // unknown distribution - should not be possible (options code should prevent this)
                metallicity = 0.0;
        }

        return metallicity;
    }


    /*
     * Draw orbital period from the distribution specified by the user
     * 
     * 
     * double SampleOrbitalPeriodDistribution(const ORBITAL_PERIOD_DISTRIBUTION p_Pdist, 
     *                                        const double                      p_PdistMax, 
     *                                        const double                      p_PdistMin)
     *
     * @param   [IN]    p_Pdist                     The distribution to use to draw orbital period
     * @param   [IN]    p_PdistMax                  Orbital period distribution maximum
     * @param   [IN]    p_PdistMin                  Orbital period distribution minimum
     * @return                                      Orbital period in days
     */
    double SampleOrbitalPeriod(const ORBITAL_PERIOD_DISTRIBUTION p_Pdist, 
                               const double                      p_PdistMax, 
                               const double                      p_PdistMin) {

        double orbitalPeriod;

        switch (p_Pdist) {                                                                                              // which distribution?

            case ORBITAL_PERIOD_DISTRIBUTION::FLATINLOG:                                                                // FLAT IN LOG

                orbitalPeriod = utils::InverseSampleFromPowerLaw(-1.0, p_PdistMax, p_PdistMin);
                break;

            default:                                                                                                    // unknown distribution
                orbitalPeriod = utils::InverseSampleFromPowerLaw(-1.0, 1000.0, 1.1);                                    // calculate orbitalPeriod using power law with default values
        }

        return orbitalPeriod;
    }


    /*
     * Draw semi-major axis from the distribution specified by the user
     * 
     * 
     * double SampleSemiMajorAxisDistribution(const SEMI_MAJOR_AXIS_DISTRIBUTION p_Adist, 
     *                                        const double                       p_AdistMax, 
     *                                        const double                       p_AdistMin, 
     *                                        const double                       p_AdistPower, 
     *                                        const double                       p_PdistMax, 
     *                                        const double                       p_PdistMin, 
     *                                        const double                       p_Mass1, 
     *                                        const double                       p_Mass2)
     *
     * @param   [IN]    p_Adist                     The distribution to use to draw semi-major axis
     * @param   [IN]    p_AdistMax                  Semi-major axis distribution maximum
     * @param   [IN]    p_AdistMin                  Semi-major axis distribution minimum
     * @param   [IN]    p_AdistPower                Semi-major axis distribution power
     * @param   [IN]    p_PdistMax                  Period distribution maximum (for SANA2012 distribution)
     * @param   [IN]    p_PdistMin                  Period distribution minimum (for SANA2012 distribution)
     * @param   [IN]    p_Mass1                     Mass of the primary
     * @param   [IN]    p_Mass2                     Mass of the secondary
     * @return                                      Semi-major axis in AU
     */
    double SampleSemiMajorAxis(const SEMI_MAJOR_AXIS_DISTRIBUTION p_Adist, 
                               const double                       p_AdistMax, 
                               const double                       p_AdistMin, 
                               const double                       p_AdistPower, 
                               const double                       p_PdistMax, 
                               const double                       p_PdistMin, 
                               const double                       p_Mass1, 
                               const double                       p_Mass2) {

        double semiMajorAxis;

        switch (p_Adist) {                                                                                              // which distribution?

            case SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG:                                                               // FLAT IN LOG

                semiMajorAxis = utils::InverseSampleFromPowerLaw(-1.0, p_AdistMax, p_AdistMin);
                break;

            case SEMI_MAJOR_AXIS_DISTRIBUTION::DUQUENNOYMAYOR1991:                                                      // Duquennoy & Mayor (1991) period distribution
                // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
                // See also the period distribution (Figure 1) of M35 in Geller+ 2013 https://arxiv.org/abs/1210.1575
                // See also the period distribution (Figure 13) of local solar type binaries from Raghavan et al 2010 https://arxiv.org/abs/1007.0414
                // They have log-normal distribution with a mean of 5.03 and a standard deviation of 2.28, with a minimum period of around 0.1 days
                // Sampling function taken from binpop.f in NBODY6

                // Make sure that the drawn semi-major axis is in the range specified by the user
                do {                                                                                                    // JR: todo: catch for non-convergence?
                    double periodInDays = PPOW(10.0, 2.3 * std::sqrt(-2.0 * log(RAND->Random())) * cos(2.0 * M_PI * RAND->Random()) + 4.8);
                    semiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, periodInDays);      // convert period in days to semi-major axis in AU
                } while (semiMajorAxis < p_AdistMin || semiMajorAxis > p_AdistMax);                                     // JR: don't use utils::Compare() here
                break;

            case SEMI_MAJOR_AXIS_DISTRIBUTION::SANA2012: {                                                              // Sana et al 2012
                // http://science.sciencemag.org/content/sci/337/6093/444.full.pdf
                // distribution of semi-major axes. Sana et al fit for the orbital period, which we sample in here, before returning the semi major axis
                // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
                // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

                double logPeriodMin = p_PdistMin > 1.0 ? log(p_PdistMin) : 0.0;                                         // smallest initial log period  JR: don't use utils::Compare() here
                double logPeriodMax = p_PdistMax > 1.0 ? log(p_PdistMax) : 0.0;                                         // largest initial log period   JR: don't use utils::Compare() here

                double periodInDays = exp(utils::InverseSampleFromPowerLaw(-0.55, logPeriodMax, logPeriodMin));         // draw a period in days from their distribution

                semiMajorAxis = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, periodInDays);          // convert period in days to semi-major axis in AU
                } break;

            default:                                                                                                    // unknown distribution
                semiMajorAxis = utils::InverseSampleFromPowerLaw(-1.0, 100.0, 0.5);                                     // calculate semiMajorAxis using power law with default values
        }

        return semiMajorAxis;
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
     *    SN_EVENT::AIC   iff AIC   bit is set
     *    SN_EVENT::SNIA  iff SNIA  bit is set
     *    SN_EVENT::HeSD  iff HeSD  bit is set
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
        if ((p_SNEvent & SN_EVENT::AIC  )                   == SN_EVENT::AIC  ) return SN_EVENT::AIC;
        if ((p_SNEvent & SN_EVENT::SNIA )                   == SN_EVENT::SNIA ) return SN_EVENT::SNIA;
        if ((p_SNEvent & SN_EVENT::HeSD   )                 == SN_EVENT::HeSD ) return SN_EVENT::HeSD;

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

        double nu = 2.0 * atan((std::sqrt((1.0 + e) / (1.0 - e))) * tan(0.5*E));                                                             // convert eccentric anomaly into true anomaly.  Equation (96) in my "A simple toy model" document

             if (utils::Compare(E, M_PI) >= 0 && utils::Compare(E, _2_PI) <= 0) nu += _2_PI;                                            // add 2PI if necessary
        else if (utils::Compare(E, 0.0)  <  0 || utils::Compare(E, _2_PI) >  0) error = ERROR::OUT_OF_BOUNDS;                           // out of bounds - set error

        return std::make_tuple(error, E, nu);
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
            double sqrtD = std::sqrt(discriminant);
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
     * Announce COMPAS
     * 
     * Constructs and returns a splash string.  Prints string to stdout if required.
     *
     *
     * std::string SplashScreen(const bool p_Print)
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

}
