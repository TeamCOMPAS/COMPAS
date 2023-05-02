#ifndef __utils_h__
#define __utils_h__

#include "constants.h"
#include "typedefs.h"

namespace utils {


    // object identifiers - all classes have these (adding here (no class) for error handling)
    inline OBJECT_ID    ObjectId()    { return static_cast<int>(OBJECT_TYPE::UTILS); }     // object id for utils - ordinal value from enum
    inline OBJECT_TYPE  ObjectType()  { return OBJECT_TYPE::UTILS; }                       // object type for utils - always "UTILS"
    inline STELLAR_TYPE StellarType() { return STELLAR_TYPE::NONE; }                       // stellar type for utils - always "NONE"


    // namespace functions - alphabetical (sort of)

    std::vector<int>                    binarySearch(const std::vector<double> p_Arr, const double p_x);

    double                              CalculateCDFKroupa(const double p_Mass, const double p_Max, const double p_Min);

    std::string                         CentreJustify(const std::string p_Str, std::size_t p_Width);

    int                                 Compare(const double p_X, const double p_Y);

    double                              ConvertPeriodInDaysToSemiMajorAxisInAU(const double p_Mass1, const double p_Mass2, const double p_Period);

    DBL_DBL                             DrawKickDirection(const KICK_DIRECTION_DISTRIBUTION p_KickDirectionDistribution, const double p_KickDirectionPower);
    
    bool                                Equals(std::string p_Str1, std::string p_Str2);

    bool                                FileExists(const std::string& p_Filename);
    bool                                FileExists(const char *p_Filename);

    /*
     * Generic function to find an element in a vector
     *
     * Determines if an element is contained within a vector.
     * Returns a tuple containing:
     *
     *   - a boolean indicating the result (true = found, false = not found)
     *   - the vector index of the element if found (-1 if not found)
     *
     *
     * std::tuple<bool, int> Find(const T &p_Elem, const std::vector<T> &p_Vector)
     *
     * @param   [IN]    p_Elem                      The element to find in the vector
     * @param   [IN]    p_Vector                    The vector in which to look for p_Elem
     * @return                                      Tuple<bool, int> indicating <found, element index>, with index = -1 if found = false
     */
    template <typename T>
    static std::tuple<bool, int> Find(const T &p_Elem, const std::vector<T> &p_Vector) {

        auto iter = std::find(p_Vector.begin(), p_Vector.end(), p_Elem);                                                        // find the element

        return iter != p_Vector.end() ? std::make_tuple(true, distance(p_Vector.begin(), iter)) : std::make_tuple(false, -1l);  // if found return index, otherwise -1
    }


    /*
     * Find a value in an unordered map and return the key if found, otherwise default value
     *
     * This function looks for the passed string value in an unordered map, and if the string
     * is found returns the key corresponding to the value found.  If the value is not found
     * the value passed as the default value is returned.
     *
     * The string comparison is case-insensitive.
     *
     * The default value passed must allow the functional return type to be deduced - so it
     * should be one of the values in the unordered_map.
     *
     *
     * template<typename M, typename E>
     * std::tuple<bool, E> GetMapKey(const std::string p_Value, const M& p_Map, const E& p_Default)
     *
     * @param   [IN]    p_Value                     The value to be located in the unordered map
     * @param   [IN]    p_Map                       The unordered map in which to locate the value
     * @param   [IN]    p_Default                   The default value to be returned if p_value is not found
     * @return                                      Tuple containing a boolean result (true if value found, else false), and the key
     *                                              corresponding to the value found, or the default if the value was not found
     */
    template<typename M, typename E>
    std::tuple<bool, E> GetMapKey(const std::string p_Value, const M& p_Map, const E& p_Default) {
        for (auto& it: p_Map)
            if (Equals(it.second, p_Value)) return std::make_tuple(true, it.first);
        return std::make_tuple(false, p_Default);
    }

    double                              intPow(const double p_Base, const int p_Exponent);

    double                              InverseSampleFromPowerLaw(const double p_Power, const double p_Xmax, const double p_Xmin);
    double                              InverseSampleFromTabulatedCDF(const double p_Y, const std::map<double, double> p_Table);

    int                                 IsBOOL(const std::string p_Str);
    bool                                IsDOUBLE(const std::string p_Str);
    bool                                IsFLOAT(const std::string p_Str);
    bool                                IsINT(const std::string p_Str);
    bool                                IsLONGDOUBLE(const std::string p_Str);
    bool                                IsLONGINT(const std::string p_Str);
    bool                                IsULONGINT(const std::string p_Str);

    bool                                IsOneOf(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE_LIST p_List);


    std::string                         PadLeadingZeros(const std::string p_Str, const std::size_t p_MaxLength);
    std::string                         PadTrailingSpaces(const std::string p_Str, const std::size_t p_MaxLength);
    
    std::string&                        ltrim(std::string& p_Str);
    std::string&                        rtrim(std::string& p_Str);
    std::string&                        trim(std::string& p_Str);

    std::string                         ToLower(std::string p_Str);
    std::string                         ToUpper(std::string p_Str);

    const std::string                   vFormat(const char* const p_zcFormat, ...);


    double                              SampleEccentricity(const ECCENTRICITY_DISTRIBUTION p_Edist, const double p_Max, const double p_Min);
    double                              SampleFromTabulatedCDF(const double p_X, const std::map<double, double> pTable);
    double                              SampleInitialMass(const INITIAL_MASS_FUNCTION p_IMF, const double p_Max, const double p_Min, const double p_Power);
    double                              SampleMassRatio(const MASS_RATIO_DISTRIBUTION p_Qdist, const double p_Max, const double p_Min);
    double                              SampleMetallicity(const METALLICITY_DISTRIBUTION p_Zdist, const double p_Max, const double p_Min);
    double                              SampleOrbitalPeriod(const ORBITAL_PERIOD_DISTRIBUTION p_Pdist, const double p_PdistMax, const double p_PdistMin);
    double                              SampleSemiMajorAxis(const SEMI_MAJOR_AXIS_DISTRIBUTION p_Adist, 
                                                            const double                       p_AdistMax, 
                                                            const double                       p_AdistMin, 
                                                            const double                       p_AdistPower, 
                                                            const double                       p_PdistMax, 
                                                            const double                       p_PdistMin, 
                                                            const double                       p_Mass1, 
                                                            const double                       p_Mass2);

    SN_EVENT                            SNEventType(const SN_EVENT p_SNEvent);

    std::tuple<ERROR, double, double>   SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity);

    std::tuple<ERROR, double>           SolveQuadratic(const double p_A, const double p_B, double p_C);

    std::string                         SplashScreen(const bool p_Print = true);

}

#endif // __utils_h__
