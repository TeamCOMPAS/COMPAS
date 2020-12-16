#ifndef __profiling_h__
#define __profiling_h__

// configure the profiling functionality here by defining or undefining these variables

#define DOPROFILING                 // comment this line out, or #undef DOPROFILING, to build production executable (i.e. no PROFILING code)
#undef DOPROFILING

#ifdef DOPROFILING                  // profiling enabled?
                                    // yes - profiling enabled

#define PROFILING_COUNTS_ONLY       // comment this line out, or #undef PROFILING_COUNTS_ONLY, to include CPU penalties for calls to profiling functions (for use with external profiling tools)
#define PROFILING_CALLER_NAME       // comment this line out, or #undef PROFILING_CALLER_NAME, to NOT include the name of the calling function in the call signature of the function being profiled


#include "constants.h"


// Put any definitions of functions to be profiled here - and ensure that their "real",
// non-profiling definitions are in the "else" section below.
// Supporting functions for the functions being profiled go inside the namespace below.

#define InitialiseProfiling         profiling::Initialise()                                         // initialise profiling functionality
#define ReportProfiling             profiling::Report()                                             // report profiling outcomes

// define profiler calls to std::pow()
#ifdef PROFILING_CALLER_NAME                                                                        // add caller name to signature?
#define PPOW(base, exponent)        profiling::pow(base, exponent, __PRETTY_FUNCTION__)             // yes
#else
#define PPOW(base, exponent)        profiling::pow(base, exponent)                                  // no
#endif


namespace profiling {


    // object identifiers - all classes have these (adding here (no class) for error handling)
    inline OBJECT_ID    ObjectId()    { return static_cast<int>(OBJECT_TYPE::PROFILING); }          // object id for profiling - ordinal value from enum
    inline OBJECT_TYPE  ObjectType()  { return OBJECT_TYPE::PROFILING; }                            // object type for profiling - always "UTILS"
    inline STELLAR_TYPE StellarType() { return STELLAR_TYPE::NONE; }                                // stellar type for profiling - always "NONE"


    // namespace functions
    void   Initialise();                                                                            // should be called once immediately after program start
    void   Report();                                                                                // should be called once immediately prior to program termination


    // put any functions here that are to be reported on by the profiling code

    // std::pow() related functions
    void   InitialisePow();                                                                         // initialisation function for std::pow() profiling
    void   ReportPow();                                                                             // reporting function for std::pow() profiling

    bool   CheckDuplicate(const std::string &p_Signature, const int p_Count);                       // checks for duplicate calls to std::pow()

    // profiling function for std::pow()
    #ifdef PROFILING_CALLER_NAME                                                                    // add caller name to signature?
    double pow(const double p_Base, const double p_Exponent, const std::string &p_CallerName);      // yes
    #else
    double pow(const double p_Base, const double p_Exponent);                                       // no
    #endif

}

#else                               // no - profiling not enabled

#define InitialiseProfiling         {}                                                              // initialise profiling functionality
#define ReportProfiling             {}                                                              // report profiling outcomes

#define PPOW(base, exponent)        std::pow(base, exponent)                                        // not profiling calls to pow()

#endif // DOPROFILING

#endif // __profiling_h__