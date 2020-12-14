// This is the profiling namespace, and exists to provide profiling of function calls
// (e.g. pow(), sqrt() etc., but not limited to system functions - user defined functions can be
// profiled here too).
//
// See the "#define PROFILING" statement in profiling.h to learn how to turn profiling on and off


#include <iostream>
#include "profiling.h"

#ifdef DOPROFILING


namespace profiling {

    std::unordered_map<std::string, int> m_CallMap;                                             // map to record call signatures and counts


    // put any variables here that are to be reported by the profiling code
    // namespace variables need to be declared and initialised here

    double m_PowBaseMin     =  __DBL_MAX__;     // minimum value of 'base' in calls to std::pow(base, exponent)
    double m_PowBaseMax     = -__DBL_MAX__;     // maximum value of 'base' in calls to std::pow(base, exponent)
    double m_PowExponentMin =  __DBL_MAX__;     // minimum value of 'exponent' in calls to std::pow(base, exponent)
    double m_PowExponentMax = -__DBL_MAX__;     // maximum value of 'exponent' in calls to std::pow(base, exponent)

    int    m_PowCallCount   = 0;                // Total (accumulated) call count for std::pow()


    // member functions

    /*
     * Initialises profiling functionality
     * 
     * Put any code here required for initialisation of profiling functionality
     * This function should be called once immediately after program start
     * 
     * 
     * void Initialise()
     * 
     */
    void Initialise() {
        InitialisePow();
    } 

    /*
     * Report profiling outcomes
     * 
     * Put any code here required for reporting of profiling functionality
     * This function should be called once immediately prior to program termination
     * 
     * 
     * void Report()
     * 
     */
    void Report() {
        
        std::cout << "\n";
        std::cout << "========== BEGIN PROFILING OUTCOMES ==========\n";
        std::cout << "\n";

        ReportPow();
        
        std::cout << "=========== END PROFILING OUTCOMES ===========\n";
    } 



    /*
     * Initialises variables used for profiling std::pow()
     * 
     * 
     * void InitialisePow()
     * 
     */
    void InitialisePow() {

        m_PowBaseMin     =  __DBL_MAX__;        // minimum value of 'base' in calls to std::pow(base, exponent)
        m_PowBaseMax     = -__DBL_MAX__;        // maximum value of 'base' in calls to std::pow(base, exponent)
        m_PowExponentMin =  __DBL_MAX__;        // minimum value of 'exponent' in calls to std::pow(base, exponent)
        m_PowExponentMax = -__DBL_MAX__;        // maximum value of 'exponent' in calls to std::pow(base, exponent)

        m_PowCallCount   = 0;                   // Total (accumulated) call count for std::pow()
    }


    /*
     * Check call map to determine if specified signature and call count exists
     * 
     * This function ostensibly checks the call map for a signature and call count, and returns
     * true if an entry exists that matches both.  The map is checked - for the same thing -
     * a number of times in a loop (number determined by call count).  This is purely to chew up
     * CPU cycles.
     * 
     * The real purpose of this function is to expend some CPU time so this call is identified 
     * by an external profiling tool so the full stack trace can be identified.  This can be 
     * disabled if all we want is the immediate calling function and call count.
     * 
     * The function could really be named 'expendCPU', but it does actually check the map for
     * duplicate calls and the return value is valid, so 'CheckDuplicate' is appropriate (even
     * though it really only does the check to prevent the code from being optimised away by
     * the compiler - a plain loop that did nothing other than spin to soak up CPU cycles
     * would be optimised out of the code by a good compiler).
     * 
     * 
     * bool CheckDuplicate(const std::string &p_Signature, const int p_Count)
     *
     * @param   [IN]    p_Signature                 The signature of the call we're profiling
     * @param   [IN]    p_Count                     The call clount we're checking for
     * @return                                      Boolean flag indicating if the signature and call count were found in the map
     *  
     */
    bool CheckDuplicate(const std::string &p_Signature, const int p_Count) {
        
        bool duplicate = false;                         // result

        #ifndef PROFILING_COUNTS_ONLY                   // counts only? (i.e. no CPU spinning...)

        int cycles     = 0;                             // measure of CPU cycles to expend

        // soak up some CPU cycles
        // the CPU time expended depends on the call count
        // p_Count should never be < 1, but just in case we include it with < 5 
        // (previously it was being included with <100)

             if (p_Count <     5) cycles = 0;
        else if (p_Count <    10) cycles = p_Count;
        else if (p_Count <   100) cycles = 10;
        else if (p_Count <   500) cycles = 20;
        else if (p_Count <  1000) cycles = 100;
        else if (p_Count <  5000) cycles = 200;
        else if (p_Count < 10000) cycles = 1000;
        else                      cycles = 2000;

        for (int i = 0; i < cycles; i++) {                      // spin some CPU cycles
            duplicate = (m_CallMap[p_Signature] == p_Count);    // and don't get optimised away...
        }
        
        #endif

        return duplicate;
    }
    

    /*
     * Profiling function to replace std::pow when profiling is enabled
     * 
     * Gathers profiling statistics for calls to std::pow(base, exponent),
     * and returns the valued of std::pow(base, exponent)
     * 
     * Statistics gathered/recored are:
     * 
     *    - minimum value of 'base'
     *    - maximum value of 'base'
     *    - minimum value of 'exponent'
     *    - maximum value of 'exponent'
     *    - call signatures
     *    - counts for call signatures
     * 
     * 
     *  double pow(const double p_Base, const double p_Exponent)                                                (if PROFILING_CALLER_NAME is not #defined)
     *  double pow(const double p_Base, const double p_Exponent, const std::string &p_CallerName)               (if PROFILING_CALLER_NAME is #defined)
     * 
     * @param   [IN]    p_Base                      Floating point number to be raised to a power
     * @param   [IN]    p_Exponent                  Floating point exponent to which p_Base is to be raised
     * @param   [IN]    p_CallerName                String containing the name of the calling function          (if PROFILING_CALLER_NAME is #defined)
     * @return                                      Value of std::pow(p_Base, p_Exponent)
     *  
     */
    #ifdef PROFILING_CALLER_NAME                                                                            // add caller name to signature?
    double pow(const double p_Base, const double p_Exponent, const std::string &p_CallerName) {             // yes
    #else
    double pow(const double p_Base, const double p_Exponent) {                                              // no
    #endif

        // update minimums and maximums for parameters
        m_PowBaseMin     = std::min(p_Base, m_PowBaseMin);                                                  // minimum value of 'base'
        m_PowBaseMax     = std::max(p_Base, m_PowBaseMax);                                                  // maximum value of 'base'
        m_PowExponentMin = std::min(p_Exponent, m_PowExponentMin);                                          // minimum value of 'exponent'
        m_PowExponentMax = std::max(p_Exponent, m_PowExponentMax);                                          // maximum value of 'exponent'

        m_PowCallCount++;                                                                                   // increment accumulated call count
        
        std::string signature = "";                                                                         // call signature

        #ifdef PROFILING_CALLER_NAME                                                                        // add caller name to signature?
        signature = p_CallerName + "->";                                                                    // yes
        #endif

        signature += "pow(" + std::to_string(p_Base) + "," + std::to_string(p_Exponent) + ")";              // call signature
            
        auto index = m_CallMap.find(signature);                                                             // look for signature map
        if (index == m_CallMap.end()) {                                                                     // found?
            m_CallMap[signature] = 1;                                                                       // no - insert it - first call
        }
        else {                                                                                              // yes - found
            int count = m_CallMap[signature]++;                                                             // increment call count
            (void)CheckDuplicate(signature, count);                                                         // check for duplicates - actually spin some CPU cycles
        }

        return std::pow(p_Base, p_Exponent);                                                                 // and return the actual value
    }


    /*
     * Report profiling statistics for std::pow()
     * 
     * Minimamly formatted for machine parsing.  If the record contains a ":", then
     * the ": " is the delimiter, the LHS is the identifier, and the RHS is the value 
     * associated with the identifier
     *      
     * 
     * void ReportPow()
     * 
     */
    void ReportPow() {

        std::cout << "Begin profiling outcomes for std::pow()\n";
        std::cout << "---------------------------------------\n";
        std::cout << "Total accumulated call count: " << m_PowCallCount << "\n";
        std::cout << "\n";
        std::cout << "Minimum base passed         : " << m_PowBaseMin << "\n";
        std::cout << "Maximum base passed         : " << m_PowBaseMax << "\n";
        std::cout << "Minimum exponent passed     : " << m_PowExponentMin << "\n";
        std::cout << "Maximum exponent passed     : " << m_PowExponentMax << "\n";

        std::cout << "\n";
        std::cout << "---------------------------------------\n";
        std::cout << "\n";

        // print the call map entries
        std::unordered_map<std::string, int>::iterator it = m_CallMap.begin();
        while (it != m_CallMap.end())
        {
            std::string signature = it->first;                  // call signature
            int count = it->second;                             // call count
            std::cout << signature << ": " << count << "\n";    // print it
            it++;                                               // next map entry
        }

        std::cout << "\n";
        std::cout << "End profiling outcomes for std::pow()\n";
        std::cout << "---------------------------------------\n";
        std::cout << "\n";
    }

}

#endif // DOPROFILING