#ifndef __Errors_h__
#define __Errors_h__


// Early versions of COMPAS (versions prior to v03.00.00) did not have a coherent, robust error-handling strategy.
// In those versions, errors were typically displayed as either errors or warnings (depending on the severity) as
// they occurred, and evolution of the star (SSE mode) or binary (BSE mode) continued - users were expected to
// check errors or warnings displayed and use results with appropriate caution.  This was not ideal.
//
// In COMPAS version 03.00.00 the error handling philosophy was changed, and more coherent and robust error-handling
// code was implemented.  The new error-handling philosophy is to stop evolution of a star or binary if an error
// occurs, and record in the (SSE/BSE) system parameters file the fact that an error occurred, and an error number
// identifying the error that occurred.  This way users can check the system paramers file at the completion of a
// run for the disposition of a star or binary and, if the evolution of that star or binary was stopped because an
// error occurred, the actual error that occurred.  Possible dispositions (for both stars and binaries) are given in
// the 'EVOLUTION_STATUS' enum class and associate label map in 'typedefs.h'.
//
// The error-handling code implemented in v03.00.00 allows developers to terminate evolution of a star or binary if
// they determine that a condition encountered is sufficiently severe that allowing the evolution of the star or
// binary to continue would produce inconsistent or untrustworthy results.  In those cases, the developers should
// terminate the evolution of the star or binary via the use of the THROW_ERROR* macros (defined in 'ErrorsMacros.h').
//
// Developers should use the SHOW_WARN* macros (defined in 'ErrorsMacros.h') to alert users to conditions they want
// to bring to the attention of users, but are not sufficiently severe to warrant termination of the evolution of the
// star or binary.
//
// The SHOW_ERROR* macros (defined in 'ErrorsMacros.h') should be used sparingly - generally only in catch blocks
// for errors thrown, or in the (very few) sections of the code not covered by catch blocks.
//
// The class member variable m_Error (in the BaseStar class for SSE; BaseBinaryStar for BSE) should not be set explicitly
// throughout the code - it is set in the appropriate catch blocks.  m_Error is the error value written to the log files.
// Note that it is possible that if users choose to add the STAR_PROPERTY::ERROR (SSE) or BINARY_PROPERTY::ERROR (BSE) to
// log files via the logfile-definitions file, the value of the error logged to those files may be 0 (ERROR::NONE) even for
// stars or binaries that ere eventually terminated due to an error - the error value is only set when the error occurs
// (and is thrown), so some records in some log files may already have been written prior to the error being identified and
// evolution terminated.
//
//
// Floating-point errors in C++
// ============================
//
// In C++ implementations that implement the IEEE floating-point standard, in ordinary operation, the division of
// a finite non-zero floating-point value by 0[.0] is well-defined and results in +infinity if the value is > zero,
// -infinity if the value is < zero, and NaN if the value is = 0.0, and in each case program execution continues
// uninterrupted.  Integer division by 0 is undefined and results in a floating-point exception and the process
// is halted.
//
// The GNU C++ implementation allows us to trap the following floating-point errors:
//
// DIVBYZERO : division by zero, or some other asymptotically infinite result (from finite arguments).
// INEXACT   : a value cannot be represented with exact accuracy (e.g. 0.1, 1.0/3.0, and sqrt(2.0)).
// INVALID   : at least one of the arguments to a floating-point library function is a value for which the function
//             is not defined (e.g. sqrt(-1.0))
// OVERFLOW  : the result of an operation is too large in magnitude to be represented as a value of the return type.
// UNDERFLOW : the result is too small in magnitude to be represented as a value of the return type.
//
// When an enabled floating-point trap is encountered, a SIGFPE signal is raised.  If we don't have a signal handler
// installed for SIGFPE the program is terminated with a floating-point exception.  If we do have a signal handler
// installed for SIGFPE, that signal handler is invoked.  Ordinarily, once the SIGFPE signal handler is invoked, there
// is no going back - after doing whetever we need to do to manage the signal, the only valid operations are to exit
// the program, or to longjmp to a specific location in the code.  Fortunately the GNU C++ designers have given us another
// option: if we compile with the -fnon-call-exceptions compiler flag we can raise an exception safely from the SIGFPE
// signal handler because the throw is just a non-local transfer of control (just like a longjmp), and then we can just
// catch the exception raised.
//
//
// Floating-point errors in COMPAS
// ===============================
//
// Instrumentation has been implemented in COMPAS v03.00.00 that traps DIVBYZERO, INVALID, OVERFLOW, and UNDERFLOW.
// Trapping INEXACT would mean we'd trap on just about every floating-point calculation (it is really just
// informational - we know there are many values we can't represent exactly in base-2).
//
// We have 3 modes for the floating-point error instrumentation:
//
//    0: instrumentation not active   - this is just the default behaviour of the C++ compiler, as described in the first
//                                      paragraph above.  In this mode, the program execution will not be interrupted in
//                                      the event of a floating-point error, but the error reported in the system
//                                      parameters file will be set to indicate if a floating-point error occurred during
//                                      evolution (and in this mode we can, and do, differentiate between DIVBYZERO,
//                                      INVALID, OVERFLOW, and UNDERFLOW).  Note that an integer divide-by-zero will cause
//                                      the execution of the program to halt (and, rather obtusely, will report
//                                      "Floating point exception").
//
//                                      This mode is enabled with the option '--fp-error-mode OFF'  (This is the default mode).
//                             
//    1: floating-point traps enabled - floating-point traps DIVBYZERO, INVALID, OVERFLOW, and UNDERFLOW are enabled.
//                                      In this mode, when a floating-point operation traps, a SIGFPE is raised and the
//                                      SIGFPE signal handler is called, and the signal handler raises a runtime_error
//                                      exception, with a what argument of "FPE" (and in this mode we cannot, and do not,
//                                      differentiate between DIVBYZERO, INVALID, OVERFLOW, and UNDERFLOW).  The exception
//                                      raised will cause the execution of the program to halt if it is not caught and
//                                      managed.  We catch runtime_error exceptions in Star::Evolve() for SSE mode, in 
//                                      BaseBinaryStar::Evolve() for BSE mode, and in main() for errors that might occur
//                                      outside the evolution of stars or binaries.
//
//                                      This mode is enabled with the option '--fp-error-mode ON'
//
//    2: debug mode                   - floating-point traps DIVBYZERO, INVALID, OVERFLOW, and UNDERFLOW are enabled.  As
//                                      for mode 1, in this mode, when a floating-point operation traps, a SIGFPE is raised
//                                      and the SIGFPE signal handler is called, but instead of raising a runtime_error
//                                      exception, the signal handler prints the stack trace that led to the error and
//                                      halts execution of the program.  In this way, the user can determine where (to the
//                                      function level - we don't determine line numbers) the floating-point error occurred.
//
//                                      The construction of the stack trace in debug mode happens inside the signal handler,
//                                      and the functions used to do that are generally not signal safe - but we call std::exit()
//                                      anyway, so that shouldn't be a problem.
//
//
//                                      This mode is enabled with the option '--fp-error-mode DEBUG'



#define ERRORS Errors::Instance()

#include <iostream>
#include <unordered_map>
#include <vector>


#include "constants.h"
#include "typedefs.h"

#include "Log.h"
#include "ErrorsMacros.h"
#include "Options.h"

class Log;


/*
 * Errors Singleton
 *
 * Provides global error handling
 *
 * Singletons and global variables are sometimes frowned-upon, but doing it this
 * way means the objects don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Errors {

private:

    Errors() { };
    Errors(Errors const&) = delete;
    Errors& operator = (Errors const&) = delete;

    static Errors *m_Instance;                                                      // pointer to the instance

    // dynamic error catalog
    // this unordered_map records the errors printed and thedetails (object type, object id, function name etc.)
    // there are 2 objectId vectors (non-stellar and stellar object), and similarly 2 funcName vectors
    // it's a bit of a hack, but splitting them out allows me to clean out the stellar vectors after each
    // star/binary is evolved, so the catalog isn't bloated by details for deleted objects
    COMPASUnorderedMap<
        ERROR,                              // error id
        std::tuple<                         // details for error id
            ERROR_SCOPE,                    //    scope
            bool,                           //    flag indicating if already printed
            std::vector<OBJECT_TYPE>,       //    object type
            ST_VECTOR,                      //    stellar type
            std::vector<OBJECT_ID>,         //    vector of non-stellar (main, utils, etc) object ids
            std::vector<OBJECT_ID>,         //    vector of stellar ids
            std::vector<std::string>,       //    vector of function names for non-stellar object ids
            std::vector<std::string>,       //    vector of function names for stellar ids
            std::string                     //    error text
        >
    > m_ErrorCatalog = {};


public:

    ~Errors() { delete m_Instance; }

    static Errors* Instance();


    bool ShowIt(const std::string  p_Prefix,
                const ERROR        p_Error,
                const std::string  p_QualifyingStr = "",
                const OBJECT_ID    p_ObjectId      = 0,
                const OBJECT_TYPE  p_ObjectType    = OBJECT_TYPE::NONE,
                const STELLAR_TYPE p_StellarType   = STELLAR_TYPE::NONE,
                const char*        p_FuncName      = "Not provided");

    void Clean();

    size_t CatalogSize() { return m_ErrorCatalog.size(); }
};

#endif // __Errors_h_
