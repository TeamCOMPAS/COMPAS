//
//  COMPAS main
//
#include <ctime>
#include <chrono>
#include <string>
#include <sstream>
#include <fstream>
#include <tuple>
#include <vector>
#include <csignal>
#include <iostream>

#include <memory>
#include <stdexcept> 
#include <setjmp.h>
#include <cstring>
#include <execinfo.h>

#include <fenv.h>

#include "constants.h"
#include "typedefs.h"

#include "errors.h"
#include "profiling.h"
#include "utils.h"
#include "yaml.h"
#include "vector3d.h"
#include "Options.h"
#include "Rand.h"
#include "Log.h"

#include "Star.h"
#include "BinaryStar.h"

OBJECT_ID globalObjectId = 1;                                   // used to uniquely identify objects - used primarily for error printing
OBJECT_ID m_ObjectId     = 0;                                   // object id for main - always 0


// The following variable supports the floating-point error instrumentation
//
// In C++ implentations that implement the IEEE floating-point standard, in ordinary operation, the division of a
// finite non-zero floating-point value by 0[.0] is well-defined and results in +infinity if the value is > zero,
// -infinity if the value is < zero, and NaN if the value is = 0.0, and in each case program execution continue
// uninterrupted.  Integer division by 0 is undefined and results in a floating-point exception and the process
// is halted.
//
// The GNU C++ implementation allows us to trap the following floating-pont errors:
//
// DIVBYZERO : division by zero, or some other asymptotically infinite result (from finite arguments).
// INEXACT   : a value cannot be represented with exact accuracy (e.g. 0.1, 1.0/3.0, and sqrt(2.0)).
// INVALID   : at least one of the arguments to a floating-point library function is a value for which the function
//             is not defined (e.g. sqrt(-1.0))
// OVERFLOW  : the result of an operation is too large in magnitude to be represented as a value of the return type.
// UNDERFLOW : the result is too small in magnitude to be represented as a value of the return type.
//
// The instrumentation implemented traps DIVBYZERO, INVALID, OVERFLOW, and UNDERFLOW.  Trapping INEXACT would mean
// we'd trap on just about every floating-point calculation (it's really just informational - we know there are many
// values we can't represent exactly in base-2).
//
// When an enabled floating-point trap is encountered, a SIGFPE signal is raised.  If we don't have a signal handler
// installed for SIGFPE the program is terminated with a floating-point exception.  If we do have a signal handler
// installed for SIGFPE, that signal handler is called.  Ordinarily, once the SIGFPE signal handler is called, there
// is no going back - after doing whetever we need to do to manage the signal, the only valid operation is to exit the
// program, or to longjmp to a specific location in the code.  Fortunately the GNU C++ designers have given us another
// option: if we compile with the -fnon-call-exceptions compiler flag, then we can raise an exception safely from the
// SIGFPE signal handler because the throw is just a non-local transfer of control (just like a longjmp), and then we
// can just catch the exception raised.
//
// We have 3 modes for the floating-point error instrumentation:
//
//    0: instrumentation not active   - this is just the default behaviour of the C++ comoiler, as described in the first
//                                      paragraph above.  In this mode, the program execution will not be interrupted in
//                                      the event of a floating-point error, but the error reported in the system
//                                      parameters file will be set to indicate if a floating-point error occurred during
//                                      evolution (and in this mode we can, and do, differentiate between DIVBYZERO,
//                                      INVALID, OVERFLOW, and UNDERFLOW).  Note that an integer divide-by-zero will cause
//                                      the excution of the program to halt (and, rather obtusely, will report
//                                      "Floating point exception").
//
//                                      This is the default mode.
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
//    2: debug mode                   - floating-point traps DIVBYZERO, INVALID, OVERFLOW, and UNDERFLOW are enabled.  As
//                                      for mode 1, in this mode, when a floating-point operation traps, a SIGFPE is raised
//                                      and the SIGFPE signal handler is called, but instead of raising a runtime_error
//                                      exception, the signal handler prints the stack trace that led to the error and
//                                      halts execution of the program.  In this way, the user can determine where (to the
//                                      function leve - we don't determine line numbers) the floating-point error occured.
//
//                                      Note here the cevaetas - not signal safe etc. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ FIX THIS - NO LONGE A GLOBAL 

/*
 * SIGFPE signal handler
 * 
 * Only handles SIGFPE.
 * 
 * DESCRIPTION HERE - JR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * 
 * 
 * void SIGFPEhandler(int p_Sig)
 * 
 * @param   [IN]        p_Sig                   The signal intercepted
 *  
 */
void SIGFPEhandler(int p_Sig) {
    if (p_Sig == SIGFPE) {                                      // SIGFPE?
        FP_ERROR_MODE fpErrorMode = OPTIONS->FPErrorMode();
        if (fpErrorMode == FP_ERROR_MODE::ON) throw std::runtime_error("FPE");
        else if (fpErrorMode == FP_ERROR_MODE::DEBUG) {
            std::cerr << "\nFloating-point error encountered: program terminated\n";
            utils::ShowStackTrace(); 
        }        
    }
    else std::cerr << "\nUnexpected signal encountered: program terminated\n";
    std::exit(1);                                               // catch-all in case we're here when we shouldn't be
}


class Star;
class BinaryStar;

OBJECT_ID          ObjectId()          { return m_ObjectId; }
OBJECT_TYPE        ObjectType()        { return OBJECT_TYPE::MAIN; }
OBJECT_PERSISTENCE ObjectPersistence() { return OBJECT_PERSISTENCE::PERMANENT; }
STELLAR_TYPE       StellarType()       { return STELLAR_TYPE::NONE; }


// The following global variables support the BSE Switch Log file
// Ideally, rather than be declared as globals, they would be in maybe the 
// LOGGING service singleton, but the Log class knows nothing about the 
// BinaryStar class...
// (maybe we could put them in the new CONSTANTS service singleton if we 
// implement it)

BinaryStar* evolvingBinaryStar      = NULL;             // pointer to the currently evolving Binary Star
bool        evolvingBinaryStarValid = false;            // flag to indicate whether the evolvingBinaryStar pointer is valid

/*
 * SIGUSR1 signal handler
 * 
 * Only handles SIGUSR1.
 * 
 * SIGUSR1 is a user generated signal - the system should not generate this signal,
 * though it is possible to send the signal to a process via the Un*x kill command,
 * or some other user-developed program that sends signals.  This code does some 
 * rudimentary sanity checks, but it is possible that sending a SIGUSR1 signal to a
 * running COMPAS process via the Un*x kill command, or otherwise, might cause a 
 * spurious entry in the BSE Switch Log file - c'est la vie.
 * 
 * We use SIGUSR1 in the Star class to signal when a Star object switches stellar 
 * type. We use a signal because the Star class knows nothing about binary stars, 
 * so can't call a binary star function to log binary star variables to the BSE 
 * Switch Log file. By raising a signal in the Star class and catching it here we 
 * can call the appropriate binary star class function to write the binary star 
 * variables to the log file.
 * 
 * The signal is raised in the Star::SwitchTo() function if OPTIONS->BSESwitchLog() 
 * is true, so the signal will be received here for every stellar type switch of 
 * every star.
 * 
 * We only process the signal here if the global variable evolvingBinaryStarValid 
 * is true. The global variable evolvingBinaryStarValid is only set true after a 
 * binary star has been constructed and is ready to evolve - so if the signal is 
 * raised, it will be ignored for SSE switches, and it will be ignored for switches 
 * inside the constructor of the binary star (and so its constituent stars).
 * 
 * This signal handler is installed in EvolveBinaryStars(), so it is installed only
 * if we're evolving binaries - signals will be ignored (by our code - the system
 * will still receive and handle them) if we're evolving single stars.
 *
 * 
 * void SIGUSR1handler(int p_Sig)
 * 
 * @param   [IN]        p_Sig                   The signal intercepted
 * 
 */
void SIGUSR1handler(int p_Sig) {   
    if (p_Sig == SIGUSR1) {                                         // SIGUSR1?  Just silently ignore anything else...
        if (evolvingBinaryStarValid && OPTIONS->SwitchLog()) {      // yes - do we have a valid binary star, and are we logging switches?
            (void)evolvingBinaryStar->PrintSwitchLog();             // yes - assume SIGUSR1 is a binary constituent star switching...
        }
    }
}


/*
 * Evolve single stars
 *
 *
 * std::tuple<int, int> EvolveSingleStars()
 * 
 * @return                                      Tuple: <number of stars requested, actual number of stars created>
 */
std::tuple<int, int> EvolveSingleStars() {

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;
    EVOLUTION_STATUS thisStarStatus  = EVOLUTION_STATUS::NOT_STARTED;

    auto    wallStart  = std::chrono::system_clock::now();                                                          // start wall timer
    clock_t clockStart = clock();                                                                                   // start CPU timer

    std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
    SAY("Start generating stars at " << std::ctime(&timeStart));

    // generate and evolve stars

    Star*  star      = nullptr;
    bool   usingGrid = !OPTIONS->GridFilename().empty();                                                            // using grid file?
    size_t index     = 0;                                                                                           // which star

    // The options specified by the user at the commandline are set to their initial values.
    // OPTIONS->AdvanceCmdLineOptionValues(), called at the end of the loop, advances the
    // options specified by the user at the commandline to their next variation (if necessary,
    // based on any ranges and/or sets specified by the user).

    while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                         // while all ok

        // generate and evolve stars

        int gridLineVariation   = 0;                                                                                // grid line variation number
        bool doneGridFile       = false;                                                                            // flags we're done with the grid file (for this commandline variation)
        bool processingGridLine = false;                                                                            // processing a gridfile line?
        while (!doneGridFile && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                    // for each star to be evolved

            try {
                if (OPTIONS->FPErrorMode() != FP_ERROR_MODE::OFF) {                                                 // floating-point error handling mode on/debug?
                    feclearexcept(FE_ALL_EXCEPT);                                                                   // yes - clear all FE traps
                    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);                         // enable FE traps (don't trap FE_INEXACT - would trap on almost all FP operations...)
                }

            bool doneGridLine = false;                                                                              // initially
            if (usingGrid) {                                                                                        // using grid file?
                gridLineVariation = 0;                                                                              // yes - first variation of this grid line
                int gridResult = OPTIONS->ApplyNextGridLine();                                                      // set options according to specified values in grid file              
                ERROR error;
                switch (gridResult) {                                                                               // handle result of grid file read
                    case -1:                                                                                        // error - unexpected end of grid file
                    case -2:                                                                                        // error - read error for grid file
                    case -3: {                                                                                      // error - invalid value in grid file
                        switch (gridResult) {
                            case -1: error = ERROR::UNEXPECTED_END_OF_FILE; break;  
                            case -2: error = ERROR::FILE_READ_ERROR; break;  
                            case -3: error = ERROR::INVALID_VALUE_IN_FILE; break;
                            default: error = ERROR::FILE_READ_ERROR; break;
                        }

                        SHOW_ERROR(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'");                 // show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // stop evolution
                    } break;
                    case  0: {                                                                                      // end of file
                        doneGridLine = true;                                                                        // flag we're done with this grid line
                        doneGridFile = true;                                                                        // flag we're done with the grid file
                        error = OPTIONS->RewindGridFile();                                                          // ready for next commandline options variation
                        if (error != ERROR::NONE) {                                                                 // rewind ok?
                            SHOW_ERROR(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'");             // no - show error (should never happen here - should be picked up at file open)
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                              // stop evolution
                        }
                        } break;
                    case  1:                                                                                        // grid line read
                        processingGridLine = true;                                                                  // not done yet...
                        break;
                    default:                                                                
                        SHOW_ERROR(ERROR::ERROR, "Accessing grid file '" + OPTIONS->GridFilename() + "'");          // unexpected error - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // stop evolution
                        break;
                }
            }
            else {                                                                                                  // no, not using a grid file
                doneGridFile = true;                                                                                // flag we're done with the grid file
            }


            // The options specified by the user in the grid file line are set to their initial values.
            // (by OPTIONS->ApplyNextGridLine()).
            // OPTIONS->AdvanceGridLineOptionValues(), called at the end of the loop, advances the
            // options specified by the user in the grid file line to their next variation (if necessary,
            // based on any ranges and/or sets specified by the user).
            // Note that `doneGridLine` may be a proxy for the command line here (it is set FALSE even if
            // there is no grid file - so the loop executes once to process the command line)

            while (!doneGridLine && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                // while all ok and not done

                // Single stars (in SSE) are provided with a random seed that is used to seed the random 
                // number generator.  The random number generator is re-seeded for each star.  Here we 
                // generate the seed for the star being evolved - by this point we have picked up the 
                // option value from either the commandline or the grid file.
                //
                // there are three scenarios:
                //
                // if the user did not specify a random seed, either on the commandline or in a grid file
                // record, we use a randomly chosen seed, based on the system time.
                //
                // if the user specified a random seed on the commandline, and not in the grid file for
                // the current star, the random seed specified on the commandline is used - and the offset 
                // applied (the index of the star being evolved).  The index of the star being evolved 
                // starts at 0 for the first star, and increments by 1 for each subsequent star evolved
                // (so the base random seed specified by the user is also the initial random seed - the 
                // random seed of the first star evolved)
                //
                // if the user specified a random seed in the grid file for the current star, regardless of
                // whether a random seed was specified on the commandline, the random seed from the grid
                // file is used, and an offset is added if the grid line also specified ranges or sets for
                // and options (if no rangers or sets were specified on the grid line then no offset is added
                // (i.e. the random seed specified is used as it)).  Note that in this scenario it is the 
                // user's responsibility to ensure that there is no duplication of seeds.

                unsigned long int randomSeed = 0l;                                                                  // random seed
                OPTIONS_ORIGIN    optsOrigin = processingGridLine ? OPTIONS_ORIGIN::GRIDFILE : OPTIONS_ORIGIN::CMDLINE; // indicate which set of program options we're using
                if (OPTIONS->FixedRandomSeedGridLine()) {                                                           // user specified a random seed in the grid file for this binary?
                                                                                                                    // yes - use it (indexed)
                    randomSeed = OPTIONS->RandomSeedGridLine() + (unsigned long int)gridLineVariation;              // random seed               
                    if (OPTIONS->SetRandomSeed(randomSeed, optsOrigin) < 0) {                                       // ok?
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                                       // no - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // and stop evolution
                    }
                }
                else if (OPTIONS->FixedRandomSeedCmdLine()) {                                                       // no - user specified a random seed on the commandline?
                                                                                                                    // yes - use it (indexed)
                    randomSeed = OPTIONS->RandomSeedCmdLine() + (unsigned long int)index;                           // random seed               
                    if (OPTIONS->SetRandomSeed(randomSeed, optsOrigin) < 0) {                                       // ok?
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                        // no - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // and stop evolution
                    }
                }
                else {                                                                                              // no
                                                                                                                    // use default seed (based on system time) + id (index)
                    randomSeed = RAND->DefaultSeed() + (unsigned long int)index;                                    // random seed               
                    if (OPTIONS->SetRandomSeed(randomSeed, optsOrigin) < 0) {                                       // ok?
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                        // no - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // and stop evolution
                    }
                }

                if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                // ok?
                                                                                                                    // yes - continue
                    randomSeed = RAND->CurrentSeed();                                                               // current random seed - to pass to star object

                    // the initial mass of the star is supplied - this is to allow binary stars to initialise
                    // the masses of their constituent stars (rather than have the constituent stars sample 
                    // their own mass).  Here we use the mass supplied by the user via the program options or, 
                    // if no mass was supplied by the user, sample the mass from the IMF.

                    double initialMass = OPTIONS->OptionSpecified("initial-mass") == 1                              // user specified mass?
                                            ? OPTIONS->InitialMass()                                                // yes, use it
                                            : utils::SampleInitialMass(OPTIONS->InitialMassFunction(),              // no, sample it
                                                                       OPTIONS->InitialMassFunctionMax(), 
                                                                       OPTIONS->InitialMassFunctionMin(), 
                                                                       OPTIONS->InitialMassFunctionPower());

                    // the metallicity of the star is supplied - this is to allow binary stars to initialise
                    // the metallicity of their constituent stars (rather than have the constituent stars sample 
                    // their own metallicity).  Here we use the mmetallicityass supplied by the user via the program
                    // options or, if no metallicity was supplied by the user, sample the metallicity.

                    double metallicity = OPTIONS->OptionSpecified("metallicity") == 1                               // user specified metallicity?
                                            ? OPTIONS->Metallicity()                                                // yes, use it
                                            : utils::SampleMetallicity(OPTIONS->MetallicityDistribution(), 
                                                                       OPTIONS->MetallicityDistributionMax(), 
                                                                       OPTIONS->MetallicityDistributionMin());      // no, sample it



                    // Single stars (in SSE) are provided with a kick structure that specifies the 
                    // values of the random number to be used to generate to kick magnitude, and the
                    // actual kick magnitude specified by the user via program option --kick-magnitude       
                    //
                    // See typedefs.h for the kick structure.
                    //
                    // We can't just pick up the values of the options inside Basestar.cpp because the
                    // constituents of binaries get different values, so use different options. The
                    // Basestar.cpp code doesn't know if the star is a single star (SSE) or a constituent
                    // of a binary (BSE) - it only knows that it is a star - so we have to setup the kick
                    // structure here.
                    //
                    // for SSE only need magnitudeRandom and magnitude - other values can just be ignored

                    KickParameters kickParameters;
                    kickParameters.magnitudeRandomSpecified = OPTIONS->OptionSpecified("kick-magnitude-random") == 1;
                    kickParameters.magnitudeRandom          = OPTIONS->KickMagnitudeRandom();
                    kickParameters.magnitudeSpecified       = OPTIONS->OptionSpecified("kick-magnitude") == 1;
                    kickParameters.magnitude                = OPTIONS->KickMagnitude();
                       
                    // create the star
                    delete star; star = nullptr;                                                                    // so we don't leak...
                    star = OPTIONS->OptionSpecified("rotational-frequency") == 1                                    // user specified rotational frequency?
                        ? new Star(randomSeed, initialMass, metallicity, kickParameters, OPTIONS->RotationalFrequency() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                        : new Star(randomSeed, initialMass, metallicity, kickParameters);                           // no - let it be calculated

                    thisStarStatus = EVOLUTION_STATUS::STARTED;

                    thisStarStatus = star->Evolve(index);                                                           // evolve the star

                    // announce the result
                    if (!OPTIONS->Quiet()) {                                                                        // quiet mode?
                        SAY(index                                     <<                                            // announce result of evolving the star
                            ": "                                      <<
                            EVOLUTION_STATUS_LABEL.at(thisStarStatus) <<                  
                            ": RandomSeed = "                         <<
                            star->RandomSeed()                        <<
                            ", Initial Mass = "                       <<
                            star->MZAMS()                             <<
                            ", Metallicity = "                        <<
                            star->Metallicity()                       <<
                            ", "                                      <<
                            STELLAR_TYPE_LABEL.at(star->StellarType()));
                    }

                    if (!LOGGING->CloseStandardFile(LOGFILE::SSE_DETAILED_OUTPUT)) {                                // close SSE detailed output file if necessary
                        SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                          // close failed - show warning
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // this will cause problems later - stop evolution
                    }

                    ERRORS->Clean();                                                                                // clean the dynamic error catalog

                    index++;                                                                                        // next...

                    if (usingGrid) {                                                                                // using grid file?
                        gridLineVariation++;                                                                        // yes - increment grid line variation number
                        int optionsStatus = OPTIONS->AdvanceGridLineOptionValues();                                 // apply next grid file options (ranges/sets)
                        if (optionsStatus < 0) {                                                                    // ok?
                            SHOW_ERROR(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                                   // no - show error
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                              // and stop evolution
                        }
                        else if (optionsStatus == 0) {                                                              // end of grid file options variations?
                            doneGridLine = true;                                                                    // yes - we're done
                        }
                    }
                    else doneGridLine = true;                                                                       // not using grid file - done    
                }
            }
            
            }
            catch (const std::runtime_error& e) {                                                                       // catch runtime exceptions
                evolutionStatus = EVOLUTION_STATUS::ERROR;                                                              // evolution terminated
//                if (std::string(e.what()) == "FPE") m_Star->SetError(ERROR::FLOATING_POINT_ERROR);                      // floating-point error
//                else                                m_Star->SetError(ERROR::ERROR);                                     // unspecified error
            }
            catch (int e) {
                evolutionStatus = EVOLUTION_STATUS::ERROR;                                                              // evolution terminated
//                if (e != static_cast<int>(ERROR::NONE)) m_Star->SetError(static_cast<ERROR>(e));                        // specified errpr
//                else                                    m_Star->SetError(ERROR::ERROR);                                 // unspecified error
            }
            catch (...) {
                evolutionStatus = EVOLUTION_STATUS::ERROR;                                                              // evolution terminated
//                m_Star->SetError(ERROR::ERROR);                                                                         // unspecified error
            }


        }
        delete star; star = nullptr;                                                                                // so we don't leak...
    
        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                        // ok?
            int optionsStatus = OPTIONS->AdvanceCmdLineOptionValues();                                              // yes - apply next commandline options (ranges/sets)
            if (optionsStatus < 0) {                                                                                // ok?
                SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                                // no - show error
                evolutionStatus = EVOLUTION_STATUS::ERROR;                                                          // and stop evolution
            }
            else if (optionsStatus == 0) {                                                                          // end of options variations?
                if (usingGrid || OPTIONS->CommandLineGrid() || (!usingGrid && index >= OPTIONS->nObjectsToEvolve())) { // created required number of stars?
                    evolutionStatus = EVOLUTION_STATUS::DONE;                                                       // yes - we're done
                }
            }
        }
    }

    int nStarsRequested = evolutionStatus == EVOLUTION_STATUS::DONE ? index : (usingGrid ? -1 : OPTIONS->nObjectsToEvolve());

    SAY("\nGenerated " << std::to_string(index) << " of " << (nStarsRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nStarsRequested)) << " stars requested");

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                        // shouldn't be
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN(ERROR::STELLAR_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));       // show warning
        }
    }

    // close SSE logfiles
    // don't check result here - let log system handle it
    (void)LOGGING->CloseAllStandardFiles();                                                                         // close any standard log files

    // announce timing stats
    double cpuSeconds = (clock() - clockStart) / (double) CLOCKS_PER_SEC;                                           // stop CPU timer and calculate seconds

    auto wallEnd = std::chrono::system_clock::now();                                                                // stop wall timer
    std::time_t timeEnd = std::chrono::system_clock::to_time_t(wallEnd);                                            // get end time and date

    SAY("\nEnd generating stars at " << std::ctime(&timeEnd));
    SAY("Clock time = " << cpuSeconds << " CPU seconds");

    std::chrono::duration<double> wallSeconds = wallEnd - wallStart;                                                // elapsed seconds

    int wallHH = (int)(wallSeconds.count() / 3600.0);                                                               // hours
    int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                                   // minutes
    int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));                  // seconds

    SAY("Wall time  = " << std::setfill('0') << std::setw(2) << wallHH << ":" << 
                           std::setfill('0') << std::setw(2) << wallMM << ":" << 
                           std::setfill('0') << std::setw(2) << wallSS << " (hh:mm:ss)");                           // Include 0 buffer 

    return  std::make_tuple(nStarsRequested, index);
}


/*
 * Evolve binary stars
 *
 *
 * std::tuple<int, int> EvolveBinaryStars()
 * 
 * @return                                      Tuple: <number of binaries requested, actual number of binaries created>
 */
std::tuple<int, int> EvolveBinaryStars() {

    signal(SIGUSR1, SIGUSR1handler);                                                                                // install SIGUSR1 signal handler

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    auto wallStart = std::chrono::system_clock::now();                                                              // start wall timer
    clock_t clockStart = clock();                                                                                   // start CPU timer

    std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
    SAY("Start generating binaries at " << std::ctime(&timeStart));

    BinaryStar* binary    = nullptr;
    bool        usingGrid = !OPTIONS->GridFilename().empty();                                                       // using grid file?
    size_t      index     = 0;                                                                                      // which binary

    // The options specified by the user at the commandline are set to their initial values.
    // OPTIONS->AdvanceCmdLineOptionValues(), called at the end of the loop, advances the
    // options specified by the user at the commandline to their next variation (if necessary,
    // based on any ranges and/or sets specified by the user).

    while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                         // while all ok and not done

        // generate and evolve binaries

        int  gridLineVariation  = 0;                                                                                // grid line variation number
        bool doneGridFile       = false;                                                                            // flags we're done with the grid file (for this commandline variation)
        bool processingGridLine = false;                                                                            // processing a gridfile line?
        while (!doneGridFile && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                    // for each binary to be evolved

            evolvingBinaryStar      = NULL;                                                                         // unset global pointer to evolving binary (for BSE Switch Log)
            evolvingBinaryStarValid = false;                                                                        // indicate that the global pointer is not (yet) valid (for BSE Switch log)

            bool doneGridLine = false;                                                                              // flags we're done with this grid file line (if using a grid file)
            if (usingGrid) {                                                                                        // using grid file?
                gridLineVariation = 0;                                                                              // yes - first variation of this grid line
                int gridResult = OPTIONS->ApplyNextGridLine();                                                      // yes - set options according to specified values in grid file              
                switch (gridResult) {                                                                               // handle result of grid file read
                    case -1:                                                                                        // error - unexpected end of file grid file read
                    case -2: {                                                                                      // error - read error for grid file
                        ERROR error = gridResult == -1 ? ERROR::UNEXPECTED_END_OF_FILE : ERROR::FILE_READ_ERROR;    // set error
                        SHOW_ERROR(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'");                 // show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // stop evolution
                    } break;
                    case  0: {                                                                                      // end of file
                        doneGridLine = true;                                                                        // flag we're done with this grid line
                        doneGridFile = true;                                                                        // flag we're done with the grid file
                        ERROR error = OPTIONS->RewindGridFile();                                                    // ready for next commandline options variation
                        if (error != ERROR::NONE) {                                                                 // rewind ok?
                            SHOW_ERROR(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'");             // no - show error (should never happen here - should be picked up at file open)
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                              // stop evolution
                        }
                    } break;
                    case  1:                                                                                        // grid line read
                        processingGridLine = true;                                                                  // not done yet...
                        break;
                    default:                                                                
                        SHOW_ERROR(ERROR::ERROR, "Accessing grid file '" + OPTIONS->GridFilename() + "'");          // unexpected error - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // stop evolution
                        break;
                }
            }
            else {                                                                                                  // no, not using a grid file
                doneGridFile = true;                                                                                // flag we're done with the grid file
            }

            // The options specified by the user in the grid file line are set to their initial values.
            // (by OPTIONS->ApplyNextGridLine()).
            // OPTIONS->AdvanceGridLineOptionValues(), called at the end of the loop, advances the
            // options specified by the user in the grid file line to their next variation (if necessary,
            // based on any ranges and/or sets specified by the user).
            // Note that `doneGridLine` may be a proxy for the command line here (it is set FALSE even if
            // there is no grid file - so the loop executes once to process the command line)

            while (!doneGridLine && (evolutionStatus == EVOLUTION_STATUS::CONTINUE)) {                              // while all ok and not done

                // we only need to pass the index number to the binary - we let the BinaryStar class do the work 
                // wrt setting the parameters for each of the constituent stars
                // (The index is really only needed for legacy comparison, so can probably be removed at any time)

                // create the binary

                // Binary stars (in BSE) are provided with a random seed that is used to seed the random 
                // number generator.  The random number generator is re-seeded for each binary.  Here we 
                // generate the seed for the binary being evolved - by this point we have picked up the 
                // option value from either the commandline or the grid file.
                //
                // there are three scenarios:
                //
                // if the user did not specify a random seed, either on the commandline or in a grid file
                // record, we use a randomly chosen seed, based on the system time.
                //
                // if the user specified a random seed on the commandline, and not in the grid file for
                // the current binary, the random seed specified on the commandline is used - and the offset 
                // applied (the index of the binary being evolved).  The index of the binary being evolved 
                // starts at 0 for the first binary, and increments by 1 for each subsequent binary evolved
                // (so the base random seed specified by the user is also the initial random seed - the 
                // random seed of the first binary evolved)
                //
                // if the user specified a random seed in the grid file for the current binary, regardless of
                // whether a random seed was specified on the commandline, the random seed from the grid
                // file is used, and an offset is added if the grid line also specified ranges or sets for
                // and options (if no rangers or sets were specified on the grid line then no offset is added
                // (i.e. the random seed specified is used as it)).  Note that in this scenario it is the 
                // user's responsibility to ensure that there is no duplication of seeds.
             
                unsigned long int thisId = index + gridLineVariation;                                               // set the id for the binary

                unsigned long int randomSeed = 0l;                                                                  // random seed
                OPTIONS_ORIGIN    optsOrigin = processingGridLine ? OPTIONS_ORIGIN::GRIDFILE : OPTIONS_ORIGIN::CMDLINE; // indicate which set of program options we're using
                if (OPTIONS->FixedRandomSeedGridLine()) {                                                           // user specified a random seed in the grid file for this binary?
                                                                                                                    // yes - use it (indexed)
                    randomSeed = OPTIONS->RandomSeedGridLine() + (unsigned long int)gridLineVariation;              // random seed               
                    if (OPTIONS->SetRandomSeed(randomSeed, optsOrigin) < 0) {                                       // ok?
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                                       // no - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // and stop evolution
                    }
                }
                else if (OPTIONS->FixedRandomSeedCmdLine()) {                                                       // no - user specified a random seed on the commandline?
                                                                                                                    // yes - use it (indexed)
                    randomSeed = OPTIONS->RandomSeedCmdLine() + (unsigned long int)index + (unsigned long int)gridLineVariation; // random seed               
                    if (OPTIONS->SetRandomSeed(randomSeed, optsOrigin) < 0) {                                       // ok?
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                        // no - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // and stop evolution
                    }
                }
                else {                                                                                              // no
                                                                                                                    // use default seed (based on system time) + id (index)
                    randomSeed = RAND->DefaultSeed() + (unsigned long int)index + (unsigned long int)gridLineVariation; // random seed               
                    if (OPTIONS->SetRandomSeed(randomSeed, optsOrigin) < 0) {                                       // ok?
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                        // no - show error
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // and stop evolution
                    }
                }

                if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                // ok?
                                                                                                                    // yes - continue
                    randomSeed = RAND->CurrentSeed();                                                               // current random seed - to pass to binary object


                    delete binary; binary = nullptr;                                                                // so we don't leak
                    binary = new BinaryStar(randomSeed, thisId);                                                    // generate binary according to the user options

                    evolvingBinaryStar      = binary;                                                               // set global pointer to evolving binary (for BSE Switch Log)
                    evolvingBinaryStarValid = true;                                                                 // indicate that the global pointer is now valid (for BSE Switch Log)

                    EVOLUTION_STATUS binaryStatus = binary->Evolve();                                               // evolve the binary

                    if (binaryStatus == EVOLUTION_STATUS::ERROR || binaryStatus == EVOLUTION_STATUS::SSE_ERROR) {   // ok?
                        SHOW_ERROR(ERROR::BINARY_EVOLUTION_STOPPED, EVOLUTION_STATUS_LABEL.at(binaryStatus));       // no - show error
                    }
                
                    // announce result of evolving the binary
                    if (!OPTIONS->Quiet()) {                                                                        // quiet mode?
                                                                                                                    // no - announce result of evolving the binary
                        if (OPTIONS->CHEMode() == CHE_MODE::NONE) {                                                 // CHE enabled?
                            SAY(thisId                                     << ": "  <<                              // no - CHE not enabled - don't need initial stellar type
                                EVOLUTION_STATUS_LABEL.at(binaryStatus)    << ": "  <<
                                STELLAR_TYPE_LABEL.at(binary->Star1Type()) << " + " <<
                                STELLAR_TYPE_LABEL.at(binary->Star2Type())
                            );
                        }
                        else {                                                                                      // CHE enabled - show initial stellar type
                            SAY(thisId                                            << ": "    <<
                                EVOLUTION_STATUS_LABEL.at(binaryStatus)           << ": ("   <<
                                STELLAR_TYPE_LABEL.at(binary->Star1InitialType()) << " -> "  <<
                                STELLAR_TYPE_LABEL.at(binary->Star1Type())        << ") + (" <<
                                STELLAR_TYPE_LABEL.at(binary->Star2InitialType()) << " -> "  <<
                                STELLAR_TYPE_LABEL.at(binary->Star2Type())        <<  ")"
                            );
                        }
                    }

                    if (!LOGGING->CloseStandardFile(LOGFILE::BSE_DETAILED_OUTPUT)) {                                // close detailed output file if necessary
                        SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                          // close failed - show warning
                        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // this will cause problems later - stop evolution
                    }

                    ERRORS->Clean();                                                                                // clean the dynamic error catalog

                    if (usingGrid) {                                                                                // using grid file?
                        gridLineVariation++;                                                                        // yes - increment grid line variation number
                        int optionsStatus = OPTIONS->AdvanceGridLineOptionValues();                                 // apply next grid file options (ranges/sets)
                        if (optionsStatus < 0) {                                                                    // ok?
                            SHOW_ERROR(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                                   // no - show error
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                              // and stop evolution
                        }
                        else if (optionsStatus == 0) {                                                              // end of grid file options variations?
                            doneGridLine = true;                                                                    // yes - we're done
                        }
                    }
                    else doneGridLine = true;                                                                       // not using grid file - done

                    if (doneGridLine) index = thisId + 1;                                                           // increment index
                }
            }
        }
        delete binary; binary = nullptr;

        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                        // ok?
            int optionsStatus = OPTIONS->AdvanceCmdLineOptionValues();                                              // apply next commandline options (ranges/sets)
            if (optionsStatus < 0) {                                                                                // ok?
                SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                                // no - show error
                evolutionStatus = EVOLUTION_STATUS::ERROR;                                                          // and stop evolution
            }
            else if (optionsStatus == 0) {                                                                          // end of options variations?
                if (usingGrid || OPTIONS->CommandLineGrid() || (!usingGrid && index >= OPTIONS->nObjectsToEvolve())) { // created required number of stars?
                    evolutionStatus = EVOLUTION_STATUS::DONE;                                                       // yes - we're done
                }
            }
        }
    }
    
    int nBinariesRequested = evolutionStatus == EVOLUTION_STATUS::DONE ? index : (usingGrid ? -1 : OPTIONS->nObjectsToEvolve());

    SAY("\nGenerated " << std::to_string(index) << " of " << (nBinariesRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nBinariesRequested)) << " binaries requested");

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                        // shouldn't be...
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN(ERROR::BINARY_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));        // show warning
        }
    }

    // close BSE logfiles
    // don't check result here - let log system handle it
    (void)LOGGING->CloseAllStandardFiles();

    double cpuSeconds = (clock() - clockStart) / (double) CLOCKS_PER_SEC;                                           // stop CPU timer and calculate seconds

    auto wallEnd = std::chrono::system_clock::now();                                                                // stop wall timer
    std::time_t timeEnd = std::chrono::system_clock::to_time_t(wallEnd);                                            // get end time and date

    SAY("\nEnd generating binaries at " << std::ctime(&timeEnd));
    SAY("Clock time = " << cpuSeconds << " CPU seconds");


    std::chrono::duration<double> wallSeconds = wallEnd - wallStart;                                                // elapsed seconds

    int wallHH = (int)(wallSeconds.count() / 3600.0);                                                               // hours
    int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                                   // minutes
    int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));                  // seconds

    SAY("Wall time  = " << std::setfill('0') << std::setw(4) << wallHH << ":" << 
                           std::setfill('0') << std::setw(2) << wallMM << ":" << 
                           std::setfill('0') << std::setw(2) << wallSS << " (hhhh:mm:ss)");                         // Include 0 buffer 

    return std::make_tuple(nBinariesRequested, index);
}


/*
 * COMPAS main program
 *
 * Does some housekeeping:
 *
 * - enables floating-point exception handling if required
 * - starts the Options service (program options)
 * - starts the Log service (for logging and debugging)
 * - starts the Rand service (random number generator)
 *
 * Then evolves either a single or binary star
 *
 */
int main(int argc, char * argv[]) {

    PROGRAM_STATUS programStatus = PROGRAM_STATUS::CONTINUE;                                        // status - initially ok

    RAND->Initialise();                                                                             // initialise the random number service
    RAND->Seed(0l);                                                                                 // set seed to 0 - ensures repeatable results

    bool ok = OPTIONS->Initialise(argc, argv);                                                      // get the program options from the commandline
    if (!ok) {                                                                                      // have commandline options ok?
        programStatus = PROGRAM_STATUS::ERROR_IN_COMMAND_LINE;                                      // no - set status
    }
    else {                                                                                          // yes - have commandline options
        if (OPTIONS->RequestedHelp()) {                                                             // user requested help?
            (void)utils::SplashScreen();                                                            // yes - show splash screen
            OPTIONS->ShowHelp();                                                                    // show help
            programStatus = PROGRAM_STATUS::SUCCESS;                                                // don't evolve anything
        }
        else if (OPTIONS->RequestedVersion()) {                                                     // user requested version?
            (void)utils::SplashScreen();                                                            // yes - show splash screen
            programStatus = PROGRAM_STATUS::SUCCESS;                                                // don't evolve anything
        }
        else if (!OPTIONS->YAMLfilename().empty()) {                                                // user requested YAML file creation?
            (void)utils::SplashScreen();                                                            // yes - show splash screen
            yaml::MakeYAMLfile(OPTIONS->YAMLfilename(), OPTIONS->YAMLtemplate());                   // create YAML file
            programStatus = PROGRAM_STATUS::SUCCESS;                                                // don't evolve anything
        }

        if (programStatus == PROGRAM_STATUS::CONTINUE) {

            if (OPTIONS->FPErrorMode() != FP_ERROR_MODE::OFF) {                                     // floating-point error handling mode on/debug?
                                                                                                    // yes
                struct sigaction sigAct;
                memset(&sigAct, 0, sizeof(sigAct));
                sigAct.sa_handler = SIGFPEhandler;                                                  // set the SIGFPE signal handler
                sigAct.sa_flags   = SA_NODEFER;                                                     // don't defer further signals
                sigaction(SIGFPE, &sigAct, NULL);                                                   // enable the signal handler

                feclearexcept(FE_ALL_EXCEPT);                                                       // clear all FE traps
                feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);             // enable FE traps (don't trap FE_INEXACT - would trap on almost all FP operations...)
            }

            InitialiseProfiling;                                                                    // initialise profiling functionality

            int objectsRequested = 0;                                                               // for logging
            int objectsCreated   = 0;                                                               // for logging

            // start the logging service
            LOGGING->Start(OPTIONS->OutputPathString(),                                             // location of logfiles
                           OPTIONS->OutputContainerName(),                                          // directory to be created for logfiles
                           OPTIONS->LogfileNamePrefix(),                                            // prefix for logfile names
                           OPTIONS->LogLevel(),                                                     // log level - determines (in part) what is written to log file
                           OPTIONS->LogClasses(),                                                   // log classes - determines (in part) what is written to log file
                           OPTIONS->DebugLevel(),                                                   // debug level - determines (in part) what debug information is displayed
                           OPTIONS->DebugClasses(),                                                 // debug classes - determines (in part) what debug information is displayed
                           OPTIONS->DebugToFile(),                                                  // should debug statements also be written to logfile?
                           OPTIONS->ErrorsToFile(),                                                 // should error messages also be written to logfile?
                           OPTIONS->LogfileType());                                                 // log file type

            (void)utils::SplashScreen();                                                            // announce ourselves

            if (!LOGGING->Enabled()) programStatus = PROGRAM_STATUS::LOGGING_FAILED;                // logging failed to start
            else {   

                if (!OPTIONS->GridFilename().empty()) {                                             // have grid filename?
                    ERROR error = OPTIONS->OpenGridFile(OPTIONS->GridFilename());                   // yes - open grid file
                    if (error != ERROR::NONE) {                                                     // open ok?
                        SHOW_ERROR(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'"); // no - show error
                        programStatus = PROGRAM_STATUS::STOPPED;                                    // set status
                    }
                }

                if (programStatus == PROGRAM_STATUS::CONTINUE) {                                    // all ok?

                    if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE) {                          // SSE?
                        std::tie(objectsRequested, objectsCreated) = EvolveSingleStars();           // yes - evolve single stars
                    }
                    else {                                                                          // no - BSE
                        std::tie(objectsRequested, objectsCreated) = EvolveBinaryStars();           // evolve binary stars
                    }

                    if (!OPTIONS->GridFilename().empty()) {                                         // have grid filename?
                        OPTIONS->CloseGridFile();                                                   // yes - close it if it's open
                    }

                    programStatus = PROGRAM_STATUS::SUCCESS;                                        // set program status, and...
                }                
            }

            LOGGING->Stop(std::make_tuple(objectsRequested, objectsCreated));                       // stop the logging service if necessary and cleanup

            ReportProfiling;                                                                        // report profiling statistics
        }
    }

    RAND->Free();                                                                                   // release gsl dynamically allocated memory

    return static_cast<int>(programStatus);                                                         // we're done
}

