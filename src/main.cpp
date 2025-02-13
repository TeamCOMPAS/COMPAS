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
#include <cstring>
#include <execinfo.h>

#include <fenv.h>

#include "constants.h"
#include "typedefs.h"

#include "Errors.h"
#include "profiling.h"
#include "utils.h"
#include "yaml.h"
#include "vector3d.h"
#include "Options.h"
#include "Rand.h"
#include "Log.h"

#include "Star.h"
#include "BinaryStar.h"


// disable floating-point error traps for MacOS for now - until I can work out how
// to do it properly for both INTEL and ARM architectures (Linux works ok).
#ifdef __APPLE__
    int feclearexcept(int p_Exceptions) { return 0; };  // no-op for macos for now
    int feenableexcept(int p_Exceptions) { return 0;};  // no-op for macos for now
    int fedisableexcept(int p_Exceptions) { return 0;}; // no-op for macos for now
#endif


OBJECT_ID globalObjectId = 1;                                   // used to uniquely identify objects - used primarily for error printing
OBJECT_ID m_ObjectId     = 0;                                   // object id for main - always 0


/*
 * SIGFPE signal handler
 * 
 * Only handles SIGFPE.
 * 
 * SIGFPE is a system generated signal raised when floating-point errors occur and
 * floating-point traps are enabled.
 * 
 * We implement two active modes for floating-point error handling: ON and DEBUG (the
 * third mode is OFF):
 * 
 *     OFF  : we just ignore the signal.
 * 
 *     ON   : we throw a runtime_error exception, with the what string set to "FPE",
 *            which is the caught by catch blocks throughout the code and handled there.
 *          
 *     DEBUG: we construct and display a stack trace, then halt the program by calling
 *            std::exit(1)
 * 
 * 
 * void SIGFPEhandler(int p_Sig)
 * 
 * @param   [IN]        p_Sig                   The signal intercepted
 *  
 */
void SIGFPEhandler(int p_Sig) {
    if (p_Sig == SIGFPE) {                                                              // SIGFPE?
                                                                                        // yes  
        FP_ERROR_MODE fpErrorMode = OPTIONS->FPErrorMode();                             // mode = ON or DEBUG?
        if (fpErrorMode == FP_ERROR_MODE::ON) throw std::runtime_error("FPE");          // mode = ON
        else if (fpErrorMode == FP_ERROR_MODE::DEBUG) {                                 // mode = DEBUG?
            std::cerr << "\nFloating-point error encountered: program terminated\n";    // announce error
            utils::ShowStackTrace();                                                    // construct and show stack trace
        }        
    }
    else {                                                                              // not SIGFPE - can't jump back from here, so...
        std::cerr << "\nUnexpected signal encountered: program terminated\n";           // announce error
        utils::ShowStackTrace();                                                        // construct and show stack trace
    }
    std::exit(1);                                                                       // catch-all in case we're here when we shouldn't be
}


class Star;
class BinaryStar;


// The following global variables support the BSE Switch Log file
// Ideally, rather than be declared as globals, they would be in maybe the 
// LOGGING service singleton, but the Log class knows nothing about the 
// BinaryStar class...
// (maybe we could put them in the new CONSTANTS service singleton if we 
// implement it)

BinaryStar* evolvingBinaryStar      = NULL;             // pointer to the currently evolving Binary Star
bool        evolvingBinaryStarValid = false;            // flag to indicate whether the evolving BinaryStar pointer is valid

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
 * SIGUSR2 signal handler
 * 
 * Only handles SIGUSR2.
 * 
 * SIGUSR2 is a user generated signal - the system should not generate this signal,
 * though it is possible to send the signal to a process via the Un*x kill command,
 * or some other user-developed program that sends signals.  This code does some 
 * rudimentary sanity checks, but it is possible that sending a SIGUSR2 signal to
 * running COMPAS process via the Un*x kill command, or otherwise, might cause the
 * program to dump a stack trace and halt - c'est la vie.
 * 
 * This signal handler, upon receipt of a SIGUSR2 signal, will construct and display
 * a stack trace, then halt the program by calling std::exit(1).  This is useful for
 * debugging code and the code path to a particular location in the code is not
 * obvious.  Inserting a
 * 
 *     raise(SIGUSR2);
 * 
 * statement at the location will, when that statement is executed, cause a SIGUSR2
 * signal to be raised, and that will invoke this signal handler which will then
 * construct and display a stack trace, and halt the program.
 *
 * void SIGUSR2handler(int p_Sig)
 * 
 * @param   [IN]        p_Sig                   The signal intercepted
 *  
 */
void SIGUSR2handler(int p_Sig) {
    if (p_Sig == SIGUSR2) {                                      // SIGUSR2?  Just silently ignore anything else...
        std::cerr << "\nSIGUSR2 raised: program terminated\n";
        utils::ShowStackTrace(); 
        std::exit(1);
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

    auto    wallStart  = std::chrono::system_clock::now();                                                                  // start wall timer
    clock_t clockStart = clock();                                                                                           // start CPU timer

    std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
    SAY("Start generating stars at " << std::ctime(&timeStart));

    // generate and evolve stars

    Star*  star      = nullptr;
    bool   usingGrid = !OPTIONS->GridFilename().empty();                                                                    // using grid file?
    size_t index     = 0;                                                                                                   // which star

    // The options specified by the user at the commandline are set to their initial values.
    // OPTIONS->AdvanceCmdLineOptionValues(), called at the end of the loop, advances the
    // options specified by the user at the commandline to their next variation (if necessary,
    // based on any ranges and/or sets specified by the user).

    while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                 // while all ok

        try {

            // generate and evolve stars

            int gridLineVariation   = 0;                                                                                    // grid line variation number
            bool doneGridFile       = false;                                                                                // flags we're done with the grid file (for this commandline variation)
            bool processingGridLine = false;                                                                                // processing a gridfile line?
            while (!doneGridFile && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                        // for each star to be evolved

                if (OPTIONS->FPErrorMode() != FP_ERROR_MODE::OFF) {                                                         // floating-point error handling mode on/debug?
                    (void)feclearexcept(FE_ALL_EXCEPT);                                                                     // yes - clear all FE traps
                    (void)feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);                           // enable FE traps (don't trap FE_INEXACT - would trap on almost all FP operations...)
                }
                else (void)fedisableexcept(FE_ALL_EXCEPT);                                                                  // disable all FE traps
                
                bool doneGridLine = false;                                                                                  // initially
                if (usingGrid) {                                                                                            // using grid file?
                    gridLineVariation = 0;                                                                                  // yes - first variation of this grid line
                    int gridResult = OPTIONS->ApplyNextGridLine();                                                          // set options according to specified values in grid file              
                    ERROR error;
                    switch (gridResult) {                                                                                   // handle result of grid file read
                        case -1:                                                                                            // error - unexpected end of grid file
                        case -2:                                                                                            // error - read error for grid file
                        case -3: {                                                                                          // error - invalid value in grid file
                            switch (gridResult) {
                                case -1: error = ERROR::UNEXPECTED_END_OF_FILE; break;  
                                case -2: error = ERROR::FILE_READ_ERROR; break;  
                                case -3: error = ERROR::INVALID_VALUE_IN_FILE; break;
                                default: error = ERROR::FILE_READ_ERROR; break;
                            }

                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                      // set status
                            THROW_ERROR_STATIC(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'");             // throw error
                        } break;

                        case  0: {                                                                                          // end of file
                            doneGridLine = true;                                                                            // flag we're done with this grid line
                            doneGridFile = true;                                                                            // flag we're done with the grid file
                            error = OPTIONS->RewindGridFile();                                                              // ready for next commandline options variation
                            if (error != ERROR::NONE) {                                                                     // rewind ok?
                                evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // no - should never happen here - should be picked up at file open
                                THROW_ERROR_STATIC(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'");         // throw error
                            }
                            } break;

                        case  1:                                                                                            // grid line read
                            processingGridLine = true;                                                                      // not done yet...
                            break;

                        default:                                                                                            // unexpected error
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                      // set status
                            THROW_ERROR_STATIC(ERROR::ERROR, "Accessing grid file '" + OPTIONS->GridFilename() + "'");      // throw error
                            break;
                    }
                }
                else {                                                                                                      // no, not using a grid file
                    doneGridFile = true;                                                                                    // flag we're done with the grid file
                }


                // The options specified by the user in the grid file line are set to their initial values.
                // (by OPTIONS->ApplyNextGridLine()).
                // OPTIONS->AdvanceGridLineOptionValues(), called at the end of the loop, advances the
                // options specified by the user in the grid file line to their next variation (if necessary,
                // based on any ranges and/or sets specified by the user).
                // Note that `doneGridLine` may be a proxy for the command line here (it is set FALSE even if
                // there is no grid file - so the loop executes once to process the command line)

                while (!doneGridLine && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                    // while all ok and not done

                    // Single stars (in SSE) are provided with a random seed that is used to seed the random 
                    // number generator.  The random number generator is re-seeded for each star.  Here we 
                    // generate the seed for the star being evolved - by this point we have picked up the 
                    // option value from either the commandline or the grid file.
                    //
                    // There are three scenarios:
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
                    // and options (if no ranges or sets were specified on the grid line then no offset is added
                    // (i.e. the random seed specified is used as is)).  Note that in this scenario it is the 
                    // user's responsibility to ensure that there is no duplication of seeds.

                    std::string       errorStr;                                                                             // error string
                    unsigned long int randomSeed = 0l;                                                                      // random seed
                    OPTIONS_ORIGIN    optsOrigin = processingGridLine ? OPTIONS_ORIGIN::GRIDFILE : OPTIONS_ORIGIN::CMDLINE; // indicate which set of program options we're using
                    if (OPTIONS->FixedRandomSeedGridLine()) {                                                               // user specified a random seed in the grid file for this binary?
                                                                                                                            // yes - use it (indexed)
                        randomSeed = OPTIONS->RandomSeedGridLine() + (unsigned long int)gridLineVariation;                  // random seed               
                        errorStr   = OPTIONS->SetRandomSeed(randomSeed, optsOrigin);                                        // set it
                        if (!errorStr.empty()) {                                                                            // ok?
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                      // no - set status
                            THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS, errorStr);                         // throw error
                        }
                    }
                    else if (OPTIONS->FixedRandomSeedCmdLine()) {                                                           // no - user specified a random seed on the commandline?
                                                                                                                            // yes - use it (indexed)
                        randomSeed = OPTIONS->RandomSeedCmdLine() + (unsigned long int)index;                               // random seed               
                        errorStr   = OPTIONS->SetRandomSeed(randomSeed, optsOrigin);                                        // set it
                        if (!errorStr.empty()) {                                                                            // ok?
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                      // no - set status
                            THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS, errorStr);                          // throw error
                        }
                    }
                    else {                                                                                                  // no
                                                                                                                            // use default seed (based on system time) + id (index)
                        randomSeed = RAND->DefaultSeed() + (unsigned long int)index;                                        // random seed               
                        errorStr   = OPTIONS->SetRandomSeed(randomSeed, optsOrigin);                                        // set it
                        if (!errorStr.empty()) {                                                                            // ok?
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                      // no - set status
                            THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS, errorStr);                          // throw error
                        }
                    }

                    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                    // ok?
                                                                                                                            // yes - continue
                        randomSeed = RAND->CurrentSeed();                                                                   // current random seed - to pass to star object

                        // the initial mass of the star is supplied - this is to allow binary stars to initialise
                        // the masses of their constituent stars (rather than have the constituent stars sample 
                        // their own mass).  Here we use the mass supplied by the user via the program options or, 
                        // if no mass was supplied by the user, sample the mass from the IMF.

                        double initialMass = OPTIONS->OptionSpecified("initial-mass")                                       // user specified mass?
                                                ? OPTIONS->InitialMass()                                                    // yes, use it
                                                : utils::SampleInitialMass(OPTIONS->InitialMassFunction(),                  // no, sample it
                                                                           OPTIONS->InitialMassFunctionMax(), 
                                                                           OPTIONS->InitialMassFunctionMin(), 
                                                                           OPTIONS->InitialMassFunctionPower());

                        // the metallicity of the star is supplied - this is to allow binary stars to initialise
                        // the metallicity of their constituent stars (rather than have the constituent stars sample 
                        // their own metallicity).  Here we use the mmetallicityass supplied by the user via the program
                        // options or, if no metallicity was supplied by the user, sample the metallicity.

                        double metallicity = OPTIONS->OptionSpecified("metallicity")                                        // user specified metallicity?
                                                ? OPTIONS->Metallicity()                                                    // yes, use it
                                                : utils::SampleMetallicity(OPTIONS->MetallicityDistribution(), 
                                                                           OPTIONS->MetallicityDistributionMax(), 
                                                                           OPTIONS->MetallicityDistributionMin());          // no, sample it

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
                        kickParameters.magnitudeRandomSpecified = OPTIONS->OptionSpecified("kick-magnitude-random");
                        kickParameters.magnitudeRandom          = OPTIONS->KickMagnitudeRandom();
                        kickParameters.magnitudeSpecified       = OPTIONS->OptionSpecified("kick-magnitude");
                        kickParameters.magnitude                = OPTIONS->KickMagnitude();
                       
                        // create the star
                        delete star; star = nullptr;                                                                        // so we don't leak...
                        star = OPTIONS->OptionSpecified("rotational-frequency")                                             // user specified rotational frequency?
                                ? new Star(randomSeed, initialMass, metallicity, kickParameters, OPTIONS->RotationalFrequency() * SECONDS_IN_YEAR) // yes - use it (convert from Hz to cycles per year - see BaseStar::CalculateZAMSAngularFrequency())
                                : new Star(randomSeed, initialMass, metallicity, kickParameters);                           // no - let it be calculated

                        thisStarStatus = EVOLUTION_STATUS::STARTED;

                        thisStarStatus = star->Evolve(index);                                                               // evolve the star

                        // announce the result
                        if (!OPTIONS->Quiet()) {                                                                            // quiet mode?
                            SAY(index                                     <<                                                // announce result of evolving the star
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

                        if (!LOGGING->CloseStandardFile(LOGFILE::SSE_DETAILED_OUTPUT)) {                                    // close SSE detailed output file if necessary
                            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                      // close failed - this will cause problems later - set status
                            THROW_ERROR_STATIC(ERROR::FILE_NOT_CLOSED);                                                     // throw error
                        }

                        ERRORS->Clean();                                                                                    // clean the dynamic error catalog

                        index++;                                                                                            // next...

                        if (usingGrid) {                                                                                    // using grid file?
                            gridLineVariation++;                                                                            // yes - increment grid line variation number
                            int optionsStatus = OPTIONS->AdvanceGridLineOptionValues();                                     // apply next grid file options (ranges/sets)
                            if (optionsStatus < 0) {                                                                        // ok?
                                evolutionStatus = EVOLUTION_STATUS::ERROR;                                                  // no - set status
                                THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                               // throw error
                            }
                            else if (optionsStatus == 0) {                                                                  // end of grid file options variations?
                                doneGridLine = true;                                                                        // yes - we're done
                            }
                        }
                        else doneGridLine = true;                                                                           // not using grid file - done    
                    }
                }
            }
            delete star; star = nullptr;                                                                                    // so we don't leak...
    
            if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                            // ok?
                int optionsStatus = OPTIONS->AdvanceCmdLineOptionValues();                                                  // yes - apply next commandline options (ranges/sets)
                if (optionsStatus < 0) {                                                                                    // ok?
                    evolutionStatus = EVOLUTION_STATUS::ERROR;                                                              // no - set status
                    THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                            // throw error
                }
                else if (optionsStatus == 0) {                                                                              // end of options variations?
                    if (usingGrid || OPTIONS->CommandLineGrid() || (!usingGrid && index >= OPTIONS->nObjectsToEvolve())) {  // created required number of stars?
                        evolutionStatus = EVOLUTION_STATUS::DONE;                                                           // yes - we're done
                    }
                }
            }
        }

        // if we catch an error here it happened outside the evolution of a binary - meaning that there is a problem
        // in recording the results of the evolution of the last binary evolved, or in setting up the next binary to be
        // evolved - either way we should halt the program here because this could result in undefined behaviour and
        // results that can't be trusted. So here we just report the error, report what we got up to (how many binaries
        // were evolved before this happened), and terminate the program.

        catch (const std::runtime_error& e) {                                                                               // catch runtime exceptions
            // anything we catch here should not already have been displayed to the user,
            // so display the error and flag termination (do not rethrow the error)
            if (std::string(e.what()) == "FPE") SHOW_ERROR_STATIC(ERROR::FLOATING_POINT_ERROR)                              // floating-point error
            else                                SHOW_ERROR_STATIC(ERROR::ERROR)                                             // unspecified error
            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                      // set status
        }
        catch (int e) {                                                                                                     // catch errors thrown
            // anything we catch here should already have been displayed to the user,
            // so just flag termination
            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                      // evolution terminated
        }
        catch (...) {                                                                                                       // catchall
            // anything we catch here should not already have been displayed to the user,
            // so display the error and flag termination (do not rethrow the error)
            SHOW_ERROR_STATIC(ERROR::ERROR);                                                                                // unspecified error
            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                      // evolution terminated
        }
    }

    int nStarsRequested = evolutionStatus == EVOLUTION_STATUS::DONE ? index : (usingGrid ? -1 : OPTIONS->nObjectsToEvolve());

    SAY("\nGenerated " << std::to_string(index) << " of " << (nStarsRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nStarsRequested)) << " stars requested");

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                                // shouldn't be
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN_STATIC(ERROR::STELLAR_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));        // show warning
        }
    }

    // close SSE logfiles
    // don't check result here - let log system handle it
    (void)LOGGING->CloseAllStandardFiles();                                                                                 // close any standard log files

    // announce timing stats
    double cpuSeconds = (clock() - clockStart) / (double) CLOCKS_PER_SEC;                                                   // stop CPU timer and calculate seconds

    auto wallEnd = std::chrono::system_clock::now();                                                                        // stop wall timer
    std::time_t timeEnd = std::chrono::system_clock::to_time_t(wallEnd);                                                    // get end time and date

    SAY("\nEnd generating stars at " << std::ctime(&timeEnd));
    SAY("Clock time = " << cpuSeconds << " CPU seconds");

    std::chrono::duration<double> wallSeconds = wallEnd - wallStart;                                                        // elapsed seconds

    int wallHH = (int)(wallSeconds.count() / 3600.0);                                                                       // hours
    int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                                           // minutes
    int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));                          // seconds

    SAY("Wall time  = " << std::setfill('0') << std::setw(2) << wallHH << ":" << 
                           std::setfill('0') << std::setw(2) << wallMM << ":" << 
                           std::setfill('0') << std::setw(2) << wallSS << " (hh:mm:ss)");

    return std::make_tuple(nStarsRequested, index);
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

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    auto wallStart = std::chrono::system_clock::now();                                                                  // start wall timer
    clock_t clockStart = clock();                                                                                       // start CPU timer

    std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
    SAY("Start generating binaries at " << std::ctime(&timeStart));

    BinaryStar* binary    = nullptr;
    bool        usingGrid = !OPTIONS->GridFilename().empty();                                                           // using grid file?
    size_t      index     = 0;                                                                                          // which binary

    signal(SIGUSR1, SIGUSR1handler);                                                                                    // install SIGUSR1 signal handler

    // The options specified by the user at the commandline are set to their initial values.
    // OPTIONS->AdvanceCmdLineOptionValues(), called at the end of the loop, advances the
    // options specified by the user at the commandline to their next variation (if necessary,
    // based on any ranges and/or sets specified by the user).

    while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                             // while all ok and not done

        try {

            // generate and evolve binaries

            int  gridLineVariation  = 0;                                                                                // grid line variation number
            bool doneGridFile       = false;                                                                            // flags we're done with the grid file (for this commandline variation)
            bool processingGridLine = false;                                                                            // processing a gridfile line?
            while (!doneGridFile && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                    // for each binary to be evolved

                if (OPTIONS->FPErrorMode() != FP_ERROR_MODE::OFF) {                                                     // floating-point error handling mode on/debug?
                    (void)feclearexcept(FE_ALL_EXCEPT);                                                                 // yes - clear all FE traps
                    (void)feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);                       // enable FE traps (don't trap FE_INEXACT - would trap on almost all FP operations...)
                }
                else (void)fedisableexcept(FE_ALL_EXCEPT);                                                              // disable all FE traps

                evolvingBinaryStar      = NULL;                                                                         // unset global pointer to evolving binary (for BSE Switch Log)
                evolvingBinaryStarValid = false;                                                                        // indicate that the global pointer is not (yet) valid (for BSE Switch log)

                bool doneGridLine = false;                                                                              // flags we're done with this grid file line (if using a grid file)
                if (usingGrid) {                                                                                        // using grid file?
                    gridLineVariation  = 0;                                                                             // yes - first variation of this grid line
                    int gridResult     = OPTIONS->ApplyNextGridLine();                                                  // set options according to specified values in grid file
                    std::string errStr = "Accessing grid file '" + OPTIONS->GridFilename() + "'";                       // common error string
                    switch (gridResult) {                                                                               // handle result of grid file read
                        case -1:                                                                                        // error - unexpected end of file grid file read
                        case -2: {                                                                                      // error - read error for grid file
                            ERROR error = gridResult == -1 ? ERROR::UNEXPECTED_END_OF_FILE : ERROR::FILE_READ_ERROR;    // set error
                            THROW_ERROR_STATIC(error, errStr);                                                          // throw error
                        } break;
                        case  0: {                                                                                      // end of file
                            doneGridLine = true;                                                                        // flag we're done with this grid line
                            doneGridFile = true;                                                                        // flag we're done with the grid file
                            ERROR error = OPTIONS->RewindGridFile();                                                    // ready for next commandline options variation
                            if (error != ERROR::NONE) {                                                                 // rewind ok?
                                THROW_ERROR_STATIC(error, errStr);                                                      // throw error
                            }
                        } break;
                        case  1:                                                                                        // grid line read
                            processingGridLine = true;                                                                  // not done yet...
                            break;
                        default:                                                                                        // unexpected error
                            THROW_ERROR_STATIC(ERROR::ERROR, errStr);                                                   // throw error
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
                    // and options (if no ranges or sets were specified on the grid line then no offset is added
                    // (i.e. the random seed specified is used as is)).  Note that in this scenario it is the 
                    // user's responsibility to ensure that there is no duplication of seeds.
             
                    unsigned long int thisId = index + gridLineVariation;                                               // set the id for the binary

                    std::string       errorStr;                                                                         // error string
                    unsigned long int randomSeed = 0l;                                                                  // random seed
                    OPTIONS_ORIGIN    optsOrigin = processingGridLine ? OPTIONS_ORIGIN::GRIDFILE : OPTIONS_ORIGIN::CMDLINE; // indicate which set of program options we're using
                    if (OPTIONS->FixedRandomSeedGridLine()) {                                                           // user specified a random seed in the grid file for this binary?
                                                                                                                        // yes - use it (indexed)
                        randomSeed = OPTIONS->RandomSeedGridLine() + (unsigned long int)gridLineVariation;              // random seed               
                        errorStr   = OPTIONS->SetRandomSeed(randomSeed, optsOrigin);                                    // set it
                        if (!errorStr.empty()) {                                                                        // ok?
                            THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS, errorStr);                     // throw error
                        }
                    }
                    else if (OPTIONS->FixedRandomSeedCmdLine()) {                                                       // no - user specified a random seed on the commandline?
                                                                                                                        // yes - use it (indexed)
                        randomSeed = OPTIONS->RandomSeedCmdLine() + (unsigned long int)index + (unsigned long int)gridLineVariation; // random seed               
                        errorStr   = OPTIONS->SetRandomSeed(randomSeed, optsOrigin);                                    // set it
                        if (!errorStr.empty()) {                                                                        // ok?
                            THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS, errorStr);                      // throw error
                        }
                    }
                    else {                                                                                              // no
                                                                                                                        // use default seed (based on system time) + id (index)
                        randomSeed = RAND->DefaultSeed() + (unsigned long int)index + (unsigned long int)gridLineVariation; // random seed               
                        errorStr   = OPTIONS->SetRandomSeed(randomSeed, optsOrigin);                                    // set it
                        if (!errorStr.empty()) {                                                                        // ok?
                            THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS, errorStr);                      // throw error
                        }
                    }

                    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                // ok?
                                                                                                                        // yes - continue
                        randomSeed = RAND->CurrentSeed();                                                               // current random seed - to pass to binary object

                        EVOLUTION_STATUS thisBinaryStatus;                                                              // status of this binary
                        bool haveBinary = true;                                                                         // binary constructed ok?
                        try {
                            delete binary; binary = nullptr;                                                            // so we don't leak
                            binary = new BinaryStar(randomSeed, thisId);                                                // generate binary according to the user options
                        }
                        catch (...) {                                                                                   // catchall
                            // if we catch an error here it happened while creating the binary - in the class constructor
                            // we don't need to halt the program here - just prevent evolution of the binary
                            // any error caught here has already been displayed to the user, so we just catch all
                            haveBinary = false;                                                                         // binary not constructed
                            thisBinaryStatus = EVOLUTION_STATUS::BINARY_ERROR;                                          // error with binary
                                        
                            // announce result
                            if (!OPTIONS->Quiet()) {                                                                    // quiet mode?
                                SAY(thisId << ": " << EVOLUTION_STATUS_LABEL.at(thisBinaryStatus) << ": not evolved");  // no - announce result of (not) evolving the binary
                            }
                        }
                        
                        if (haveBinary) {                                                                               // binary constructed ok?
                            evolvingBinaryStar      = binary;                                                           // yes - set global pointer to evolving binary (for BSE Switch Log)
                            evolvingBinaryStarValid = true;                                                             // indicate that the global pointer is now valid (for BSE Switch Log)

                            thisBinaryStatus = binary->Evolve();                                                        // evolve the binary
                                        
                            // announce result
                            if (!OPTIONS->Quiet()) {                                                                    // quiet mode?
                                                                                                                        // no - announce result of evolving the binary
                                SAY(thisId                                            << ": "    <<
                                    EVOLUTION_STATUS_LABEL.at(thisBinaryStatus)       << ": ("   <<
                                    STELLAR_TYPE_LABEL.at(binary->Star1InitialType()) << " -> "  <<
                                    STELLAR_TYPE_LABEL.at(binary->Star1Type())        << ") + (" <<
                                    STELLAR_TYPE_LABEL.at(binary->Star2InitialType()) << " -> "  <<
                                    STELLAR_TYPE_LABEL.at(binary->Star2Type())        <<  ")"
                                );
                            }

                            if (!LOGGING->CloseStandardFile(LOGFILE::BSE_DETAILED_OUTPUT)) {                            // close detailed output file if necessary
                                THROW_ERROR_STATIC(ERROR::FILE_NOT_CLOSED);                                             // close failed - this will cause problems later - throw error
                            }
                        }

                        ERRORS->Clean();                                                                                // clean the dynamic error catalog

                        if (usingGrid) {                                                                                // using grid file?
                            gridLineVariation++;                                                                        // yes - increment grid line variation number
                            int optionsStatus = OPTIONS->AdvanceGridLineOptionValues();                                 // apply next grid file options (ranges/sets)
                            if (optionsStatus < 0) {                                                                    // ok?
                                THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                           // no - throw error
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
                    evolutionStatus = EVOLUTION_STATUS::ERROR;                                                          // no - set status
                    THROW_ERROR_STATIC(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                        // throw error
                }
                else if (optionsStatus == 0) {                                                                          // end of options variations?
                    if (usingGrid || OPTIONS->CommandLineGrid() || (!usingGrid && index >= OPTIONS->nObjectsToEvolve())) { // created required number of stars?
                        evolutionStatus = EVOLUTION_STATUS::DONE;                                                       // yes - we're done
                    }
                }
            }
        }

        // if we catch an error here it happened outside the evolution of a star - meaning that there is a problem
        // in recording the results of the evolution of the last star evolved, or in setting up the next star to be
        // evolved - either way we should halt the program here because this could result in undefined behaviour and
        // results that can't be trusted. So here we just report the error, report what we got up to (how many stars
        // were evolved before this happened), and terminate the program.

        catch (const std::runtime_error& e) {                                                                           // catch runtime exceptions
            // anything we catch here should not already have been displayed to the user,
            // so display the error and flag termination (do not rethrow the error)
            if (std::string(e.what()) == "FPE") SHOW_ERROR_STATIC(ERROR::FLOATING_POINT_ERROR)                          // floating-point error
            else                                SHOW_ERROR_STATIC(ERROR::ERROR)                                         // unspecified error
            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                  // set status
        }
        catch (int e) {                                                                                                 // catch errors thrown
            // anything we catch here should already have been displayed to the user,
            // so just flag termination
            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                  // evolution terminated
        }
        catch (...) {                                                                                                   // catchall
            // anything we catch here should not already have been displayed to the user,
            // so display the error and flag termination (do not rethrow the error)
            SHOW_ERROR_STATIC(ERROR::ERROR);                                                                            // unspecified error
            evolutionStatus = EVOLUTION_STATUS::ERROR;                                                                  // evolution terminated
        }
    }
    
    int nBinariesRequested = evolutionStatus == EVOLUTION_STATUS::DONE ? index : (usingGrid ? -1 : OPTIONS->nObjectsToEvolve());

    SAY("\nGenerated " << std::to_string(index) << " of " << (nBinariesRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nBinariesRequested)) << " binaries requested");

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                            // shouldn't be...
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN_STATIC(ERROR::BINARY_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));     // show warning
        }
    }

    // close BSE logfiles
    // don't check result here - let log system handle it
    (void)LOGGING->CloseAllStandardFiles();

    double cpuSeconds = (clock() - clockStart) / (double) CLOCKS_PER_SEC;                                               // stop CPU timer and calculate seconds

    auto wallEnd = std::chrono::system_clock::now();                                                                    // stop wall timer
    std::time_t timeEnd = std::chrono::system_clock::to_time_t(wallEnd);                                                // get end time and date

    SAY("\nEnd generating binaries at " << std::ctime(&timeEnd));
    SAY("Clock time = " << cpuSeconds << " CPU seconds");


    std::chrono::duration<double> wallSeconds = wallEnd - wallStart;                                                    // elapsed seconds

    int wallHH = (int)(wallSeconds.count() / 3600.0);                                                                   // hours
    int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                                       // minutes
    int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));                      // seconds

    SAY("Wall time  = " << std::setfill('0') << std::setw(4) << wallHH << ":" << 
                           std::setfill('0') << std::setw(2) << wallMM << ":" << 
                           std::setfill('0') << std::setw(2) << wallSS << " (hhhh:mm:ss)");

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

    signal(SIGUSR2, SIGUSR2handler);                                                                // install SIGUSR2 signal handler

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

                (void)feclearexcept(FE_ALL_EXCEPT);                                                 // clear all FE traps
                (void)feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);       // enable FE traps (don't trap FE_INEXACT - would trap on almost all FP operations...)
            }
            else (void)fedisableexcept(FE_ALL_EXCEPT);                                              // disable all FE traps

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
                        programStatus = PROGRAM_STATUS::STOPPED;                                    // no - set status - stop evolution
                        SHOW_ERROR_STATIC(error, "Accessing grid file '" + OPTIONS->GridFilename() + "'"); // show error (don't throw here - nothing to catch it)
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

