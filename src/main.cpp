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

#include "constants.h"
#include "typedefs.h"

#include "profiling.h"
#include "utils.h"
#include "Options.h"
#include "Rand.h"
#include "Log.h"

#include "AIS.h"
#include "Star.h"
#include "BinaryStar.h"

OBJECT_ID globalObjectId = 1;                                   // used to uniquely identify objects - used primarily for error printing
OBJECT_ID m_ObjectId     = 0;                                   // object id for main - always 0


class AIS;
class Star;
class BinaryStar;

OBJECT_ID    ObjectId()    { return m_ObjectId; }
OBJECT_TYPE  ObjectType()  { return OBJECT_TYPE::MAIN; }
STELLAR_TYPE StellarType() { return STELLAR_TYPE::NONE; }


// The following global variables support the BSE Switch Log file
// Ideally, rather than be declared as globals, they would be in maybe the 
// LOGGING service singleton, but the Log class knows nothing about the 
// BinaryStar class...
// (maybe we could put them in the new CONSTANTS service singleton if we 
// implement it)

BinaryStar* evolvingBinaryStar      = NULL;             // pointer to the currently evolving Binary Star
bool        evolvingBinaryStarValid = false;            // flag to indicate whether the evolvingBinaryStar pointer is valid

/*
 * Signal handler
 * 
 * Only handles SIGUSR1; all other signals are left to the system to handle.
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
 * void sigHandler(int p_Sig)
 * 
 * @param   [IN]        p_Sig                   The signal intercepted
 * 
 */
void sigHandler(int p_Sig) {   
    if (p_Sig == SIGUSR1) {                                         // SIGUSR1?  Just silently ignore anything else...
        if (evolvingBinaryStarValid && OPTIONS->SwitchLog()) {      // yes - do we have a valid binary star, and are we logging switches?
            evolvingBinaryStar->PrintSwitchLog();                   // yes - assume SIGUSR1 is a binary constituent star switching...
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

    auto wallStart = std::chrono::system_clock::now();                                                                  // start wall timer
    clock_t clockStart = clock();                                                                                       // start CPU timer

    std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
    SAY("Start generating stars at " << std::ctime(&timeStart));

    // generate and evolve stars

    bool usingGrid     = !OPTIONS->GridFilename().empty();                                                              // using grid file?
    int  nStars        = usingGrid ? 1 : OPTIONS->nObjectsToEvolve();                                                   // how many stars? (grid file is 1 per record...)
    int  nStarsCreated = 0;                                                                                             // number of stars actually created
    int  index         = 0;                                                                                             // which star

    // The options specified by the user at the commandline are set to their initial values.
    // OPTIONS->AdvanceCmdLineOptionValues(), called at the end of the loop, advances the
    // options specified by the user at the commandline to their next variation (if necessary,
    // based on any ranges and/or sets specified by the user).

    while (evolutionStatus ==EVOLUTION_STATUS::CONTINUE) {                                                              // while all ok

        // generate and evolve stars

        int nStarCount = 0;                                                                                             // count of stars created for this options variation
        Star* star = nullptr;
        while (evolutionStatus == EVOLUTION_STATUS::CONTINUE && nStarCount < nStars) {                                  // for each star to be evolved

            double initialMass = 0.0;                                                                                   // initial mass of star

            bool doneThisGridLine = false;                                                                              // flags we're done with this grid file line (if using a grid file)
            if (usingGrid) {                                                                                            // using grid file?
                int gridResult = OPTIONS->ApplyNextGridLine();                                                          // yes - set options according to specified values in grid file              
                switch (gridResult) {                                                                                   // handle result of grid file read
                    case -1: evolutionStatus = EVOLUTION_STATUS::STOPPED; break;                                        // read error - stop evolution
                    case  0:                                                                                            // end of file
                        doneThisGridLine = true;                                                                        // flag we're done with this grid line
                        nStars--;                                                                                       // decrement count required
                        OPTIONS->RewindGridFile();                                                                      // ready for next commandline options variation
                        break;

                    case  1:                                                                                            // grid record read - not done yet...
                        nStars++;                                                                                       // increment count of stars requested to be evolved
                        initialMass = OPTIONS->InitialMass();                                                           // set initial mass for the star being evolved
                        break;             

                    default: evolutionStatus = EVOLUTION_STATUS::STOPPED; break;                                        // problem - stop evolution
                }
            }
            else {                                                                                                      // no, not using a grid file
                initialMass = OPTIONS->OptionSpecified("initial-mass") == 1                                             // user specified mass?
                                ? OPTIONS->InitialMass()                                                                // yes, use it
                                : utils::SampleInitialMassDistribution(OPTIONS->InitialMassFunction(),                  // no, sample it
                                                                       OPTIONS->InitialMassFunctionMax(), 
                                                                       OPTIONS->InitialMassFunctionMin(), 
                                                                       OPTIONS->InitialMassFunctionPower());
            }


            // The options specified by the user in the grid file line are set to their initial values.
            // (by OPTIONS->ApplyNextGridLine()).
            // OPTIONS->AdvanceGridLineOptionValues(), called at the end of the loop, advances the
            // options specified by the user in the grid file line to their next variation (if necessary,
            // based on any ranges and/or sets specified by the user).

            while (!doneThisGridLine && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                // while all ok and not done

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
                // if the user specified a random seed in the grid file for the current str, regardless of
                // whether a random seed was specified on the commandline, the random seed from the grid
                // file is used, and no offset is added (i.e. the random seed specified is used as it).
                // note that in this scenario it is the user's responsibility to ensure that there is no
                // duplication of seeds.

                unsigned long int randomSeed = 0l;
                if (OPTIONS->FixedRandomSeedGridLine()) {                                                               // user specified a random seed in the grid file for this star?
                    randomSeed = OPTIONS->RandomSeedGridLine();                                                         // yes - use it as is
                }
                else if (OPTIONS->FixedRandomSeedCmdLine()) {                                                           // no - user supplied seed for the random number generator?
                    randomSeed = RAND->Seed(OPTIONS->RandomSeedCmdLine() + (long int)index);                            // yes - this allows the user to reproduce results for each star
                }
                else {                                                                                                  // no
                    randomSeed = RAND->Seed(RAND->DefaultSeed() + (long int)index);                                     // use default seed (based on system time) + id (index)
                }

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
                delete star;                                                                                            // so we don't leak...
                star = new Star(randomSeed, initialMass, kickParameters);                                               // create star according to the user-specified options

                // evolve the star
                star->Evolve(index);

                // announce results if required
                if (!OPTIONS->Quiet()) {                                                                                // quiet mode?
                    SAY(index               <<                                                                          // no - announce result of evolving the star
                        ": RandomSeed = "   <<
                        randomSeed          <<
                        ", Initial Mass = " <<
                        initialMass         <<
                        ", Metallicity = "  <<
                        star->Metallicity() <<
                        ", "                <<
                        STELLAR_TYPE_LABEL.at(star->StellarType()));
                }
                nStarsCreated++;                                                                                        // increment the number of stars created

                if (!LOGGING->CloseStandardFile(LOGFILE::SSE_DETAILED_OUTPUT)) {                                        // close SSE detailed output file
                    SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                  // close failed - show warning
                    evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                        // this will cause problems later - stop evolution
                }

                if (!LOGGING->CloseStandardFile(LOGFILE::SSE_SWITCH_LOG)) {                                             // close SSE switch log file if necessary
                    SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                  // close failed - show warning
                    evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                        // this will cause problems later - stop evolution
                }
                ERRORS->Clean();                                                                                        // clean the dynamic error catalog

                nStarCount++;
                index++;                                                                                                // next...

                if (usingGrid) {                                                                                        // using grid file?
                    int optionsStatus = OPTIONS->AdvanceGridLineOptionValues();                                         // apply next grid file options (ranges/sets)
                    if (optionsStatus < 0) {                                                                            // ok?
                        evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                    // no - stop evolution
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                                           // show error
                    }
                    else if (optionsStatus == 0) {                                                                      // end of grid file options variations?
                        doneThisGridLine = true;                                                                        // yes - we're done
                    }
                }
                else doneThisGridLine = true;                                                                           // not using grid file - done    

            }
        }
        delete star;
    
        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                            // ok?
            int optionsStatus = OPTIONS->AdvanceCmdLineOptionValues();                                                  // yes - apply next commandline options (ranges/sets)
            if (optionsStatus < 0) {                                                                                    // ok?
                evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                            // no - stop evolution
                SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                                    // show error
            }
            else if (optionsStatus == 0) {                                                                              // end of options variations?
                evolutionStatus = EVOLUTION_STATUS::DONE;                                                               // yes - we're done
            }
        }
    }

    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE && index >= nStars) evolutionStatus = EVOLUTION_STATUS::DONE;     // set done

    int nStarsRequested = !usingGrid ? OPTIONS->nObjectsToEvolve() : (evolutionStatus == EVOLUTION_STATUS::DONE ? nStarsCreated : -1);

    SAY("\nGenerated " << std::to_string(nStarsCreated) << " of " << (nStarsRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nStarsRequested)) << " stars requested");

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                            // shouldn't be
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN(ERROR::STELLAR_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));           // show warning
        }
    }

    // close SSE logfiles
    // don't check result here - let log system handle it
    (void)LOGGING->CloseAllStandardFiles();                                                                             // close any standard log files

    // announce timing stats
    double cpuSeconds = (clock() - clockStart) / (double) CLOCKS_PER_SEC;                                               // stop CPU timer and calculate seconds

    auto wallEnd = std::chrono::system_clock::now();                                                                    // stop wall timer
    std::time_t timeEnd = std::chrono::system_clock::to_time_t(wallEnd);                                                // get end time and date

    SAY("\nEnd generating stars at " << std::ctime(&timeEnd));
    SAY("Clock time = " << cpuSeconds << " CPU seconds");

    std::chrono::duration<double> wallSeconds = wallEnd - wallStart;                                                    // elapsed seconds

    int wallHH = (int)(wallSeconds.count() / 3600.0);                                                                   // hours
    int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                                       // minutes
    int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));                      // seconds

    SAY("Wall time  = " << wallHH << ":" << wallMM << ":" << wallSS << " (hh:mm:ss)");

    return  std::make_tuple(nStarsRequested, nStarsCreated);
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

    signal(SIGUSR1, sigHandler);                                                                                        // install signal handler

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    auto wallStart = std::chrono::system_clock::now();                                                                  // start wall timer
    clock_t clockStart = clock();                                                                                       // start CPU timer

    std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
    SAY("Start generating binaries at " << std::ctime(&timeStart));

    AIS ais;                                                                                                            // Adaptive Importance Sampling (AIS)

    if (OPTIONS->AIS_ExploratoryPhase()) ais.PrintExploratorySettings();                                                // print the selected options for AIS Exploratory phase in the beginning of the run
    if (OPTIONS->AIS_RefinementPhase() ) ais.DefineGaussians();                                                         // if we are sampling using AIS (step 2):read in gaussians

    bool usingGrid        = !OPTIONS->GridFilename().empty();                                                           // using grid file?
    int  nBinaries        = usingGrid ? 1 : OPTIONS->nObjectsToEvolve();                                                // how many binaries? (grid file is 1 per record...)
    int  nBinariesCreated = 0;                                                                                          // number of binaries actually created
    int  index            = 0;                                                                                          // which binary

    // The options specified by the user at the commandline are set to their initial values.
    // OPTIONS->AdvanceCmdLineOptionValues(), called at the end of the loop, advances the
    // options specified by the user at the commandline to their next variation (if necessary,
    // based on any ranges and/or sets specified by the user).

    while (evolutionStatus ==EVOLUTION_STATUS::CONTINUE) {                                                              // while all ok and not done

        // generate and evolve binaries

        int nBinaryCount = 0;                                                                                           // count of binaries created for this commandline options variation
        BinaryStar *binary = nullptr;
        while (evolutionStatus == EVOLUTION_STATUS::CONTINUE && nBinaryCount < nBinaries) {                             // for each binary to be evolved

            evolvingBinaryStar      = NULL;                                                                             // unset global pointer to evolving binary (for BSE Switch Log)
            evolvingBinaryStarValid = false;                                                                            // indicate that the global pointer is not (yet) valid (for BSE Switch log)

            bool doneThisGridLine = false;                                                                              // flags we're done with this grid file line (if using a grid file)
            if (usingGrid) {                                                                                            // using grid file?
                int gridResult = OPTIONS->ApplyNextGridLine();                                                          // yes - set options according to specified values in grid file              
                switch (gridResult) {                                                                                   // handle result of grid file read
                    case -1: evolutionStatus = EVOLUTION_STATUS::STOPPED; break;                                        // read error - stop evolution
                    case  0:                                                                                            // end of file
                        doneThisGridLine = true;                                                                        // flag we're done with this grid line
                        nBinaries--;                                                                                    // decrement count required
                        OPTIONS->RewindGridFile();                                                                      // ready for next commandline options variation
                        break;

                    case  1:                                                                                            // grid record read - not done yet...
                        nBinaries++;                                                                                    // increment count of binaries requested to be evolved
                        break;             

                    default: evolutionStatus = EVOLUTION_STATUS::STOPPED; break;                                        // problem - stop evolution
                }
            }

            // The options specified by the user in the grid file line are set to their initial values.
            // (by OPTIONS->ApplyNextGridLine()).
            // OPTIONS->AdvanceGridLineOptionValues(), called at the end of the loop, advances the
            // options specified by the user in the grid file line to their next variation (if necessary,
            // based on any ranges and/or sets specified by the user).

            while (!doneThisGridLine && evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                // while all ok and not done

                // we only need to pass the AIS structure and the index number to the binary - we let the 
                // BinaryStar class do the work wrt setting the parameters for each of the constituent stars
                // (The index is really only needed for legacy comparison, so can probably be removed at any time)
                //
                // Note: the AIS structure will probably go away when Stroopwafel is completeley moved to  outside COMPAS.

                delete binary;
                binary = new BinaryStar(ais, (long int)index);                                                          // generate binary according to the user options
std::cout << "JRPRINT EVOLVING BINARY, metallicity = " << OPTIONS->Metallicity() << 
                                      ", WR factor = " << OPTIONS->WolfRayetFactor() << 
                                      ", LBV factor = " << OPTIONS->LuminousBlueVariableFactor() << 
                                      ", Common envelpe alpha = " << OPTIONS->CommonEnvelopeAlpha() << 
                                      ", Mass1 = " << OPTIONS->InitialMass1() << 
                                      ", Mass2 = " << OPTIONS->InitialMass2() << "\n";

                evolvingBinaryStar      = binary;                                                                       // set global pointer to evolving binary (for BSE Switch Log)
                evolvingBinaryStarValid = true;                                                                         // indicate that the global pointer is now valid (for BSE Switch Log)

                EVOLUTION_STATUS binaryStatus = binary->Evolve();                                                       // evolve the binary

                // announce result of evolving the binary
                if (!OPTIONS->Quiet()) {                                                                                // quiet mode?
                                                                                                                        // no - announce result of evolving the binary
                    if (OPTIONS->CHE_Option() == CHE_OPTION::NONE) {                                                    // CHE enabled?
                        SAY(index                                      << ": "  <<                                      // no - CHE not enabled - don't need initial stellar type
                            EVOLUTION_STATUS_LABEL.at(binaryStatus)    << ": "  <<
                            STELLAR_TYPE_LABEL.at(binary->Star1Type()) << " + " <<
                            STELLAR_TYPE_LABEL.at(binary->Star2Type())
                        );
                    }
                    else {                                                                                              // CHE enabled - show initial stellar type
                        SAY(index                                             << ": "    <<
                            EVOLUTION_STATUS_LABEL.at(binaryStatus)           << ": ("   <<
                            STELLAR_TYPE_LABEL.at(binary->Star1InitialType()) << " -> "  <<
                            STELLAR_TYPE_LABEL.at(binary->Star1Type())        << ") + (" <<
                            STELLAR_TYPE_LABEL.at(binary->Star2InitialType()) << " -> "  <<
                            STELLAR_TYPE_LABEL.at(binary->Star2Type())        <<  ")"
                        );
                    }
                }

                nBinariesCreated++;                                                                                     // increment the number of binaries created

                if (OPTIONS->AIS_ExploratoryPhase() && ais.ShouldStopExploratoryPhase(index)) {                         // AIS says should stop simulation?
                    evolutionStatus = EVOLUTION_STATUS::AIS_EXPLORATORY;                                                // ... and stop
                }

                if (!LOGGING->CloseStandardFile(LOGFILE::BSE_DETAILED_OUTPUT)) {                                        // close detailed output file if necessary
                    SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                  // close failed - show warning
                    evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                        // this will cause problems later - stop evolution
                }

                if (!LOGGING->CloseStandardFile(LOGFILE::BSE_SWITCH_LOG)) {                                             // close BSE switch log file if necessary
                    SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                  // close failed - show warning
                    evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                        // this will cause problems later - stop evolution
                }

                ERRORS->Clean();                                                                                        // clean the dynamic error catalog

                nBinaryCount++;
                index++;                                                                                                // next...

                if (usingGrid) {                                                                                        // using grid file?
                    int optionsStatus = OPTIONS->AdvanceGridLineOptionValues();                                         // apply next grid file options (ranges/sets)
                    if (optionsStatus < 0) {                                                                            // ok?
                        evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                    // no - stop evolution
                        SHOW_ERROR(ERROR::ERROR_PROCESSING_GRIDLINE_OPTIONS);                                           // show error
                    }
                    else if (optionsStatus == 0) {                                                                      // end of grid file options variations?
                        doneThisGridLine = true;                                                                        // yes - we're done
                    }
                }
                else doneThisGridLine = true;                                                                           // not using grid file - done    
            }
        }
        delete binary;
    
        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                            // ok?
            int optionsStatus = OPTIONS->AdvanceCmdLineOptionValues();                                                  // apply next commandline options (ranges/sets)
            if (optionsStatus < 0) {                                                                                    // ok?
                evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                            // no - stop evolution
                SHOW_ERROR(ERROR::ERROR_PROCESSING_CMDLINE_OPTIONS);                                                    // show error
            }
            else if (optionsStatus == 0) {                                                                              // end of commandline options variations?
                evolutionStatus = EVOLUTION_STATUS::DONE;                                                               // yes - we're done
            }
        }
    }

    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE && index >= nBinaries) evolutionStatus = EVOLUTION_STATUS::DONE;  // set done

    int nBinariesRequested = usingGrid ? (evolutionStatus == EVOLUTION_STATUS::DONE ? nBinariesCreated : -1) : OPTIONS->nObjectsToEvolve();

    SAY("\nGenerated " << std::to_string(nBinariesCreated) << " of " << (nBinariesRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nBinariesRequested)) << " binaries requested");

    if (evolutionStatus == EVOLUTION_STATUS::AIS_EXPLORATORY) {                                                         // AIS said stop?
        SHOW_WARN(ERROR::BINARY_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(evolutionStatus));                        // yes - show warning
        evolutionStatus = EVOLUTION_STATUS::DONE;                                                                       // set done
    }

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                            // shouldn't be...
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN(ERROR::BINARY_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));            // show warning
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

    SAY("Wall time  = " << wallHH << ":" << wallMM << ":" << wallSS << " (hh:mm:ss)");

    return std::make_tuple(nBinariesRequested, nBinariesCreated);
}


/*
 * COMPAS main program
 *
 * Does some housekeeping:
 *
 * - starts the Options service (program options)
 * - starts the Log service (for logging and debugging)
 * - starts the Rand service (random number generator)
 *
 * Then evolves either a single or binary star (single star only at the moment...)
 *
 */
int main(int argc, char * argv[]) {

    PROGRAM_STATUS programStatus = PROGRAM_STATUS::CONTINUE;                                        // status - initially ok

    RAND->Initialise();                                                                             // initialise the random number service

    bool ok = OPTIONS->Initialise(argc, argv);                                                      // get the program options from the commandline
    if (!ok) {                                                                                      // have commandline options ok?
                                                                                                    // no - commandline options not ok
//            SAY(cmdline_options);                                                                 // and help                       JRFIX 
        programStatus = PROGRAM_STATUS::ERROR_IN_COMMAND_LINE;                                      // set status
    }
    else {                                                                                          // yes - have commandline options
        if (OPTIONS->RequestedHelp()) {                                                             // user requested help?
            utils::SplashScreen();                                                                  // yes - show splash screen
//            SAY(cmdline_options);                                                                 // and help                       JRFIX 
            programStatus = PROGRAM_STATUS::SUCCESS;                                                // don't evolve anything
        }
        else if (OPTIONS->RequestedVersion()) {                                                     // user requested version?
            SAY("COMPAS v" << VERSION_STRING);                                                      // yes, show version string
            programStatus = PROGRAM_STATUS::SUCCESS;                                                // don't evolve anything
        }
    

        if (programStatus == PROGRAM_STATUS::CONTINUE) {

            InitialiseProfiling;                                                                    // initialise profiling functionality

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
                           DELIMITERValue.at(OPTIONS->LogfileDelimiter()));                         // log record field delimiter

            (void)utils::SplashScreen();                                                            // announce ourselves

            if (!LOGGING->Enabled()) programStatus = PROGRAM_STATUS::LOGGING_FAILED;                // logging failed to start
            else {   
                if (!OPTIONS->GridFilename().empty()) {                                             // have grid filename?
                    ERROR error = OPTIONS->OpenGridFile(OPTIONS->GridFilename());                   // yes - open grid file
                    if (error != ERROR::NONE) {                                                     // open ok?
                        SHOW_ERROR(error, "Opening grid file '" + OPTIONS->GridFilename() + "'");   // no - show error
                        programStatus = PROGRAM_STATUS::STOPPED;                                    // set status
                    }
                }

                int objectsRequested = 0;                                                           // for logging
                int objectsCreated   = 0;                                                           // for logging

                if (programStatus == PROGRAM_STATUS::CONTINUE) {                                    // all ok?

                    if(OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE) {                           // SSE?
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
                
                LOGGING->Stop(std::make_tuple(objectsRequested, objectsCreated));                   // stop the logging service
            }

            ReportProfiling;                                                                        // report profiling statistics
        }
    }

    RAND->Free();                                                                                   // release gsl dynamically allocated memory

    return static_cast<int>(programStatus);                                                         // we're done
}

