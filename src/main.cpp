//
//  COMPAS main
//
#include <ctime>

#include "constants.h"
#include "typedefs.h"

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


/*
 * Evolve single stars
 *
 *
 * void EvolveSingleStars()
 */
void EvolveSingleStars() {

    clock_t clockStart = clock();

    int    nStepsMass = OPTIONS->SingleStarMassSteps();
    double massMin    = OPTIONS->SingleStarMassMin();
    double massMax    = OPTIONS->SingleStarMassMax();
    double stepMass   = (massMax - massMin) / nStepsMass;

    // Loop over masses to evolve - create and evolve star for each mass

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;
    int stepNum = 0;
    do {

        if (!OPTIONS->Quiet()) SAY("Evolving step " << stepNum << ", mass = " << massMin + (stepNum * stepMass));

        // single stars are provided with a random seed
        // when they are constituent stars of a binary the binary provides the seed
        // here we generate the seed for the single star

        unsigned long int randomSeed = 0l;
        if (OPTIONS->FixedRandomSeed()) {                                                           // user supplied seed for the random number generator?
            randomSeed = RAND->Seed(OPTIONS->RandomSeed() + stepNum);                               // yes - this allows the user to reproduce results for each binary
            if (!OPTIONS->Quiet()) SAY("Using supplied random seed " << randomSeed);
        }
        else {                                                                                      // no
            randomSeed = RAND->Seed(RAND->DefaultSeed() + stepNum);                                 // use default seed (based on system time) + binary id
            if (!OPTIONS->Quiet()) SAY("Using default random seed " << randomSeed);
        }

        Star* star = new Star(randomSeed, massMin + (stepNum * stepMass), OPTIONS->Metallicity());  // create a star with required mass and metallicity, and...
        star->Evolve(stepNum);                                                                      // ... evolve it
        delete star;                                                                                // ... delete it

        if (!LOGGING->CloseStandardFile(LOGFILE::SSE_PARAMETERS)) {                                 // close single star output file
            SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                      // close failed - show warning
            evolutionStatus = EVOLUTION_STATUS::STOPPED;                                            // this will cause problems later - stop evolution
        }

    } while (evolutionStatus == EVOLUTION_STATUS::CONTINUE && ++stepNum <= nStepsMass);             // really should be just <, but <= was how the original code did it...

    if (!OPTIONS->Quiet()) SAY("\nDone. Clock time = " << (clock() - clockStart) / (double) CLOCKS_PER_SEC << " seconds" << "\n");
}


/*
 * Evolve binary stars
 *
 *
 * void EvolveBinaryStars()
 */
void EvolveBinaryStars() {

    clock_t clockStart = clock();

    if (!OPTIONS->Quiet()) {
        SAY("Now generating binaries.")
    }

    AIS ais;                                                                                                            // Adaptive Importance Sampling (AIS)

    if (OPTIONS->AIS_ExploratoryPhase()) ais.PrintExploratorySettings();                                                // print the selected options for AIS Exploratory phase in the beginning of the run
    if (OPTIONS->AIS_RefinementPhase() ) ais.DefineGaussians();                                                         // if we are sampling using AIS (step 2):read in gaussians

    // generate and evolve binaries
    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;
    int index = 0;

    do {
        BinaryStar binary = BinaryStar(ais, index);                                                                     // generate a random binary according to the user options    Floor : added means and covariances as variables in case using Adaptive Importance Sampling

        if (!OPTIONS->Quiet()) SAY("Evolve binary " << index);
        binary.Evolve(index); // JR todo: - fold this into constructor!

        if (OPTIONS->AIS_ExploratoryPhase() && ais.ShouldStopExploratoryPhase(index)) {                                 // AIS says should stop simulation?
            SHOW_WARN(ERROR::BINARY_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::AIS_EXPLORATORY));  // yes - show warning
            break;                                                                                                      // ... and stop
        }

        if (!LOGGING->CloseStandardFile(LOGFILE::BSE_DETAILED_OUTPUT)) {                                                // close detailed output file
            SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                          // close failed - show warning
            evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                // this will cause problems later - stop evolution
        }
    } while (evolutionStatus == EVOLUTION_STATUS::CONTINUE && ++index < OPTIONS->NBinaries());                          // assumption is OPTIONS->nBinaries > 0 JR: todo: should check options...

    // JR: todo: should check/display status here

    // close BSE logfiles
    // don't check result here - let log system handle it

    (void)LOGGING->CloseAllStandardFiles();

    if (!OPTIONS->Quiet()) SAY("\nDone. Clock time = " << (clock() - clockStart) / (double) CLOCKS_PER_SEC << " seconds" << "\n");
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

    COMMANDLINE_STATUS programStatus = OPTIONS->Initialise(argc, argv);                 // Get the program options from the commandline

    if (programStatus == COMMANDLINE_STATUS::CONTINUE) {

        // start the logging service
        LOGGING->Start(OPTIONS->OutputPathString(),                                     // location of logfiles
                       OPTIONS->LogfileNamePrefix(),                                    // prefix for logfile names
                       OPTIONS->LogLevel(),                                             // log level - determines (in part) what is written to log file
                       OPTIONS->LogClasses(),                                           // log classes - determines (in part) what is written to log file
                       OPTIONS->DebugLevel(),                                           // debug level - determines (in part) what debug information is displayed
                       OPTIONS->DebugClasses(),                                         // debug classes - determines (in part) what debug information is displayed
                       OPTIONS->DebugToFile(),                                          // should debug statements also be written to logfile?
                       OPTIONS->ErrorsToFile(),                                         // should error messages also be written to logfile?
                       DELIMITERValue.at(OPTIONS->LogfileDelimiter()));                 // log record field delimiter

        utils::SplashScreen();                                                          // announce ourselves

        if (!LOGGING->Enabled()) programStatus = COMMANDLINE_STATUS::LOGGING_FAILED;    // logging failed to start
        else {

            RAND->Initialise();                                                         // initialise the random number service

            if(OPTIONS->SingleStar()) {                                                 // Single star?
                EvolveSingleStars();                                                    // yes - evolve single stars
            }
            else {                                                                      // no - binary
                EvolveBinaryStars();                                                    // evolve binary stars
            }

            if (!OPTIONS->Quiet()) {                                                    // verbose?
                SAY("\nSuccess!");                                                      // yes - gloat
            }

            RAND->Free();                                                               // release gsl dynamically allocated memory

            LOGGING->Stop();                                                            // stop the logging service

            programStatus = COMMANDLINE_STATUS::SUCCESS;                                // set program status, and...
        }
    }

    return static_cast<int>(programStatus);                                             // we're done
}

