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
 *   PUT THESE IN LOG - IF LOGGING NOT ENABLED WE WON'T BE USING THEM ANYWAY...
 * 
 * 
 * 
 * void sigusr1Handler(int p_Sig)
 * 
 * @param   [IN]        p_Sig                   The signal intercepted
 * 
 */
void sigHandler(int p_Sig) {   
    if (p_Sig == SIGUSR1) {                                         // SIGUSR1?  Just silently ignore anything else...
        if (evolvingBinaryStarValid && OPTIONS->BSESwitchLog()) {   // yes - do we have a valid binary star, and are we logging BSE switches?
            evolvingBinaryStar->PrintSwitchLog();                   // yes - assume SIGUSR1 is a binary constituent star switching...
        }
    }
}


/*
 * Open the SSE grid file, read and parse the header record
 *
 *
 * std::tuple<int, std::vector<std::string>, std::string> OpenSSEGridFile(std::ifstream &p_Grid, const std::string p_Filename)
 *
 * @param   [IN/OUT]    p_Grid                  The std::ifstream  to which the file should be opened
 * @param   [IN]        p_Filename              The filename of the Grid file
 * @return                                      tuple containing the line number of the next line to be read and a vector of header strings
 */
std::tuple<int, std::vector<std::string>> OpenSSEGridFile(std::ifstream &p_Grid, const std::string p_Filename) {

    std::vector<std::string> gridHeaders = {};                                                                                  // grid file headers

    int mass               = 0;                                                                                                 // count 'Mass" occurrences
    int metallicity        = 0;                                                                                                 // count 'Metallicity" occurrences

    int kickMagnitudeRandom = 0;                                                                                                 // count 'Kick_Magnitude_Random' occurrences
    int kickMagnitude       = 0;                                                                                                 // count 'Kick_Magnitude" occurrences

    int unknown            = 0;                                                                                                 // count unknown occurrences

    int tokenCount         = 0;                                                                                                 // token count - to check if header should be present
    int currentPos         = p_Grid.tellg();                                                                                    // current file position - to rewind if no header

    int lineNo = 1;                                                                                                             // line number - for error messages

    if (!OPTIONS->GridFilename().empty()) {                                                                                     // have grid filename?

        p_Grid.open(OPTIONS->GridFilename());                                                                                   // yes - open the file
        if (p_Grid.fail()) {                                                                                                    // open ok?
            SAY(ERR_MSG(ERROR::FILE_OPEN_ERROR) << " " << OPTIONS->GridFilename());                                             // no - show error
        }
        else {                                                                                                                  // file open ok

            bool readFailed  = false;                                                                                           // record read ok?
            bool emptyRecord = true;                                                                                            // record empty?
            std::string record;
            while (emptyRecord && !readFailed) {                                                                                // get the header - first non-empty record

                currentPos = p_Grid.tellg();                                                                                    // current file position
                if (std::getline(p_Grid, record)) {                                                                             // read next record
                    readFailed = false;                                                                                         // ok

                    size_t hashPos = record.find("#");                                                                          // find first occurrence of "#"
                    if (hashPos != std::string::npos) record.erase(hashPos, record.size() - hashPos);                           // if "#" found, prune it and everything after it (ignore comments)

                    while (record.size() > 0 && (record[0] == ' ' || record[0] == '\t')) record.erase(0, 1);                    // strip any leading ' ' or TAB ('\t') characters from the record

                    while (record.size() > 0               &&
                          (record[record.size()-1] == '\r' ||
                           record[record.size()-1] == '\n' ||
                           record[record.size()-1] == '\t' ||
                           record[record.size()-1] == ' '  )) record.erase(record.size()-1, 1);                                 // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the record

                    if (record.empty()) {
                        lineNo++;                                                                                               // increment line number
                        continue;                                                                                               // skip empty record
                    }

                    emptyRecord = false;                                                                                        // have non-empty record

                    record = utils::ToUpper(record);                                                                            // upshift

                    std::stringstream recordSS(record);                                                                         // open stream on record
                    std::string token;
                    while (std::getline(recordSS, token, ',')) {                                                                // get token from record read

                        tokenCount++;                                                                                           // increment token count - even if it's empty, it's a token

                        while (token.size() > 0 && (token[0] == ' ' || token[0] == '\t')) token.erase(0, 1);                    // strip any leading ' ' or TAB ('\t') characters from the token

                        while (token.size() > 0              &&
                              (token[token.size()-1] == '\r' ||
                               token[token.size()-1] == '\n' ||
                               token[token.size()-1] == '\t' ||
                               token[token.size()-1] == ' ')) token.erase(token.size()-1, 1);                                   // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the token

                        if (token.empty()) {                                                                                    // empty token?
                            SAY(ERR_MSG(ERROR::GRID_FILE_EMPTY_HEADER));                                                        // show error
                            continue;                                                                                           // next token
                        }

                        switch (_(token.c_str())) {                                                                             // which column header?

                            case _("MASS")                : mass++;                 gridHeaders.push_back(token); break;        // Mass

                            case _("METALLICITY")         : metallicity++;          gridHeaders.push_back(token); break;        // Metallicity

                            case _("KICK_MAGNITUDE_RANDOM"): kickMagnitudeRandom++;   gridHeaders.push_back(token); break;        // Kick magnitude random number

                            case _("KICK_MAGNITUDE")       : kickMagnitude++;         gridHeaders.push_back(token); break;        // Kick magnitude

                            default                       : unknown++;                                                          // unknown - deal with this later
                        }
                    }

                    lineNo++;                                                                                                   // increment line number
                }
                else readFailed = true;                                                                                         // read failed - may be EOF
            }
        }
    }

    if (mass != 1 || metallicity > 1 || kickMagnitude > 1 || kickMagnitudeRandom > 1 || unknown > 0) {                            // check we have all the headers we need, and in the right numbers, and no extraneous headers
                                                                                                                                // we don't, but maybe this wasn't a header record
        if (tokenCount > 1 || mass >= 1 || metallicity >= 1 || kickMagnitude >= 1 || kickMagnitudeRandom >= 1) {                  // more than 1 column, or we got some header strings, so should have been a header
            bool error = true;                                                                                                  // error

            if (mass < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Mass")                                             // no 'Mass'
            else if (mass > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Mass");                                     // duplicate 'Mass'

            if (tokenCount > 1 && metallicity < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Metallicity")             // no 'Metallicity'
            else if (metallicity > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Metallicity");                       // duplicate 'Metallicity'

            if (kickMagnitude > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Magnitude");                         // duplicate 'Kick_Magnitude'

            if (kickMagnitudeRandom > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Magnitude_Random");            // duplicate 'Kick_Magnitude_Random'

            if (unknown > 0) SAY(ERR_MSG(ERROR::GRID_FILE_UNKNOWN_HEADER));                                                     // unknown header string

            if (error) {                                                                                                        // problem?
                if (p_Grid.is_open()) p_Grid.close();                                                                           // yes - close the grid file
            }
        }
        else {                                                                                                                  // otherwise leave the file open and ...
            gridHeaders = { "MASS" };                                                                                           // ... assume we have a header for "Mass" and ...
            p_Grid.seekg(currentPos);                                                                                           // ... rewind to the record just read
        }
    }

    return std::make_tuple(lineNo, gridHeaders);
}


/*
 * Read and parse the next record from the SSE grid file
 *
 * Returns the following data values from the Grid file in the following order:
 *
 *     Mass
 *     Metallicity
 *     Kick_Magnitude_Random 
 *     Kick_Magnitude 
 * 
 * Plus boolean flags (see definition of KickParameters in typedefs.h):
 * 
 *     supplied{1,2}          - true if kick values were supplied in the grid file
 *     useMagnitudeRandom{1,2} - true if the user supplied the kick magnitude magnitude random number
 * 
 * Missing values are treated as zero (0.0) - a warning will be issued, and reading of the Grid file continues
 * (A value is considered missing only if there is a header for the column, but no data value in the column)
 * Invalid values are treated as errors - an error will be issued and reading of the Grid file stopped
 * Negative values are treated as errors - an error will be issued and reading of the Grid file stopped
 *
 * If the SSE Grid file contains Mass values only, Metallicity will be set to the value returned by OPTIONS->Metallicity()
 *
 * std::tuple<int, std::vector<double>> ReadSSEGridRecord(std::ifstream &p_Grid, const std::vector<std::string> p_GridHeaders, const int p_LineNo)
 *
 * @param   [IN/OUT]    p_Grid                  The Grid file (std::ifstream)
 * @param   [IN]        p_GridHeaders           Vector of header strings returned by OpenSSEGridFile()
 * @param   [IN]        p_LineNo                The line number of the next line to be read
 * @return                                      tuple containing the error status, line number of the next line to be read, and a vector of data values
 */
std::tuple<bool, int, SSEGridParameters> ReadSSEGridRecord(std::ifstream &p_Grid, const std::vector<std::string> p_GridHeaders, const int p_LineNo) {

    bool error = false;

    SSEGridParameters gridValues;                                                                                       // grid record values

    gridValues.mass                             = 0.0;
    gridValues.metallicity                      = 0.0;
    gridValues.kickParameters.supplied          = false;
    gridValues.kickParameters.useMagnitudeRandom = false;
    gridValues.kickParameters.magnitudeRandom    = 0.0;
    gridValues.kickParameters.magnitude          = 0.0;

    int lineNo = p_LineNo;                                                                                              // line number - for error messages

    bool emptyRecord = true;
    std::string record;
    while (!error && emptyRecord && std::getline(p_Grid, record)) {                                                     // read the next record of the file

        size_t hashPos = record.find("#");                                                                              // find first occurrence of "#"
        if (hashPos != std::string::npos) record.erase(hashPos, record.size() - hashPos);                               // if "#" found, prune it and everything after it (ignore comments)

        for (size_t pos = 0; pos < record.size(); pos++) if (record[pos] == '\t') record[pos] = ' ';                    // replace tab with space

        while (record.size() > 0 && (record[0] == ' ' || record[0] == '\t')) record.erase(0, 1);                        // strip any leading ' ' or TAB ('\t') characters from the record

        while (record.size() > 0               &&
              (record[record.size()-1] == '\r' ||
               record[record.size()-1] == '\n' ||
               record[record.size()-1] == '\t' ||
               record[record.size()-1] == ' '  )) record.erase(record.size()-1, 1);                                     // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the record

        if (record.empty()) {
            lineNo++;                                                                                                   // increment line number
            continue;                                                                                                   // skip empty record
        }
        emptyRecord = false;                                                                                            // have non-empty record

        std::istringstream recordSS(record);                                                                            // open input stream on record

        size_t column = 0;                                                                                              // start at column 0
        std::string token;                                                                                              // token from record
        while (!error && std::getline(recordSS, token, ',')) {                                                          // get next token from record read

            while (token.size() > 0 && (token[0] == ' ' || token[0] == '\t')) token.erase(0, 1);                        // strip any leading ' ' or TAB ('\t') characters from the token

            while (token.size() > 0              &&
                  (token[token.size()-1] == '\r' ||
                   token[token.size()-1] == '\n' ||
                   token[token.size()-1] == '\t' ||
                   token[token.size()-1] == ' '  )) token.erase(token.size()-1, 1);                                     // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the token

            if (column >= p_GridHeaders.size()) {                                                                       // extraneous value?
                SAY(ERR_MSG(ERROR::GRID_FILE_EXTRA_COLUMN) << " ignored");                                              // yes - show error
            }
            else {                                                                                                      // no - expected

                if (token.empty()) {                                                                                    // empty token?
                    if (utils::Equals(p_GridHeaders[column], "METALLICITY")) {                                          // yes metallicity column?
                        token = std::to_string(OPTIONS->Metallicity());                                                 // yes - use default value from program option
                        SAY(ERR_MSG(ERROR::GRID_FILE_DEFAULT_METALLICITY) << " at line " << lineNo);                    // let the user know we're using default metallicity
                    }
                    else {                                                                                              // no - something else
                        token = "0.0";                                                                                  // use 0.0
                        SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_DATA) << " at line " << lineNo << ": 0.0 used");           // show warning
                    }
                }

                double value;
                std::istringstream tokenSS(token);                                                                      // open input stream on token
                tokenSS >> value;                                                                                       // read double from token
                if (tokenSS.fail()) {                                                                                   // check for valid conversion
                    error = true;                                                                                       // set error flag
                    SAY(ERR_MSG(ERROR::GRID_FILE_INVALID_DATA) << " at line " << lineNo << ": " << token);              // show error
                }

                // don't use utils::Compare() here for bounds/range checks

                if (!error) {
                    std::string columnName = p_GridHeaders[column];                                                     // get column name for value
                    switch (_(columnName.c_str())) {                                                                    // which column?

                        case _("MASS"):                                                                                 // Mass
                            if (value < 0.0) {                                                                          // value < 0?
                                error = true;                                                                           // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token); // show error
                            }
                            else gridValues.mass = value;                                                               // no - proceed
                            break;

                        case _("METALLICITY"):                                                                          // Metallicity
                            if (value < 0.0) {                                                                          // value < 0?
                                error = true;                                                                           // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token); // show error
                            }
                            else gridValues.metallicity = value;                                                        // no - proceed
                            break;

                        case _("KICK_MAGNITUDE_RANDOM"):                                                                 // Kick magnitude random number
                            if (value < 0.0 || value >= 1.0) {                                                          // in the range [0.0, 1.0)? 
                                error = true;                                                                           // no - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_INVALID_DATA) << " at line " << lineNo << ": " << token);  // show error
                            }
                            else {                                                                                      // yes - proceed
                                gridValues.kickParameters.supplied          = true;                                     // kick parameters supplied
                                gridValues.kickParameters.magnitudeRandom    = value;                                    // Kick magnitude random number
                                gridValues.kickParameters.useMagnitudeRandom = true;                                     // use this in preference to actual kick value
                            }
                            break;

                        case _("KICK_MAGNITUDE"):                                                                        // Kick magnitude (magnitude only, so must be +ve - probably technically "speed" rather than "velocity")                        
                            if (value < 0.0) {                                                                          // value < 0?
                                error = true;                                                                           // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token); // show error
                            }
                            else {                                                                                      // no - proceed
                                gridValues.kickParameters.supplied = true;                                              // kick parameters supplied
                                gridValues.kickParameters.magnitude = value;                                             // Kick magnitude
                            }
                            break;

                        default:                                                                                        // unknown - this shouldn't happen
                            SAY(ERR_MSG(ERROR::GRID_FILE_UNKNOWN_HEADER) << ": " << columnName);                        // show error
                    }
                }
            }
            column++;                                                                                                   // next column
        }

        if (!error) {
            if (column < p_GridHeaders.size()) {                                                                        // missing value(s)?
                for (size_t c = column; c < p_GridHeaders.size(); c++) {                                                // yes - use default value if metallicity, otherwise 0.0
                    if (utils::Equals(p_GridHeaders[c], "METALLICITY")) {                                               // metallicity column?
                        SAY(ERR_MSG(ERROR::GRID_FILE_DEFAULT_METALLICITY) << " at line " << lineNo);                    // yes - let the user know we're using default metallicity
                        gridValues.metallicity = OPTIONS->Metallicity();                                                // use program option for metallicity
                    }
                    else {                                                                                              // no - something else
                        SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_DATA) << " at line " << lineNo << ": 0.0 used");           // show warning
                    }
                }
            }

            if (p_GridHeaders.size() == 1 && utils::Equals(p_GridHeaders[0], "Mass")) {                                 // grid file has only Mass values?
                SAY(ERR_MSG(ERROR::GRID_FILE_DEFAULT_METALLICITY));                                                     // yes - let the user know we're using default metallicity
                gridValues.metallicity = OPTIONS->Metallicity();                                                        // use program option for metallicity
            }
        }
        lineNo++;                                                                                                       // increment line number
    }

    return std::make_tuple(error, lineNo, gridValues);
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

    int nStarsCreated = 0;                                                                                                      // number of stars actually created

    auto wallStart = std::chrono::system_clock::now();                                                                          // start wall timer
    clock_t clockStart = clock();                                                                                               // start CPU timer

    if (!OPTIONS->Quiet()) {
        std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
        SAY("Start generating stars at " << std::ctime(&timeStart));
    }

    std::ifstream grid;                                                                                                         // grid file
    std::vector<std::string> gridHeaders;                                                                                       // grid file headers

    int nStars;                                                                                                                 // the number of stars to evolve

    double massIncrement = (OPTIONS->SingleStarMassMax() - OPTIONS->SingleStarMassMin()) / OPTIONS->SingleStarMassSteps();      // use user-specified values if no grid file

    int lineNo;                                                                                                                 // grid file line number - for error messages
    if (!OPTIONS->GridFilename().empty()) {                                                                                     // have grid filename?
        std::tie(lineNo, gridHeaders) = OpenSSEGridFile(grid, OPTIONS->GridFilename());                                         // yes - read grid file headers
        if (!grid.is_open()) evolutionStatus = EVOLUTION_STATUS::ERROR;                                                         // failed to open/read grid file
        nStars = 1;                                                                                                             // ... for each data record in the grid file
    }
    else {
        nStars = OPTIONS->SingleStarMassSteps();                                                                                // user specified
    }

    // Loop over stars to evolve

    evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    long int index  = 0;
    Star* star = nullptr;
    while (evolutionStatus == EVOLUTION_STATUS::CONTINUE && index < nStars) {

        // single stars are provided with a random seed
        // when they are constituent stars of a binary the binary provides the seed
        // here we generate the seed for the single star

        unsigned long int randomSeed = 0l;
        if (OPTIONS->FixedRandomSeed()) {                                                                                       // user supplied seed for the random number generator?
            randomSeed = RAND->Seed(OPTIONS->RandomSeed() + index);                                                             // yes - this allows the user to reproduce results for each binary
        }
        else {                                                                                                                  // no
            randomSeed = RAND->Seed(RAND->DefaultSeed() + index);                                                               // use default seed (based on system time) + binary id
        }

        // create the star

        double initialMass = 0.0;
        if (!OPTIONS->GridFilename().empty()) {                                                                                 // have grid filename?
            if (grid.is_open()) {                                                                                               // yes - grid file open?

                bool error;
                SSEGridParameters gV;                                                                                           // yes - read the next record of grid values
                std::tie(error, lineNo, gV) = ReadSSEGridRecord(grid, gridHeaders, lineNo);                                     // read grid file record
                if (grid.fail() || error) {                                                                                     // read failed?
                    grid.close();                                                                                               // yes - EOF or error - stop reading, close file

                    if (error) evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                     // if error, set error
                    else       evolutionStatus = EVOLUTION_STATUS::DONE;                                                        // otherwise, set done
                }
                else {                                                                                                          // no
                    initialMass = gV.mass;                                                                                      // for status message
                    nStars++;                                                                                                   // not done yet...
                    delete star;
                    star = new Star(randomSeed, gV.mass, gV.metallicity, gV.kickParameters, OPTIONS->LuminousBlueVariableFactor(), OPTIONS->WolfRayetFactor()); // create a star with required attricbutes
                }
            }
            else evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                   // must have been a problem
        }
        else {                                                                                                                  // no grid file - use user-specified values
            initialMass = OPTIONS->SingleStarMassMin() + (index * massIncrement);                                               // for status message

            double fLBV = OPTIONS->SampleLuminousBlueVariableMultiplier()                                                       // LBV winds factor
                            ? RAND->Random(OPTIONS->SampleLuminousBlueVariableMultiplierMin(), OPTIONS->SampleLuminousBlueVariableMultiplierMax())
                            : OPTIONS->LuminousBlueVariableFactor();

            double fWR  = OPTIONS->SampleWolfRayetMultiplier()                                                                  // Wolf Rayet winds factor
                            ? RAND->Random(OPTIONS->SampleWolfRayetMultiplierMin(), OPTIONS->SampleWolfRayetMultiplierMax())
                            : OPTIONS->WolfRayetFactor();

            delete star;
            star = new Star(randomSeed, initialMass, OPTIONS->Metallicity(), {}, fLBV, fWR);                                    // create a star with required attributes
        }

        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                    // still good?

            star->Evolve(index);                                                                                                // yes  - evolve the star

            if (!OPTIONS->Quiet()) {                                                                                            // announce result of evolving the star
                SAY(index               <<
                    ": RandomSeed = "   <<
                    randomSeed          <<
                    ", Initial Mass = " <<
                    initialMass         <<
                    ", Metallicity = "  <<
                    star->Metallicity() <<
                    ", "                <<
                    STELLAR_TYPE_LABEL.at(star->StellarType()));
            }

            nStarsCreated++;                                                                                                    // increment the number of stars created

        }

        if (!LOGGING->CloseStandardFile(LOGFILE::SSE_PARAMETERS)) {                                                             // close single star output file
            SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                                  // close failed - show warning
            evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                        // this will cause problems later - stop evolution
        }

        if (!LOGGING->CloseStandardFile(LOGFILE::SSE_SWITCH_LOG)) {                                                             // close SSE switch log file if necessary
            SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                                  // close failed - show warning
            evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                        // this will cause problems later - stop evolution
        }

        ERRORS->Clean();                                                                                                        // clean the dynamic error catalog

        index++;                                                                                                                // next...
    }
    delete star;

    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE && index >= nStars) evolutionStatus = EVOLUTION_STATUS::DONE;             // set done

    int nStarsRequested = OPTIONS->GridFilename().empty() ? OPTIONS->SingleStarMassSteps() : (evolutionStatus == EVOLUTION_STATUS::DONE ? nStarsCreated : -1);

    SAY("\nGenerated " << std::to_string(nStarsCreated) << " of " << (nStarsRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nStarsRequested)) << " stars requested");

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                                    // shouldn't be
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN(ERROR::STELLAR_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));                   // show warning
        }
    }


    double cpuSeconds = (clock() - clockStart) / (double) CLOCKS_PER_SEC;                                                       // stop CPU timer and calculate seconds

    auto wallEnd = std::chrono::system_clock::now();                                                                            // stop wall timer
    std::time_t timeEnd = std::chrono::system_clock::to_time_t(wallEnd);                                                        // get end time and date

    SAY("\nEnd generating stars at " << std::ctime(&timeEnd));
    SAY("Clock time = " << cpuSeconds << " CPU seconds");


    std::chrono::duration<double> wallSeconds = wallEnd - wallStart;                                                            // elapsed seconds

    int wallHH = (int)(wallSeconds.count() / 3600.0);                                                                           // hours
    int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                                               // minutes
    int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));                              // seconds

    SAY("Wall time  = " << wallHH << ":" << wallMM << ":" << wallSS << " (hh:mm:ss)");

    return  std::make_tuple(nStarsRequested, nStarsCreated);
}


/*
 * Open the BSE grid file, read and parse the header record
 *
 *
 * std::tuple<int, std::vector<std::string>> OpenBSEGridFile(std::ifstream &p_Grid, const std::string p_Filename)
 *
 * @param   [IN/OUT]    p_Grid                  The std::ifstream  to which the file should be opened
 * @param   [IN]        p_Filename              The filename of the Grid file
 * @return                                      tuple containing the line number of the next line to be read, and a vector of header strings
 */
std::tuple<int, std::vector<std::string>> OpenBSEGridFile(std::ifstream &p_Grid, const std::string p_Filename) {

    std::vector<std::string> gridHeaders = {};                                                                                      // grid file headers

    int mass1               = 0;                                                                                                    // count 'Mass_1" occurrences
    int mass2               = 0;                                                                                                    // count 'Mass_2" occurrences
    int metallicity1        = 0;                                                                                                    // count 'Metallicity_1" occurrences
    int metallicity2        = 0;                                                                                                    // count 'Metallicity_2" occurrences

    int separation          = 0;                                                                                                    // count 'Separation" occurrences
    int eccentricity        = 0;                                                                                                    // count 'Eccentricity" occurrences
    int period              = 0;                                                                                                    // count 'Period" occurrences

    int kickMagnitudeRandom1 = 0;                                                                                                    // count 'Kick_Magnitude_Random_1' occurrences
    int kickMagnitude1       = 0;                                                                                                    // count 'Kick_Magnitude_1" occurrences
    int kickTheta1          = 0;                                                                                                    // count 'Kick_Theta_1" occurrences
    int kickPhi1            = 0;                                                                                                    // count 'Kick_Phi_1" occurrences
    int kickMeanAnomaly1    = 0;                                                                                                    // count 'Kick_Mean_Anomaly_1" occurrences

    int kickMagnitudeRandom2 = 0;                                                                                                    // count 'Kick_Magnitude_Random_1' occurrences
    int kickMagnitude2       = 0;                                                                                                    // count 'Kick_Magnitude_2" occurrences
    int kickTheta2          = 0;                                                                                                    // count 'Kick_Theta_2" occurrences
    int kickPhi2            = 0;                                                                                                    // count 'Kick_Phi_2" occurrences
    int kickMeanAnomaly2    = 0;                                                                                                    // count 'Kick_Mean_Anomaly_2" occurrences

    int lineNo = 1;                                                                                                                 // line number - for error messages

    if (!OPTIONS->GridFilename().empty()) {                                                                                         // have grid filename?

        p_Grid.open(OPTIONS->GridFilename());                                                                                       // yes - open the file
        if (p_Grid.fail()) {                                                                                                        // open ok?
            SAY(ERR_MSG(ERROR::FILE_OPEN_ERROR) << " " << OPTIONS->GridFilename());                                                 // no - show error
        }
        else {                                                                                                                      // file open ok

            bool emptyRecord = true;
            std::string record;
            while (emptyRecord && std::getline(p_Grid, record)) {                                                                   // read the header - first non-empty record
                size_t hashPos = record.find("#");                                                                                  // find first occurrence of "#"
                if (hashPos != std::string::npos) record.erase(hashPos, record.size() - hashPos);                                   // if "#" found, prune it and everything after it (ignore comments)

                while (record.size() > 0 && (record[0] == ' ' || record[0] == '\t')) record.erase(0, 1);                            // strip any leading ' ' or TAB ('\t') characters from the record

                while (record.size() > 0               &&
                      (record[record.size()-1] == '\r' ||
                       record[record.size()-1] == '\n' ||
                       record[record.size()-1] == '\t' ||
                       record[record.size()-1] == ' '  )) record.erase(record.size()-1, 1);                                         // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the record

                if (record.empty()) {
                    lineNo++;                                                                                                       // increment line number
                    continue;                                                                                                       // skip empty record
                }

                emptyRecord = false;                                                                                                // have non-empty record

                record = utils::ToUpper(record);                                                                                    // upshift

                std::stringstream recordSS(record);                                                                                 // open stream on record
                std::string token;
                while (std::getline(recordSS, token, ',')) {                                                                        // get token from record read

                    while (token.size() > 0 && (token[0] == ' ' || token[0] == '\t')) token.erase(0, 1);                            // strip any leading ' ' or TAB ('\t') characters from the token

                    while (token.size() > 0              &&
                          (token[token.size()-1] == '\r' ||
                           token[token.size()-1] == '\n' ||
                           token[token.size()-1] == '\t' ||
                           token[token.size()-1] == ' ')) token.erase(token.size()-1, 1);                                           // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the token

                    if (token.empty()) {                                                                                            // empty token?
                        SAY(ERR_MSG(ERROR::GRID_FILE_EMPTY_HEADER));                                                                // show error
                        continue;                                                                                                   // next token
                    }

                    switch (_(token.c_str())) {                                                                                     // which column header?

                        case _("MASS_1")                : mass1++;                  gridHeaders.push_back(token); break;            // Star 1 Mass
                        case _("MASS_2")                : mass2++;                  gridHeaders.push_back(token); break;            // Star 2 Mass

                        case _("METALLICITY_1")         : metallicity1++;           gridHeaders.push_back(token); break;            // Star 1 Metallicity
                        case _("METALLICITY_2")         : metallicity2++;           gridHeaders.push_back(token); break;            // Star 2 Metallicity

                        case _("SEPARATION")            : separation++;             gridHeaders.push_back(token); break;            // Separation

                        case _("ECCENTRICITY")          : eccentricity++;           gridHeaders.push_back(token); break;            // Eccentricity

                        case _("PERIOD")                : period++;                 gridHeaders.push_back(token); break;            // Period

                        case _("KICK_MAGNITUDE_RANDOM_1"): kickMagnitudeRandom1++;    gridHeaders.push_back(token); break;            // Star 1 Kick magnitude random number

                        case _("KICK_MAGNITUDE_1")       : kickMagnitude1++;          gridHeaders.push_back(token); break;            // Star 1 Kick magnitude

                        case _("KICK_THETA_1")          : kickTheta1++;             gridHeaders.push_back(token); break;            // Star 1 Kick theta

                        case _("KICK_PHI_1")            : kickPhi1++;               gridHeaders.push_back(token); break;            // Star 1 Kick phi

                        case _("KICK_MEAN_ANOMALY_1")   : kickMeanAnomaly1++;       gridHeaders.push_back(token); break;            // Star 1 Kick mean anomaly

                        case _("KICK_MAGNITUDE_RANDOM_2"): kickMagnitudeRandom2++;    gridHeaders.push_back(token); break;            // Star 2 Kick magnitude random number

                        case _("KICK_MAGNITUDE_2")       : kickMagnitude2++;          gridHeaders.push_back(token); break;            // Star 2 Kick magnitude

                        case _("KICK_THETA_2")          : kickTheta2++;             gridHeaders.push_back(token); break;            // Star 2 Kick theta

                        case _("KICK_PHI_2")            : kickPhi2++;               gridHeaders.push_back(token); break;            // Star 2 Kick phi

                        case _("KICK_MEAN_ANOMALY_2")   : kickMeanAnomaly2++;       gridHeaders.push_back(token); break;            // Star 2 Kick mean anomaly

                        default:                                                                                                    // unknown column header
                            SAY(ERR_MSG(ERROR::GRID_FILE_UNKNOWN_HEADER) << ": " << token);                                         // show error
                    }
                }

                if (mass1 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Mass_1")                                          // no 'Mass_1'
                else if (mass1 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Mass_1");                                  // duplicate 'Mass_1'

                if (mass2 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Mass_2")                                          // no 'Mass_2'
                else if (mass2 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Mass_2");                                  // duplicate 'Mass_2'

                if (metallicity1 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Metallicity_1")                            // no 'Metallicity_1'
                else if (metallicity1 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Metallicity_1");                    // duplicate 'Metallicity_1'

                if (metallicity2 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Metallicity_2")                            // no 'Metallicity_2'
                else if (metallicity2 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Metallicity_2");                    // duplicate 'Metallicity_2'

                if (eccentricity < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Eccentricity")                             // no 'Eccentricity'
                else if (eccentricity > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Eccentricity");                     // duplicate 'Eccentricity'

                if (separation < 1 && period < 1) {                                                                                 // neither 'Separation' nor 'Period'
                    SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": One of {Separation, Period}");
                }
                else {
                    if (separation > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Separation");                          // duplicate 'Separation'
                    if (period > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Period");                                  // duplicate 'Period'
                }


                if ((kickMagnitudeRandom1 + kickMagnitude1 + kickTheta1 + kickPhi1 + kickMeanAnomaly1 + 
                     kickMagnitudeRandom2 + kickMagnitude2 + kickTheta2 + kickPhi2 + kickMeanAnomaly2) > 0) {                         // at least one kick* header present, so all are required

                    if (kickMagnitudeRandom1 < 1 && kickMagnitude1 < 1) {                                                             // neither 'Kick_Magnitude_Random_1' nor 'Kick_Magnitude_1'
                        SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": One of {Kick_Magnitude_Random_1, Kick_Magnitude_1}");
                    }
                    else {
                        if (kickMagnitudeRandom1 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Magnitude_Random_1"); // duplicate 'Kick_Magnitude_Random_1'
                        if (kickMagnitude1 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Magnitude_1");              // duplicate 'Kick_Magnitude_1'
                    }
                    
                    if (kickTheta1 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Kick_Theta_1")                           // no 'Kick_Theta_1'
                    else if (kickTheta1 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Theta_1");                   // duplicate 'Kick_Theta_1'

                    if (kickPhi1 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Kick_Phi_1")                               // no 'Kick_Phi_1'
                    else if (kickPhi1 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Phi_1");                       // duplicate 'Kick_Phi_1'

                    if (kickMeanAnomaly1 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Kick_Mean_Anomaly_1")              // no 'Kick_Mean_Anomaly_1'
                    else if (kickMeanAnomaly1 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Mean_Anomaly_1");      // duplicate 'Kick_Mean_Anomaly_1'


                    if (kickMagnitudeRandom2 < 1 && kickMagnitude2 < 1) {                                                             // neither 'Kick_Magnitude_Random_2' nor 'Kick_Magnitude_2'
                        SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": One of {Kick_Magnitude_Random_2, Kick_Magnitude_2}");
                    }
                    else {
                        if (kickMagnitudeRandom2 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Magnitude_Random_2"); // duplicate 'Kick_Magnitude_Random_2'
                        if (kickMagnitude2 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Magnitude_2");              // duplicate 'Kick_Magnitude_2'
                    }

                    if (kickTheta2 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Kick_Theta_2")                           // no 'Kick_Theta_2'
                    else if (kickTheta2 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Theta_2");                   // duplicate 'Kick_Theta_2'

                    if (kickPhi2 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Kick_Phi_2")                               // no 'Kick_Phi_2'
                    else if (kickPhi2 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Phi_2");                       // duplicate 'Kick_Phi_2'

                    if (kickMeanAnomaly2 < 1) SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_HEADER) << ": Kick_Mean_Anomaly_2")              // no 'Kick_Mean_Anomaly_2'
                    else if (kickMeanAnomaly2 > 1) SAY(ERR_MSG(ERROR::GRID_FILE_DUPLICATE_HEADER) << ": Kick_Mean_Anomaly_2");      // duplicate 'Kick_Mean_Anomaly_2'
                }

                lineNo++;                                                                                                           // increment line number
            }
        }
    }

    // check we have all the headers we need, and in the right numbers

    bool error = false;                                                                                                             // assume no error initially

    error = !(mass1                == 1 && mass2        == 1  &&                                                                    // must have exactly one each of mass1 and mass2
              metallicity1         == 1 && metallicity2 == 1  &&                                                                    // must have exactly one each of metallicity1 and metallicity2
              separation           <= 1 && period       <= 1  &&                                                                    // must have at most one each of separation and period
             (separation + period) >  0 &&                                                                                          // must have at least one of separation and period
              eccentricity         == 1);                                                                                           // must have exactly one eccentricity

    if ((kickMagnitudeRandom1 + kickMagnitude1 + kickTheta1 + kickPhi1 + kickMeanAnomaly1 +                                           // if any kick parameter is present
         kickMagnitudeRandom2 + kickMagnitude2 + kickTheta2 + kickPhi2 + kickMeanAnomaly2) > 0) {     

        error = error || !(
        
                kickMagnitudeRandom1 <= 1 && kickMagnitude1 <= 1 &&                                                                   // must have at most one each of kickMagnitudeRandom1 and kickMagnitude1
               (kickMagnitudeRandom1 + kickMagnitude1)      >  0 &&                                                                   // must have at least one of kickMagnitudeRandom1 and kickMagnitude1
                kickTheta1                                == 1 &&                                                                   // must have exactly one kickTheta1
                kickPhi1                                  == 1 &&                                                                   // must have exactly one kickPhi1
                kickMeanAnomaly1                          == 1 &&                                                                   // must have exactly one kickMeanAnomaly1

                kickMagnitudeRandom2 <= 1 && kickMagnitude2 <= 1 &&                                                                   // must have at most one each of kickMagnitudeRandom2 and kickMagnitude2
               (kickMagnitudeRandom2 + kickMagnitude2)      >  0 &&                                                                   // must have at least one of kickMagnitudeRandom2 and kickMagnitude2
                kickTheta2                                == 1 &&                                                                   // must have exactly one kickTheta2
                kickPhi2                                  == 1 &&                                                                   // must have exactly one kickPhi2
                kickMeanAnomaly2                          == 1                                                                      // must have exactly one kickMeanAnomaly2
            );
    }
    if (error && p_Grid.is_open()) p_Grid.close();                                                                                  // we don't - must have been a problem - close the grid file

    return std::make_tuple(lineNo, gridHeaders);
}


/*
 * Read and parse the next record from the BSE grid file
 *
 * Expects the following units:
 * 
 *    Mass      : Msol
 *    Separation: AU
 *    Period    : Days
 *    Random    : floating point number uniformly distributed in the range [0.0, 1.0) 
 *    Velocity  : kms^-1
 *    Theta, Phi: radians
 *    Anomaly   : radians (?)
 * 
 * 
 * Returns the following data values from the Grid file:
 *
 *     Mass_1,
 *     Mass_2,
 *     Metallicity_1, 
 *     Metallicity_2, 
 *     Separation, 
 *     Eccentricity,
 *     Kick_Magnitude_Random_1, 
 *     Kick_Magnitude_1, 
 *     Kick_Theta_1, 
 *     Kick_Phi_1,
 *     Kick_Mean_Anomaly_1, 
 *     Kick_Magnitude_Random_2, 
 *     Kick_Magnitude_2, 
 *     Kick_Theta_2, 
 *     Kick_Phi_2,
 *     Kick_Mean_Anomaly_2
 * 
 * Plus boolean flags (see definition of KickParameters in typedefs.h):
 * 
 *     supplied{1,2}          - true if kick values were supplied in the grid file
 *     useMagnitudeRandom{1,2} - true if the user supplied the kick magnitude magnitude random number
 *
 * If the user specifies Period rather than Separation, the separation is calculated using the masses and the orbital period
 * If the user specifies both Separation and Period, Separation is used in preference to Period
 * 
 * If the user specifies the kick magnitude magnitude random number, the appropriate flag is set (per star)
 * If the use specifies both the kick magnitude magnitude random number and the kick magnitude, the random number will be used in preference to the supplied velocity
 * 
 * Missing values are treated as zero (0.0) - a warning will be issued, and reading of the Grid file continues
 * (A value is considered missing only if there is a header for the column, but no data value in the column)
 * Invalid values are treated as errors - an error will be issued and reading of the Grid file stopped
 * Negative values for Mass, Metallicity, Separation and Eccentricity are treated as errors - an error will be issued and reading of the Grid file stopped
 *
 * JR todo: there are some hard-coded numbers here (e.g. size of p_Headers vector) that should be made constants one day
 * 
 * 
 * std::tuple<int, std::vector<double>> ReadBSEGridRecord(std::ifstream &p_Grid, const std::vector<std::string> p_GridHeaders, const int p_LineNo)
 *
 * @param   [IN/OUT]    p_Grid                  The Grid file (std::ifstream)
 * @param   [IN]        p_GridHeaders           Vector of header strings returned by OpenBSEGridFile()
 * @param   [IN]        p_LineNo                The line number of the next line to be read
 * @return                                      tuple containing the error status, line number of the next line to be read, and the data values
 */
std::tuple<bool, int, BSEGridParameters> ReadBSEGridRecord(std::ifstream &p_Grid, const std::vector<std::string> p_GridHeaders, const int p_LineNo) {

    bool error = false;

    BSEGridParameters gridValues;

    // initialise grid values

    gridValues.mass1                                 = 0.0;
    gridValues.mass2                                 = 0.0;
    gridValues.metallicity1                          = 0.0; 
    gridValues.metallicity2                          = 0.0;
    gridValues.separation                            = 0.0;
    gridValues.eccentricity                          = 0.0;

    gridValues.star1KickParameters.supplied          = false;
    gridValues.star1KickParameters.useMagnitudeRandom = false;
    gridValues.star1KickParameters.magnitudeRandom    = 0.0;
    gridValues.star1KickParameters.magnitude          = 0.0;
    gridValues.star1KickParameters.theta             = 0.0;
    gridValues.star1KickParameters.phi               = 0.0;
    gridValues.star1KickParameters.meanAnomaly       = 0.0;

    gridValues.star2KickParameters.supplied          = false;
    gridValues.star2KickParameters.useMagnitudeRandom = false;
    gridValues.star2KickParameters.magnitudeRandom    = 0.0;
    gridValues.star2KickParameters.magnitude          = 0.0;
    gridValues.star2KickParameters.theta             = 0.0;
    gridValues.star2KickParameters.phi               = 0.0;
    gridValues.star2KickParameters.meanAnomaly       = 0.0;

    int lineNo = p_LineNo;                                                                                                  // line number - for error messages

    bool emptyRecord = true;
    std::string record;
    while (!error && emptyRecord && std::getline(p_Grid, record)) {                                                         // read the next record of the file

        size_t hashPos = record.find("#");                                                                                  // find first occurrence of "#"
        if (hashPos != std::string::npos) record.erase(hashPos, record.size() - hashPos);                                   // if "#" found, prune it and everything after it (ignore comments)

        for (size_t pos = 0; pos < record.size(); pos++) if (record[pos] == '\t') record[pos] = ' ';                        // replace tab with space

        while (record.size() > 0 && (record[0] == ' ' || record[0] == '\t')) record.erase(0, 1);                            // strip any leading ' ' or TAB ('\t') characters from the record

        while (record.size() > 0               &&
              (record[record.size()-1] == '\r' ||
               record[record.size()-1] == '\n' ||
               record[record.size()-1] == '\t' ||
               record[record.size()-1] == ' '  )) record.erase(record.size()-1, 1);                                         // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the record

        if (record.empty()) {
            lineNo++;                                                                                                       // increment line number
            continue;                                                                                                       // skip empty record
        }
        emptyRecord = false;                                                                                                // have non-empty record

        std::istringstream recordSS(record);                                                                                // open input stream on record

        double period = 0.0;

        size_t column = 0;                                                                                                  // start at column 0
        std::string token;                                                                                                  // token from record
        while (!error && std::getline(recordSS, token, ',')) {                                                              // get next token from record read

            while (token.size() > 0 && (token[0] == ' ' || token[0] == '\t')) token.erase(0, 1);                            // strip any leading ' ' or TAB ('\t') characters from the token

            while (token.size() > 0              &&
                  (token[token.size()-1] == '\r' ||
                   token[token.size()-1] == '\n' ||
                   token[token.size()-1] == '\t' ||
                   token[token.size()-1] == ' '  )) token.erase(token.size()-1, 1);                                         // strip any trailing '\r', '\n', TAB ('\t'), and ' ' characters from the token

            if (column >= p_GridHeaders.size()) {
                SAY(ERR_MSG(ERROR::GRID_FILE_EXTRA_COLUMN) << " ignored");                                                  // show error
            }
            else {

                if (token.empty()) {                                                                                        // empty token?
                    token = "0.0";                                                                                          // yes - use 0.0
                    SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_DATA) << " at line " << lineNo << ": 0.0 used");                   // show warning
                }

                double value;
                std::istringstream tokenSS(token);                                                                          // open input stream on token
                tokenSS >> value;                                                                                           // read double from token
                if (tokenSS.fail()) {                                                                                       // check for valid conversion
                    error = true;                                                                                           // set error flag
                    SAY(ERR_MSG(ERROR::GRID_FILE_INVALID_DATA) << " at line " << lineNo << ": " << token);                  // show error
                }

                // don't use utils::Compare() here for bounds/range checks

                if (!error) {
                    std::string columnName = p_GridHeaders[column];                                                         // get column name for value
                    switch (_(columnName.c_str())) {                                                                        // which column?

                        case _("MASS_1"):                                                                                   // Star 1 Mass
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else gridValues.mass1 = value;                                                                  // no - proceed
                            break;

                        case _("MASS_2"):                                                                                   // Star 2 Mass
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else gridValues.mass2 = value;                                                                  // no - proceed
                            break;

                        case _("METALLICITY_1"):                                                                            // Star 1 Metallicity
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else gridValues.metallicity1 = value;                                                           // no - proceed
                            break;

                        case _("METALLICITY_2"):                                                                            // Star 2 Metallicity
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else gridValues.metallicity2 = value;                                                           // no - proceed
                            break;

                        case _("SEPARATION"):                                                                               // Separation
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else gridValues.separation = value;                                                             // no - proceed
                            break;

                        case _("ECCENTRICITY"):                                                                             // Eccentricity
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else gridValues.eccentricity = value;                                                           // no - proceed
                            break;

                        case _("PERIOD"):                                                                                   // Period
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else period = value;                                                                            // no - proceed
                            break;

                        case _("KICK_MAGNITUDE_RANDOM_1"):                                                                   // Star 1 Kick magnitude random number
                            if (value < 0.0 || value >= 1.0) {                                                              // in the range [0.0, 1.0)? 
                                error = true;                                                                               // no - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_INVALID_DATA) << " at line " << lineNo << ": " << token);      // show error
                            }
                            else {                                                                                          // yes - proceed
                                gridValues.star1KickParameters.supplied          = true;                                    // Star 1 kick parameters supplied
                                gridValues.star1KickParameters.magnitudeRandom    = value;                                   // Star 1 Kick magnitude random number
                                gridValues.star1KickParameters.useMagnitudeRandom = true;                                    // use this in preference to actual kick value
                            }
                            break;

                        case _("KICK_MAGNITUDE_1"):                                                                          // Star 1 Kick magnitude (magnitude only, so must be +ve - probably technically "speed" rather than "velocity")                        
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else {                                                                                          // no - proceed
                                gridValues.star1KickParameters.supplied = true;                                             // Star 1 kick parameters supplied
                                gridValues.star1KickParameters.magnitude = value;                                            // Star 1 Kick magnitude
                            }
                            break;

                        case _("KICK_THETA_1"):                                                                             // Star 1 Kick theta
                            gridValues.star1KickParameters.supplied = true;                                                 // Star 1 kick parameters supplied
                            gridValues.star1KickParameters.theta    = value;                                                // Star 1 Kick theta
                            break;  

                        case _("KICK_PHI_1"):                                                                               // Star 1 Kick phi
                            gridValues.star1KickParameters.supplied = true;                                                 // Star 1 kick parameters supplied
                            gridValues.star1KickParameters.phi      = value;                                                // Star 1 Kick phi
                            break;  

                        case _("KICK_MEAN_ANOMALY_1"):                                                                      // Star 1 Kick mean anomaly
                            if (value < 0.0 || value >= _2_PI) {                                                            // in the range [0.0, _2_PI)? 
                                error = true;                                                                               // no - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_INVALID_DATA) << " at line " << lineNo << ": " << token);      // show error
                            }
                            else {                                                                                          // yes - proceed
                                gridValues.star1KickParameters.supplied    = true;                                          // Star 1 kick parameters supplied
                                gridValues.star1KickParameters.meanAnomaly = value;                                         // Star 1 Kick mean anomaly
                            }
                            break;  

                        case _("KICK_MAGNITUDE_RANDOM_2"):                                                                   // Star 1 Kick magnitude random number
                            if (value < 0.0 || value >= 1.0) {                                                              // in the range [0.0, 1.0)? 
                                error = true;                                                                               // no - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_INVALID_DATA) << " at line " << lineNo << ": " << token);      // show error
                            }
                            else {                                                                                          // yes - proceed                        
                                gridValues.star2KickParameters.supplied          = true;                                    // Star 2 kick parameters supplied
                                gridValues.star2KickParameters.magnitudeRandom    = value;                                   // Star 2 Kick magnitude random number
                                gridValues.star2KickParameters.useMagnitudeRandom = true;                                    // use this in preference to actual kick value
                            }
                            break;

                        case _("KICK_MAGNITUDE_2"):                                                                          // Star 2 Kick magnitude (magnitude only, so must be +ve - probably technically "speed" rather than "velocity")                        
                            if (value < 0.0) {                                                                              // value < 0?
                                error = true;                                                                               // yes - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_NEGATIVE_DATA) << " at line " << lineNo << ": " << token);     // show error
                            }
                            else {                                                                                          // no - proceed
                                gridValues.star2KickParameters.supplied = true;                                             // Star 2 kick parameters supplied
                                gridValues.star2KickParameters.magnitude = value;                                            // Star 2 Kick magnitude
                            }
                            break;

                        case _("KICK_THETA_2"):                                                                             // Star 2 Kick theta
                            gridValues.star2KickParameters.supplied = true;                                                 // Star 2 kick parameters supplied
                            gridValues.star2KickParameters.theta    = value;                                                // Star 2 Kick theta
                            break;

                        case _("KICK_PHI_2"):                                                                               // Star 2 Kick phi
                            gridValues.star2KickParameters.supplied = true;                                                 // Star 2 kick parameters supplied
                            gridValues.star2KickParameters.phi      = value;                                                // Star 2 Kick phi
                            break;

                        case _("KICK_MEAN_ANOMALY_2"):                                                                      // Star 2 Kick mean anomaly
                            if (value < 0.0 || value >= _2_PI) {                                                            // in the range [0.0, _2_PI)? 
                                error = true;                                                                               // no - set error flag
                                SAY(ERR_MSG(ERROR::GRID_FILE_INVALID_DATA) << " at line " << lineNo << ": " << token);      // show error
                            }
                            else {                                                                                          // yes - proceed
                                gridValues.star2KickParameters.supplied    = true;                                          // Star 2 kick parameters supplied
                                gridValues.star2KickParameters.meanAnomaly = value;                                         // Star 2 Kick mean anomaly
                            }
                            break;

                        default:                                                                                            // unknown - this shouldn't happen
                            SAY(ERR_MSG(ERROR::GRID_FILE_UNKNOWN_HEADER) << ": " << columnName);                            // show error
                    }
                }
            }
            column++;                                                                                                       // next column
        }

        if (!error) {
            if (column < p_GridHeaders.size()) {
                for (size_t c = column; c < p_GridHeaders.size(); c++) {
                    SAY(ERR_MSG(ERROR::GRID_FILE_MISSING_DATA) << " at line " << lineNo << ": 0.0 used");                   // show warning
                }
            }

            if (gridValues.separation <= 0.0 && period > 0.0 && gridValues.mass1 > 0.0 && gridValues.mass2 > 0.0) {         // already have separation? Have required values to calculate from period?  Don't use utils::Compare() here
// JR: change the next line to a warning that can be suppressed            
//                SAY(ERR_MSG(ERROR::GRID_FILE_USING_PERIOD));                                                                // let the user know we're using Period
                gridValues.separation = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(gridValues.mass1, gridValues.mass2, period);  // calculate separation from period
            }
        }
        lineNo++;                                                                                                           // increment line number
    }

    return std::make_tuple(error, lineNo, gridValues);
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

    signal(SIGUSR1, sigHandler);                                                                                            // install signal handler

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;

    int nBinariesCreated = 0;                                                                                               // number of binaries actually created

    auto wallStart = std::chrono::system_clock::now();                                                                      // start wall timer
    clock_t clockStart = clock();                                                                                           // start CPU timer

    if (!OPTIONS->Quiet()) {
        std::time_t timeStart = std::chrono::system_clock::to_time_t(wallStart);
        SAY("Start generating binaries at " << std::ctime(&timeStart));
    }

    AIS ais;                                                                                                                // Adaptive Importance Sampling (AIS)

    if (OPTIONS->AIS_ExploratoryPhase()) ais.PrintExploratorySettings();                                                    // print the selected options for AIS Exploratory phase in the beginning of the run
    if (OPTIONS->AIS_RefinementPhase() ) ais.DefineGaussians();                                                             // if we are sampling using AIS (step 2):read in gaussians

    std::ifstream grid;                                                                                                     // grid file
    std::vector<std::string> gridHeaders;                                                                                   // grid file headers

    int lineNo;                                                                                                             // grid file line number - for error messages
    if (!OPTIONS->GridFilename().empty()) {                                                                                 // have grid filename?
        std::tie(lineNo, gridHeaders) = OpenBSEGridFile(grid, OPTIONS->GridFilename());                                     // yes - read grid file headers
        if (!grid.is_open()) evolutionStatus = EVOLUTION_STATUS::ERROR;                                                     // failed to open/read grid file
    }

    // generate and evolve binaries
    int index    = 0;                                                                                                       // which binary
    int nBinaries = grid.is_open() ? 1 : OPTIONS->NBinaries();                                                              // assumption is OPTIONS->nBinaries > 0 JR: todo: should check options...

    BinaryStar *binary = nullptr;
    while (evolutionStatus == EVOLUTION_STATUS::CONTINUE && index < nBinaries) {

        evolvingBinaryStar      = NULL;                                                                                     // unset global pointer to evolving binary (for BSE Switch Log)
        evolvingBinaryStarValid = false;                                                                                    // indicate that the global pointer is not (yet) valid (for BSE Switch log)

        if (!OPTIONS->GridFilename().empty()) {                                                                             // have grid filename?
            if (grid.is_open()) {                                                                                           // yes - grid file open?

                bool error;
                BSEGridParameters gV;                                                                                       // yes - read the next record of grid values
                std::tie(error, lineNo, gV) = ReadBSEGridRecord(grid, gridHeaders, lineNo);                                 // read grid file record
                if (grid.fail() || error) {                                                                                 // read failed?
                    grid.close();                                                                                           // yes - EOF or error - stop reading, close file

                    if (error) evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                 // if error, set error
                    else       evolutionStatus = EVOLUTION_STATUS::DONE;                                                    // otherwise, set done
                }
                else {                                                                                                      // no

                    // generate a binary according to specified values in grid file
                    
                    nBinaries++;                                                                                            // not done yet...
                    delete binary;
                    binary = new BinaryStar(ais, 
                                            gV.mass1, gV.mass2, gV.metallicity1, gV.metallicity2, gV.separation, gV.eccentricity,
                                            gV.star1KickParameters, 
                                            gV.star2KickParameters, 
                                            index); 
                }
            }
            else evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                               // must have been a problem
        }
        else {                                                                                                              // no grid file
            delete binary;
            binary = new BinaryStar(ais, index);                                                                            // generate a random binary according to the user options
        }

        if (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {                                                                // still good?

            evolvingBinaryStar      = binary;                                                                               // set global pointer to evolving binary (for BSE Switch Log)
            evolvingBinaryStarValid = true;                                                                                 // indicate that the global pointer is now valid (for BSE Switch Log)

            EVOLUTION_STATUS binaryStatus = binary->Evolve();                                                               // evolve the binary

            // announce result of evolving the binary
            if (!OPTIONS->Quiet()) {
                if (OPTIONS->CHE_Option() == CHE_OPTION::NONE) {
                    SAY(index                                      << ": "  <<
                        EVOLUTION_STATUS_LABEL.at(binaryStatus)    << ": "  <<
                        STELLAR_TYPE_LABEL.at(binary->Star1Type()) << " + " <<
                        STELLAR_TYPE_LABEL.at(binary->Star2Type())
                    );
                }
                else {
                    SAY(index                                             << ": "    <<
                        EVOLUTION_STATUS_LABEL.at(binaryStatus)           << ": ("   <<
                        STELLAR_TYPE_LABEL.at(binary->Star1InitialType()) << " -> "  <<
                        STELLAR_TYPE_LABEL.at(binary->Star1Type())        << ") + (" <<
                        STELLAR_TYPE_LABEL.at(binary->Star2InitialType()) << " -> "  <<
                        STELLAR_TYPE_LABEL.at(binary->Star2Type())        <<  ")"
                    );
                }
            }

            nBinariesCreated++;                                                                                             // increment the number of binaries created

            if (OPTIONS->AIS_ExploratoryPhase() && ais.ShouldStopExploratoryPhase(index)) {                                 // AIS says should stop simulation?
                evolutionStatus = EVOLUTION_STATUS::AIS_EXPLORATORY;                                                        // ... and stop
            }

            if (!LOGGING->CloseStandardFile(LOGFILE::BSE_DETAILED_OUTPUT)) {                                                // close detailed output file if necessary
                SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                          // close failed - show warning
                evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                // this will cause problems later - stop evolution
            }

            if (!LOGGING->CloseStandardFile(LOGFILE::BSE_SWITCH_LOG)) {                                                     // close BSE switch log file if necessary
                SHOW_WARN(ERROR::FILE_NOT_CLOSED);                                                                          // close failed - show warning
                evolutionStatus = EVOLUTION_STATUS::STOPPED;                                                                // this will cause problems later - stop evolution
            }
        }

        ERRORS->Clean();                                                                                                    // clean the dynamic error catalog

        index++;                                                                                                            // next...
    }
    delete binary;

    if (evolutionStatus == EVOLUTION_STATUS::CONTINUE && index >= nBinaries) evolutionStatus = EVOLUTION_STATUS::DONE;      // set done

    int nBinariesRequested = OPTIONS->GridFilename().empty() ? OPTIONS->NBinaries() : (evolutionStatus == EVOLUTION_STATUS::DONE ? nBinariesCreated : -1);

    SAY("\nGenerated " << std::to_string(nBinariesCreated) << " of " << (nBinariesRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(nBinariesRequested)) << " binaries requested");

    if (evolutionStatus == EVOLUTION_STATUS::AIS_EXPLORATORY) {                                                             // AIS said stop?
        SHOW_WARN(ERROR::BINARY_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(evolutionStatus));                            // yes - show warning
        evolutionStatus = EVOLUTION_STATUS::DONE;                                                                           // set done
    }

    // announce result
    if (!OPTIONS->Quiet()) {
        if (evolutionStatus != EVOLUTION_STATUS::CONTINUE) {                                                                // shouldn't be...
            SAY("\n" << EVOLUTION_STATUS_LABEL.at(evolutionStatus));
        }
        else {
            SHOW_WARN(ERROR::BINARY_SIMULATION_STOPPED, EVOLUTION_STATUS_LABEL.at(EVOLUTION_STATUS::ERROR));                // show warning
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

    PROGRAM_STATUS programStatus = OPTIONS->Initialise(argc, argv);                     // Get the program options from the commandline

    if (programStatus == PROGRAM_STATUS::CONTINUE) {

        // start the logging service
        LOGGING->Start(OPTIONS->OutputPathString(),                                     // location of logfiles
                       OPTIONS->OutputContainerName(),                                  // directory to be created for logfiles
                       OPTIONS->LogfileNamePrefix(),                                    // prefix for logfile names
                       OPTIONS->LogLevel(),                                             // log level - determines (in part) what is written to log file
                       OPTIONS->LogClasses(),                                           // log classes - determines (in part) what is written to log file
                       OPTIONS->DebugLevel(),                                           // debug level - determines (in part) what debug information is displayed
                       OPTIONS->DebugClasses(),                                         // debug classes - determines (in part) what debug information is displayed
                       OPTIONS->DebugToFile(),                                          // should debug statements also be written to logfile?
                       OPTIONS->ErrorsToFile(),                                         // should error messages also be written to logfile?
                       DELIMITERValue.at(OPTIONS->LogfileDelimiter()));                 // log record field delimiter

        (void)utils::SplashScreen();                                                    // announce ourselves

        if (!LOGGING->Enabled()) programStatus = PROGRAM_STATUS::LOGGING_FAILED;        // logging failed to start
        else {

            RAND->Initialise();                                                         // initialise the random number service

            int objectsRequested = 0;
            int objectsCreated   = 0;

            if(OPTIONS->SingleStar()) {                                                 // Single star?
                std::tie(objectsRequested, objectsCreated) = EvolveSingleStars();       // yes - evolve single stars
            }
            else {                                                                      // no - binary
                std::tie(objectsRequested, objectsCreated) = EvolveBinaryStars();       // evolve binary stars
            }

            RAND->Free();                                                               // release gsl dynamically allocated memory

            LOGGING->Stop(std::make_tuple(objectsRequested, objectsCreated));           // stop the logging service

            programStatus = PROGRAM_STATUS::SUCCESS;                                    // set program status, and...
        }
    }

    return static_cast<int>(programStatus);                                             // we're done
}

