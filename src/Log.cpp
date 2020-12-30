// Class Log
//
// This is where all logging and debugging is performed.
//
// JR: todo: complete this documentation

#include "Log.h"

Log* Log::m_Instance = nullptr;


Log* Log::Instance() {
    if (!m_Instance) {
        m_Instance = new Log();
    }
    return m_Instance;
}


/*
 * Start logging.
 *
 * Set application logging parameters, including debug classes and levels.
 *
 * Calls to Say(), Put() and Debug() (especially via the macros provided) specify loggging and debug levels respectively.
 * Only calls which specify classes == enabled classes and levels <= enabled level will be actioned.  That is, calls to
 * Say(), Put() or Debug() that specify classes that are not enabled or levels above the respective application level will
 * be ignored.  In this way logging and, in particular, debugging statements can be made throughout the code and activated
 * or deactivated at runtime by specifying the enabled classes and application logging and debugging levels (via commandline
 * parameters).
 *
 *
 * Start(const string              p_LogBasePath,
 *       const string              p_LogContainerName,
 *       const string              p_LogNamePrefix,
 *       const int                 p_LogLevel,
 *       const std::vector<string> p_LogClasses,
 *       const int                 p_DbgLevel,
 *       const std::vector<string> p_DbgClasses,
 *       const bool                p_DbgToFile,
 *       const bool                p_ErrToFile,
 *       const string              p_LogfileType)
 *
 * @param   [IN]    p_LogBasePath               The path at which log files should be created
 * @param   [IN]    p_LogContainerName          The name of the directory that should be created at p_LogBasePath to hold all log files
 * @param   [IN]    p_LogNamePrefix             String to be prepended to logfile names - can be blank
 * @param   [IN]    p_LogLevel                  The application logging level
 * @param   [IN]    p_LogClasses                List of classes enabled for logging (vector<string>)
 * @param   [IN]    p_DbgLevel                  The application debug level
 * @param   [IN]    p_DbgClasses                List of classes enabled for debugging (vector<string>)
 * @param   [IN]    p_DbgToFile                 Boolean indicating whether debug records should also be written to a log file
 * @param   [IN]    p_ErrorsToFile              Boolean indicating whether error records should also be written to a log file
 * @param   [IN]    p_LogfileType               Log file type
 */
void Log::Start(const string              p_LogBasePath,
                const string              p_LogContainerName,
                const string              p_LogNamePrefix,
                const int                 p_LogLevel,
                const std::vector<string> p_LogClasses,
                const int                 p_DbgLevel,
                const std::vector<string> p_DbgClasses,
                const bool                p_DbgToLogfile,
                const bool                p_ErrorsToLogfile,
                const LOGFILETYPE         p_LogfileType) {

    H5Eset_auto (0, NULL, NULL);

    if (!m_Enabled) {

        m_Enabled       = true;                                                                                     // logging enabled;
        m_LogBasePath   = p_LogBasePath;                                                                            // set base path
        m_LogNamePrefix = p_LogNamePrefix;                                                                          // set log file name prefix
        m_LogLevel      = p_LogLevel;                                                                               // set log level
        m_LogClasses    = p_LogClasses;                                                                             // set enabled log classes
        m_DbgLevel      = p_DbgLevel;                                                                               // set debug level
        m_DbgClasses    = p_DbgClasses;                                                                             // set enagled debug classes
        m_DbgToLogfile  = p_DbgToLogfile;                                                                           // write debug records to logfile?
        m_ErrToLogfile  = p_ErrorsToLogfile;                                                                        // write error records to logfile?
        m_LogfileType   = p_LogfileType;                                                                            // set log file type

        m_Logfiles.clear();                                                                                         // clear all entries

        m_Enabled = UpdateAllLogfileRecordSpecs();                                                                  // update all logfile record specifications - disable logging upon failure

        if (m_Enabled) {

            // first create the container at p_LogBasePath
            // use boost filesystem here - easier...
        
            string containerName = p_LogContainerName;                                                              // container name with path ("/" works on Uni*x and Windows)
            string dirName       = containerName;                                                                   // directory name to create

            int version = 0;                                                                                        // container version number if required - start at 1
            while (boost::filesystem::exists(m_LogBasePath + "/" + dirName)) {                                      // container already exists?
                dirName = containerName + "_" + std::to_string(++version);                                          // yes - add a version number and generate new container name
            }
            m_LogContainerName = dirName;                                                                           // record actual container directory name

            boost::system::error_code err;
            try {
                boost::filesystem::create_directory(m_LogBasePath + "/" + m_LogContainerName, err);                 // create container - let boost throw an exception if it fails
                if (err.value() == 0) {                                                                             // ok?

                    if (m_DbgToLogfile) {                                                                           // write dubug output to a logfile?
                        string filename = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::DEBUG_LOG));                   // extract filename from descriptor
                        int id = Open(filename, false, true, false);                                                // open the log file - new file, timestamps, no record labels, space delimited
                        if (id >= 0) {                                                                              // success
                            m_DbgLogfileId = id;                                                                    // record the file id
                        }
                        else {                                                                                      // failure
                            Squawk("ERROR: Unable to create log file for debug output with file name " + filename); // announce error
                            Squawk("Debug output logging disabled");                                                // show disabled warning
                        }
                    }

                    if (m_ErrToLogfile) {                                                                           // write dubug output to a logfile?
                        string filename = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::ERROR_LOG));                   // extract filename from descriptor
                        int id = Open(filename, false, true, false);                                                // open the log file - new file, timestamps, no record labels, space delimited
                        if (id >= 0) {                                                                              // success
                            m_ErrLogfileId = id;                                                                    // record the file id
                        }
                        else {                                                                                      // failure
                            Squawk("ERROR: Unable to create log file for error output with file name " + filename); // announce error
                            Squawk("Error output logging disabled");                                                // show disabled warning
                        }
                    }
                }
                else  {                                                                                             // not ok...
                    Squawk("ERROR: Unable to create log file container with name " + dirName);                      // announce error
                    Squawk("Boost filesystem error = " + err.message());                                            // plus details
                    Squawk("Logging disabled");                                                                     // show disabled warning
                    m_Enabled = false;                                                                              // disable
                }
            
            }
            catch (...) {                                                                                           // unhandled problem...
                Squawk("ERROR: Unable to create log file container with name " + dirName);                          // announce error
                Squawk("Logging disabled");                                                                         // show disabled warning
                m_Enabled = false;                                                                                  // disable
            }

            // create run details file if all ok
            if (m_Enabled) {                                                                                        // ok?
                                                                                                                    // yes
                m_WallStart = std::chrono::system_clock::now();                                                     // start wall timer
                m_ClockStart = clock();                                                                             // start CPU timer

                string filename = m_LogBasePath + "/" + m_LogContainerName + "/" + RUN_DETAILS_FILE_NAME;           // run details filename with container name
                try {
                    m_RunDetailsFile.open(filename, std::ios::out);                                                 // create run details file
                    m_RunDetailsFile.exceptions(std::ofstream::failbit | std::ofstream::badbit);                    // enable exceptions on run details file

                    // file should be open - write the run details

                    try {                                                                                             

                        m_RunDetailsFile << utils::SplashScreen(false) << std::endl;                                // write splash string with version number to file

                        std::time_t timeStart = std::chrono::system_clock::to_time_t(m_WallStart);                  // record start time

                        // record start time and whether evolving single stars or binaries   
                        if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE)
                            m_RunDetailsFile << "Start generating stars at " << std::ctime(&timeStart) << std::endl;
                        else
                            m_RunDetailsFile << "Start generating binaries at " << std::ctime(&timeStart) << std::endl;

                        // run details file will be updated and closed in Log::Stop()
                    }
                    catch (const std::ofstream::failure &e) {                                                       // problem...
                        Squawk("ERROR: Unable to write to run details file with name " + filename);                 // announce error
                        Squawk(e.what());                                                                           // plus details
                    }
                }
                catch (const std::ofstream::failure &e) {                                                           // fs problem...
                    Squawk("ERROR: Unable to create run details file with file name " + filename);                  // announce error
                    Squawk(e.what());                                                                               // plus details
                    Squawk("Run details will not be recorded");                                                     // show warning
                }
                catch (...) {                                                                                       // unhandled problem...
                    Squawk("ERROR: Unable to create log file with file name " + filename);                          // announce error
                    Squawk("Run details will not be recorded");                                                     // show warning
                }
            }
        }
    }
}


/*
 * Stop logging
 *
 * Closes any open logfiles
 *
 *
 * void Stop()
 *
 */
void Log::Stop(std::tuple<int, int> p_ObjectStats) {
    if (m_Enabled) {
        CloseAllStandardFiles();                                                                                // first close all standard log files
        for(unsigned int index = 0; index < m_Logfiles.size(); index++) {                                       // check for open logfiles (even if not active)
            if (IsActiveId(index)) {                                                                            // logfile active?
                if (m_Logfiles[index].file.is_open()) {                                                         // open file?
                    try {                                                                                       // yes
                        m_Logfiles[index].file.flush();                                                         // flush output and
                        m_Logfiles[index].file.close();                                                         // close it
                    }
                    catch (const std::ofstream::failure &e) {                                                   // problem...
                        Squawk("ERROR: Unable to close log file with file name " + m_Logfiles[index].name);     // announce error
                        Squawk(e.what());                                                                       // plus details
                    }
                }
            }
        }

        // update run details file
        string filename = m_LogBasePath + "/" + m_LogContainerName + "/" + RUN_DETAILS_FILE_NAME;               // run details filename with container name
        try {  
            double cpuSeconds = (clock() - m_ClockStart) / (double) CLOCKS_PER_SEC;                             // stop CPU timer and calculate seconds

            m_WallEnd = std::chrono::system_clock::now();                                                       // stop wall timer
            std::time_t timeEnd = std::chrono::system_clock::to_time_t(m_WallEnd);                              // get end time and date
        
            // record end time and whether evolving single stars or binaries   
            if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE) {
                m_RunDetailsFile << "Generated " << std::to_string(std::get<1>(p_ObjectStats)) << " of " << (std::get<0>(p_ObjectStats) < 0 ? "<INCOMPLETE GRID>" : std::to_string(std::get<0>(p_ObjectStats))) << " stars requested" << std::endl;
                m_RunDetailsFile << "\nEnd generating stars at " << std::ctime(&timeEnd) << std::endl;
            }
            else {
                m_RunDetailsFile << "Generated " << std::to_string(std::get<1>(p_ObjectStats)) << " of " << (std::get<0>(p_ObjectStats) < 0 ? "<INCOMPLETE GRID>" : std::to_string(std::get<0>(p_ObjectStats))) << " binaries requested" << std::endl;
                m_RunDetailsFile << "\nEnd generating binaries at " << std::ctime(&timeEnd) << std::endl; 
            }

            m_RunDetailsFile << "Clock time = " << cpuSeconds << " CPU seconds" << std::endl;                   // record cpu second

            std::chrono::duration<double> wallSeconds = m_WallEnd - m_WallStart;                                // elapsed seconds

            int wallHH = (int)(wallSeconds.count() / 3600.0);                                                   // hours
            int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                       // minutes
            int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));      // seconds

            m_RunDetailsFile << "Wall time  = " << 
                                std::setfill('0') << std::setw(2) << wallHH << ":" << 
                                std::setfill('0') << std::setw(2) << wallMM << ":" << 
                                std::setfill('0') << std::setw(2) << wallSS << " (hh:mm:ss)" << std::endl;      // Include 0 buffer 

            m_RunDetailsFile << "\n\n" << OPTIONS->CmdLineOptionsDetails();                                     // record the commandline options details string
            m_RunDetailsFile << "Actual random seed = " << (OPTIONS->FixedRandomSeedCmdLine() ? OPTIONS->RandomSeedCmdLine() : RAND->DefaultSeed()) << ", CALCULATED, UNSIGNED_LONG" << std::endl;    // actual random seed
        }
        catch (const std::ofstream::failure &e) {                                                               // problem...
            Squawk("ERROR: Unable to update run details file with file name " + filename);                      // announce error
            Squawk(e.what());                                                                                   // plus details
        }

        // flush and close the file
        try {
            m_RunDetailsFile.flush();
            m_RunDetailsFile.close();
        }
        catch (const std::ofstream::failure &e) {                                                               // problem...
            Squawk("ERROR: Unable to close run details file with file name " + filename);                       // announce error
            Squawk(e.what());                                                                                   // plus details
        }
    }

    m_Logfiles.clear();                                                                                         // clear all entries
    m_Enabled = false;                                                                                          // set not enabled
}


/*
 * Create and open new log file
 *
 * New log file is created at path m_LogBasePath
 *
 * The file extension is based on the parameter p_Delimiter:
 *     - if the delimiter is TAB   the file extension is ".tsv" (Tab Separated Variables)
 *     - if the delimiter is COMMA the file extension is ".csv" (Comma Separated Variables)
 *     - if the delimiter is SPACE the file extension is ".txt"
 * JR FIX DESCRIOTION ABOVE <<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *
 * int Open(const string p_LogFileName, const bool p_Append, const bool p_TimeStamp, const bool p_Label, const LOGFILE p_StandardLogfile)
 *
 * @param   [IN]    p_LogFileName               The name of the logfile to be created and opened - filename only - path, prefix and extension are added
 * @param   [IN]    p_Append                    Boolean indicating whether an existing file of the same name should be opened and appended to
 *                                              (or whether a new (versioned) file should be opened).
 * @param   [IN]    p_Timestamp                 Boolean indicating whether a timestamp should be written with each log record
 * @param   [IN]    p_Label                     Boolean indicating whether a record label should be written with each log record
 * @param   [IN]    p_StandardLogfile           If Standard logfile, which (optional, default = LOGFILE::NONE)
 * @return                                      Logfile id (integer index into m_Logfiles vector).  A value of -1 indicates log file not opened successfully.
 */
int Log::Open(const string p_LogFileName, const bool p_Append, const bool p_Timestamp, const bool p_Label, const LOGFILE p_StandardLogfile) {

    bool ok = true;
    int id  = -1;  

    if (m_Enabled) {                                                                                            // logging enabled?   

        string basename = m_LogBasePath + "/" + m_LogContainerName + "/" + m_LogNamePrefix + p_LogFileName;     // base filename with path and container ("/" works on Uni*x and Windows)
        string fileext  = LOGFILETYPEFileExt.at(OPTIONS->LogfileType());                                        // file extension
        string filename = basename + "." + fileext;                                                             // full filename

        int version = 0;                                                                                        // logfile version number if required - start at 1
        while (utils::FileExists(filename) && !p_Append) {                                                      // file already exists - and we don't want to append?
            filename = basename + "_" + std::to_string(++version) + "." + fileext;                              // yes - add a version number and generate new filename
        }

        if (m_LogfileType == LOGFILETYPE::HDF5) {                                                               // HDF5 file?
                                                                                                                // yes
            // if we're logging to HDF5 files we should have a containing HDF5 file open.
            //
            // for log files other than detailed output files (SSE and BSE) that is a container
            // file to which all regular log files are written as groups in the container file.
            // The datasets (columns) are then written to each group.  The container file should
            // not already contain a group corresponding to this logfile - that is created here.
            //
            // for detailed output files (SSE and BSE) the containing HDF5 file is a separate
            // HDF5 file for each detailed output file.  Detailed output HDF5 files do not 
            // contain groups - the datasets (columns) are written directly to the file.

            hid_t  h5FileId = -1;                                                                                   // HDF5 file id
            hid_t  h5GroupId = -1;                                                                                  // HDF5 file group id
            string h5GroupName = "";                                                                                // HDF5 group name

            if (p_StandardLogfile == LOGFILE::SSE_DETAILED_OUTPUT || p_StandardLogfile == LOGFILE::BSE_DETAILED_OUTPUT) { // detailed output file?
                h5FileId  = m_HDF5DetailedId;                                                                       // yes - use detailed file id
                h5GroupId = h5FileId;                                                                               // no group for detailed file - just use the file id
            }
            else {                                                                                                  // no, not detailed ouput file
                h5FileId    = m_HDF5ContainerId;                                                                    // container file id
                h5GroupName = m_LogNamePrefix + p_LogFileName;                                                      // HDF5 file group name

                if (m_HDF5ContainerId >= 0) {                                                                       // yes - have HDF5 container file?
                    h5GroupId = H5Gopen(m_HDF5ContainerId, h5GroupName.c_str(), H5P_DEFAULT);                       // yes - open the group
                    if (h5GroupId < 0) {                                                                            // group open (and therefore already exists)?
                        h5GroupId = H5Gcreate(m_HDF5ContainerId, h5GroupName.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT); // no - create the group
                        if (h5GroupId < 0) {                                                                        // group created ok?
                            Squawk("ERROR: Error creating HDF5 group with name " + h5GroupName);                    // announce error
                            ok = false;                                                                             // fail
                        }
                    }
                    else {                                                                                          // yes - group exists and is open               
                        if (p_Append) {                                                                             // that's ok if we're appending - are we appending?
                            if (H5Gclose(h5GroupId) < 0) {                                                          // group closed ok?
                                Squawk("ERROR: Error closing HDF5 group with name " + h5GroupName);                 // no - announce error
                                ok = false;                                                                         // fail
                            }
                        }
                        else {                                                                                      // not appending - that's not ok...
                            Squawk("ERROR: HDF5 group with name " + h5GroupName + " already exists");               // announce error
                            (void)H5Gclose(h5GroupId);                                                              // close the group
                            ok = false;                                                                             // fail
                        }
                    }
                }
                else {                                                                                              // no - don't have HDF5 container file
                    Squawk("ERROR: HDF5 container file does not exist");                                            // that's an issue... announce error
                    ok = false;                                                                                     // fail
                }
            }

            if (ok) {                                                                                           // still ok?

                // record attributes
                // find an empty slot in m_Logfiles vector if there is one
                // this way is a bit slower for opening logfiles, but faster for writing to them
                id = -1;  
                for(unsigned int index = 0; index < m_Logfiles.size(); index++) {
                    if (!m_Logfiles[index].active) {                                                            // empty slot?
                        id = index;                                                                             // yes - use it
                        break;                                                                                  // and stop looking
                    }
                }

                if (id < 0) {                                                                                   // have empty slot?
                    logfileAttrT attr;                                                                          // no - create new attributes struct
                    id = m_Logfiles.size();                                                                     // set new id
                    m_Logfiles.push_back(std::move(attr));                                                      // append to vector
                }

                m_Logfiles[id].active           = true;                                                         // this entry now active
                m_Logfiles[id].logfiletype      = p_StandardLogfile;                                            // standard logfile type
                m_Logfiles[id].filetype         = m_LogfileType;                                                // set filetype for this log file
                m_Logfiles[id].name             = h5GroupName;                                                  // log file name (HDF5 group name)
                m_Logfiles[id].timestamp        = p_Timestamp;                                                  // set timestamp flag for this log file
                m_Logfiles[id].label            = p_Label;                                                      // set label flag for this log file
                m_Logfiles[id].h5file.fileId    = h5FileId;                                                     // HDF5 file id
                m_Logfiles[id].h5file.groupId   = h5GroupId;                                                    // HDF5 group id
                m_Logfiles[id].h5file.dataSets  = {};                                                           // HDF5 data sets
                m_Logfiles[id].h5file.dataTypes = {};                                                           // HDF5 data types
            }
        }
        else {                                                                                                  // no - not HDF5
            // record attributes

            // find an empty slot in m_Logfiles vector if there is one
            // this way is a bit slower for opening logfiles, but faster for writing to them
            for(unsigned int index = 0; index < m_Logfiles.size(); index++) {
                if (!m_Logfiles[index].active) {                                                                // empty slot?
                    id = index;                                                                                 // yes - use it
                    break;                                                                                      // and stop looking
                }
            }

            if (id < 0) {                                                                                       // have empty slot?
                id = m_Logfiles.size();                                                                         // no - set new id
                logfileAttrT attr;                                                                              // create new attributes struct
                attr.active = false;                                                                            // for now...
                m_Logfiles.push_back(std::move(attr));                                                          // append to vector
            }
               
            try {
                m_Logfiles[id].file.open(filename, std::ios::out | std::ios::app);                              // create fs log file
                m_Logfiles[id].file.exceptions(std::ofstream::failbit | std::ofstream::badbit);                 // enable exceptions on log file

                m_Logfiles[id].active         = true;                                                           // this entry now active
                m_Logfiles[id].logfiletype    = p_StandardLogfile;                                              // standard logfile type
                m_Logfiles[id].filetype       = m_LogfileType;                                                  // set filetype for this log file
                m_Logfiles[id].name           = filename;                                                       // log file name
                m_Logfiles[id].timestamp      = p_Timestamp;                                                    // set timestamp flag for this log file
                m_Logfiles[id].label          = p_Label;                                                        // set label flag for this log file
                m_Logfiles[id].h5file.fileId  = -1;                                                             // not HDF5 file
                m_Logfiles[id].h5file.groupId = -1;                                                             // not HDF5 file
            }
            catch (const std::ofstream::failure &e) {                                                           // fs problem...
                Squawk("ERROR: Unable to create fs log file with file name " + m_Logfiles[id].name);            // announce error
                Squawk(e.what());                                                                               // plus details

                // it's possible this m_Logfiles entry was just appended - it can be used next time
                if (id >= 0) ClearEntry(id);                                                                    // clear entry
                ok = false;                                                                                     // fail
            }
        }
    }

    return ok ? id : -1;                                                                                        // return log file id
}


/*
 * Close specified log file
 *
 *
 * bool Close(const int p_LogfileId)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to be closed
 * @return                                      Boolean indicating whether log file was closed successfully
 */
bool Log::Close(const int p_LogfileId) {

    bool result = true;

    if (m_Enabled) {                                                                                            // logging enabled?
                                                                                                                // yes
        int id = p_LogfileId;                                                                                   // file id to close
                                                                                                                // check if the file is an open standard file
        bool standardFile;
        LOGFILE logfile;
        std::tie(standardFile, logfile) = GetStandardLogfileKey(id);                                            // look in open standard file map
        if (standardFile) {                                                                                     // file is an open standard file
            COMPASUnorderedMap<LOGFILE, LogfileDetailsT>::const_iterator iter;                                  // iterator
            iter = m_OpenStandardLogFileIds.find(logfile);                                                      // get the details
            if (iter != m_OpenStandardLogFileIds.end()) {                                                       // found
                LogfileDetailsT fileDetails = iter->second;                                                     // existing file details
                id = fileDetails.id;                                                                            // standard file id
            }
        }

        result = Close_(id);                                                                                    // close the file

        if (result && standardFile) {                                                                           // if it was a standard file...
            m_OpenStandardLogFileIds.erase(logfile);                                                            // ...remove it from the map
        }
    }

    return result;
}


/*
 * Close specified log file - without first checking if it is a standard log file
 *
 *
 * bool Close_(const int p_LogfileId)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to be closed
 * @return                                      Boolean indicating whether log file was closed successfully
 */
bool Log::Close_(const int p_LogfileId) {

    bool result = true;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging enabled and logfile active?
        if (m_Logfiles[p_LogfileId].filetype == LOGFILETYPE::HDF5) {                                                // yes - HDF5 logfile?

            // first, close all open datasets
            for (auto &dataSet : m_Logfiles[p_LogfileId].h5file.dataSets) {                                         // for each open dataset
                if (H5Dclose(dataSet) < 0) {                                                                        // closed ok?
                    Squawk("ERROR: Unable to close HDF5 dataset for log file " + m_Logfiles[p_LogfileId].name);     // no - announce error
                    result = false;                                                                                 // fail
                }
            }

            // try closing the group - even if closing datasets failed
            if (m_Logfiles[p_LogfileId].logfiletype != LOGFILE::SSE_DETAILED_OUTPUT && 
                m_Logfiles[p_LogfileId].logfiletype != LOGFILE::BSE_DETAILED_OUTPUT) {                              // detailed output file?

                if (m_Logfiles[p_LogfileId].h5file.groupId >= 0) {                                                  // no - HDF5 group open?
                    if (H5Gclose(m_Logfiles[p_LogfileId].h5file.groupId) < 0) {                                     // yes - closed ok?
                        Squawk("ERROR: Unable to close HDF5 group for log file " + m_Logfiles[p_LogfileId].name);   // no - announce error
                        result = false;                                                                             // fail
                    }
                }
            }
            else {                                                                                                  // yes - detailed output file
                if (H5Fclose(m_Logfiles[p_LogfileId].h5file.fileId) < 0) {                                          // HDF5 file closed ok?
                    Squawk("ERROR: Unable to close HDF5 file " + m_Logfiles[p_LogfileId].name);                     // no - announce error
                    result = false;                                                                                 // fail
                }
                m_HDF5DetailedId = -1;                                                                              // (should have) no open detailed output file
            }
        }
        else {                                                                                                      // no, FS logfile
            if (m_Logfiles[p_LogfileId].file.is_open()) {                                                           // log file open?
                try {                                                                                               // yes
                    m_Logfiles[p_LogfileId].file.flush();                                                           // flush output and
                    m_Logfiles[p_LogfileId].file.close();                                                           // close it
                }
                catch (const std::ofstream::failure &e) {                                                           // problem...
                    Squawk("ERROR: Unable to close log file with file name " + m_Logfiles[p_LogfileId].name);       // announce error
                    Squawk(e.what());                                                                               // plus details
                    result = false;                                                                                 // fail
                }
            }
        }
        ClearEntry(p_LogfileId);                                                                                    // clear entry whether the close succeeded or not
    }

    return result;
}


/*
 * Writes string to stderr
 *
 *
 * void Squawk(const string p_SquawkStr)
 *
 * @param   [IN]    p_SquawkStr                 String to be written to stderr
 */
void Log::Squawk(const string p_SquawkStr) {
    std::cerr << p_SquawkStr << std::endl;                                                                      // write to stderr
}


/*
 * Determine whether record should be logged/debug string should be written based on log/debug classes and log/debug level.
 *
 *
 * bool DoIt(const string p_Class, const int p_Level, const std::vector<string> p_EnabledClasses, const int p_EnabledLevel)
 *
 * @param   [IN]    p_Class                     The class (logging or debug) of the record being evaluated
 * @param   [IN]    p_LogLevel                  The level (logging or debug) of the record being evaluated
 * @param   [IN]    p_EnabledClasses            List of classes enabled for logging or debugging (vector<string>)
 * @param   [IN]    p_EnabledLevel              The application logging or debug level
 * @return                                      Boolean indicating whether record should be logged/debug string should be written
 */
bool Log::DoIt(const string p_Class, const int p_Level, const std::vector<string> p_EnabledClasses, const int p_EnabledLevel) {

    bool doIt = (p_Level <= p_EnabledLevel);                                                                    // first check level

    if (doIt) {                                                                                                 // now check class
        if (!(p_Class.empty() || p_EnabledClasses.empty())) {                                                   // log class restricted?
            doIt = false;                                                                                       // until we know better...
            for (auto iter = p_EnabledClasses.begin(); iter != p_EnabledClasses.end(); ++iter) {                // check enabled classes
                if (utils::Equals(p_Class, *iter)) {                                                            // matches an enabled class?
                    doIt = true;                                                                                // yes, and...
                    break;                                                                                      // we're done
                }
            }
        }
    }

    return doIt;
}


/*
 * Write to stdout
 *
 * If class and level are enabled (see DoIt()), write string to stdout
 *
 *
 * void Say(const string p_SayClass, const int p_SayLevel, const string p_SayStr)
 *
 * @param   [IN]    p_SayClass                  Class to determine if string should be written
 * @param   [IN]    p_SayLevel                  Level to determine if string should be written
 * @param   [IN]    p_SayStr                    The string to be written
 */
void Log::Say(const string p_SayClass, const int p_SayLevel, const string p_SayStr) {
    if (m_Enabled) {                                                                                                // logging service enabled?
        if (DoIt(p_SayClass, p_SayLevel, m_LogClasses, m_LogLevel)) {                                               // logging this class and level?
            Say_(p_SayStr);                                                                                         // say it
        }
    }
}


/*
 * Say() with no class or level check - internal use only
 *
 *
 * void Say_(const string p_SayStr)
 *
 * @param   [IN]    p_SayStr                    The string to be written
 */
void Log::Say_(const string p_SayStr) {
    std::cout << p_SayStr << std::endl;                                                                             // write to stdout
}


/*
 * Write to specified log file
 *
 * If logging is enabled and the specified log file is active, and class and level are enabled (see DoIt()),
 * write string to log file
 *
 * Disable the specified log file if errors occur.
 *
 *
 * bool Write(const id p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log string should be written
 * @param   [IN]    p_LogClass                  Class to determine if string should be written
 * @param   [IN]    p_LogLevel                  Level to determine if string should be written
 * @param   [IN]    p_LogStr                    The string to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Write(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging service enabled and specified log file active?
        if (DoIt(p_LogClass, p_LogLevel, m_LogClasses, m_LogLevel)) {                                               // yes - logging this class and level?
            result = Write_(p_LogfileId, p_LogStr);                                                                 // yes - log it
        }
    }

    return result;
}


/*
 * Write a string record to specified log file with no class or level check - internal use only
 * Used for CSV, TSV, and TXT files
 * 
 * Disable the specified log file if errors occur.
 *
 *
 * bool Write_(const id p_LogfileId, const string p_LogStr)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log string should be written
 * @param   [IN]    p_LogStr                    The string to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Write_(const int p_LogfileId, const string p_LogStr) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging service enabled and specified log file active?
        try {
            m_Logfiles[p_LogfileId].file << p_LogStr << std::endl;                                                  // write string to log file
            m_Logfiles[p_LogfileId].file.flush();                                                                   // flush data to log file

            result = true;                                                                                          // set result
        }
        catch (const std::ofstream::failure &e) {                                                                   // problem...
            Squawk("ERROR: Unable to write to log file with file name " + m_Logfiles[p_LogfileId].name);            // announce error
            Squawk(e.what());                                                                                       // plus details
            Squawk("LOG RECORD: " + p_LogStr);                                                                      // show log record

            ClearEntry(p_LogfileId);                                                                                // clear entry
        }
    }
    else {                                                                                                          // not enabled or not active
        Squawk(p_LogStr);                                                                                           // show log record on stderr
    }

    return result;
}


/*
 * Write a multi-value record to specified log file with no class or level check - internal use only
 * Used for HDF5 files
 * 
 * This is where the work is done
 *
 * Disable the specified log file if errors occur.
 *
 *
 * bool Write_(const id p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const std::vector<TYPENAME> p_LogRecordTypes)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the multi-value record should be written
 * @param   [IN]    p_LogRecordValues           Vector of COMPAS_VARIABLE_TYPE values to be written
 * @param   [IN]    p_LogRecordTypes            Vector of TYPENAME values indicating the COMPAS data types of the values to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Write_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const std::vector<TYPENAME> p_LogRecordTypes) {

    herr_t ok = -1;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging service enabled and specified log file active?


// JR CHECK ALL POINTS OF FAILURE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // CHECK IF HDF5 FILE EXISTS
    // HOW TO HANDLE DETAILED OUTPUT FILES

// m_HDF5ContainerId


        if (m_Logfiles[p_LogfileId].h5file.groupId >= 0) {                                                      // HDF5 group open?
            for (size_t idx = 0; idx < p_LogRecordValues.size(); idx++) {                 // for each property

                hid_t dSet  = m_Logfiles[p_LogfileId].h5file.dataSets[idx];

                hid_t dType = m_Logfiles[p_LogfileId].h5file.dataTypes[idx];




                // create a 1-d HDF5 dataspace

                
//hsize_t size = H5Dget_storage_size(dSet) + H5Tget_size(dType);


                hsize_t currentSize = H5Dget_storage_size(dSet) / H5Tget_size(dType);
                hsize_t h5Dims[1] = {1};
//                h5Dims[0] = currentSize + 1;   // chunk size is 1
                hid_t   h5Dspace  = H5Screate_simple(1, h5Dims, NULL);                                                  // create the dataspace
        if (h5Dspace < 0) {
Squawk("Create simple failed!");
ok = -1;
        }


    // Extend dataset
    h5Dims[0] = currentSize + 1;   // chunk size is 1
    ok = H5Dset_extent(dSet, h5Dims);
    if (ok < 0) {
Squawk("Extend failed!");
    }


        // create file space
        hid_t   h5FSpace = H5Dget_space(dSet);
        if (h5FSpace < 0) {
Squawk("Get space failed!");
ok = -1;
        }

    // Select hyperslab on file space
    hsize_t h5Start[1] = {0};
    h5Start[0] = currentSize;
    hsize_t h5Count[1] = {1};
    ok = H5Sselect_hyperslab(h5FSpace, H5S_SELECT_SET, h5Start, NULL, h5Count, NULL);
    if (ok < 0) {
Squawk("Hyperslab failed!");
    }



    if (ok >= 0) {

                // can't use switch here...

                if (dType == H5T_NATIVE_UCHAR) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<bool>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_SHORT) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<short int>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_INT) {

                    int v = 0;
                    switch (p_LogRecordTypes[idx]) {
                        case TYPENAME::INT: v = static_cast<int>(boost::get<int>(p_LogRecordValues[idx])); break;
                        case TYPENAME::ERROR: v = static_cast<int>(boost::get<ERROR>(p_LogRecordValues[idx])); break;
                        case TYPENAME::STELLAR_TYPE: v = static_cast<int>(boost::get<STELLAR_TYPE>(p_LogRecordValues[idx])); break;
                        case TYPENAME::MT_CASE: v = static_cast<int>(boost::get<MT_CASE>(p_LogRecordValues[idx])); break;
                        case TYPENAME::MT_TRACKING: v = static_cast<int>(boost::get<MT_TRACKING>(p_LogRecordValues[idx])); break;
                        case TYPENAME::SN_EVENT: v = static_cast<int>(boost::get<SN_EVENT>(p_LogRecordValues[idx])); break;
                        case TYPENAME::SN_STATE: v = static_cast<int>(boost::get<SN_STATE>(p_LogRecordValues[idx])); break;
                        default: 
                            Squawk("ERROR: Unable to format data to write to HDF5 group for log file " + m_Logfiles[p_LogfileId].name);            // announce error
                            v = 0;
                    }

                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&v);
                }
                else if (dType == H5T_NATIVE_LONG) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<long int>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_USHORT) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<unsigned short int>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_UINT) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<unsigned int>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_ULONG) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<unsigned long int>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_FLOAT) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<float>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_DOUBLE) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<double>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_NATIVE_LDOUBLE) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&boost::get<long double>(p_LogRecordValues[idx]));
                }
                else if (dType == H5T_C_S1) {
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)boost::get<string>(p_LogRecordValues[idx]).c_str());
                }
                else { // assume custom string type
                    ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)boost::get<string>(p_LogRecordValues[idx]).c_str());
                }

    }

    (void)H5Sclose(h5FSpace);                                                                                       // close filespace
    (void)H5Sclose(h5Dspace);                                                                                       // close dataspace
//    (void)H5DClose(h5Dset);                                                                                     // close dataset

//H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dd);


//             herr_t H5Dwrite (hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buf)

////                            m_Logfiles[id].h5file.dataTypes[idx] = fileDetails.typeStrings[idx] == "BOOL"  ? H5T_NATIVE_UCHAR  :                // COMPAS BOOL datatype
////                                                                  (fileDetails.typeStrings[idx] == "INT"   ? H5T_NATIVE_INT    :                // COMPAS INT datatype
////                                                                  (fileDetails.typeStrings[idx] == "FLOAT" ? H5T_NATIVE_DOUBLE :                // COMPAS FLOAT datatype
////                                                                                                             H5T_C_S1));    

if (ok < 0) {
    break;
}
            }

        }
        else {                                                                                                  // no, HDF5 group not open
            Squawk("ERROR: Unable to write to HDF5 group for log file " + m_Logfiles[p_LogfileId].name);            // announce error
            ClearEntry(p_LogfileId);                                                                                // clear entry
        }
    }
    else {                                                                                                          // not enabled or not active
//        Squawk(p_LogStr);                                                                                           // show log record on stderr

//                            boost::variant<string> fmtStr(fileDetails.fmtStrings[index++]);                                             // get format string
//                            valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);                                       // format value
//                            logRecord += valueStr + delimiter;                                                                          // add value string to log record - with delimiter
        ok = 0;
    }

    return (ok >= 0);
}


/*
 * Write a minimally formatted record to the specified log file.
 *
 * If logging is enabled and the specified log file is active, and class and level are enabled (see DoIt()),
 * write string to log file
 *
 * Disable the specified log file if errors occur.
 *
 * Takes parameter p_LogStr, prepends a timestamp and label if required and writes it to the log file.
 *
 * Timestamp format is 'yyyymmdd hh:mm:ss'
 * Label is the logging/debug class associated with the log/debug string
 *
 *
 * bool Put(const id p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log string should be written
 * @param   [IN]    p_LogClass                  Class to determine if string sould be written
 * @param   [IN]    p_LogLevel                  Level to determine if string should be written
 * @param   [IN]    p_LogStr                    The string to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Put(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging service enabled and specified log file active?
        if (DoIt(p_LogClass, p_LogLevel, m_LogClasses, m_LogLevel)) {                                               // yes - logging this class and level?
            result = Put_(p_LogfileId, p_LogStr, p_LogClass);                                                       // yes - log it
        }
    }

    return result;
}


/*
 * Put() a string record to specified log file with no class or level check - internal use only
 * Used for CSV, TSV, and TXT files
 *
 * Timestamps and labels are added here if required.
 *
 * Disable the specified log file if errors occur.
 *
 *
 * bool Put_(const int p_LogfileId, const string p_LogStr, const string p_Label)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log string should be written
 * @param   [IN]    p_LogStr                    The string to be written
 * @param   [IN]    p_Label                     The record label to be written (if required)
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Put_(const int p_LogfileId, const string p_LogStr, const string p_Label) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging service enabled and specified log file active?

        string delimiter = "";                                                                                      // field delimiter
        switch (m_Logfiles[p_LogfileId].filetype) {
            case LOGFILETYPE::HDF5: delimiter = ""; break;                                                          // HDF5
            case LOGFILETYPE::CSV : delimiter = DELIMITERValue.at(DELIMITER::COMMA); break;                         // CSV
            case LOGFILETYPE::TSV : delimiter = DELIMITERValue.at(DELIMITER::TAB); break;                           // TSV
            case LOGFILETYPE::TXT : delimiter = DELIMITERValue.at(DELIMITER::SPACE); break;                         // TXT
            default               : delimiter = ""; break;                                                          // default
        }

        string timestamp;
        if (m_Logfiles[p_LogfileId].timestamp) {                                                                    // timestamp enabled?
            time_t currentTime;                                                                                     // yes - add it
            currentTime = time(NULL);
            tm *now     = localtime(&currentTime);

            timestamp = utils::PadLeadingZeros(std::to_string(1900 + now->tm_year), 4) +                            // year
                        utils::PadLeadingZeros(std::to_string(now->tm_mon + 1    ), 2) +                            // month
                        utils::PadLeadingZeros(std::to_string(now->tm_mday       ), 2) + delimiter +                // day
                        utils::PadLeadingZeros(std::to_string(now->tm_hour       ), 2) + ":" +                      // hour
                        utils::PadLeadingZeros(std::to_string(now->tm_min        ), 2) + ":" +                      // minute
                        utils::PadLeadingZeros(std::to_string(now->tm_sec        ), 2);                             // second
        }

        string logStr = "";                                                                                         // initialise the output string
               logStr += m_Logfiles[p_LogfileId].timestamp ? timestamp + delimiter : "";                            // add timestamp if required
               logStr += m_Logfiles[p_LogfileId].label && p_Label.length() > 0 ? p_Label + delimiter : "";          // add record label if required (may be blank)
               logStr += p_LogStr;                                                                                  // add the log string

        return Write_(p_LogfileId, logStr);                                                                         // log it
    }

    return result;
}


/*
 * Put() a multi-value record to specified log file with no class or level check - internal use only
 * Used for HDF5 files
 *
 * Timestamps and labels should be passed in the p_LogRecordValues vector - they will not be added here.
 *
 * Disable the specified log file if errors occur.
 *
 *
 * bool Put_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const std::vector<TYPENAME> p_LogRecordTypes)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log string should be written
 * @param   [IN]    p_LogRecordValues           Vector of COMPAS_VARIABLE_TYPE values to be written
 * @param   [IN]    p_LogRecordTypes            Vector of TYPENAME values indicating the COMPAS data types of the values to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Put_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const std::vector<TYPENAME> p_LogRecordTypes) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                             // logging service enabled and specified log file active?
        return Write_(p_LogfileId, p_LogRecordValues, p_LogRecordTypes);    // yes - log it
    }

    return result;
}


/*
 * Writes debug string to stdout
 * Also writes debug string to log file if logging is active and so configured (m_DbgToFile via Program Options)
 *
 *
 * bool Debug(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr)
 *
 * @param   [IN]    p_DbgClass                  Class to determine if string should be written
 * @param   [IN]    p_DbgLevel                  Level to determine if string should be written
 * @param   [IN]    p_DbgStr                    The string to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Debug(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr) {

    bool result = false;

    if (m_Enabled) {                                                                                                // logging service enabled?
        if (DoIt(p_DbgClass, p_DbgLevel, m_DbgClasses, m_DbgLevel)) {                                               // debugging this class and level?
            result = Debug_(p_DbgStr);                                                                              // debug it
        }
    }

    return result;
}


/*
 * Debug() with no class or level check - internal use only
 *
 * This is where the work is done
 *
 *
 * bool Debug_(const string p_DbgStr)
 *
 * @param   [IN]    p_DbgStr                    The string to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Debug_(const string p_DbgStr) {

    bool result = false;

    if (m_Enabled) {                                                                                                // logging service enabled?
        Say_("DEBUG: " + p_DbgStr);                                                                                 // to stdout
        result = true;                                                                                              // interim result
        if (m_DbgToLogfile) {                                                                                       // also to log file?
            result = Put_(m_DbgLogfileId, p_DbgStr);                                                                // log it
        }
    }

    return result;
}


/*
 * Writes debug string to stdout and waits for user input before proceeding
 * Also writes debug string to log file if logging is active and so configured (m_DbgToFile via Program Options)
 *
 *
 * bool DebugWait(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr)
 *
 * @param   [IN]    p_DbgClass                  Class to determine if string should be written
 * @param   [IN]    p_DbgLevel                  Level to determine if string should be written
 * @param   [IN]    p_DbgStr                    The string to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::DebugWait(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr) {

    bool result = false;

    if (m_Enabled) {                                                                                                // logging service enabled?
        if (DoIt(p_DbgClass, p_DbgLevel, m_DbgClasses, m_DbgLevel)) {                                               // debugging this class and level?
            result = Debug_(p_DbgStr);                                                                              // debug it
            std::cout << "DEBUG: Press any key to continue...";                                                     // announce
            string tmp; std::cin >> tmp;                                                                            // and wait for input
        }
    }

    return result;
}


/*
 * Writes error string to stderr
 * Also writes error string to log file if logging is active and so configured (m_ErrToFile via Program Options)
 *
 *
 * bool Error(const string p_ErrStr)
 *
 * @param   [IN]    p_DbgStr                    The string to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Error(const string p_ErrStr) {

    bool result = true;

    Squawk(p_ErrStr);                                                                                               // don't need logging enabled to squawk error

    if (m_ErrToLogfile && m_Enabled) {                                                                              // but do need it enabled to write error to logfile
        result = Put_(m_ErrLogfileId, p_ErrStr);                                                                    // log it
    }

    return result;
}


/*
 * The following functions are effectively a wrapper around the base logging functions.
 * These functions are higher-level functions for logging the various SSE and BSE records
 * to COMPAS "standrad log files" (pre-defined logfiles with pre-defined default record definitions)
 *
 * JR: todo: complete this documentation
 */


/*
 * Get property details for STELLAR (e.g. individual rather than binary stars) properties.
 *
 * These are any of:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *
 * Property details are:
 *
 * { TYPENAME, Header string, Units string, fields width, precision }
 *
 *
 * PROPERTY_DETAILS StellarPropertyDetails(ANY_STAR_PROPERTY p_Property)
 *
 * @param   [IN]    p_Property                  The property for which the details are required
 * @return                                      Tuple containing the properties (default properties if p_Property not found)
 */
PROPERTY_DETAILS Log::StellarPropertyDetails(ANY_STAR_PROPERTY p_Property) {

    PROPERTY_DETAILS details;
    try { details = ANY_STAR_PROPERTY_DETAIL.at(p_Property); }          // get stellar property details
    catch (const std::exception& e) {                                   // unknown property
        details = std::make_tuple(TYPENAME::NONE, "", "", 0, 0);        // empty details
        DBG_WARN(ERR_MSG(ERROR::UNKNOWN_STELLAR_PROPERTY));             // show warning
    }

    return details;
}


/*
 * Get property details for BINARY (e.g. binary rather than individual stars) properties.
 *
 * Property details are:
 *
 * { TYPENAME, Header string, Units string, fields width, precision }
 *
 *
 * PROPERTY_DETAILS BinaryPropertyDetails(BINARY_PROPERTY p_Property)
 *
 * @param   [IN]    p_Property                  The property for which the details are required
 * @return                                      Tuple containing the properties (default properties if p_Property not found)
 */
PROPERTY_DETAILS Log::BinaryPropertyDetails(BINARY_PROPERTY p_Property) {

    PROPERTY_DETAILS details;

    try { details = BINARY_PROPERTY_DETAIL.at(p_Property); }            // get binary property details
    catch (const std::exception& e) {                                   // unknown property
        details = std::make_tuple(TYPENAME::NONE, "", "", 0, 0);        // empty details
        DBG_WARN(ERR_MSG(ERROR::UNKNOWN_BINARY_PROPERTY));              // show warning
    }

    return details;
}


/*
 * Get property details for PROGRAM OPTIONS
 *
 * Property details are:
 *
 * { TYPENAME, Header string, Units string, fields width, precision }
 *
 *
 * PROPERTY_DETAILS ProgramOptionDetails(ProgramOptionDetails p_Property)
 *
 * @param   [IN]    p_Property                  The property for which the details are required
 * @return                                      Tuple containing the properties (default properties if p_Property not found)
 */
PROPERTY_DETAILS Log::ProgramOptionDetails(PROGRAM_OPTION p_Property) {

    PROPERTY_DETAILS details;

    try { details = PROGRAM_OPTION_DETAIL.at(p_Property); }             // get program option details
    catch (const std::exception& e) {                                   // unknown property
        details = std::make_tuple(TYPENAME::NONE, "", "", 0, 0);        // empty details
        DBG_WARN(ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION));               // show warning
    }

    return details;
}


/*
 * Format field header strings (header, units, type, format)
 *
 * This function takes the property details, and a suffix string and assembles strings to be printed
 * as the headers for the property described by the property details.
 *
 * The content of the header strings is a combination of the property details and the suffix string supplied.
 * The final field width for the property is determined here, and is calculated based upon the widths of the
 * header deatils and the configured field withd of the property (part of the property details passed as a
 * parameter).
 *
 * The strings formatted are:
 *
 *    Header - The header string is the title of the field: it describes the contents (e.g. "Metallicity")
 *             The header string will have the header suffix string appended to it.  The suffix string
 *             provides differentiation for the constituent stars of a binary.  For example, the suffix
 *             could be "(1)" or "(2)" to indicate the primary and secondary stars, or "(SN)" or "(CN)" to
 *             indicated the supernova or companion stars in a supernova event.
 *
 *    Units  - the units string indicates the units of the property.  This string is taken directly from
 *             the property details passed as a parameter.
 *
 *    Type   - the type string indicates the data type of the property.  This will be the "short name" of
 *             the property data type, retrieved from the TYPENAME_LABEL enum defined in constants.h for the
 *             property.
 *
 *    Format - the format string is used by logging code to format the property value.  The format string is
 *             constructed he because the final field width is determined here.
 *
 *
 * STR_STR_STR_STR FormatFieldHeaders(PROPERTY_DETAILS p_PropertyDetails, string p_HeaderSuffix)
 *
 * @param   [IN]    p_PropertyDetails           The property details for the property for which the headers are to be formatted
 * @param   [IN]    p_HeaderSuffix              The suffix string to be appended to the header string
 * @return                                      Tuple containing formatted strings for the property requested: <header, units, type, format>
 *                                              If the property details are not valid (e.g. unknown data type), error strings will be returned
 */
STR_STR_STR_STR Log::FormatFieldHeaders(PROPERTY_DETAILS p_PropertyDetails, string p_HeaderSuffix) {

    TYPENAME typeName = std::get<0>(p_PropertyDetails);                                                             // data type
    if (typeName == TYPENAME::NONE) {                                                                               // valid data type?
        return std::make_tuple("ERROR!", "ERROR!", "ERROR!", "ERROR!");                                             // return error values
    }

    string headerStr = std::get<1>(p_PropertyDetails) + p_HeaderSuffix;                                             // header string
    string unitsStr  = std::get<2>(p_PropertyDetails);                                                              // units string
    string typeStr   = std::get<1>(TYPENAME_LABEL.at(typeName));                                                    // type will be one of "BOOL", "INT", "FLOAT" and "STRING" (non-primitive types coerced to INT)

    int fieldWidth     = std::get<3>(p_PropertyDetails);
    int fieldPrecision = std::get<4>(p_PropertyDetails);                                                            // field precision (for double and int)

    fieldWidth = std::max(fieldWidth, std::max((int)headerStr.length(), std::max((int)unitsStr.length(), 6)));      // field width - maximum of requested width, header width, units width and type width ("STRING" is max type)

    headerStr = utils::CentreJustify(headerStr, fieldWidth);                                                        // centre-justify header string
    unitsStr  = utils::CentreJustify(unitsStr, fieldWidth);                                                         // centre-justify units string
    typeStr   = utils::CentreJustify(typeStr, fieldWidth);                                                          // centre-justify type string

    string fmtStr = std::to_string(fieldWidth) + "." + std::to_string(fieldPrecision);                              // field width and precision spcifiers

    return std::make_tuple(headerStr, unitsStr, typeStr, fmtStr);
}


/*
 * Find a value in the LOGFILE_DESCRIPTOR map and return the key if found, otherwise defaut value
 *
 * This function looks for the passed string value in the LOGFILE_DESCRIPTOR map, and if the string
 * is found returns the key correspoding to the value found.  If the value is not found the default
 * value LOGFILE::NONE is returned.
 *
 * The string comparison is case-insensitive.
 *
 *
 * static std::tuple<bool, LOGFILE> GetLogfileDescriptorKey(const string p_Value)
 *
 * @param   [IN]    p_Value                     The value to be located in the LOGFILE_DESCRIPTOR map
 * @return                                      Tuple containing a boolean result (true if value found, else false), and the key
 *                                              corresponding to the value found, or LOGFILE::NONE if the value was not found
 */
std::tuple<bool, LOGFILE> Log::GetLogfileDescriptorKey(const string p_Value) {
    for (auto& it: LOGFILE_DESCRIPTOR)
        if (utils::Equals(std::get<3>(it.second), p_Value)) return std::make_tuple(true, it.first);
    return std::make_tuple(false, LOGFILE::NONE);
}


/*
 * Find a value in the m_OpenStandardLogFileIds map and return the key if found, otherwise defaut value
 *
 * This function looks for the file id value in the m_OpenStandardLogFileIds map, and if the id is found
 * returns the key correspoding to the id found.  If the value is not found the default value LOGFILE::NONE
 * is returned.
 *
 *
 * static std::tuple<bool, LOGFILE> GetStandardLogfileKey(const int p_FileId)
 *
 * @param   [IN]    p_Value                     The value to be located in the LOGFILE_DESCRIPTOR map
 * @return                                      Tuple containing a boolean result (true if id found, else false), and the key
 *                                              corresponding to the id found, or LOGFILE::NONE if the value was not found
 */
std::tuple<bool, LOGFILE> Log::GetStandardLogfileKey(const int p_FileId) {
    for (auto& it: m_OpenStandardLogFileIds)
        if (it.second.id == p_FileId) return std::make_tuple(true, it.first);
    return std::make_tuple(false, LOGFILE::NONE);
}


/*
 * Get standard log file record properties and format vector from the logfile record specifier
 *
 * This function is a very reduced version of Log::StandardLogFileDetails(), and exists mainly
 * to support writing to the (new) SSE Supernova logfile (it was written specifically for that
 * purpose, but was left general enough to retrieve the properties and format vector for any
 * of the logfiles).
 * 
 * The reason this function is needed is that because we (currently) save the state of a
 * single star in SSE and revert to the previous state if we find we've evolved too far and
 * possibly missed something interesting, we can't write a record to the SSE supernova file
 * at the time it occurs because if we revert to the previous state after writing the record
 * we can't (easily) unwind the write - and the logging code is written to only create a log
 * file on the first write to the file (so that we don't have files created that have no
 * records other than the header if we never write to them), so even if we could unwind the
 * write, we might have a file created that may never have data records written to it (we
 * could specifically check for that and delete the file, but that's a little inelegant -
 * better to not create the file in the first place).  So, this function will enable me to 
 * format a SSE Supernova record at the right time, but delay writing it to after we decide 
 * that we'll accept the current state and not revert.
 * 
 * With hindsight, Log::StandardLogFileDetails() should probably have been written with this
 * part separated out.  Ideally Log::StandardLogFileDetails() would call this function to get
 * these details - that way if we ever need to change how this is done we don't need to change
 * it in two places.  But Log::StandardLogFileDetails() is, the way it was initially written,
 * a bit too complex (this whole flexible printing code is a complex beast - unfortunately it 
 * has to be to get it to work) and this code is a bit too intertwined to easily and quickly
 * disentangle it from Log::StandardLogFileDetails() - that's probably a good code cleanup to 
 * do some time in the future, but for now this will have to suffice.
 *
 *
 * std::tuple<ANY_PROPERTY_VECTOR, std::vector<string>> StandardLogFileRecordDetails(const LOGFILE p_Logfile)
 *
 * @param   [IN]    p_Logfile                   Logfile for which details are to be retrieved (see enum class LOGFILE in constants.h)
 * @return                                      Tuple containing:
 *                                                 - a vector of logfile record properties
 *                                                 - a vector of format strings
 */
std::tuple<ANY_PROPERTY_VECTOR, std::vector<string>> Log::GetStandardLogFileRecordDetails(const LOGFILE p_Logfile) {

    ANY_PROPERTY_VECTOR  recordProperties = {};                                                                                     // default is empty
    std::vector<string>  fmtVector = {};                                                                                            // default is empty

    try {
        // get record properties for this file

        switch (p_Logfile) {                                                                                                        // which logfile?

            case LOGFILE::BSE_BE_BINARIES:                                                                                          // BSE_BE_BINARIES
                recordProperties = m_BSE_BE_Binaries_Rec;                                                                           // record properties
                break;

            case LOGFILE::BSE_COMMON_ENVELOPES:                                                                                     // BSE_COMMON_ENVELOPES
                recordProperties = m_BSE_CEE_Rec;                                                                                   // record properties
                break;

            case LOGFILE::BSE_DETAILED_OUTPUT:                                                                                      // BSE_DETAILED_OUTPUT
                recordProperties = m_BSE_Detailed_Rec;                                                                              // record properties
                break;

            case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS:                                                                               // BSE_DOUBLE_COMPACT_OBJECTS
                recordProperties = m_BSE_DCO_Rec;                                                                                   // record properties
                break;

            case LOGFILE::BSE_PULSAR_EVOLUTION:                                                                                     // BSE_PULSAR_EVOLUTION
                recordProperties = m_BSE_Pulsars_Rec;                                                                               // record properties
                break;

            case LOGFILE::BSE_RLOF_PARAMETERS:                                                                                      // BSE_RLOF_PARAMETERS
                recordProperties = m_BSE_RLOF_Rec;                                                                                  // record properties
                break;

            case LOGFILE::BSE_SUPERNOVAE:                                                                                           // BSE_SUPERNOVAE
                recordProperties = m_BSE_SNE_Rec;                                                                                   // record properties
                break;

            case LOGFILE::BSE_SWITCH_LOG:                                                                                           // BSE_SWITCH_LOG
                recordProperties = m_BSE_Switch_Rec;                                                                                // record properties
                break;

            case LOGFILE::BSE_SYSTEM_PARAMETERS:                                                                                    // BSE_SYSTEM_PARAMETERS
                recordProperties = m_BSE_SysParms_Rec;                                                                              // record properties
                break;

            case LOGFILE::SSE_DETAILED_OUTPUT:                                                                                      // SSE_DETAILED_OUTPUT
                recordProperties = m_SSE_Detailed_Rec;                                                                              // record properties
                break;

            case LOGFILE::SSE_SUPERNOVAE:                                                                                           // SSE_SUPERNOVAE
                recordProperties = m_SSE_SNE_Rec;                                                                                   // record properties
                break;

            case LOGFILE::SSE_SWITCH_LOG:                                                                                           // SSE_SWITCH_LOG
                recordProperties = m_SSE_Switch_Rec;                                                                                // record properties
                break;

            case LOGFILE::SSE_SYSTEM_PARAMETERS:                                                                                    // SSE_SYSTEM_PARAMETERS
                recordProperties = m_SSE_SysParms_Rec;                                                                              // record properties
                break;

            default:                                                                                                                // unknown logfile
                recordProperties = {};                                                                                              // no record properties
        }

        if (!recordProperties.empty()) {                                                                                            // have properties?

            // get field format strings

            bool ok = true;                                                                                                         // ok so far...

            for (auto &property : recordProperties) {                                                                               // for each property to be included in the log record

                ANY_PROPERTY_TYPE propertyType = boost::apply_visitor(VariantPropertyType(), property);                             // property type
                            
                string fmtStr = "";

                switch (propertyType) {                                                                                             // which property type?

                    case ANY_PROPERTY_TYPE::T_STAR_PROPERTY: {                                                                      // single star
                        ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<STAR_PROPERTY>(property));        // property
                        PROPERTY_DETAILS details = StellarPropertyDetails(anyStarProp);                                             // property details
                        std::tie(std::ignore, std::ignore, std::ignore, fmtStr) = FormatFieldHeaders(details);                      // get format string
                        } break;

                    case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY: {                                                                    // star 1 of binary
                        ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<STAR_1_PROPERTY>(property));      // property
                        PROPERTY_DETAILS details = StellarPropertyDetails(anyStarProp);                                             // property details
                        std::tie(std::ignore, std::ignore, std::ignore, fmtStr) = FormatFieldHeaders(details, "(1)");               // get format string
                        } break;

                    case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY: {                                                                    // star 2 of binary
                        ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<STAR_2_PROPERTY>(property));      // property
                        PROPERTY_DETAILS details = StellarPropertyDetails(anyStarProp);                                             // property details
                        std::tie(std::ignore, std::ignore, std::ignore, fmtStr) = FormatFieldHeaders(details, "(2)");               // get format string
                        } break;

                    case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: {                                                                 // supernova star of binary that contains a supernova
                        ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<SUPERNOVA_PROPERTY>(property));   // property
                        PROPERTY_DETAILS details = StellarPropertyDetails(anyStarProp);                                             // property details
                        std::tie(std::ignore, std::ignore, std::ignore, fmtStr) = FormatFieldHeaders(details, "(SN)");              // get format string
                        } break;

                    case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: {                                                                 // companion star of binary that contains a supernova
                        ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<COMPANION_PROPERTY>(property));   // property
                        PROPERTY_DETAILS details = StellarPropertyDetails(anyStarProp);                                             // property details
                        std::tie(std::ignore, std::ignore, std::ignore, fmtStr) = FormatFieldHeaders(details, "(CP)");              // get format string
                        } break;

                    case ANY_PROPERTY_TYPE::T_BINARY_PROPERTY: {                                                                    // binary
                        BINARY_PROPERTY binaryProp = boost::get<BINARY_PROPERTY>(property);                                         // property
                        PROPERTY_DETAILS details = BinaryPropertyDetails(binaryProp);                                               // property details
                        std::tie(std::ignore, std::ignore, std::ignore, fmtStr) = FormatFieldHeaders(details);                      // get format string
                        } break;

                    case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION: {                                                                     // program option
                        PROGRAM_OPTION programOption = boost::get<PROGRAM_OPTION>(property);                                        // property
                        PROPERTY_DETAILS details = ProgramOptionDetails(programOption);                                             // property details
                        std::tie(std::ignore, std::ignore, std::ignore, fmtStr) = FormatFieldHeaders(details);                      // get format string
                        } break;

                    default:                                                                                                        // unknown property type
                        ok = false;                                                                                                 // that's not ok...
                }

                if (ok) {
                    fmtVector.push_back(fmtStr);                                                                                    // record format string for field
                }
            }

            if (!ok) {                                                                                                              // have format vectr ok?
                fmtVector = {};                                                                                                     // no format vector
            }
        }
    }
    catch (const std::exception& e) {                                                                                               // oops...
        recordProperties = {};                                                                                                      // no record properties
        fmtVector = {};                                                                                                             // no format vector
    }

    return std::make_tuple(recordProperties, fmtVector);
}


/*
 * Get standard log file details and open file if necessary
 *
 * This function will retrieve the details for the logfile specified, and open the file if it
 * is not already open.
 *
 * The function first checks map of currently open standard logfiles (m_OpenStandardLogFileIds),
 * and if the logfile is found just returns the details recorded there.
 *
 * If the logfile is not found in m_OpenStandardLogFileIds, the logfile details are assembled
 * by retrieving the record specifier for the logfile (default or user-specified) and formatting
 * the header strings according to the contents of the record specifier.  A new log file is then
 * opened (with version number appended if a file of the same name already exists), and the header
 * strings are written to the file.  The logfile name is constructed using the LOGFILE_DESCRIPTOR
 * enum from constants.h, and the p_FileSuffix parameter.  The file suffix (p_FileSuffix) can be
 * used to indicate, for example, the ordinal number of the star for which information is being
 * logged (for example, when loopong through a population of binary star, p_FileSuffix would indicate
 * the loop index of the star for which information is being logged). The file remains open until
 * explicitly closed (by calling CloseStandardFile(), possibly via CloseAllStandardFiles().
 *
 * The logfile details are returned.
 *
 *
 * LogfileDetailsT StandardLogFileDetails(const LOGFILE p_Logfile, const string p_FileSuffix)
 *
 * @param   [IN]    p_Logfile                   Logfile for which details are to be retrieved (see enum class LOGFILE in constants.h)
 * @param   [IN]    p_Suffix                    String suffix to be appended to the logfile name
 * @return                                      Struct with logfile details - see typedefs.h
 */
LogfileDetailsT Log::StandardLogFileDetails(const LOGFILE p_Logfile, const string p_FileSuffix) {

    bool                 ok = true;
    LogfileDetailsT      retVal = {-1, "", {}, {}, {}, {}, {}, {}};                                                                             // default return value

    LogfileDetailsT      fileDetails = retVal;                                                                                                  // logfile details
    LOGFILE_DESCRIPTOR_T fileDescriptor;                                                                                                        // logfile descriptor

    COMPASUnorderedMap<LOGFILE, LogfileDetailsT>::const_iterator logfile;                                                                       // iterator
    logfile = m_OpenStandardLogFileIds.find(p_Logfile);                                                                                         // look for open logfile
    if (logfile == m_OpenStandardLogFileIds.end()) {                                                                                            // doesn't exist
        try {                                                                                                                                   // get record properties for this file
            switch (p_Logfile) {                                                                                                                // which logfile?

                case LOGFILE::BSE_BE_BINARIES:                                                                                                  // BSE_BE_BINARIES
                    fileDetails.filename         = OPTIONS->LogfileBeBinaries();
                    fileDetails.recordProperties = m_BSE_BE_Binaries_Rec;
                    break;

                case LOGFILE::BSE_COMMON_ENVELOPES:                                                                                             // BSE_COMMON_ENVELOPES
                    fileDetails.filename         = OPTIONS->LogfileCommonEnvelopes();
                    fileDetails.recordProperties = m_BSE_CEE_Rec;
                    break;

                case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS:                                                                                       // BSE_DOUBLE_COMPACT_OBJECTS
                    fileDetails.filename         = OPTIONS->LogfileDoubleCompactObjects();
                    fileDetails.recordProperties = m_BSE_DCO_Rec;
                    break;
               
                case LOGFILE::BSE_PULSAR_EVOLUTION:                                                                                             // BSE_PULSAR_EVOLUTION
                    fileDetails.filename         = OPTIONS->LogfilePulsarEvolution();
                    fileDetails.recordProperties = m_BSE_Pulsars_Rec;
                    break;

                case LOGFILE::BSE_RLOF_PARAMETERS:                                                                                              // BSE_RLOF_PARAMETERS
                    fileDetails.filename         = OPTIONS->LogfileRLOFParameters();
                    fileDetails.recordProperties = m_BSE_RLOF_Rec;
                    break;

                case LOGFILE::BSE_SUPERNOVAE:                                                                                                   // BSE_SUPERNOVAE
                    fileDetails.filename         = OPTIONS->LogfileSupernovae();
                    fileDetails.recordProperties = m_BSE_SNE_Rec;
                    break;

                case LOGFILE::BSE_SWITCH_LOG:                                                                                                   // BSE_SWITCH_LOG
                    fileDetails.filename         = OPTIONS->LogfileSwitchLog();
                    fileDetails.recordProperties = m_BSE_Switch_Rec;
                    break;

                case LOGFILE::BSE_SYSTEM_PARAMETERS:                                                                                            // BSE_SYSTEM_PARAMETERS
                    fileDetails.filename         = OPTIONS->LogfileSystemParameters();
                    fileDetails.recordProperties = m_BSE_SysParms_Rec;
                    break;

                case LOGFILE::SSE_SUPERNOVAE:                                                                                                   // SSE_SUPERNOVAE
                    fileDetails.filename         = OPTIONS->LogfileSupernovae();
                    fileDetails.recordProperties = m_SSE_SNE_Rec;
                    break;

                case LOGFILE::SSE_SWITCH_LOG:                                                                                                   // SSE_SWITCH_LOG
                    fileDetails.filename         = OPTIONS->LogfileSwitchLog();
                    fileDetails.recordProperties = m_SSE_Switch_Rec;
                    break;

                case LOGFILE::SSE_SYSTEM_PARAMETERS:                                                                                            // SSE_SYSTEM_PARAMETERS
                    fileDetails.filename         = OPTIONS->LogfileSystemParameters();
                    fileDetails.recordProperties = m_SSE_SysParms_Rec;
                    break;

                case LOGFILE::SSE_DETAILED_OUTPUT:                                                                                              // SSE_DETAILED_OUTPUT
                case LOGFILE::BSE_DETAILED_OUTPUT: {                                                                                            // BSE_DETAILED_OUTPUT

                    // first check if the detailed output directory exists - if not, create it
                    // use boost filesystem here - easier...

                    bool detailedOutputDirectoryExists = false;                                                                                 // detailed output directory exists?  Start with no

                    string detailedDirName = m_LogBasePath + "/" + m_LogContainerName + "/" + DETAILED_OUTPUT_DIRECTORY_NAME;                   // directory name with path ("/" works on Uni*x and Windows)

                    if (boost::filesystem::exists(detailedDirName)) {                                                                           // directory already exists?
                        detailedOutputDirectoryExists = true;                                                                                   // yes
                    }
                    else {                                                                                                                      // no - create directory
                        boost::system::error_code err;
                        try {
                            boost::filesystem::create_directory(detailedDirName, err);                                                          // create container - let boost throw an exception if it fails
                            if (err.value() == 0) {                                                                                             // ok?
                                detailedOutputDirectoryExists = true;                                                                           // yes
                            }
                            else  {                                                                                                             // not ok...
                                Squawk("ERROR: Unable to create detailed output directory " + detailedDirName);                                 // announce error
                                Squawk("Boost filesystem error = " + err.message());                                                            // plus details
                                Squawk("Detailed Output logging disabled");                                                                     // show disabled warning
                            }
                        }
                        catch (...) {                                                                                                           // unhandled problem...
                                Squawk("ERROR: Unable to create detailed output directory " + detailedDirName);                                 // announce error
                                Squawk("Detailed Output logging disabled");                                                                     // show disabled warning
                        }
                    }

                    if (detailedOutputDirectoryExists) {                                                                                        // detailed output directory exists?
                                                                                                                                                // yes - add path to filename
                        switch (p_Logfile) {                                                                                                    // which logfile?

                            case LOGFILE::SSE_DETAILED_OUTPUT:                                                                                  // SSE_DETAILED_OUTPUT
                                fileDetails.filename         = DETAILED_OUTPUT_DIRECTORY_NAME + "/" + OPTIONS->LogfileDetailedOutput();         // logfile filename with directory
                                fileDetails.recordProperties = m_SSE_Detailed_Rec;                                                              // record properties
                                break;

                            case LOGFILE::BSE_DETAILED_OUTPUT:                                                                                  // BSE_DETAILED_OUTPUT
                                fileDetails.filename         = DETAILED_OUTPUT_DIRECTORY_NAME + "/" + OPTIONS->LogfileDetailedOutput();         // logfile filename with directory
                                fileDetails.recordProperties = m_BSE_Detailed_Rec;                                                              // record properties
                                break;

                            default: break;
                       }
                    }
                    } break;

                default:                                                                                                                        // unknown logfile      
                    fileDetails.filename = ""; 
                    fileDetails.recordProperties = {};
                    Squawk(ERR_MSG(ERROR::UNKNOWN_LOGFILE) + ": Logging disabled for this file");                                               // announce error
            }

            if (!fileDetails.filename.empty() && !fileDetails.recordProperties.empty()) {                                                       // have filename and properties?

                fileDetails.filename += p_FileSuffix;                                                                                           // add suffix to filename

                // if we're logging to HDF5 files:
                //    - all logfiles except detailed output files are included in a single HDF5 file
                //    - detailed output files are created individually as HDF5 files inside their containing directory
                if (m_LogfileType == LOGFILETYPE::HDF5) {                                                                                       // logging to HDF5 files?
                                                                                                                                                // yes
                    string fileExt = "." + LOGFILETYPEFileExt.at(OPTIONS->LogfileType());                                                       // file extension for HDF5 files

                    if (p_Logfile == LOGFILE::SSE_DETAILED_OUTPUT || p_Logfile == LOGFILE::BSE_DETAILED_OUTPUT) {                               // yes - detailed output file (SSE or BSE)?
                        if (m_HDF5DetailedId < 0) {                                                                                             // have HDF5 detailed file?
                                                                                                                                                // no - create it
                            string h5Filename = m_LogBasePath + "/" + m_LogContainerName + "/" + fileDetails.filename + fileExt;                // full filename with path, container, and extension ("/" works on Uni*x and Windows)
                            m_HDF5DetailedId = H5Fcreate(h5Filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);                           // create HDF5 detailed file
                            if (m_HDF5DetailedId < 0) {                                                                                         // created ok?                        
                                Squawk("ERROR: Unable to create HDF5 detailed file with file name " + h5Filename);                              // no - announce error
                                Squawk("Logging disabled");                                                                                     // show disabled warning
                                m_Enabled = false;                                                                                              // disable logging
                                ok = false;                                                                                                     // fail
                            }
                        }
                    }
                    else {                                                                                                                      // no - not detailed output file
                        if (m_HDF5ContainerId < 0) {                                                                                            // have HDF5 container?
                                                                                                                                                // no - create it
                            string h5Filename = m_LogBasePath + "/" + m_LogContainerName + "/" + m_LogContainerName + fileExt;                  // full filename with path, container, and extension ("/" works on Uni*x and Windows)
                            m_HDF5ContainerId = H5Fcreate(h5Filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);                          // create HDF5 container file
                            if (m_HDF5ContainerId < 0) {                                                                                        // created ok?                        
                                Squawk("ERROR: Unable to create HDF5 container file with file name " + h5Filename);                             // no - announce error
                                Squawk("Logging disabled");                                                                                     // show disabled warning
                                m_Enabled = false;                                                                                              // disable logging
                                ok = false;                                                                                                     // fail
                            }
                        }
                    }
                }

                // open the specific logfile - if logging is enabled
                // (may have been disabled by HDF5 container open fail)
                if (m_Enabled && ok) {                                                                                                          // logging enabled?
                    if ((fileDetails.id = Open(fileDetails.filename, false, false, false, p_Logfile)) < 0) {                                    // yes - open the log file - new file, no timestamps, no record labels (all same type here)
                        Squawk(ERR_MSG(ERROR::FILE_OPEN_ERROR) + ": Logging disabled for this file");                                           // announce error if open failed
                        ok = false;                                                                                                             // fail
                    }
                }

                if (ok) {                                                                                                                       // ok?
                                                                                                                                                // yes
                    string delimiter = "";                                                                                                      // field delimiter
                    switch (m_Logfiles[fileDetails.id].filetype) {
                        case LOGFILETYPE::HDF5: delimiter = ""; break;                                                                          // HDF5
                        case LOGFILETYPE::CSV : delimiter = DELIMITERValue.at(DELIMITER::COMMA); break;                                         // CSV
                        case LOGFILETYPE::TSV : delimiter = DELIMITERValue.at(DELIMITER::TAB); break;                                           // TSV
                        case LOGFILETYPE::TXT : delimiter = DELIMITERValue.at(DELIMITER::SPACE); break;                                         // TXT
                        default               : delimiter = ""; break;                                                                          // default
                    }

                    // get and format field headers for printing; get field format strings
                    if (!fileDetails.recordProperties.empty()) {

                        for (auto &property : fileDetails.recordProperties) {                                                                   // for each property to be included in the log record

                            string headerStr = "";
                            string unitsStr  = "";
                            string typeStr   = "";
                            string fmtStr    = "";

                            bool ok = true;

                            ANY_PROPERTY_TYPE propertyType = boost::apply_visitor(VariantPropertyType(), property);                             // property type
                            PROPERTY_DETAILS  details;                                                                                          // property details

                            switch (propertyType) {                                                                                             // which property type?

                                case ANY_PROPERTY_TYPE::T_STAR_PROPERTY: {                                                                      // single star
                                    ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<STAR_PROPERTY>(property));        // property
                                    details = StellarPropertyDetails(anyStarProp);                                                              // property details
                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details);                               // format the headers
                                    } break;

                                case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY: {                                                                    // star 1 of binary
                                    ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<STAR_1_PROPERTY>(property));      // property
                                    details = StellarPropertyDetails(anyStarProp);                                                              // property details
                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details, "(1)");                        // format the headers
                                    } break;

                                case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY: {                                                                    // star 2 of binary
                                    ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<STAR_2_PROPERTY>(property));      // property
                                    details = StellarPropertyDetails(anyStarProp);                                                              // property details
                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details, "(2)");                        // format the headers
                                    } break;

                                case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: {                                                                 // supernova star of binary that contains a supernova
                                    ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<SUPERNOVA_PROPERTY>(property));   // property
                                    details = StellarPropertyDetails(anyStarProp);                                                              // property details
                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details, "(SN)");                       // format the headers
                                    } break;

                                case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: {                                                                 // companion star of binary that contains a supernova
                                    ANY_STAR_PROPERTY anyStarProp = static_cast<ANY_STAR_PROPERTY>(boost::get<COMPANION_PROPERTY>(property));   // property
                                    details = StellarPropertyDetails(anyStarProp);                                                              // property details
                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details, "(CP)");                       // format the headers
                                    } break;

                                case ANY_PROPERTY_TYPE::T_BINARY_PROPERTY: {                                                                    // binary
                                    BINARY_PROPERTY binaryProp = boost::get<BINARY_PROPERTY>(property);                                         // property
                                    details = BinaryPropertyDetails(binaryProp);                                                                // property details
                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details);                               // format the headers
                                    } break;

                                case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION: {                                                                     // program option
                                    PROGRAM_OPTION programOption = boost::get<PROGRAM_OPTION>(property);                                        // property
                                    details = ProgramOptionDetails(programOption);                                                              // property details
                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details);                               // format the headers
                                    } break;

                                default:                                                                                                        // unknown property type
                                    ok = false;                                                                                                 // that's not ok...
                                    Squawk(ERR_MSG(ERROR::UNKNOWN_PROPERTY_TYPE));                                                              // show warning
                            }

                            if (ok) {
                                fileDetails.propertyTypes.push_back(std::get<0>(details));                                                      // append property typename
                                fileDetails.hdrStrings.push_back(headerStr);                                                                    // append header string for field
                                fileDetails.unitsStrings.push_back(unitsStr);                                                                   // append units string for field
                                fileDetails.typeStrings.push_back(typeStr);                                                                     // append type string for field
                                fileDetails.fmtStrings.push_back(fmtStr);                                                                       // append format string for field
                            }
                        }

                        // if we are writing to the SSE Switch file we add two pre-defined columns
                        // to the end of the log record.  These are:
                        //
                        // ( i) the steller type from which the star is switching
                        // (ii) the stellar type to which the star is switching
                        //
                        // if we are writing to the BSE Switch file we add three pre-defined columns
                        // to the end of the log record.  These are:
                        //
                        // (  i) the star switching - 1 = primary, 2 = secondary
                        // ( ii) the steller type from which the star is switching
                        // (iii) the stellar type to which the star is switching

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// JR FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        if (p_Logfile == LOGFILE::BSE_SWITCH_LOG) {                                                                             // BSE Switch Log
                            fileDetails.hdrStrings.push_back("STAR_SWITCHING");                                                                 // append header string for field
                            fileDetails.unitsStrings.push_back("-");                                                                            // append units string for field
                            fileDetails.typeStrings.push_back("INT");                                                                           // append type string for field
                            fileDetails.fmtStrings.push_back("4.1");                                                                            // append fromat string for field
                        }

                        if (p_Logfile == LOGFILE::BSE_SWITCH_LOG || p_Logfile == LOGFILE::SSE_SWITCH_LOG) {                                     // BSE Switch Log or SSE Switch Log
                            fileDetails.hdrStrings.push_back("SWITCHING_FROM");                                                                 // append header string for field
                            fileDetails.hdrStrings.push_back("SWITCHING_TO");                                                                   // append header string for field

                            fileDetails.unitsStrings.push_back("-");                                                                            // append units string for field
                            fileDetails.unitsStrings.push_back("-");                                                                            // append units string for field

                            fileDetails.typeStrings.push_back("INT");                                                                           // append type string for field
                            fileDetails.typeStrings.push_back("INT");                                                                           // append type string for field

                            fileDetails.fmtStrings.push_back("4.1");                                                                            // append fromat string for field
                            fileDetails.fmtStrings.push_back("4.1");                                                                            // append format string for field
                        }
                    }
                    
                    // record new open file details
                    m_OpenStandardLogFileIds.insert({p_Logfile, fileDetails});                                                                  // record the new file details and format strings

                    // initialise files:
                    //    - create datasets for HDF5 files
                    //    - write header/units/types strings for CSV/TSV/TXT files
                    if (OPTIONS->LogfileType() == LOGFILETYPE::HDF5) {                                                                          // logging to HDF5 files?
                        for (size_t idx = 0; idx < fileDetails.hdrStrings.size(); idx++) {                                                      // for each property
                            // set HDF5 datatype
                            switch (fileDetails.propertyTypes[idx]) {                                                                           // which type?
                                case TYPENAME::BOOL        : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_UCHAR); break;
                                case TYPENAME::SHORTINT    : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_SHORT); break;
                                case TYPENAME::INT         : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_INT); break;
                                case TYPENAME::LONGINT     : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_LONG); break;
                                case TYPENAME::USHORTINT   : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_USHORT); break;
                                case TYPENAME::UINT        : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_UINT); break;
                                case TYPENAME::ULONGINT    : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_ULONG); break;
                                case TYPENAME::FLOAT       : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_FLOAT); break;
                                case TYPENAME::DOUBLE      : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_DOUBLE); break;
                                case TYPENAME::LONGDOUBLE  : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_LDOUBLE); break;
                                case TYPENAME::OBJECT_ID   : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_ULONG); break;
                                case TYPENAME::ERROR       : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_INT); break;
                                case TYPENAME::STELLAR_TYPE: m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_INT); break;
                                case TYPENAME::MT_CASE     : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_INT); break;
                                case TYPENAME::MT_TRACKING : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_INT); break;
                                case TYPENAME::SN_EVENT    : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_INT); break;
                                case TYPENAME::SN_STATE    : m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(H5T_NATIVE_INT); break;
                                case TYPENAME::STRING      : {
                                        hid_t h5DType = H5Tcopy(H5T_C_S1);                                                                      // HDF5 c-string datatype
                                        (void)H5Tset_size(h5DType, (int)std::stod(fileDetails.fmtStrings[idx]) + 1);                            // size is field width + 1 (for NULL terminator)
                                        (void)H5Tset_cset(h5DType, H5T_CSET_ASCII);                                                             // ASCII (rather than UTF-8)
                                        m_Logfiles[fileDetails.id].h5file.dataTypes.push_back(h5DType);
                                    } break;
                                default:                                                                                                        // unknown property type
Squawk(ERR_MSG(ERROR::UNKNOWN_DATA_TYPE));                                                                  // announce error
ok = false; // fail
                                    ok = false;                                                                                                 // that's not ok...
                            }

                            if (ok) {
                                // create HDF5 dataset

                                herr_t h5Result;

                                // create a 1-d HDF5 dataspace
                                hsize_t h5Dims[1]    = {0};                                                                                     // initially 0, but...
                                hsize_t h5MaxDims[1] = {H5S_UNLIMITED};                                                                         // ... unlimited
                                hid_t   h5Dspace     = H5Screate_simple(1, h5Dims, h5MaxDims);                                                  // create the dataspace
                                hid_t   h5CPlist     = H5Pcreate(H5P_DATASET_CREATE);                                                           // create the dataset creation property list
                                (void)H5Pset_alloc_time(h5CPlist, H5D_ALLOC_TIME_INCR);                                                         // allocate space on disk incrementally
                                (void)H5Pset_layout(h5CPlist, H5D_CHUNKED);                                                                     // must be chunked when using unlimited dimensions

                                hsize_t h5ChunkDims[1] = {1};                                                                                   // chunk size affects performance - here just 1 value, 1-d
                                h5Result = H5Pset_chunk(h5CPlist, 1, h5ChunkDims);                                                              // set chunk size
                                if (h5Result < 0) {                                                                                             // ok?
Squawk("set chunk size failed");   // no - announce error
ok = false; // fail
                                }
                                else {                                                                                                          // yes - chunk size set ok
                                    // create HDF5 dataset
                                    string h5DsetName = utils::trim(fileDetails.hdrStrings[idx]);                                               // dataset name is just the header string
                                    hid_t  h5Dset     = H5Dcreate(m_Logfiles[fileDetails.id].h5file.groupId,                                    // create the dataset
                                                                  h5DsetName.c_str(),                                                           // dataset name
                                                                  m_Logfiles[fileDetails.id].h5file.dataTypes[idx],                             // datatype
                                                                  h5Dspace,                                                                     // dataspace
                                                                  H5P_DEFAULT,                                                                  // dataset link property list                                                                     
                                                                  h5CPlist,                                                                     // dataset creation property list
                                                                  H5P_DEFAULT);                                                                 // dataset access property list
                                    if (h5Dset < 0) {                                                                                           // dataset created ok?
Squawk("ERROR: Unable to create HDF5 dataSet " + h5DsetName);   // no - announce error
ok = false; // fail
                                    }
                                    else {                                                                                                      // yes - dataset created ok

                                        m_Logfiles[fileDetails.id].h5file.dataSets.push_back(h5Dset);                                           // dataset id

                                        // create attribute for units
                                        hid_t h5Dspace = H5Screate(H5S_SCALAR);                                                                 // HDF5 scalar dataspace
                                        hid_t h5DType  = H5Tcopy(H5T_C_S1);                                                                     // HDF5 c-string datatype

                                        string h5UnitsStr = utils::trim(fileDetails.unitsStrings[idx]);                                         // units string
                                        (void)H5Tset_size(h5DType, h5UnitsStr.length() + 1);                                                    // size is strlen + 1 (for NULL terminator)
                                        (void)H5Tset_cset(h5DType, H5T_CSET_ASCII);                                                             // ASCII (rather than UTF-8)
                                        hid_t h5Attr = H5Acreate(h5Dset, "units", h5DType, h5Dspace, H5P_DEFAULT, H5P_DEFAULT);                 // create attribute for units
                                        if (h5Attr < 0) {                                                                                       // attribute created ok?
Squawk("ERROR: Unable to create HDF5 attribute " + h5UnitsStr + " for HDF5 dataSet " + h5DsetName);   // no - announce error
ok = false; // fail
                                        }
                                        else {
                                            if (H5Awrite(h5Attr, h5DType, (const void *)h5UnitsStr.c_str()) < 0) {                              // write units attributes to file - ok?
Squawk("ERROR: Unable to write HDF5 attribute " + h5UnitsStr + " for HDF5 dataSet " + h5DsetName);   // no - announce error
ok = false; // fail
                                            }
                                        }
                                        (void)H5Aclose(h5Attr);                                                                                 // close attribute 
                                    }
//                                    (void)H5DClose(h5Dset);                                                                                     // close dataset
                                }
                                (void)H5Sclose(h5CPlist);                                                                                       // close creation property list
                                (void)H5Sclose(h5Dspace);                                                                                       // close scalar dataspace
                            }
                        }
                    }
                    else {                                                                                                                      // no - FS file
                        string fullHdrsStr  = "";                                                                                               // initialise full headers string
                        string fullUnitsStr = "";                                                                                               // initialise full units string
                        string fullTypesStr = "";                                                                                               // initialise full types string
                        for (size_t idx = 0; idx < fileDetails.hdrStrings.size(); idx++) {                                                      // for each property
                            fullHdrsStr += fileDetails.hdrStrings[idx] + delimiter;                                                             // append field header string to full header string
                            fullUnitsStr  += fileDetails.unitsStrings[idx] + delimiter;                                                         // append field units string to full units string
                            fullTypesStr   += fileDetails.typeStrings[idx] + delimiter;                                                         // append field type string to full type string
                        }

                        if (!fullHdrsStr.empty())  fullHdrsStr.pop_back();                                                                      // remove the trailing delimiter from the header string
                        if (!fullUnitsStr.empty()) fullUnitsStr.pop_back();                                                                     // remove the trailing delimiter from the units string
                        if (!fullTypesStr.empty()) fullTypesStr.pop_back();                                                                     // remove the trailing delimiter from the type string

                        // write the headers to file
                        if (!(ok = Put_(fileDetails.id, fullTypesStr))) {                                                                       // types string written ok?
                            Squawk(ERR_MSG(ERROR::FILE_WRITE_ERROR) + ": Type String");                                                         // no - show warning
                        }
                        else {                                                                                                                  // yes - ok
                            if (!(ok = Put_(fileDetails.id, fullUnitsStr))) {                                                                   // units string written ok?
                                Squawk(ERR_MSG(ERROR::FILE_WRITE_ERROR) + ": Units String");                                                    // no - show warning
                            }
                            else {                        
                                if (!(ok = Put_(fileDetails.id, fullHdrsStr))) {                                                                // header string written ok?
                                    Squawk(ERR_MSG(ERROR::FILE_WRITE_ERROR) + ": Header String");                                               // no - show warning
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (const std::exception& e) {                                                                                                       // unknown logfile
            fileDetails.recordProperties = {};                                                                                                  // no record properties
            Squawk(ERR_MSG(ERROR::UNKNOWN_LOGFILE) + ": Logging disabled for this file");                                                       // show warning
        }
    }
    else {                                                                                                                                      // already exists
        fileDetails = logfile->second;                                                                                                          // get existing file details
    }

    return ok ? fileDetails : retVal;
}


/*
 * Close a (currently open) standard logfile
 *
 * The logfile indicated by the p_Logfile parameter is closed and its details removed
 * from the map of currently open standard logfiles
 *
 *
 * bool CloseStandardFile(const LOGFILE p_Logfile, const bool p_Erase)
 *
 * @param   [IN]    p_Logfile                   The logfile to be closed
 * @param   [IN]    p_Erase                     Indicates whether the logfile entry should be erased
 * @return                                      True if the logfile was closed successfully, false if not
 */
bool Log::CloseStandardFile(const LOGFILE p_Logfile, const bool p_Erase) {

    bool result = true;                                                                                             // default is success

    COMPASUnorderedMap<LOGFILE, LogfileDetailsT>::const_iterator logfile;                                           // iterator
    logfile = m_OpenStandardLogFileIds.find(p_Logfile);                                                             // look for open logfile
    if (logfile != m_OpenStandardLogFileIds.end()) {                                                                // found
        result = Close_(logfile->second.id);                                                                        // close the file
        if (result && p_Erase) {                                                                                    // closed ok and erase required?
            m_OpenStandardLogFileIds.erase(logfile);                                                                // yes - remove from map
        }
    }

    return result;
}


/*
 * Close all currently open standard logfiles
 *
 * All currently open standard logfiles are closed and their details removed
 * from the map of currently open standard logfiles
 *
 *
 * bool CloseAllStandardFiles()
 *
 * @return                                      True if all currently open logfiles were closed successfully, false if not
 */
bool Log::CloseAllStandardFiles() {

    bool result = true;                                                                                             // default = success
    for (auto& iter: m_OpenStandardLogFileIds) {                                                                    // for each open standard log file
        if (!CloseStandardFile(iter.first, false)) result = false;                                                  // close it - flag if fail
    }
    if (result) m_OpenStandardLogFileIds.clear();                                                                   // remove all entries

    // close any open detailed output HDF5 container
    if (m_HDF5DetailedId >= 0) {                                                                                    // have open HDF5 detailed output file?
        if (H5Fclose(m_HDF5DetailedId) < 0) {                                                                       // yes - closed ok?
            result = false;                                                                                         // no - fail
        }
        m_HDF5DetailedId = -1;                                                                                      // (should have) no open HDF5 detailed output file
    }

    // close any open HDF5 container
    if (m_HDF5ContainerId >= 0) {                                                                                   // have open HDF5 detailed output file?
        if (H5Fclose(m_HDF5ContainerId) < 0) {                                                                      // yes - closed ok?
            result = false;                                                                                         // no - fail
        }
        m_HDF5ContainerId = -1;                                                                                     // (should have) no open HDF5 container file
    }

    return result;
}


/*
 * Prints text representation of the specification of a logfile record
 *
 * Written for testing the parse - left here for convenience in case it's needed later...
 *
 *
 * void PrintLogfileRecordDetails(const ANY_PROPERTY_VECTOR& p_LogfileRecord, const string p_LogfileRecordName)
 *
 * @param   [IN]    p_LogfileRecord             The logfile record for which details are to be printed
 * @param   [IN]    p_LogfileRecordName         The name of the logfile record for which details are to be printed
 */
void Log::PrintLogfileRecordDetails(const ANY_PROPERTY_VECTOR& p_LogfileRecord, const string p_LogfileRecordName) {

    SAY(p_LogfileRecordName << ": {");                                                                              // announce logfile record name
    for (auto const& property: p_LogfileRecord) {                                                                   // for each property in the record definition

        ANY_PROPERTY_TYPE propertyType = boost::apply_visitor(VariantPropertyType(), property);                     // get property type
        switch (propertyType) {                                                                                     // which property type?

            case ANY_PROPERTY_TYPE::T_STAR_PROPERTY:                                                                // single star
                SAY("    STAR_PROPERTY::" << STAR_PROPERTY_LABEL.at(boost::get<STAR_PROPERTY>(property)));
                break;

            case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY:                                                              // star 1 of binary
                SAY("    STAR_1_PROPERTY::" << STAR_PROPERTY_LABEL.at(static_cast<STAR_PROPERTY>(boost::get<STAR_1_PROPERTY>(property))));
                break;

            case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY:                                                              // star 2 of binary
                SAY("    STAR_2_PROPERTY::" << STAR_PROPERTY_LABEL.at(static_cast<STAR_PROPERTY>(boost::get<STAR_2_PROPERTY>(property))));
                break;

            case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY:                                                           // supernova star of binary
                SAY("    SUPERNOVA_PROPERTY::" << STAR_PROPERTY_LABEL.at(static_cast<STAR_PROPERTY>(boost::get<SUPERNOVA_PROPERTY>(property))));
                break;

            case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY:                                                           // companion star of binary
                SAY("    COMPANION_PROPERTY::" << STAR_PROPERTY_LABEL.at(static_cast<STAR_PROPERTY>(boost::get<COMPANION_PROPERTY>(property))));
                break;

            case ANY_PROPERTY_TYPE::T_BINARY_PROPERTY:                                                              // binary star
                SAY("    BINARY_PROPERTY::" << BINARY_PROPERTY_LABEL.at(boost::get<BINARY_PROPERTY>(property)));
                break;

            case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION:                                                               // program option
                SAY("    PROGRAM_OPTION::" << PROGRAM_OPTION_LABEL.at(boost::get<PROGRAM_OPTION>(property)));
                break;

            default:
                SAY("    UNKNOWN PROPERTY TYPE!");                                                                  // unknown - property type - show warning
        }
    }
    SAY("}");                                                                                                       // close
}


/*
 * Update logfile record specifier for given logfile
 *
 * The logfile record specifier for the given logfile is updated based on the parameters
 * passed (see descriptions below)
 *
 * The class member variable relevant to the given logfile is modified to reflect the
 * updated set of properties.
 *
 * Note the order in which the updates are applied.  Usually only one of p_AddProps and
 * p_SubtractProps will be populated, so order isn't relevant in that case.  Order is
 * relevant of both p_AddProps and p_SubtractProps are non-empty.  The order is:
 *
 *    1. The base set of parameters is set empty, or populated with the default set, according
 *       to p_UseDefaultProps
 *    2. Properties which are present in the base set and not present in p_SubtractProps, are
 *       copied to the new set of properties
 *    3. Properties which are present in p_AddProps and not already present in the new set of
 *       properties are copied to the new set
 *
 * The corollaries of the order of application are (for a single execution of this function):
 *
 *    1. Properties in p_SubtractProps that are not present in the base (or new) set are ignored
 *    2. Properties in p_AddProps that are already present in the base (or new) set are ignored
 *    3. Properties in p_AddProps will always be present in the new set of properties (i.e. they
 *       will not be subtracted)
 *
 *    Note that 1 & 2 don't generate parse errors because in each case the user's intent is satisfied
 *
 *
 * void UpdateLogfileRecordSpecs(const LOGFILE             p_Logfile,
 *                               bool                      p_UseDefaultProps,
 *                               const ANY_PROPERTY_VECTOR p_AddProps,
 *                               const ANY_PROPERTY_VECTOR p_SubtractProps)
 *
 * @param   [IN]    p_UseDefaultProps           indicates whether the default properties of the given logfile should be
 *                                              be used as the base set of properties.  If p_UseDefaultProps is true,
 *                                              the base set of properties is set to thedefault set (from constants.h).
 *                                              If p_UseDefaultProps is false, the base set of of properties is set
 *                                              empty
 * @param   [IN]    p_AddProps                  vector containing the properties to be added to the given logfile properties
 * @param   [IN]    p_SubtractProps             vector containing the properties to be subtracted from the given logfile properties
 */
void Log::UpdateLogfileRecordSpecs(const LOGFILE             p_Logfile,
                                   bool                      p_UseDefaultProps,
                                   const ANY_PROPERTY_VECTOR p_AddProps,
                                   const ANY_PROPERTY_VECTOR p_SubtractProps) {

    ANY_PROPERTY_VECTOR baseProps = {};                                                                                 // base props for the given logfile - deault is {}
    ANY_PROPERTY_VECTOR newProps = {};                                                                                  // new props for the given logfile - deault is {}

    if (p_UseDefaultProps) {                                                                                            // use logfile default props as base?
                                                                                                                        // yes - get existing props for given logfile
        switch (p_Logfile) {
            case LOGFILE::BSE_BE_BINARIES           : baseProps = m_BSE_BE_Binaries_Rec; break;
            case LOGFILE::BSE_COMMON_ENVELOPES      : baseProps = m_BSE_CEE_Rec;         break;
            case LOGFILE::BSE_DETAILED_OUTPUT       : baseProps = m_BSE_Detailed_Rec;    break;
            case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS: baseProps = m_BSE_DCO_Rec;         break;
            case LOGFILE::BSE_PULSAR_EVOLUTION      : baseProps = m_BSE_Pulsars_Rec;     break;
            case LOGFILE::BSE_RLOF_PARAMETERS       : baseProps = m_BSE_RLOF_Rec;        break;
            case LOGFILE::BSE_SUPERNOVAE            : baseProps = m_BSE_SNE_Rec;         break;
            case LOGFILE::BSE_SWITCH_LOG            : baseProps = m_BSE_Switch_Rec;      break;
            case LOGFILE::BSE_SYSTEM_PARAMETERS     : baseProps = m_BSE_SysParms_Rec;    break;
            case LOGFILE::SSE_DETAILED_OUTPUT       : baseProps = m_SSE_Detailed_Rec;    break;
            case LOGFILE::SSE_SUPERNOVAE            : baseProps = m_SSE_SNE_Rec;         break;
            case LOGFILE::SSE_SWITCH_LOG            : baseProps = m_SSE_Switch_Rec;      break;
            case LOGFILE::SSE_SYSTEM_PARAMETERS     : baseProps = m_SSE_SysParms_Rec;    break;
            default: break;                                                                                             // avoids compiler warning
        }
    }


    // copy all base props that are not to be subtracted
    if (!baseProps.empty()) {                                                                                           // nothing to copy if baseProps is empty
        for (auto const& baseProperty: baseProps) {                                                                     // for each property in baseProps
            ANY_PROPERTY_TYPE basePropertyType = boost::apply_visitor(VariantPropertyType(), baseProperty);             // base property type

            bool isSubtract = false;
            for (auto const& subtractProperty: p_SubtractProps) {                                                       // for each property in p_SubtractProps
                ANY_PROPERTY_TYPE subtractPropertyType = boost::apply_visitor(VariantPropertyType(), subtractProperty); // subtract property type
                if (basePropertyType == subtractPropertyType) {                                                         // props same property type?
                    switch (basePropertyType) {                                                                         // yes - which property type
                        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : isSubtract = boost::get<STAR_PROPERTY>(baseProperty)      == boost::get<STAR_PROPERTY>(subtractProperty);      break; // STAR_PROPERTY
                        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : isSubtract = boost::get<STAR_1_PROPERTY>(baseProperty)    == boost::get<STAR_1_PROPERTY>(subtractProperty);    break; // STAR_1_PROPERTY
                        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : isSubtract = boost::get<STAR_2_PROPERTY>(baseProperty)    == boost::get<STAR_2_PROPERTY>(subtractProperty);    break; // STAR_2_PROPERTY
                        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: isSubtract = boost::get<SUPERNOVA_PROPERTY>(baseProperty) == boost::get<SUPERNOVA_PROPERTY>(subtractProperty); break; // SUPERNOVA_PROPERTY
                        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: isSubtract = boost::get<COMPANION_PROPERTY>(baseProperty) == boost::get<COMPANION_PROPERTY>(subtractProperty); break; // STAR_1_PROPERTY
                        case ANY_PROPERTY_TYPE::T_BINARY_PROPERTY   : isSubtract = boost::get<BINARY_PROPERTY>(baseProperty)    == boost::get<BINARY_PROPERTY>(subtractProperty);    break; // BINARY_PROPERTY
                        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION    : isSubtract = boost::get<PROGRAM_OPTION>(baseProperty)     == boost::get<PROGRAM_OPTION>(subtractProperty);     break; // PROGRAM_OPTION
                    }
                    if (isSubtract) break;
                }
            }

            if (!isSubtract) newProps.push_back(baseProperty);                                                          // add property to newProps if necessary
        }
    }


    // copy all props to be added that don't already exist
    if (!p_AddProps.empty()) {                                                                                          // nothing to copy if p_AddProps is empty
        for (auto const& addProperty: p_AddProps) {                                                                     // for each property in p_AddProps
            ANY_PROPERTY_TYPE addPropertyType = boost::apply_visitor(VariantPropertyType(), addProperty);                             // add property type

            bool isAlready = false;
            for (auto const& newProperty: newProps) {                                                                   // for each property in newProps
                ANY_PROPERTY_TYPE newPropertyType = boost::apply_visitor(VariantPropertyType(), newProperty);           // new property type
                if (addPropertyType == newPropertyType) {                                                               // props same property type?
                    switch (addPropertyType) {                                                                          // yes - which propert type?
                        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : isAlready = boost::get<STAR_PROPERTY>(addProperty)      == boost::get<STAR_PROPERTY>(newProperty);      break; // STAR_PROPERTY
                        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : isAlready = boost::get<STAR_1_PROPERTY>(addProperty)    == boost::get<STAR_1_PROPERTY>(newProperty);    break; // STAR_1_PROPERTY
                        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : isAlready = boost::get<STAR_2_PROPERTY>(addProperty)    == boost::get<STAR_2_PROPERTY>(newProperty);    break; // STAR_2_PROPERTY
                        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: isAlready = boost::get<SUPERNOVA_PROPERTY>(addProperty) == boost::get<SUPERNOVA_PROPERTY>(newProperty); break; // SUPERNOVA_PROPERTY
                        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: isAlready = boost::get<COMPANION_PROPERTY>(addProperty) == boost::get<COMPANION_PROPERTY>(newProperty); break; // STAR_1_PROPERTY
                        case ANY_PROPERTY_TYPE::T_BINARY_PROPERTY   : isAlready = boost::get<BINARY_PROPERTY>(addProperty)    == boost::get<BINARY_PROPERTY>(newProperty);    break; // BINARY_PROPERTY
                        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION    : isAlready = boost::get<PROGRAM_OPTION>(addProperty)     == boost::get<PROGRAM_OPTION>(newProperty);     break; // PROGRAM_OPTION
                    }
                    if (isAlready) break;
                }
            }

            if (!isAlready) newProps.push_back(addProperty);                                                            // add property to newProps if necessary
        }
    }

    // replace  existing props for given logfile
    switch (p_Logfile) {
        case LOGFILE::BSE_BE_BINARIES           : m_BSE_BE_Binaries_Rec = newProps; break;
        case LOGFILE::BSE_COMMON_ENVELOPES      : m_BSE_CEE_Rec         = newProps; break;
        case LOGFILE::BSE_DETAILED_OUTPUT       : m_BSE_Detailed_Rec    = newProps; break;
        case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS: m_BSE_DCO_Rec         = newProps; break;
        case LOGFILE::BSE_PULSAR_EVOLUTION      : m_BSE_Pulsars_Rec     = newProps; break;
        case LOGFILE::BSE_RLOF_PARAMETERS       : m_BSE_RLOF_Rec        = newProps; break;
        case LOGFILE::BSE_SUPERNOVAE            : m_BSE_SNE_Rec         = newProps; break;
        case LOGFILE::BSE_SWITCH_LOG            : m_BSE_Switch_Rec      = newProps; break;
        case LOGFILE::BSE_SYSTEM_PARAMETERS     : m_BSE_SysParms_Rec    = newProps; break;
        case LOGFILE::SSE_DETAILED_OUTPUT       : m_SSE_Detailed_Rec    = newProps; break;
        case LOGFILE::SSE_SUPERNOVAE            : m_SSE_SNE_Rec         = newProps; break;
        case LOGFILE::SSE_SWITCH_LOG            : m_SSE_Switch_Rec      = newProps; break;
        case LOGFILE::SSE_SYSTEM_PARAMETERS     : m_SSE_SysParms_Rec    = newProps; break;
        default: break;                                                                                                 // avoids compiler warning...
    }
}


/*
 * Construct logfile record definitions based on user specifications
 *
 * Parses logfile definitions file and constructs logfile record definitions based
 * on the user specifications in the file.  Default record definitions are defined
 * in constants.h - the user specifications in the definitions file can use the
 * defaults as a base and add or subtract from the defaults, or the user specifications
 * can specify completley new definitions.
 *
 * The parser here is not very sophisticated - the syntax of the definitions file is
 * fairly simple.  The definitions file is expected to contain zero or more logfile
 * record specifications, as explained below.
 *
 * For the following specification:
 *
 *      ::=   means "expands to" or "is defined as"
 *     { x }  means (possible) repetition: x may appear zero or more times
 *     [ x ]  means x is optional: x may appear, or not
 *     <name> is a term (expression)
 *     "abc"  means literal string "abc"
 *       |    means "or"
 *       #    indicates the start of a comment
 *
 *
 * Logfile Definitions File specification:
 *
 * <def_file>   ::= { <rec_spec> }
 *
 * <rec_spec>   ::= <rec_name> <op> "{" { [ <props_list> ] } "}" <spec_delim>
 *
 * <rec_name>   ::= "SSE_SYSPARMS_REC"       |				# SSE only
 *                  "SSE_SWITCH_REC"         |				# SSE only
 *                  "SSE_SNE_REC"            |				# SSE only
 *                  "SSE_DETAILED_REC"       |				# SSE only
 *                  "BSE_SYSPARMS_REC"       |				# BSE only
 *                  "BSE_DCO_REC"            |				# BSE only
 *                  "BSE_SNE_REC"            |				# BSE only
 *                  "BSE_CEE_REC"            |				# BSE only
 *                  "BSE_RLOF_REC"           |				# BSE only
 *                  "BSE_BE_BINARIES_REC"    |				# BSE only
 *                  "BSE_PULSARS_REC"        |				# BSE only
 *                  "BSE_DETAILED_REC"	     |				# BSE only
 *                  "BSE_SWITCH_REC"		   			    # BSE only
 *
 * <op>         ::= "=" | "+=" | "-="
 *
 * <props_list> ::= <prop_spec> [ <prop_delim> <props_list> ]
 *
 * <prop_spec>  ::= <prop_type> "::" <prop_name> <prop_delim>
 *
 * <spec_delim> ::= " " | EOL
 *
 * <prop_delim> ::= "," | <spec_delim>
 *
 * <prop_type>  ::= "STAR_PROPERTY"      |				# SSE only
 *                  "STAR_1_PROPERTY"    |				# BSE only
 *                  "STAR_2_PROPERTY"    |				# BSE only
 *                  "SUPERNOVA_PROPERTY" |				# BSE only
 *                  "COMPANION_PROPERTY" |				# BSE only
 *                  "BINARY_PROPERTY"    |				# BSE only
 *                  "PROGRAM_OPTION"      				# SSE or BSE
 *
 * <prop_name>  ::= valid property name for specified property type
 *                  (see definitions in constants.h)
 *
 *
 * The file may contain comments.  Comments are denoted by the hash/pound character ('#').
 * The hash character and any text following it on the line in which the hash character appears
 * is ignored by the parser.  The hash character can appear anywhere on a line - if it is the
 * first character then the entire line is a comment and ignored by the parser, or it can follow
 * valid symbols on a line, in which case the symbols before the hash character are parsed and
 * interpreted by the parser.
 *
 * A logfile record is initially set to its default value (the default logfile record specifications
 * are defined in constants.h).  The definitions file informs the code as to the modifications to the
 * default values the user wants.  This means that the definitions logfile is not mandatory, and if
 * the definitions file is not present, or empty, the logfile record definitions will remain at their
 * default values.
 *
 * The assignment operator given in a record specification (<op> in the file specification above) can
 * be one of “=”, “+=”, and “-=”.  The meanings of these are:
 *
 *    “=”  means that the record specifier should be assigned the list of properties specified in the
 *         braced-list following the “=” operator.  The value of the record specifier prior to the
 *         assignment is discarded, and the new value set as described.
 *
 *    “+=” means that the list of properties specified in the braced-list following the “+=” operator
 *         should be appended to the existing value of the record specifier.  Note that the new
 *         properties are appended to the existing list, so will appear at the end of the list
 *         (properties are printed in the order they appear in the list).
 *
 *    “-=” means that the list of properties specified in the braced-list following the “-=” operator
 *         should be subtracted from the existing value of the record specifier.
 *
 *
 * bool UpdateAllLogfileRecordSpecs()
 *
 * @return                                      boolean indicating if if lohfile record specifications were updated successfully:
 *                                              true = yes - updated successfully, false = no - an error occurred
 */
bool Log::UpdateAllLogfileRecordSpecs() {

    string filename = OPTIONS->LogfileDefinitionsFilename();                                                                    // get user-specified definitions file

    // do some sanity checks in the definitions file

    if (filename.empty()) return true;                                                                                          // user did not specify a file - nothing to do

	if (!utils::FileExists(filename)) {                                                                                         // check definitions file exists
        SAY("");
        SAY(ERR_MSG(ERROR::BAD_LOGFILE_RECORD_SPECIFICATIONS));                                                                 // announce error
        SAY(ERR_MSG(ERROR::FILE_DOES_NOT_EXIST) + ": " + filename);                                                             // file does not exist - show warning, and ...
        return false;                                                                                                           // ... bail/bale out
	}

    std::ifstream defFile;
    defFile.open(filename);                                                                                                     // open the definitions file
	if (defFile.fail()) {                                                                                                       // check open succeeded
        SAY("");
        SAY(ERR_MSG(ERROR::BAD_LOGFILE_RECORD_SPECIFICATIONS));                                                                 // announce error
        SAY(ERR_MSG(ERROR::FILE_OPEN_ERROR) + ": " + filename);                                                                 // failed - show warning, and ...
        return false;                                                                                                           // ... bail/bale out
    }

    ERROR       error = ERROR::NONE;                                                                                            // initially no error
    std::size_t errorPos;                                                                                                       // position of error in input record

    enum class TOKEN_TYPE: int { LOGFILE_RECORD_NAME, ASSIGN, COMMA, OPEN_BRACE, CLOSE_BRACE, PROPERTY_SPECIFIER };             // token types

    std::vector<std::tuple<string, std::size_t>> strTokens = {};                                                                // parsed tokens - token value and column position

    bool useDefaultProps;                                                                                                       // indicates whether the default props for a logfile should be the base set of props
    ANY_PROPERTY_VECTOR addProps;                                                                                               // properties user wants added to the base properties
    ANY_PROPERTY_VECTOR subtractProps;                                                                                          // properties user wants subtracted from the base properties

    // read and parse the file records

    TOKEN_TYPE expecting            = TOKEN_TYPE::LOGFILE_RECORD_NAME;                                                          // token type we're expecting to see - initally a logfile record name
    LOGFILE currentLogfile          = LOGFILE::NONE;                                                                            // the logfile definition being modified
    LOGFILE_TYPE currentLogfileType = LOGFILE_TYPE::NONE;                                                                       // the type of the logfile definition being modified (STELLAR or BINARY)

    string recIn;                                                                                                               // input record
    string parseRec;                                                                                                            // record to be parsed - input record after stripping spaces and comments
    string recParsed;                                                                                                           // the most recent record parsed (will be last non-empty record) - used for error handling
    while (std::getline(defFile, recIn)) {                                                                                      // read the next record from the file

        parseRec = recIn;                                                                                                       // copy the record just read
        size_t hashPos = parseRec.find("#");                                                                                    // find first occurrence of "#"
        if (hashPos != string::npos) parseRec.erase(hashPos, parseRec.size() - hashPos);                                        // if "#" found, prune it and everything after it (ignore comments)

        if (parseRec.empty()) continue;                                                                                         // ignore empty records

        recParsed = recIn;                                                                                                      // for error handling
        errorPos  = recParsed.size();                                                                                           // initially

        // tokenise the input record

        strTokens.clear();                                                                                                      // clear the vector of tokens

        std::size_t prev = 0;                                                                                                   // previous position in the input record (token start)
        std::size_t pos  = 0;                                                                                                   // current position in the input record (delimiter position)
        while ((pos = parseRec.find_first_of(" ,+-={}", prev)) != string::npos) {                                               // find the next delimiter

            if (pos > prev) {                                                                                                   // found - token string before delimiter?
                string tokStr = parseRec.substr(prev, pos - prev);                                                              // yes - extract token string
                tokStr.erase(remove_if(tokStr.begin(), tokStr.end(), ::isspace), tokStr.end());                                 // remove whitespace from token
                if (!tokStr.empty()) strTokens.push_back(std::make_tuple(tokStr, prev));                                        // stash the token string and position
            }

            string delimStr = parseRec.substr(pos, 1);                                                                          // extract the delimiter string
            if (delimStr != " ") {                                                                                              // discard whitespace

                // delimiter is also a token
                // if delimiter found is one of {"-", "+"} check for 2-character delimiter

                std::size_t delimPos = pos;                                                                                     // position of the delimiter
                if (delimStr == "-" || delimStr =="+") {                                                                        // possible 2-character delimiter?
                    if (pos + 1 < parseRec.size()) {                                                                            // yes - at end of input record?
                        if (parseRec[pos + 1] == '=') {                                                                         // no - next character "="?
                            delimStr += "=";                                                                                    // yes - 2-character delimiter
                            pos++;                                                                                              // skip next character
                        }
                    }
                }
                strTokens.push_back(std::make_tuple(delimStr, delimPos));                                                       // stash the delimiter string and position
            }
            prev = pos + 1;                                                                                                     // advance past delimiter just found
        }

        if (pos > prev) {                                                                                                       // found - non-empty token string before end of line?
            string tokStr = parseRec.substr(prev, pos - prev);                                                                  // yes - extract token string
            tokStr.erase(remove_if(tokStr.begin(), tokStr.end(), ::isspace), tokStr.end());                                     // remove whitespace from token
            if (!tokStr.empty()) strTokens.push_back(std::make_tuple(tokStr, prev));                                            // stash the token string and position
        }

        // process the tokenised input record

        bool addAssign;                                                                                                         // assignment type - adding or subtracting to the base properties

        for (auto const& strTok: strTokens) {                                                                                   // for each token

            string      tokStr = std::get<0>(strTok);                                                                           // token string
            std::size_t tokPos = std::get<1>(strTok);                                                                           // token position in input record
            errorPos           = tokPos;                                                                                        // in most instances

            switch (expecting) {                                                                                                // what token type are we expecting?

                case TOKEN_TYPE::LOGFILE_RECORD_NAME:                                                                           // ...logfile record name

                    bool found;
                    LOGFILE logfile;                                                                                            // logfile name
                    std::tie(found, logfile) = GetLogfileDescriptorKey(tokStr);                                                 // check token against known logfile record names
                    if (found) {                                                                                                // found - token is a logfile record name
                        currentLogfile     = logfile;                                                                           // record the logfile name
                        currentLogfileType = std::get<4>(LOGFILE_DESCRIPTOR.at(currentLogfile));                                // and type (STELLAR or BINARY)
                        addProps      = {};                                                                                     // start with empty set of properties to be added
                        subtractProps = {};                                                                                     // start with empty set of properties to be subtracted

                        expecting = TOKEN_TYPE::ASSIGN;                                                                         // now expecting assignment operator {"=", "-=", "+="}
                    }
                    else {                                                                                                      // token is not a known logfile record name - error
                        error = ERROR::EXPECTED_LOGFILE_RECORD_NAME;                                                            // set error
                    }
                    break;                                                                                                      // end expecting LOGFILE_RECORD_NAME::ASSIGN

                case TOKEN_TYPE::ASSIGN:                                                                                        // ... assignment operator

                    switch (_(tokStr.c_str())) {                                                                                // token is...?

                        case _("="):                                                                                            // plain assign - empty base and add
                            useDefaultProps = false;                                                                            // start with empty base properties vector - don't use default props
                            addAssign = true;                                                                                   // add new properties to base properties
                            expecting = TOKEN_TYPE::OPEN_BRACE;                                                                 // now expecting open brace "{" to start properties list
                            break;

                        case _("-="):                                                                                           // subtract from base properties
                        case _("+="):                                                                                           // add to base properties
                            useDefaultProps = true;                                                                             // start with base properties = default properties for current logfile
                            addAssign = (tokStr == "+=");                                                                       // add or subtract new properties to/from props?
                            expecting = TOKEN_TYPE::OPEN_BRACE;                                                                 // now expecting open brace "{" to start properties list
                            break;

                        default:                                                                                                // didn't get assign token - error
                            error = ERROR::EXPECTED_ASSIGNMENT_OPERATOR;                                                        // set error
                    }
                    break;                                                                                                      // end expecting TOKEN_TYPE::ASSIGN

                case TOKEN_TYPE::OPEN_BRACE:                                                                                    // ... open brace (begin properties list)

                    if (tokStr == "{")                                                                                          // token is open brace "{"?
                        expecting = TOKEN_TYPE::PROPERTY_SPECIFIER;                                                             // yes - now expecting the name of a property to add or subtract ( or close brace)
                    else                                                                                                        // no - didn't get open brace - error
                        error = ERROR::EXPECTED_OPEN_BRACE;                                                                     // set error
                    break;                                                                                                      // end expecting TOKEN_TYPE::OPEN_BRACE

                case TOKEN_TYPE::COMMA:                                                                                         // ... comma (or close brace)

                    if (tokStr == "}") {                                                                                        // token is open brace "{"?
                        UpdateLogfileRecordSpecs(currentLogfile, useDefaultProps, addProps, subtractProps);                     // update current logfile record specifications

                        expecting = TOKEN_TYPE::LOGFILE_RECORD_NAME;                                                            // now expecting new logfile record name
                    }
                    else {                                                                                                      // not close brace - try comma
                        if (tokStr == ",")                                                                                      // token is comma ","?
                            expecting = TOKEN_TYPE::PROPERTY_SPECIFIER;                                                         // yes - now expecting the name of a property to add or subtract (or close brace)
                        else                                                                                                    // no - didn't get comma or close brace - error
                            error = ERROR::EXPECTED_COMMA_OR_CLOSE_BRACE;                                                       // set error
                    }
                    break;                                                                                                      // end expecting TOKEN_TYPE::OPEN_BRACE

                case TOKEN_TYPE::PROPERTY_SPECIFIER:                                                                            // ... property specifier

                    // check for close brace to end property list
                    if (tokStr == "}") {                                                                                        // have close brace - end of list
                        UpdateLogfileRecordSpecs(currentLogfile, useDefaultProps, addProps, subtractProps);                     // update current logfile record specifications

                        expecting = TOKEN_TYPE::LOGFILE_RECORD_NAME;                                                            // now expecting new logfile record name
                    }
                    else {                                                                                                      // no close brace - look for property specifier
                        // property specifier must be of the form PROPERTY_TYPE::PROPERTY_NAME
                        string      propTypeStr;                                                                                // first part of property specifier - property type
                        string      propNameStr;                                                                                // second part of property specifier - property name
                        std::size_t propTypeLen;                                                                                // length of the property type string

                        if ((propTypeLen = tokStr.find("::")) != string::npos) {                                                // find :: separator
                            if (propTypeLen > 0) {                                                                              // :: separator found - have property type?
                                propTypeStr = tokStr.substr(0, propTypeLen);                                                    // yes - extract property type from token

                                if (tokStr.size() - propTypeStr.size() - 2 > 0) {                                               // have property name?
                                    std::size_t namePos = propTypeStr.size() + 2;                                               // yes - start position of property name in token
                                    std::size_t nameLen = tokStr.size() - propTypeStr.size() - 2;                               // length of property name in token
                                    propNameStr = tokStr.substr(namePos, nameLen);                                              // extract property name from token
                                }
                                else {                                                                                          // didn't get property name - error
                                    error    = ERROR::EXPECTED_PROPERTY_SPECIFIER;                                              // set error
                                    errorPos = tokPos + propTypeLen + 2;                                                        // position of error in input record
                                }
                            }
                            else {                                                                                              // didn't get property type - error
                                error = ERROR::EXPECTED_PROPERTY_SPECIFIER;                                                     // set error
                            }
                        }
                        else {                                                                                                  // didn't get :: separator - error
                            error = ERROR::EXPECTED_PROPERTY_SPECIFIER;                                                         // set error
                        }

                        // if no error at this stage, we have the property type and name
                        // check that property type and name are a valid, known property type and name

                        if (error == ERROR::NONE) {                                                                             // parse error?
                                                                                                                                // no - have property type and name
                            bool found;
                            PROPERTY_TYPE propertyType;                                                                         // check property type against known property types
                            std::tie(found, propertyType) = utils::GetMapKey(propTypeStr, PROPERTY_TYPE_LABEL, PROPERTY_TYPE::NONE);
                            if (!found) {                                                                                       // found known property type?
                                error = ERROR::UNKNOWN_PROPERTY_TYPE;                                                           // no - set error
                            }
                            else {                                                                                              // yes - property type is a known property type
                                switch (propertyType) {                                                                         // which (known) property type?

                                    case PROPERTY_TYPE::STAR_PROPERTY : {                                                       // STAR_PROPERTY

                                        bool found;
                                        STAR_PROPERTY property;                                                                 // lookup property name
                                        std::tie(found, property) = utils::GetMapKey(propNameStr, STAR_PROPERTY_LABEL, STAR_PROPERTY::ID);
                                        if (!found) {                                                                           // property name found?
                                            error = ERROR::UNKNOWN_STELLAR_PROPERTY;                                            // no - set error
                                        }
                                        else {                                                                                  // found known property name
                                            if (currentLogfileType == LOGFILE_TYPE::STELLAR) {                                  // current logfile type STELLAR?
                                                                                                                                // yes - ok
                                                if (addAssign)                                                                  // add property?
                                                    addProps.push_back(property);                                               // yes - add it to 'add' list
                                                else                                                                            // no - subtract property
                                                    subtractProps.push_back(property);                                          // add it to 'subtract' list

                                                expecting = TOKEN_TYPE::COMMA;                                                  // now expecting a comma or close brace
                                            }
                                            else {                                                                              // current logfile type is not STELLAR - not ok
                                                error = ERROR::EXPECTED_BINARY_PROPERTY;                                        // set error - assumption is BINARY if not STELLAR
                                            }
                                        }
                                        } break;

                                    case PROPERTY_TYPE::STAR_1_PROPERTY   :                                                     // STAR_1_PROPERTY, or
                                    case PROPERTY_TYPE::STAR_2_PROPERTY   :                                                     // STAR_2_PROPERTY, or
                                    case PROPERTY_TYPE::SUPERNOVA_PROPERTY:                                                     // SUPERNOVA_PROPERTY, or
                                    case PROPERTY_TYPE::COMPANION_PROPERTY: {                                                   // COMPANION_PROPERTY

                                        bool found;
                                        STAR_PROPERTY property;                                                                 // lookup property name
                                        std::tie(found, property) = utils::GetMapKey(propNameStr, STAR_PROPERTY_LABEL, STAR_PROPERTY::ID);
                                        if (!found) {                                                                           // property name found?
                                            error = ERROR::UNKNOWN_BINARY_PROPERTY;                                             // no - set error
                                        }
                                        else {                                                                                  // found known property name
                                            if (currentLogfileType == LOGFILE_TYPE::BINARY) {                                   // current logfile type BINARY?
                                                                                                                                // yes - ok
                                                if (addAssign) {                                                                // add property?
                                                    switch (propertyType) {                                                     // yes - add it to 'add' list
                                                        case PROPERTY_TYPE::STAR_PROPERTY     : addProps.push_back(static_cast<STAR_PROPERTY>(property));      break;
                                                        case PROPERTY_TYPE::STAR_1_PROPERTY   : addProps.push_back(static_cast<STAR_1_PROPERTY>(property));    break;
                                                        case PROPERTY_TYPE::STAR_2_PROPERTY   : addProps.push_back(static_cast<STAR_2_PROPERTY>(property));    break;
                                                        case PROPERTY_TYPE::SUPERNOVA_PROPERTY: addProps.push_back(static_cast<SUPERNOVA_PROPERTY>(property)); break;
                                                        case PROPERTY_TYPE::COMPANION_PROPERTY: addProps.push_back(static_cast<COMPANION_PROPERTY>(property)); break;
                                                        default: break;                                                         // avoids compiler warning
                                                    }
                                                    expecting = TOKEN_TYPE::COMMA;                                              // now expecting a comma or close brace
                                                }
                                                else {                                                                          // no - subtract property
                                                    switch (propertyType) {                                                     // add it to 'subtract' list
                                                        case PROPERTY_TYPE::STAR_PROPERTY     : subtractProps.push_back(static_cast<STAR_PROPERTY>(property));      break;
                                                        case PROPERTY_TYPE::STAR_1_PROPERTY   : subtractProps.push_back(static_cast<STAR_1_PROPERTY>(property));    break;
                                                        case PROPERTY_TYPE::STAR_2_PROPERTY   : subtractProps.push_back(static_cast<STAR_2_PROPERTY>(property));    break;
                                                        case PROPERTY_TYPE::SUPERNOVA_PROPERTY: subtractProps.push_back(static_cast<SUPERNOVA_PROPERTY>(property)); break;
                                                        case PROPERTY_TYPE::COMPANION_PROPERTY: subtractProps.push_back(static_cast<COMPANION_PROPERTY>(property)); break;
                                                        default: break;                                                         // avoids compiler warning
                                                    }
                                                    expecting = TOKEN_TYPE::COMMA;                                              // now expecting a comma or close brace
                                                }
                                            }
                                            else {                                                                              // current logfile type is not BINARY - not ok
                                                error = ERROR::EXPECTED_STELLAR_PROPERTY;                                       // set error - assumption is STELLAR if not BINARY
                                            }
                                        }
                                        } break;

                                    case PROPERTY_TYPE::BINARY_PROPERTY: {                                                      // BINARY_PROPERTY

                                        bool found;
                                        BINARY_PROPERTY property;                                                               // lookup property name
                                        std::tie(found, property) = utils::GetMapKey(propNameStr, BINARY_PROPERTY_LABEL, BINARY_PROPERTY::ID);
                                        if (!found) {                                                                           // property name found?
                                            error = ERROR::UNKNOWN_BINARY_PROPERTY;                                             // no - set error
                                        }
                                        else {                                                                                  // found known property name
                                            if (currentLogfileType == LOGFILE_TYPE::BINARY) {                                   // current logfile type BINARY?
                                                                                                                                // yes - ok
                                                if (addAssign)                                                                  // add property?
                                                    addProps.push_back(property);                                               // yes - add it to 'add' list
                                                else                                                                            // no - subtract property
                                                    subtractProps.push_back(property);                                          // add it to 'subtract' list

                                                expecting = TOKEN_TYPE::COMMA;                                                  // now expecting a comma or close brace
                                            }
                                            else {                                                                              // current logfile type is not BINARY - not ok
                                                error = ERROR::EXPECTED_STELLAR_PROPERTY;                                       // set error - assumption is STELLAR if not BINARY
                                            }
                                        }
                                        } break;

                                    case PROPERTY_TYPE::PROGRAM_OPTION: {                                                       // PROGRAM_OPTION

                                        bool found;
                                        PROGRAM_OPTION property;                                                                // lookup property name
                                        std::tie(found, property) = utils::GetMapKey(propNameStr, PROGRAM_OPTION_LABEL, PROGRAM_OPTION::RANDOM_SEED);
                                        if (!found) {                                                                           // property name found?
                                            error = ERROR::UNKNOWN_PROGRAM_OPTION;                                              // no - set error
                                        }
                                        else {                                                                                  // found known property name
                                            if (addAssign)                                                                      // add property?
                                                addProps.push_back(property);                                                   // yes - add it to 'add' list
                                            else                                                                                // no - subtract property
                                                subtractProps.push_back(property);                                              // add it to 'subtract' list
                                        }
                                        expecting = TOKEN_TYPE::COMMA;                                                          // now expecting a comma or close brace
                                        } break;

                                    default:                                                                                    // unknown property type
                                        error = ERROR::UNKNOWN_PROPERTY_TYPE;                                                   // set error
                                }
                            }
                        }
                    }                                                                                                           // end expecting TOKEN_TYPE::PROPERTY_NAME
                    break;

                case TOKEN_TYPE::CLOSE_BRACE:                                                                                   // ... close brace (end properties list)
                    break;                                                                                                      // nothing to do but wait for next token

            }                                                                                                                   // end expecting switch
            if (error != ERROR::NONE) break;                                                                                    // stop processing tokens on error
        }                                                                                                                       // end for all tokens loop
        if (error != ERROR::NONE) break;                                                                                        // stop processing records on error
    }                                                                                                                           // end reading file

    // show error details if an error occurred
    // parsing is considered incomplete if after parsing the complete file we are left
    // expecting anything other than a logfile record name - that is considered an error
    // case an will generate an "Unexpected end of file" error

    if (error != ERROR::NONE || expecting != TOKEN_TYPE::LOGFILE_RECORD_NAME) {                                                 // error?

        SAY("");
        SAY(ERR_MSG(ERROR::BAD_LOGFILE_RECORD_SPECIFICATIONS) << " in file: " << filename);                                     // announce error

        if (error == ERROR::NONE) {                                                                                             // must be unexpected end of file
            error = ERROR::UNEXPECTED_END_OF_FILE;                                                                              // set error
            SAY(ERR_MSG(error));                                                                                                // announce error
            size_t hashPos = recParsed.find("#");                                                                               // find first occurrence of "#" in the last record parsed
            errorPos = hashPos == string::npos ? recParsed.size() : hashPos;                                                    // set location for caret indicator ("^")
            error = ERROR::EXPECTED_LOGFILE_RECORD_NAME;                                                                        // set error
        }

        SAY(recParsed);                                                                                                         // show the last record parsed

        string loc(errorPos, ' ');                                                                                              // leading spaces for caret indicator
        loc += "^";                                                                                                             // add the caret indicator
        SAY(loc);                                                                                                               // show the caret indicator

        SAY(ERR_MSG(error));                                                                                                    // announce error
    }

	return (error == ERROR::NONE);
}
