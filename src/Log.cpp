// Class Log
//
// This is where all logging and debugging is performed.
//
// JR: todo: complete this documentation
// JR: todo: move error/warning strings to error catalogue in constants.h
// JR: todo: clean up use of Squawk() vs SAY() etc

#include "Log.h"

Log* Log::m_Instance = nullptr;


Log* Log::Instance() {
    if (!m_Instance) {
        m_Instance = new Log();
    }
    return m_Instance;
}


/*
 * Open the run details file inside the HDF5 container (if logging to HDF5 files)
 * 
 * Creates the file (group) inside the HDF5 container, and creates the columns
 * (datasets) required.  Ciluns (datasets) are created for the preamble/stats
 * information written to the run details file: 
 * 
 *  - COMPAS version (STRING: 'xx.yy.zz')
 *  - run start time (STRING: formatted system time)
 *  - run end time (STRING: formatted system time)
 *  - number of objects (stars/binaries) requested (INT)
 *  - number of objects (stars/binaries) created (INT)
 *  - CPU (clock) time (DOUBLE: seconds)
 *  - Wall (elapsed) time (STRING: 'hhhh:mm:ss')
 *  - Actual random seed used (UNSIGNED LONG INT) * 
 * 
 * as well as columns (datasets) for each of the program options (whether the
 * use specified them on the commandline or not).
 *
 * Additionally, all columns (datasets) in the HDF5 copy of the run details file have
 * a "shadow" column (dataset) (the "derivation" columns) that has '-Derivation' appended
 * to its name (this shadow columns is not strictly necessary for the preamble/stats columns
 * (datasets), but included for consistency).
 *
 * The "derivation" column indicates how the data was derived, and will be one of the following
 * strings (particularly relevant for program options):
 *
 *  - 'USER_SUPPLIED' : indicates the user supplied a value, and the user-supplied value was used
 *  - 'DEFAULT_USED'  : indicates the user did not supply a value, and the COMPAS C++ default value was used
 *  - 'CALCULATED'    : the value was calculated by COMPAS
 * 
 * None of the run details columns have assocatied "units" (we could add them for the preamble/stats
 * datasets, but too hard to work them out for the program options - we could with a bit more coding,
 * but I don't really think we need them here (we've not had them in the run details file in the past
 * and it hasn't caused a problem))
 * 
 * 
 * bool Log::OpenHDF5RunDetailsFile(const string p_Filename)
 * 
 * @param   [IN]    p_Filename                  The run details filename (group name for the HDF5 file)
 * @return                                      Boolean status - true = created ok; false = open failed
 */
bool Log::OpenHDF5RunDetailsFile(const string p_Filename) {

    if (!m_Enabled) return m_Enabled;                                                                                   // logging not enabled - no business being here

    bool ok = true;                                                                                                     // return value

    // The run details file is not treated as a standard logfile, so we need to 
    // open it manually rather than have Log::StandardLogFileDetails() open it.
    //
    // Note that the HDF5 run details file columns (datasets) have the "units" attribute
    // set to '-' for all columns (datasets) - the contents of the run details file are
    // just a reflection of the values supplied by the user (or calculated by COMPAS), 
    // so no units.

    m_Run_Details_H5_File.fileId = m_HDF5ContainerId;                                                                   // record run details HDF5 fileid - just the HDF5 container id

    // open the run details file inside the HDF5 container
    string h5GroupName = p_Filename;                                                                                    // HDF5 group name for run details file
    h5GroupName        = utils::trim(h5GroupName);                                                                      // remove leading and trailing blanks
    hid_t h5GroupId = H5Gopen(m_Run_Details_H5_File.fileId, h5GroupName.c_str(), H5P_DEFAULT);                          // open the group
    if (h5GroupId >= 0) {                                                                                               // group open (and therefore already exists)?
        Squawk("ERROR: HDF5 group with name " + h5GroupName + " already exists");                                       // that's not ok - announce error
        (void)H5Gclose(h5GroupId);                                                                                      // close the group
        ok = false;                                                                                                     // fail
    }
    else {                                                                                                              // group does not exist/is not open               
        h5GroupId = H5Gcreate(m_HDF5ContainerId, h5GroupName.c_str(), 0, H5P_DEFAULT, H5P_DEFAULT);                     // create the group
        if (h5GroupId < 0) {                                                                                            // group created ok?
            Squawk("ERROR: Error creating HDF5 group with name " + h5GroupName);                                        // no - announce error
            ok = false;                                                                                                 // fail
        }
        else {                                                                                                          // group created ok
            m_Run_Details_H5_File.groupId = h5GroupId;                                                                  // record group id for run details file

            // We now have the run details file in the HDF5 container, so we create the required
            // columns (datasets) here (all columns, not just the preamble/stats).
            //
            // Preamble/stats columns (datasets) are:
            //      - COMPAS version (STRING: 'xx.yy.zz')
            //      - run start time (STRING: formatted system time)
            //      - run end time (STRING: formatted system time)
            //      - number of objects (stars/binaries) requested (INT)
            //      - number of objects (stars/binaries) created (INT)
            //      - CPU (clock) time (DOUBLE: seconds)
            //      - Wall (elapsed) time (STRING: 'hhhh:mm:ss')
            //      - Actual random seed used (UNSIGNED LONG INT)
            //
            // We also create columns (datasets) for each of the program options.
            //
            // All columns (datasets) in the HDF5 copy of the run details file have a 
            // "shadow" column (dataset) that has '-Derivation' appended to its name.
            // (Not strictly necessary for the preamble/stats columns (datasets), but
            // included for consistency).
            //
            // The "derivation" column indicates how the data was derived, and will be
            // one of the following strings (particularly relevant for program options):
            //
            //      - 'USER_SUPPLIED' : indicates the user supplied a value, and the user-supplied value was used
            //      - 'DEFAULT_USED'  : indicates the user did not supply a value, and the COMPAS C++ default value was used
            //      - 'CALCULATED'    : the value was calculated by COMPAS

            if (ok) {                                                                                                   // still ok?

                string h5DatasetName;
                hid_t  h5DataType;
                hid_t  h5Dset;
                hid_t  h5String13DataType = GetHDF5DataType(TYPENAME::STRING, 13);                                      // HDF5 data type for 13-character string (derivation columns)
                  
                string h5Filename = OPTIONS->OutputContainerName();                                                     // HDF5 container file name

                size_t chunkSize = HDF5_MINIMUM_CHUNK_SIZE;                                                             // chunk size
                size_t IOBufSize = OPTIONS->HDF5BufferSize() * chunkSize;                                               // IO buffer size

                m_Run_Details_H5_File.chunkSize = chunkSize;                                                            // record chunk size for file
                m_Run_Details_H5_File.IOBufSize = IOBufSize;                                                            // record IO buf size for file

                // preamble/stats datasets

                for (int dSetIdx = static_cast<int>(RUN_DETAILS_COLUMNS::COMPAS_VERSION); dSetIdx < static_cast<int>(RUN_DETAILS_COLUMNS::SENTINEL); dSetIdx++ ) {

                    std::tuple<std::string, TYPENAME, std::size_t> runDetails;
                    try { runDetails = RUN_DETAILS_DETAIL.at(static_cast<RUN_DETAILS_COLUMNS>(dSetIdx)); }              // get run details details
                    catch (const std::exception& e) {                                                                   // unknown property
                        Squawk("ERROR: Unknown property for HDF5 file with name " + h5Filename);                        // announce error
                        ok = false;                                                                                     // fail
                    }
                    
                    if (ok) {                                                                                           // have valid property
                        h5DatasetName = std::get<0>(runDetails);                                                        // dataset name
                        TYPENAME compasType = std::get<1>(runDetails);                                                  // COMPAS data type
                        h5DataType = GetHDF5DataType(compasType, std::get<2>(runDetails));                              // HDF5 data type
                        h5Dset = CreateHDF5Dataset(h5Filename, h5GroupId, h5DatasetName, h5DataType, "-", chunkSize);   // create dataset
                        if (h5Dset < 0) {                                                                               // dataset not created
                            Squawk("ERROR: Error creating HDF5 dataset with name " + h5DatasetName);                    // announce error
                            ok = false;                                                                                 // fail
                        }
                        else {                                                                                          // dataset created ok

                            m_Run_Details_H5_File.dataSets.push_back({h5Dset, h5DataType, compasType, STRING_QUALIFIER::FIXED_LENGTH, {}}); // record dataset details

                            // derivation
                            h5DatasetName += "-Derivation";                                                             // derivation
                            h5Dset = CreateHDF5Dataset(h5Filename, h5GroupId, h5DatasetName, h5String13DataType, "-", chunkSize); // create dataset
                            if (h5Dset < 0) {                                                                           // dataset not created
                                Squawk("ERROR: Error creating HDF5 dataset with name " + h5DatasetName);                // announce error
                                ok = false;                                                                             // fail
                            }
                            else
                                m_Run_Details_H5_File.dataSets.push_back({h5Dset, h5String13DataType, TYPENAME::STRING, STRING_QUALIFIER::FIXED_LENGTH, {}}); // dataset created ok - record details
                        }
                    }
                    if (!ok) break;                                                                                     // something went wrong - fail
                }

                // program options datasets

                for (std::size_t idx = 0; idx < m_OptionDetails.size(); idx++) {                                        // for each program option
                    // option
                    TYPENAME compasType = std::get<4>(m_OptionDetails[idx]);                                            // COMPAS data type
                    h5DataType = GetHDF5DataType(compasType, (std::get<1>(m_OptionDetails[idx])).length());             // HDF5 data type for COMPAS data type
                    h5DatasetName = std::get<0>(m_OptionDetails[idx]);                                                  // dataset (option name)
                    h5Dset = CreateHDF5Dataset(h5Filename, h5GroupId, h5DatasetName, h5DataType, "-", chunkSize);       // create dataset
                    if (h5Dset < 0) {                                                                                   // dataset not created
                        Squawk("ERROR: Error creating HDF5 dataset with name " + h5DatasetName);                        // announce error
                        ok = false;                                                                                     // fail
                    }
                    else {                                                                                              // dataset created ok

                        m_Run_Details_H5_File.dataSets.push_back({h5Dset, h5DataType, compasType, STRING_QUALIFIER::FIXED_LENGTH, {}}); // record dataset details

                        // derivation
                        h5DatasetName += "-Derivation";                                                                 // derivation
                        h5Dset = CreateHDF5Dataset(h5Filename, h5GroupId, h5DatasetName, h5String13DataType, "-", chunkSize); // create dataset
                        if (h5Dset < 0) {                                                                               // dataset not created
                            Squawk("ERROR: Error creating HDF5 dataset with name " + h5DatasetName);                    // announce error
                            ok = false;                                                                                 // fail
                        }
                        else
                            m_Run_Details_H5_File.dataSets.push_back({h5Dset, h5String13DataType, TYPENAME::STRING, STRING_QUALIFIER::FIXED_LENGTH, {}}); // dataset created ok - record details
                    }
                    if (!ok) break;                                                                                     // something went wrong - fail
                }
            }
        }   
    }
    return ok;
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
 * Any error here disables logging.
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

    if (!m_Enabled) {                                                                                                       // logging enabled?
                                                                                                                            // no...
        // start timers etc.
        m_WallStart           = std::chrono::system_clock::now();                                                           // start wall timer
        m_ClockStart          = clock();                                                                                    // start CPU timer

        // enable logging
        m_Enabled       = true;                                                                                             // enabled logging
        m_LogBasePath   = p_LogBasePath;                                                                                    // set base path
        m_LogNamePrefix = p_LogNamePrefix;                                                                                  // set log file name prefix
        m_LogLevel      = p_LogLevel;                                                                                       // set log level
        m_LogClasses    = p_LogClasses;                                                                                     // set enabled log classes
        m_DbgLevel      = p_DbgLevel;                                                                                       // set debug level
        m_DbgClasses    = p_DbgClasses;                                                                                     // set enagled debug classes
        m_DbgToLogfile  = p_DbgToLogfile;                                                                                   // write debug records to logfile?
        m_ErrToLogfile  = p_ErrorsToLogfile;                                                                                // write error records to logfile?
        m_LogfileType   = p_LogfileType;                                                                                    // set log file type

        m_Logfiles.clear();                                                                                                 // clear all entries

        m_OptionDetails = OPTIONS->CmdLineOptionsDetails();                                                                 // get commandline option details

        // if the PROGRAM_OPTION::NOTES property is present in the record specification for a logfile at this stage,
        // we want all notes present in the logfile, so set the default for the annotations to true.  Note that this
        // may be changed if a logfile definitions file is present and processed.

        // BSE
        if (NotesPropertyPresent(m_BSE_BE_Binaries_Rec)) m_BSE_BE_Binaries_Notes = std::vector<bool>(m_BSE_BE_Binaries_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_CEE_Rec        )) m_BSE_CEE_Notes         = std::vector<bool>(m_BSE_CEE_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_DCO_Rec        )) m_BSE_DCO_Notes         = std::vector<bool>(m_BSE_DCO_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_Detailed_Rec   )) m_BSE_Detailed_Notes    = std::vector<bool>(m_BSE_Detailed_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_Pulsars_Rec    )) m_BSE_Pulsars_Notes     = std::vector<bool>(m_BSE_Pulsars_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_RLOF_Rec       )) m_BSE_RLOF_Notes        = std::vector<bool>(m_BSE_RLOF_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_SNE_Rec        )) m_BSE_SNE_Notes         = std::vector<bool>(m_BSE_SNE_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_Switch_Rec     )) m_BSE_Switch_Notes      = std::vector<bool>(m_BSE_Switch_Notes.size(), true);
        if (NotesPropertyPresent(m_BSE_SysParms_Rec   )) m_BSE_SysParms_Notes    = std::vector<bool>(m_BSE_SysParms_Notes.size(), true);

        // SSE
        if (NotesPropertyPresent(m_SSE_Detailed_Rec   )) m_SSE_Detailed_Notes    = std::vector<bool>(m_SSE_Detailed_Notes.size(), true);
        if (NotesPropertyPresent(m_SSE_SNE_Rec        )) m_SSE_SNE_Notes         = std::vector<bool>(m_SSE_SNE_Notes.size(), true);
        if (NotesPropertyPresent(m_SSE_Switch_Rec     )) m_SSE_Switch_Notes      = std::vector<bool>(m_SSE_Switch_Notes.size(), true);
        if (NotesPropertyPresent(m_SSE_SysParms_Rec   )) m_SSE_SysParms_Notes    = std::vector<bool>(m_SSE_SysParms_Notes.size(), true);

        // process the logfile definitions file if specified
        m_Enabled = UpdateAllLogfileRecordSpecs();                                                                          // update all logfile record specifications - disable logging upon failure

        if (m_Enabled) {                                                                                                    // still ok?
                                                                                                                            // yes
            // first create the container folder at p_LogBasePath
            // use boost filesystem here - easier...
        
            string containerName = p_LogContainerName;                                                                      // container name
            m_HDF5ContainerName  = p_LogContainerName;                                                                      // HDF5 container name
            string dirName       = containerName;                                                                           // directory name to create

            int version = 0;                                                                                                // container version number if required - start at 1
            while (boost::filesystem::exists(m_LogBasePath + "/" + dirName)) {                                              // container already exists?
                dirName = containerName + "_" + std::to_string(++version);                                                  // yes - add a version number and generate new container name
            }
            m_LogContainerName = dirName;                                                                                   // record actual container directory name

            boost::system::error_code err;
            try {
                boost::filesystem::create_directory(m_LogBasePath + "/" + m_LogContainerName, err);                         // create container - let boost throw an exception if it fails
                if (err.value() == 0) {                                                                                     // ok?

                    if (m_DbgToLogfile) {                                                                                   // write dubug output to a logfile?
                        string filename = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::DEBUG_LOG));                           // extract filename from descriptor
                        int id = Open(filename, false, true, false);                                                        // open the log file - new file, timestamps, no record labels, space delimited
                        if (id >= 0) {                                                                                      // success
                            m_DbgLogfileId = id;                                                                            // record the file id
                        }
                        else {                                                                                              // failure
                            Squawk("ERROR: Unable to create log file for debug output with file name " + filename);         // announce error
                            Squawk("Debug output logging disabled");                                                        // show disabled warning
                        }
                    }

                    if (m_ErrToLogfile) {                                                                                   // write dubug output to a logfile?
                        string filename = std::get<0>(LOGFILE_DESCRIPTOR.at(LOGFILE::ERROR_LOG));                           // extract filename from descriptor
                        int id = Open(filename, false, true, false);                                                        // open the log file - new file, timestamps, no record labels, space delimited
                        if (id >= 0) {                                                                                      // success
                            m_ErrLogfileId = id;                                                                            // record the file id
                        }
                        else {                                                                                              // failure
                            Squawk("ERROR: Unable to create log file for error output with file name " + filename);         // announce error
                            Squawk("Error output logging disabled");                                                        // show disabled warning
                        }
                    }
                }
                else  {                                                                                                     // not ok...
                    Squawk("ERROR: Unable to create log file container with name " + dirName);                              // announce error
                    Squawk("Boost filesystem error = " + err.message());                                                    // plus details
                    Squawk("Logging disabled");                                                                             // show disabled warning
                    m_Enabled = false;                                                                                      // disable
                }
            
            }
            catch (...) {                                                                                                   // unhandled problem...
                Squawk("ERROR: Unable to create log file container with name " + dirName);                                  // announce error
                Squawk("Logging disabled");                                                                                 // show disabled warning
                m_Enabled = false;                                                                                          // disable
            }
        }

        if (m_Enabled) {                                                                                                    // still ok?
                                                                                                                            // yes
            // containing folder now exists
            // now create the run details file

            // if we're logging to HDF5 files We put a copy of the run details file in the HDF5 container file,
            // so if we are in fact logging to HDF5 files, we need to create the HDF5 container file here.

            if (m_LogfileType == LOGFILETYPE::HDF5) {                                                                       // logging to HDF5 files?
                                                                                                                            // yes
                string fileExt = "." + LOGFILETYPEFileExt.at(OPTIONS->LogfileType());                                       // file extension for HDF5 files
                string h5Filename = m_LogBasePath + "/" + m_LogContainerName + "/" + m_HDF5ContainerName + fileExt;         // full filename with path, container, and extension ("/" works on Uni*x and Windows)
                m_HDF5ContainerId = H5Fcreate(h5Filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);                  // create HDF5 container file
                if (m_HDF5ContainerId < 0) {                                                                                // created ok?                        
                    Squawk("ERROR: Unable to create HDF5 container file with file name " + h5Filename);                     // no - announce error
                    Squawk("Logging disabled");                                                                             // show disabled warning
                    m_Enabled = false;                                                                                      // disable logging
                }
                else {                                                                                                      // HDF5 container file now exists and open
                    m_Enabled = OpenHDF5RunDetailsFile(RUN_DETAILS_FILE_NAME);                                              // open HDF5 run details file inside HDF5 container
                    if (!m_Enabled) Squawk("Logging disabled");                                                             // show disabled warning
                }
            }

            // If logging is still enabled then if we are logging to HDF5 files we
            // have an open HDF5 container with an open group for the run details file
            //
            // We still need to create the run details text file
            // We only create the file here - the run details file is populated in Log::Stop()

            if (m_Enabled) {                                                                                                // still ok?
                                                                                                                            // yes
                string filename = m_LogBasePath + "/" + m_LogContainerName + "/" + RUN_DETAILS_FILE_NAME;                   // run details (text) filename with container name
                try {
                    m_RunDetailsFile.open(filename, std::ios::out);                                                         // create run details (text) file
                    m_RunDetailsFile.exceptions(std::ofstream::failbit | std::ofstream::badbit);                            // enable exceptions on run details file
                }
                catch (const std::ofstream::failure &e) {                                                                   // fs problem...
                    Squawk("ERROR: Unable to create run details file with file name " + filename);                          // announce error
                    Squawk(e.what());                                                                                       // plus details
                    m_Enabled = false;                                                                                      // fail
                }
                catch (...) {                                                                                               // unhandled problem...
                    Squawk("ERROR: Unable to create log file with file name " + filename);                                  // announce error
                    m_Enabled = false;                                                                                      // fail
                }
            }

            // store input files if required
            // use Boost to do the copy - copy_file() is available in standard c++17
            if (OPTIONS->StoreInputFiles()) {                                                                               // user wants input files stored in output container?
                                                                                                                            // yes
                string dstPath = m_LogBasePath + "/" + m_LogContainerName + "/";                                            // destination path (output container)
                if (!OPTIONS->GridFilename().empty()) {                                                                     // user specified a grid file?
                    try {                                                                                                   // yes - copy it
                        boost::filesystem::path srcPath(OPTIONS->GridFilename());                                           // grid file fully-qualified name
                        string dstFn = dstPath + srcPath.filename().string();                                               // fully-qualified grid filename (inside container)
                        boost::filesystem::copy_file(OPTIONS->GridFilename(), dstFn, boost::filesystem::copy_option::overwrite_if_exists); // copy grid file - overwrite any existing file (shouldn't be one, but just in case we want this one)
                    } catch(const boost::filesystem::filesystem_error& e) {
                        Squawk("ERROR: Unable to copy grid file " + OPTIONS->GridFilename() + " to output container " + dstPath); // announce error
                        m_Enabled = false;                                                                                  // fail
                    }
                }

                // if the user specified a logfile-definitions file, copy it to the output container

                if (m_Enabled && !OPTIONS->LogfileDefinitionsFilename().empty()) {                                          // user specified a logfile-definitions file?
                    try {                                                                                                   // yes - copy it
                        boost::filesystem::path srcPath(OPTIONS->LogfileDefinitionsFilename());                             // logfile-definitions file fully-qualified name
                        string dstFn = dstPath + srcPath.filename().string();                                               // fully-qualified logfile-definitions filename (inside container)
                        boost::filesystem::copy_file(OPTIONS->LogfileDefinitionsFilename(), dstFn, boost::filesystem::copy_option::overwrite_if_exists); // copy logfile-definitions file - overwrite any existing file (shouldn't be one, but just in case we want this one)
                    } catch(const boost::filesystem::filesystem_error& e) {
                        Squawk("ERROR: Unable to copy logfile-definitions file " + OPTIONS->LogfileDefinitionsFilename() + " to output container " + dstPath); // announce error
                        m_Enabled = false;                                                                                  // fail
                    }
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
 * @param   [IN]    p_ObjectStats               Tuple containg the number of objects requested and the number created
 *                                                 - the number requested is a calculated number: it could just be the number the user requested,
 *                                                   but if a grid file or ranges/sets are used, the number will be calculated.  Furthermore,
 *                                                   the number will be -1 if the simulation was stopped before all grid file entries (or
 *                                                   ranges or sets) were completed - indication we don't really know how many were requested...
 *                                                 - the number created is the actual number created (which may be short of the number requested...)
 */
void Log::Stop(std::tuple<int, int> p_ObjectStats) {

    if (m_Enabled) {                                                                                                                    // only need to do most of this if logging is enabled 

        // get some run stats
     
        double cpuSeconds = (clock() - m_ClockStart) / (double) CLOCKS_PER_SEC;                                                         // stop CPU timer and calculate seconds

        m_WallEnd = std::chrono::system_clock::now();                                                                                   // stop wall timer
        std::time_t timeEnd = std::chrono::system_clock::to_time_t(m_WallEnd);                                                          // get end time and date

        std::chrono::duration<double> wallSeconds = m_WallEnd - m_WallStart;                                                            // elapsed seconds

        int wallHH = (int)(wallSeconds.count() / 3600.0);                                                                               // hours
        int wallMM = (int)((wallSeconds.count() - ((double)wallHH * 3600.0)) / 60.0);                                                   // minutes
        int wallSS = (int)(wallSeconds.count() - ((double)wallHH * 3600.0) - ((double)wallMM * 60.0));                                  // seconds

        std::ostringstream wallTimeSS;
        wallTimeSS << std::setfill('0') << std::setw(4) << wallHH << ":"                                                                // hours, padded with leading '0' if necessary
                   << std::setfill('0') << std::setw(2) << wallMM << ":"                                                                // minutes, padded with leading '0' if necessary
                   << std::setfill('0') << std::setw(2) << wallSS;                                                                      // seconds, padded with leading '0' if necessary
        string wallTime = wallTimeSS.str();

        std::time_t timeStart = std::chrono::system_clock::to_time_t(m_WallStart);                                                      // convert start time

        int objectsRequested = std::get<0>(p_ObjectStats);                                                                              // objects requested (may be -1)
        int objectsCreated   = std::get<1>(p_ObjectStats);                                                                              // objects created

        unsigned long int actualRandomSeed = OPTIONS->FixedRandomSeedCmdLine() ? OPTIONS->RandomSeedCmdLine() : RAND->DefaultSeed();    // actual random seed used

        // update run details file

        if (m_LogfileType == LOGFILETYPE::HDF5) {                                                                                       // logging to HDF5 files?
              
            bool ok = true;                                                                                                             // status
                                                                                                                                        // yes - write run details data to HDF5 output file
            // update run HDF5 details file
            
            string h5DatasetName;
            string derivation;
            int    dSetIdx;

            // preamble/stats datasets

            std::ostringstream ss;

            for ( int idx = static_cast<int>(RUN_DETAILS_COLUMNS::COMPAS_VERSION); idx < static_cast<int>(RUN_DETAILS_COLUMNS::SENTINEL); idx++ ) {

                std::tuple<std::string, TYPENAME, std::size_t> runDetails;
                try { runDetails = RUN_DETAILS_DETAIL.at(static_cast<RUN_DETAILS_COLUMNS>(idx)); }                                      // get run details details
                catch (const std::exception& e) {                                                                                       // unknown property
                    Squawk("ERROR: Unknown property for HDF5 file with name " + OPTIONS->OutputContainerName());                        // announce error
                    ok = false;                                                                                                         // fail
                }
                
                if (ok) {                                                                                                               // have valid property
                    h5DatasetName = std::get<0>(runDetails);                                                                            // dataset name
                    dSetIdx       = idx * 2;
                    derivation    = "CALCULATED";
                    switch (static_cast<RUN_DETAILS_COLUMNS>(idx)) {                                                                    // which dataset?
                        case RUN_DETAILS_COLUMNS::COMPAS_VERSION:                                                                       // COMPAS_Version
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(VERSION_STRING);                                      // add write data to buffer
                            break;

                        case RUN_DETAILS_COLUMNS::RUN_START:                                                                            // Run_Start
                            ss << std::ctime(&timeStart);                                                                               // get start time string
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(ss.str());                                            // add write data to buffer
                            ss.str(std::string());ss.clear();
                            break;
                            
                        case RUN_DETAILS_COLUMNS::RUN_END:                                                                              // Run_End
                            ss << std::ctime(&timeEnd);                                                                                 // get end time string
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(ss.str());                                            // add write data to buffer
                            ss.str(std::string());ss.clear();
                            break;
                            
                        case RUN_DETAILS_COLUMNS::OBJECTS_REQUESTED:                                                                    // Objects_Requested
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(objectsRequested);                                    // add write data to buffer
                            if (objectsRequested >= 0 && (int)OPTIONS->nObjectsToEvolve() == objectsRequested) derivation = "USER_SUPPLIED"; // should be right most of the time (not critical)
                            break;
                            
                        case RUN_DETAILS_COLUMNS::OBJECTS_CREATED:                                                                      // Objects_Created
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(objectsCreated);                                      // add write data to buffer
                            break;
                            
                        case RUN_DETAILS_COLUMNS::CLOCK_TIME:                                                                           // Clock_Time (CPU seconds)
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(cpuSeconds);                                          // add write data to buffer
                            break;
                          
                        case RUN_DETAILS_COLUMNS::WALL_TIME:                                                                            // Wall_Time (elapsed time: hhhh:mm:ss)
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(wallTime);                                            // add write data to buffer
                            break;
                            
                        case RUN_DETAILS_COLUMNS::ACTUAL_RANDOM_SEED:                                                                   // Actual_Random_Seed
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(actualRandomSeed);                                    // add write data to buffer
                            break;

                        default:                                                                                                        // unknown dataset - how did that happen?
                            Squawk("ERROR: Invalid HDF5 dataset with name " + h5DatasetName);                                           // announce error
                            ok = false;                                                                                                 // fail
                    }

                    if (ok) {
                        if (!WriteHDF5_(m_Run_Details_H5_File, RUN_DETAILS_FILE_NAME, dSetIdx)) {                                       // write to file ok?
                            Squawk("ERROR: Error writing to HDF5 dataset with name " + h5DatasetName);                                  // no - announce error
                            ok = false;                                                                                                 // fail
                        }
                        else {                                                                                                          // write succeeded
                            // Derivation
                            dSetIdx += 1;                                                                                               // increment dataset
                            m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::string("CALCULATED"));                           // add write data to buffer
                            if (!WriteHDF5_(m_Run_Details_H5_File, RUN_DETAILS_FILE_NAME, dSetIdx)) {                                   // write to file ok?
                                Squawk("ERROR: Error writing to HDF5 dataset with name " + h5DatasetName);                              // no - announce error
                                ok = false;                                                                                             // fail
                            }
                        }
                    }
                }
                if (!ok) break;                                                                                                         // something went wrong
            }

            if (ok) {
                // program options datasets

                try {
                    for (std::size_t idx = 0; idx < m_OptionDetails.size(); idx++) {                                                    // for eav program option

                        h5DatasetName = std::get<0>(m_OptionDetails[idx]);                                                              // dataset name
                        string strValue = std::get<1>(m_OptionDetails[idx]);                                                            // value formatted as string

                        dSetIdx++;                                                                                                      // incremement run details dataset
                        TYPENAME compasType = std::get<4>(m_OptionDetails[idx]);                                                        // COMPAS datatype
                        switch (compasType) {
                            case TYPENAME::INT         : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stoi(strValue));   break;
                            case TYPENAME::LONGINT     : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stol(strValue));   break;
                            case TYPENAME::LONGLONGINT : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stoll(strValue));  break;
                            case TYPENAME::ULONGINT    : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stoul(strValue));  break;                   
                            case TYPENAME::ULONGLONGINT: m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stoull(strValue)); break;                   
                            case TYPENAME::FLOAT       : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stof(strValue));   break;
                            case TYPENAME::DOUBLE      : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stod(strValue));   break;
                            case TYPENAME::LONGDOUBLE  : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::stold(strValue));  break;
                            case TYPENAME::STRING      : m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(strValue);              break;
                            case TYPENAME::BOOL        :
                                // boolean string value here is "TRUE or "FALSE"
                                // convert to 1 or 0 if necessary
                                if (OPTIONS->PrintBoolAsString())
                                    m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(strValue);
                                else
                                    m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(strValue == "TRUE" ? true : false);
                                break;
                    
                            default:                                                                                                    // invalid datatype
                                Squawk("ERROR: Invalid datatype for HDF5 dataset with name " + h5DatasetName);                          // announce error
                                ok = false;                                                                                             // fail
                        }

                        if (ok) {
                            if (!WriteHDF5_(m_Run_Details_H5_File, RUN_DETAILS_FILE_NAME, dSetIdx)) {                                   // write to file ok?
                                Squawk("ERROR: Error writing to HDF5 dataset with name " + h5DatasetName);                              // no - announce error
                                ok = false;                                                                                             // fail
                            }
                            else {                                                                                                      // write succeeded
                            // Derivation
                                dSetIdx += 1;                                                                                           // incremement dataset
                                m_Run_Details_H5_File.dataSets[dSetIdx].buf.push_back(std::string(std::get<2>(m_OptionDetails[idx])));  // add write data to buffer
                                if (!WriteHDF5_(m_Run_Details_H5_File, RUN_DETAILS_FILE_NAME, dSetIdx)) {                               // write to file ok?
                                    Squawk("ERROR: Error writing to HDF5 dataset with name " + h5DatasetName);                          // no - announce error
                                    ok = false;                                                                                         // fail
                                }
                            }
                        }
                        if (!ok) break;                                                                                                 // something went wrong
                    }
                } 
                catch (const std::out_of_range& e) {                                                                                    // type conversion failed
                    Squawk("ERROR: Error converting option value to correct datatype for HDF5 dataset with name " + h5DatasetName);     // announce error
                    ok = false;
                }
                catch (const std::invalid_argument& e) {                                                                                // type conversion failed
                    Squawk("ERROR: Error converting option value to correct datatype for HDF5 dataset with name " + h5DatasetName);     // announce error
                    ok = false;
                }
            }
        }

        // update run details text file
        string filename = m_LogBasePath + "/" + m_LogContainerName + "/" + RUN_DETAILS_FILE_NAME;                                       // run details filename with container name
        try {  
            m_RunDetailsFile << utils::SplashScreen(false) << std::endl;                                                                // write splash string with version number to file

            // record start time and whether evolving single stars or binaries   
            if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE)
                m_RunDetailsFile << "Start generating stars at " << std::ctime(&timeStart) << std::endl;
            else
                m_RunDetailsFile << "Start generating binaries at " << std::ctime(&timeStart) << std::endl;

            // record end time and whether evolving single stars or binaries   
            if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE) {
                m_RunDetailsFile << "Generated " << std::to_string(objectsCreated) << " of " << (objectsRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(objectsRequested)) << " stars requested" << std::endl;
                m_RunDetailsFile << "\nEnd generating stars at " << std::ctime(&timeEnd) << std::endl;
            }
            else {
                m_RunDetailsFile << "Generated " << std::to_string(objectsCreated) << " of " << (objectsRequested < 0 ? "<INCOMPLETE GRID>" : std::to_string(objectsRequested)) << " binaries requested" << std::endl;
                m_RunDetailsFile << "\nEnd generating binaries at " << std::ctime(&timeEnd) << std::endl; 
            }

            m_RunDetailsFile << "Clock time = " << cpuSeconds << " CPU seconds" << std::endl;                                           // record cpu seconds

            m_RunDetailsFile << "Wall time  = " << wallTime << " (hhhh:mm:ss)" << std::endl;                                            // wall time 

            // add commandline options
            // moved this code here from Options.cpp
            // have to add a small kludge here to get it to look the same (someone might be relying on format)
            // (run through twice - first time specified options, second time calculated options - only need to do once per run, so not a disaster...)

            // first, specified options

            m_RunDetailsFile << "\n\nCOMMAND LINE OPTIONS\n--------------------\n\n";                                                   // add commandline options (all of them...)
            for (std::size_t idx = 0; idx < m_OptionDetails.size(); idx++) {                                                            // and add them to the run details file

                if (utils::Equals(std::get<2>(m_OptionDetails[idx]), "CALCULATED")) continue;                                           // CALCULATED later

                m_RunDetailsFile << std::get<0>(m_OptionDetails[idx]) << " = ";                                                         // option name

                if (std::get<1>(m_OptionDetails[idx]) == "")                                                                            // empty option?
                    m_RunDetailsFile << "<EMPTY_OPTION>\n";                                                                             // yes - say so
                else                                                                                                                    // no - add option details
                    m_RunDetailsFile << std::get<1>(m_OptionDetails[idx]) + ", "                                                        // value
                                     << std::get<2>(m_OptionDetails[idx]) + ", "                                                        // defaulted
                                     << std::get<3>(m_OptionDetails[idx]) << "\n";                                                      // datatype
            }

            // next, calculated options

            m_RunDetailsFile << "\n\nOTHER PARAMETERS\n----------------\n\n";
            for (std::size_t idx = 0; idx < m_OptionDetails.size(); idx++) {                                                            // and add them to the run details file

                if (!utils::Equals(std::get<2>(m_OptionDetails[idx]), "CALCULATED")) continue;                                          // only CALCULATED here

                m_RunDetailsFile << std::get<0>(m_OptionDetails[idx]) << " = ";                                                         // option name

                if (std::get<1>(m_OptionDetails[idx]) == "")                                                                            // empty option?
                    m_RunDetailsFile << "<EMPTY_OPTION>\n";                                                                             // yes - say so
                else                                                                                                                    // no - add option details
                    m_RunDetailsFile << std::get<1>(m_OptionDetails[idx]) + ", "                                                        // value
                                     << std::get<2>(m_OptionDetails[idx]) + ", "                                                        // defaulted
                                     << std::get<3>(m_OptionDetails[idx]) << "\n";                                                      // datatype
            }

            m_RunDetailsFile << "Actual random seed = " << actualRandomSeed  << ", CALCULATED, UNSIGNED_LONG" << std::endl;             // actual random seed

            // done writing - flush and close the file
            try {
                m_RunDetailsFile.flush();
                m_RunDetailsFile.close();
            }
            catch (const std::ofstream::failure &e) {                                                                                   // problem...
                Squawk("ERROR: Unable to close run details file with file name " + filename);                                           // announce error
                Squawk(e.what());                                                                                                       // plus details
            }
        }
        catch (const std::ofstream::failure &e) {                                                                                       // problem...
            Squawk("ERROR: Unable to write to run details file with name " + filename);                                                 // announce error
            Squawk(e.what());                                                                                                           // plus details
        }

        // close standard log files

        CloseAllStandardFiles();                                                                                                        // close all standard log files
        for(unsigned int index = 0; index < m_Logfiles.size(); index++) {                                                               // check for open logfiles (even if not active)
            if (IsActiveId(index)) {                                                                                                    // logfile active?
                if (m_Logfiles[index].file.is_open()) {                                                                                 // open file?
                    try {                                                                                                               // yes
                        m_Logfiles[index].file.flush();                                                                                 // flush output and
                        m_Logfiles[index].file.close();                                                                                 // close it
                    }
                    catch (const std::ofstream::failure &e) {                                                                           // problem...
                        Squawk("ERROR: Unable to close log file with file name " + m_Logfiles[index].name);                             // announce error
                        Squawk(e.what());                                                                                               // plus details
                    }
                }
            }
        }
    }

    m_Logfiles.clear();                                                                                                                 // clear all entries
    m_Enabled = false;                                                                                                                  // set not enabled
}


/*
 * Create and open new log file
 *
 * New log file is created at path m_LogBasePath
 *
 * The file extension anmd delimier are based on the file type of the logfile:
 *     - CSV : the file extension is "csv" (Comma Separated Variables), and the delimiter is the COMMA character (",")
 *     - TSV : the file extension is "tsv" (Tab Separated Variables), and the delimiter is the TAB character ("\t")
 *     - TXT : the file extension is "txt" (Plain text file), and the delimiter is the SPACE character (" ")
 *     - HDF5: the file extension is "h5"  (Hierarchical Data Format, version 5).  HDF5 files are not delimited.
 *
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

    if (m_Enabled) {                                                                                                // logging enabled?   

        string basename = m_LogBasePath + "/" + m_LogContainerName + "/" + m_LogNamePrefix + p_LogFileName;         // base filename with path and container ("/" works on Uni*x and Windows)
        string fileext  = LOGFILETYPEFileExt.at(OPTIONS->LogfileType());                                            // file extension
        string filename = basename + "." + fileext;                                                                 // full filename

        int version = 0;                                                                                            // logfile version number if required - start at 1
        while (utils::FileExists(filename) && !p_Append) {                                                          // file already exists - and we don't want to append?
            filename = basename + "_" + std::to_string(++version) + "." + fileext;                                  // yes - add a version number and generate new filename
        }

        if (m_LogfileType == LOGFILETYPE::HDF5) {                                                                   // HDF5 file?
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

            if (ok) {                                                                                               // still ok?

                // record attributes
                // find an empty slot in m_Logfiles vector if there is one
                // this way is a bit slower for opening logfiles, but faster for writing to them
                id = -1;  
                for(unsigned int index = 0; index < m_Logfiles.size(); index++) {
                    if (!m_Logfiles[index].active) {                                                                // empty slot?
                        id = index;                                                                                 // yes - use it
                        break;                                                                                      // and stop looking
                    }
                }

                if (id < 0) {                                                                                       // have empty slot?
                    logfileAttrT attr;                                                                              // no - create new attributes struct
                    id = m_Logfiles.size();                                                                         // set new id
                    m_Logfiles.push_back(std::move(attr));                                                          // append to vector
                }

                m_Logfiles[id].active           = true;                                                             // this entry now active
                m_Logfiles[id].logfiletype      = p_StandardLogfile;                                                // standard logfile type
                m_Logfiles[id].filetype         = m_LogfileType;                                                    // set filetype for this log file
                m_Logfiles[id].name             = h5GroupName;                                                      // log file name (HDF5 group name)
                m_Logfiles[id].timestamp        = p_Timestamp;                                                      // set timestamp flag for this log file
                m_Logfiles[id].label            = p_Label;                                                          // set label flag for this log file
                m_Logfiles[id].h5File.fileId    = h5FileId;                                                         // HDF5 file id
                m_Logfiles[id].h5File.groupId   = h5GroupId;                                                        // HDF5 group id
                m_Logfiles[id].h5File.dataSets  = {};                                                               // HDF5 data sets
            }
        }
        else {                                                                                                      // no - not HDF5
            // record attributes

            // find an empty slot in m_Logfiles vector if there is one
            // this way is a bit slower for opening logfiles, but faster for writing to them
            for(unsigned int index = 0; index < m_Logfiles.size(); index++) {
                if (!m_Logfiles[index].active) {                                                                    // empty slot?
                    id = index;                                                                                     // yes - use it
                    break;                                                                                          // and stop looking
                }
            }

            if (id < 0) {                                                                                           // have empty slot?
                id = m_Logfiles.size();                                                                             // no - set new id
                logfileAttrT attr;                                                                                  // create new attributes struct
                attr.active = false;                                                                                // for now...
                m_Logfiles.push_back(std::move(attr));                                                              // append to vector
            }
               
            try {
                m_Logfiles[id].file.open(filename, std::ios::out | std::ios::app);                                  // create fs log file
                m_Logfiles[id].file.exceptions(std::ofstream::failbit | std::ofstream::badbit);                     // enable exceptions on log file

                m_Logfiles[id].active         = true;                                                               // this entry now active
                m_Logfiles[id].logfiletype    = p_StandardLogfile;                                                  // standard logfile type
                m_Logfiles[id].filetype       = m_LogfileType;                                                      // set filetype for this log file
                m_Logfiles[id].name           = filename;                                                           // log file name
                m_Logfiles[id].timestamp      = p_Timestamp;                                                        // set timestamp flag for this log file
                m_Logfiles[id].label          = p_Label;                                                            // set label flag for this log file
                m_Logfiles[id].h5File.fileId  = -1;                                                                 // not HDF5 file
                m_Logfiles[id].h5File.groupId = -1;                                                                 // not HDF5 file
            }
            catch (const std::ofstream::failure &e) {                                                               // fs problem...
                Squawk("ERROR: Unable to create fs log file with file name " + m_Logfiles[id].name);                // announce error
                Squawk(e.what());                                                                                   // plus details

                // it's possible this m_Logfiles entry was just appended - it can be used next time
                if (id >= 0) ClearEntry(id);                                                                        // clear entry
                ok = false;                                                                                         // fail
            }
        }
    }

    return ok ? id : -1;                                                                                            // return log file id
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
            Flush_(p_LogfileId);                                                                                    // first, flush any unwritten data

            // close all open datasets
            for (auto &dataSet : m_Logfiles[p_LogfileId].h5File.dataSets) {                                         // for each dataset
                if (dataSet.dataSetId >= 0) {                                                                       // dataset open?
                                                                                                                    // yes - close it
                    if (H5Dclose(dataSet.dataSetId) < 0) {                                                          // closed ok?
                        Squawk("ERROR: Unable to close HDF5 dataset for log file " + m_Logfiles[p_LogfileId].name); // no - announce error
                        result = false;                                                                             // fail
                    }
                }
            }

            // try closing the group - even if closing datasets failed
            if (m_Logfiles[p_LogfileId].logfiletype != LOGFILE::SSE_DETAILED_OUTPUT && 
                m_Logfiles[p_LogfileId].logfiletype != LOGFILE::BSE_DETAILED_OUTPUT) {                              // detailed output file?

                if (m_Logfiles[p_LogfileId].h5File.groupId >= 0) {                                                  // no - HDF5 group open?
                    if (H5Gclose(m_Logfiles[p_LogfileId].h5File.groupId) < 0) {                                     // yes - closed ok?
                        Squawk("ERROR: Unable to close HDF5 group for log file " + m_Logfiles[p_LogfileId].name);   // no - announce error
                        result = false;                                                                             // fail
                    }
                }
            }
            else {                                                                                                  // yes - detailed output file
                if (H5Fclose(m_Logfiles[p_LogfileId].h5File.fileId) < 0) {                                          // HDF5 file closed ok?
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
 * @param   [IN]    p_Level                     The level (logging or debug) of the record being evaluated
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
 * Write a string record to specified log file
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
        else result = true;                                                                                         // not logging this class and level - but ok
    }

    return result;
}


/*
 * Write a multi-value record to specified log file
 *
 * If logging is enabled and the specified log file is active, and class and level are enabled (see DoIt()),
 * write the log record values to log file
 *
 * Disable the specified log file if errors occur.
 *
 *
 * bool Write(const id                                p_LogfileId, 
 *            const string                            p_LogClass, 
 *            const int                               p_LogLevel, 
 *            const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, 
 *            const bool                              p_Flush)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log string should be written
 * @param   [IN]    p_LogClass                  Class to determine if string should be written
 * @param   [IN]    p_LogLevel                  Level to determine if string should be written
 * @param   [IN]    p_LogRecordValues           Vector of COMPAS_VARIABLE_TYPE values to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Write(const int                               p_LogfileId, 
                const string                            p_LogClass, 
                const int                               p_LogLevel, 
                const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, 
                const bool                              p_Flush) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging service enabled and specified log file active?
        if (DoIt(p_LogClass, p_LogLevel, m_LogClasses, m_LogLevel)) {                                               // yes - logging this class and level?
            result = Write_(p_LogfileId, p_LogRecordValues, p_Flush);                                               // yes - log it
        }
        else result = true;                                                                                         // not logging this class and level - but ok
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
 * Write data from buffer to HDF5 files
 * 
 * This is where the real work is done for HDF5 files
 * 
 * Note that the first parameter, p_H5file, will be modified
 * (the contents of the write buf will be cleared after writing)
 *
 *
 * bool WriteHDF5_(h5AttrT& p_H5file, const string p_H5filename, const size_t p_DataSetIdx)
 *
 * @param   [IN]    p_H5file                    Struct containing details of the HDF5 file to which the buffer should be written - contains the buffer to write
 * @param   [IN]    p_H5filename                String filename of the HDF5 file - for error logging should an error occur
 * @param   [IN]    p_DataSetIdx                Index of the dataset within the HDF5 file to which the buf should be written (assumed to exist and be open)
 * @return                                      Boolean indicating whether buffer was written successfully
 */
bool Log::WriteHDF5_(h5AttrT& p_H5file, const string p_H5filename, const size_t p_DataSetIdx) {

    herr_t ok = 0;                                                                                                          // return value

    // setup write:                                                                                 
    //    - create dataspaces
    //    - extend dataset
    //    - setup hyperslab

    size_t  bufSize         = p_H5file.dataSets[p_DataSetIdx].buf.size();                                                   // size of write buffer

    hid_t   dSet            = p_H5file.dataSets[p_DataSetIdx].dataSetId;                                                    // dataset id
    hid_t   dType           = p_H5file.dataSets[p_DataSetIdx].h5DataType;                                                   // HDF5 datatype
    hsize_t dSetCurrentSize = H5Dget_storage_size(dSet) / H5Tget_size(dType);                                               // current size (entries) of HDF5 dataset
    hsize_t h5Dims[1]       = {bufSize};                                                                                    // size of buffer to be written
    hid_t   h5Dspace        = H5Screate_simple(1, h5Dims, NULL);                                                            // create memory dataspace for write
    hid_t   h5FSpace;                                                                                                       // filespace for write - allocated later

    if (h5Dspace < 0) {                                                                                                     // created ok?
        Squawk("ERROR: Unable to allocate memory to write to HDF5 group for log file " + p_H5filename);                     // no - announce error
        ok = -1;                                                                                                            // fail
    }
    else {                                                                                                                  // yes - memory dataspace created ok                           
        h5Dims[0] = dSetCurrentSize + bufSize;                                                                              // new size for dataset
        if ((ok = H5Dset_extent(dSet, h5Dims)) < 0) {                                                                       // extend dataset - ok?
            Squawk("ERROR: Unable to extend file to write to HDF5 group for log file " + p_H5filename);                     // no - announce error
        }
        else {                                                                                                              // yes - dataset extended ok
            h5FSpace = H5Dget_space(dSet);                                                                                  // allocate filespace
            if (h5FSpace < 0) {                                                                                             // allocated ok?
                Squawk("ERROR: Unable to allocate filespace to write to HDF5 group for log file " + p_H5filename);          // no - announce error
                ok = -1;                                                                                                    // fail
            }
            else {                                                                                                          // yes - filespace allocated ok

                // setup hyperslab on file space
                // the hyperslab mirrors the dataset:
                //    - start is the start position in the hyperslab of the write
                //    - count is the number of entries to write (the chunk size)
                hsize_t h5Start[1] = {dSetCurrentSize};
                hsize_t h5Count[1] = {bufSize};
                if ((ok = H5Sselect_hyperslab(h5FSpace, H5S_SELECT_SET, h5Start, NULL, h5Count, NULL)) < 0) {               // hyperslab setup ok?
                    Squawk("ERROR: Unable to set location to write to HDF5 group for log file " + p_H5filename);            // no - announce error
                }
            }                                    
        }
    }

    // setup done - if no errors, write the data to the file

    if (ok >= 0) {                                                                                                          // good to write?
                                                                                                                            // yes
        // can't use switch here - the HDF5 constants are not really constants...
        // for each datatype:
        //   - create a buffer of correct C++
        //   - extract the boost variant values from the write buffer and populate the C++-typed buffer
        //   - clear the write buffer
        //   - write the C++-typed buffer to the file

        if (dType == H5T_NATIVE_UCHAR) {
            bool buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<bool>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);                                  // guaranteed to release memory
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_SHORT) {
            short int buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<short int>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_INT) {
            // enum classes are cast to type int - need to cast and extract here
            int buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) {
                int v = 0;
                switch (p_H5file.dataSets[p_DataSetIdx].dataType) {
                    case TYPENAME::INT         : v = static_cast<int>(boost::get<int>(p_H5file.dataSets[p_DataSetIdx].buf[i])); break;
                    case TYPENAME::ERROR       : v = static_cast<int>(boost::get<ERROR>(p_H5file.dataSets[p_DataSetIdx].buf[i])); break;
                    case TYPENAME::STELLAR_TYPE: v = static_cast<int>(boost::get<STELLAR_TYPE>(p_H5file.dataSets[p_DataSetIdx].buf[i])); break;
                    case TYPENAME::MT_CASE     : v = static_cast<int>(boost::get<MT_CASE>(p_H5file.dataSets[p_DataSetIdx].buf[i])); break;
                    case TYPENAME::MT_TRACKING : v = static_cast<int>(boost::get<MT_TRACKING>(p_H5file.dataSets[p_DataSetIdx].buf[i])); break;
                    case TYPENAME::SN_EVENT    : v = static_cast<int>(boost::get<SN_EVENT>(p_H5file.dataSets[p_DataSetIdx].buf[i])); break;
                    case TYPENAME::SN_STATE    : v = static_cast<int>(boost::get<SN_STATE>(p_H5file.dataSets[p_DataSetIdx].buf[i])); break;
                    default: 
                        Squawk("ERROR: Unable to format data to write to HDF5 group for log file " + p_H5filename);         // announce error
                        ok = -1;                                                                                            // fail
                }
                buf[i] = v;
            }
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            if (ok >=0) {                                                                                                   // data formatted ok?
                ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);                            // yes - write it
            }
        }
        else if (dType == H5T_NATIVE_LONG) {
            long int buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<long int>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_USHORT) {
            unsigned short int buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<unsigned short int>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_UINT) {
            unsigned int buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<unsigned int>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_ULONG) {
            unsigned long int buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<unsigned long int>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_FLOAT) {
            float buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<float>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_DOUBLE) {
            double buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<double>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (dType == H5T_NATIVE_LDOUBLE) {
            long double buf[bufSize];
            for (size_t i = 0; i < bufSize; i++) buf[i] = boost::get<long double>(p_H5file.dataSets[p_DataSetIdx].buf[i]);
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);
            ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)&buf);
        }
        else if (p_H5file.dataSets[p_DataSetIdx].dataType == TYPENAME::STRING) {

            // here we handle strings - both fixed length and variable length
            // not that HDF5 functions don't know about std::string - they need c-style char arrays

            bool fixedLength = p_H5file.dataSets[p_DataSetIdx].stringType == STRING_QUALIFIER::FIXED_LENGTH;

            string buf[bufSize];
            size_t elemLen = H5Tget_size(dType) - 1;                                                                        // for fixed-length strings (-1 for null terminator)

            // format each string:
            // format bool as "TRUE" or "FALSE" if required
            // for fixed-length strings, pad the string to the defined length

            for (size_t i = 0; i < bufSize; i++) {
                // if user specified "print-bool-as-string" option, need to translate bool value to "TRUE" or "FALSE"
                string v = p_H5file.dataSets[p_DataSetIdx].dataType == TYPENAME::BOOL                                       // bool variable (printing as string "TRUE" or "FALSE")?
                            ? boost::get<bool>(p_H5file.dataSets[p_DataSetIdx].buf[i]) ? string("TRUE") : string("FALSE")   // yes
                            : boost::get<string>(p_H5file.dataSets[p_DataSetIdx].buf[i]);                                   // no, regular STRING variable
                        
                buf[i] = fixedLength ? utils::PadTrailingSpaces(v, elemLen) : v;
            }
            std::vector<COMPAS_VARIABLE_TYPE>().swap(p_H5file.dataSets[p_DataSetIdx].buf);

            // write the strings - need c-style char array
            if (fixedLength) {                                                                                              // fixed-length string?
                char *cBuf = new char[bufSize * (buf[0].length() + 1)];                                                     // H5Dwrite() needs contiguous memory
                size_t pos = 0;                                                                                             // start position
                for (size_t i = 0; i < bufSize; i++) {                                                                      // for each entry in the buffer
                    strcpy(cBuf + pos, buf[i].c_str());                                                                     // copy chars + null terminator
                    pos += buf[i].length() + 1;                                                                             // next start position
                }

                ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)cBuf);                            // write the data
                delete[] cBuf;                                                                                              // release allocated memory
            }
            else {                                                                                                          // no - variable length
                char* cBuf[bufSize];                                                                                        // char array for H5Dwrite()
                for (size_t i = 0; i < bufSize; i++) {                                                                      // for each entry in the buffer
                    cBuf[i] = new char[buf[i].length() + 1];                                                                // H5Dwrite() needs contiguous memory only for each element
                    strcpy(cBuf[i], buf[i].c_str());                                                                        // copy chars + null terminator
                }

                ok = H5Dwrite(dSet, dType, h5Dspace, h5FSpace, H5P_DEFAULT, (const void *)cBuf);;                            // write the data
            
                // release allocated memory
                for (size_t i = 0; i < bufSize; i++) {
                    delete[] cBuf[i];
                }
            }
        }
        else {
            Squawk("ERROR: Unable to format data to write to HDF5 group for log file " + p_H5filename);                     // announce error
            ok = -1;                                                                                                        // fail
        }

        (void)H5Sclose(h5FSpace);                                                                                           // close filespace
        (void)H5Sclose(h5Dspace);                                                                                           // close dataspace
    }
    return (ok >= 0);
}


/*
 * Write a multi-value record to specified log file with no class or level check - internal use only
 * Used for HDF5 files
 * 
 * This is where (most of) the work is done
 *
 * Disable the specified log file if errors occur.
 *
 *
 * bool Write_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const bool p_Flush)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the multi-value record should be written
 * @param   [IN]    p_LogRecordValues           Vector of COMPAS_VARIABLE_TYPE values to be written
 * @param   [IN]    p_Flush                     Boolean indicating whether the writebuffer should be flushed regardless of chunk size (optional, default = false)
 *                                              If p_Flush is true, no data is added to the write buffer, and the entire write buffer is written to the file
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Write_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const bool p_Flush) {

    herr_t ok = 0;

    if (m_Logfiles[p_LogfileId].filetype != LOGFILETYPE::HDF5) return ok;                                                   // shouldn't be here if not HDF5 logfile

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                             // logging service enabled and specified log file active?
        if ((ok = m_Logfiles[p_LogfileId].h5File.fileId) < 0) {                                                             // HDF5 file open?
            Squawk("ERROR: Unable to write to HDF5 log file " + m_Logfiles[p_LogfileId].name);                              // no - announce error
        }
        else {
            if (m_Logfiles[p_LogfileId].h5File.groupId < 0) {                                                               // HDF5 group open?
                Squawk("ERROR: Unable to write to HDF5 group for log file " + m_Logfiles[p_LogfileId].name);                // no - announce error
            }
            else {

                for (size_t idx = 0; idx < m_Logfiles[p_LogfileId].h5File.dataSets.size(); idx++) {                         // for each dataset

                    hid_t dSet  = m_Logfiles[p_LogfileId].h5File.dataSets[idx].dataSetId;                                   // dataset id

                    if (dSet >= 0) {                                                                                        // dataset open?
                                                                                                                            // yes
                        if (!p_Flush) {                                                                                     // flush only?
                            m_Logfiles[p_LogfileId].h5File.dataSets[idx].buf.push_back(p_LogRecordValues[idx]);             // no - add write data to buffer
                        }

                        if ((m_Logfiles[p_LogfileId].h5File.dataSets[idx].buf.size() >= m_Logfiles[p_LogfileId].h5File.IOBufSize) || p_Flush) { // need to write?
                            ok = WriteHDF5_(m_Logfiles[p_LogfileId].h5File, m_Logfiles[p_LogfileId].name, idx);             // do the write 
                        }
                    }
                }
            }
        }
    }
    else {                                                                                                                  // logging not enabled or not active          

        // construct a log record and display it on stderr
        string logRecord = "";        

        for (auto &value : p_LogRecordValues) {
            string valueStr = boost::apply_visitor(FormatVariantValueDefault(), value);                                     // format value
            logRecord += valueStr + ",";                                                                                    // append to output string + delimiter
        }
        logRecord = logRecord.substr(0, logRecord.size()-1);                                                                // remove the last character - extraneous delimiter

        Squawk(logRecord);                                                                                                  // show log record on stderr
        ok = -1;                                                                                                            // fail
    }

    return (ok >= 0);
}


/*
 * Write a minimally formatted record to the specified log file.  Used for CSV, TSV, and TXT files.
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
        else result = true;                                                                                         // not logging this class and level - but ok
    }

    return result;
}


/*
 * Write a multi-value record to specified log file.
 *
 * If logging is enabled and the specified log file is active, and class and level are enabled (see DoIt()),
 * write the log record values to the log file
 *
 * Disable the specified log file if errors occur.
 *
 *
 * bool Put(const id p_LogfileId, const string p_LogClass, const int p_LogLevel, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log record values should be written
 * @param   [IN]    p_LogClass                  Class to determine if log record values sould be written
 * @param   [IN]    p_LogLevel                  Level to determine if log record values should be written
 * @param   [IN]    p_LogRecordValues           Vector of COMPAS_VARIABLE_TYPE values to be written
 * @return                                      Boolean indicating whether log record values were written successfully
 */
bool Log::Put(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {                                                                     // logging service enabled and specified log file active?
        if (DoIt(p_LogClass, p_LogLevel, m_LogClasses, m_LogLevel)) {                                               // yes - logging this class and level?
            result = Put_(p_LogfileId, p_LogRecordValues);                                                          // yes - log it
        }
        else result = true;                                                                                         // not logging this class and level - but ok
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
 * bool Put_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues)
 *
 * @param   [IN]    p_LogfileId                 The id of the log file to which the log string should be written
 * @param   [IN]    p_LogRecordValues           Vector of COMPAS_VARIABLE_TYPE values to be written
 * @return                                      Boolean indicating whether record was written successfully
 */
bool Log::Put_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues) {

    bool result = false;

    if (m_Enabled && IsActiveId(p_LogfileId)) {             // logging service enabled and specified log file active?
        return Write_(p_LogfileId, p_LogRecordValues);      // yes - log it
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
        else result = true;                                                                                         // not logging this class and level - but ok
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
        else result = true;                                                                                         // not logging this class and level - but ok
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
 * PROPERTY_DETAILS StellarPropertyDetails(const ANY_STAR_PROPERTY p_Property)
 *
 * @param   [IN]    p_Property                  The property for which the details are required
 * @return                                      Tuple containing the properties (default properties if p_Property not found)
 */
PROPERTY_DETAILS Log::StellarPropertyDetails(const ANY_STAR_PROPERTY p_Property) {

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
 * PROPERTY_DETAILS BinaryPropertyDetails(const BINARY_PROPERTY p_Property)
 *
 * @param   [IN]    p_Property                  The property for which the details are required
 * @return                                      Tuple containing the properties (default properties if p_Property not found)
 */
PROPERTY_DETAILS Log::BinaryPropertyDetails(const BINARY_PROPERTY p_Property) {

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
 * { TYPENAME, Header string, Units string, field width, precision }
 *
 *
 * PROPERTY_DETAILS ProgramOptionDetails(const ProgramOptionDetails p_Property, const size_t p_Idx)
 *
 * @param   [IN]    p_Property                  The property for which the details are required
 * @param   [IN]    p_Idx                       Index for vector properties (only PROGRAM_OPTION::NOTES currently)
 *                                              Defaults to 0
 * @return                                      Tuple containing the properties (default properties if p_Property not found)
 */
PROPERTY_DETAILS Log::ProgramOptionDetails(const PROGRAM_OPTION p_Property, const size_t p_Idx) {

    PROPERTY_DETAILS details;

    // PROGRAM_OPTION::NOTES is special
    // There are no entry for PROGRAM_OPTION::NOTES in the PROGRAM_OPTION_DETAIL map
    // (at least there shouldn't be - and if there is we just ignore it)
    // We construct the PROPERTY_DETAILS tuple for PROGRAM_OPTION::NOTES here

    if (p_Property == PROGRAM_OPTION::NOTES) {                                              // PROGRAM_OPTION::NOTES?
        details = std::make_tuple(TYPENAME::STRING, OPTIONS->NotesHdrs(p_Idx), "-", 0, 1);  // yes - construct details
    }
    else {                                                                                  // not PROGRAM_OPTION::NOTES
        try { details = PROGRAM_OPTION_DETAIL.at(p_Property); }                             // get program option details
        catch (const std::exception& e) {                                                   // unknown property
            details = std::make_tuple(TYPENAME::NONE, "", "", 0, 0);                        // empty details
            DBG_WARN(ERR_MSG(ERROR::UNKNOWN_PROGRAM_OPTION));                               // show warning
        }
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
 *             constructed here because the final field width is determined here.
 *
 *
 * STR_STR_STR_STR FormatFieldHeaders(const PROPERTY_DETAILS p_PropertyDetails, string p_HeaderSuffix)
 *
 * @param   [IN]    p_PropertyDetails           The property details for the property for which the headers are to be formatted
 * @param   [IN]    p_HeaderSuffix              The suffix string to be appended to the header string
 * @return                                      Tuple containing formatted strings for the property requested: <header, units, type, format>
 *                                              If the property details are not valid (e.g. unknown data type), error strings will be returned
 */
STR_STR_STR_STR Log::FormatFieldHeaders(const PROPERTY_DETAILS p_PropertyDetails, string p_HeaderSuffix) {

    TYPENAME typeName = std::get<0>(p_PropertyDetails);                                                                 // data type
    if (typeName == TYPENAME::NONE) {                                                                                   // valid data type?
        return std::make_tuple("ERROR!", "ERROR!", "ERROR!", "ERROR!");                                                 // return error values
    }

    string headerStr = std::get<1>(p_PropertyDetails) + p_HeaderSuffix;                                                 // header string
    string unitsStr  = std::get<2>(p_PropertyDetails);                                                                  // units string
    string typeStr   = std::get<1>(TYPENAME_LABEL.at(typeName));                                                        // type will be one of "BOOL", "INT", "FLOAT" and "STRING" (non-primitive types coerced to INT)

    int fieldWidth     = std::get<3>(p_PropertyDetails);
    int fieldPrecision = std::get<4>(p_PropertyDetails);                                                                // field precision (for double and int)

    string fmtStr;
    if (typeName == TYPENAME::STRING) {                                                                                 // type string
        fmtStr = fieldWidth == 0 ? "" : std::to_string(fieldWidth);                                                     // precision not used for strings
    }
    else {                                                                                                              // no - fieldwidth limited
        fieldWidth = std::max(fieldWidth, std::max((int)headerStr.length(), std::max((int)unitsStr.length(), 6)));      // maximum of requested width, header width, units width and type width ("STRING" is max type)

        headerStr = utils::CentreJustify(headerStr, fieldWidth);                                                        // centre-justify header string
        unitsStr  = utils::CentreJustify(unitsStr, fieldWidth);                                                         // centre-justify units string
        typeStr   = utils::CentreJustify(typeStr, fieldWidth);                                                          // centre-justify type string

        fmtStr    = std::to_string(fieldWidth) + "." + std::to_string(fieldPrecision);                                  // field width and precision spcifiers
    }

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
 * Get standard log file record properties, format vector, and annotations, from the logfile
 * record specifier
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
 * it in two places.  Although, this function needs to be able to get the properties even if
 * the file is not open...  Still, should be able to put the code in one place instead of two.
 * But Log::StandardLogFileDetails() is, the way it was initially written, a bit too complex 
 * (this whole flexible printing code is a complex beast - unfortunately it  has to be to get
 * it to work) and this code is a bit too intertwined to easily and quickly disentangle it from
 * Log::StandardLogFileDetails() - that's probably a good code cleanup to do some time in the 
 * future, but for now this will have to suffice.
 *
 *
 * std::tuple<ANY_PROPERTY_VECTOR, std::vector<string>, std::vector<bool>> GetStandardLogFileRecordDetails(const LOGFILE p_Logfile)
 *
 * @param   [IN]    p_Logfile                   Logfile for which details are to be retrieved (see enum class LOGFILE in constants.h)
 * @return                                      Tuple containing:
 *                                                 - vector of logfile record properties
 *                                                 - vector of format strings
 *                                                 - vector of annotations
 */
std::tuple<ANY_PROPERTY_VECTOR, std::vector<string>, std::vector<bool>> Log::GetStandardLogFileRecordDetails(const LOGFILE p_Logfile) {

    ANY_PROPERTY_VECTOR  recordProperties = {};                                                                                     // default is empty
    std::vector<string>  fmtVector = {};                                                                                            // default is empty
    std::vector<bool>  annotations = {};                                                                                            // default is empty

    try {
        // get record properties for this file

        switch (p_Logfile) {                                                                                                        // which logfile?

            case LOGFILE::BSE_BE_BINARIES:                                                                                          // BSE_BE_BINARIES
                recordProperties = m_BSE_BE_Binaries_Rec;                                                                           // record properties
                annotations      = m_BSE_BE_Binaries_Notes;                                                                         // logfile annotations
                break;

            case LOGFILE::BSE_COMMON_ENVELOPES:                                                                                     // BSE_COMMON_ENVELOPES
                recordProperties = m_BSE_CEE_Rec;                                                                                   // record properties
                annotations      = m_BSE_CEE_Notes;                                                                                 // logfile annotations
                break;

            case LOGFILE::BSE_DETAILED_OUTPUT:                                                                                      // BSE_DETAILED_OUTPUT
                recordProperties = m_BSE_Detailed_Rec;                                                                              // record properties
                annotations      = m_BSE_Detailed_Notes;                                                                            // logfile annotations
                break;

            case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS:                                                                               // BSE_DOUBLE_COMPACT_OBJECTS
                recordProperties = m_BSE_DCO_Rec;                                                                                   // record properties
                annotations      = m_BSE_DCO_Notes;                                                                                 // logfile annotations
                break;

            case LOGFILE::BSE_PULSAR_EVOLUTION:                                                                                     // BSE_PULSAR_EVOLUTION
                recordProperties = m_BSE_Pulsars_Rec;                                                                               // record properties
                annotations      = m_BSE_Pulsars_Notes;                                                                             // logfile annotations
                break;

            case LOGFILE::BSE_RLOF_PARAMETERS:                                                                                      // BSE_RLOF_PARAMETERS
                recordProperties = m_BSE_RLOF_Rec;                                                                                  // record properties
                annotations      = m_BSE_RLOF_Notes;                                                                                // logfile annotations
                break;

            case LOGFILE::BSE_SUPERNOVAE:                                                                                           // BSE_SUPERNOVAE
                recordProperties = m_BSE_SNE_Rec;                                                                                   // record properties
                annotations      = m_BSE_SNE_Notes;                                                                                 // logfile annotations
                break;

            case LOGFILE::BSE_SWITCH_LOG:                                                                                           // BSE_SWITCH_LOG
                recordProperties = m_BSE_Switch_Rec;                                                                                // record properties
                annotations      = m_BSE_Switch_Notes;                                                                              // logfile annotations
                break;

            case LOGFILE::BSE_SYSTEM_PARAMETERS:                                                                                    // BSE_SYSTEM_PARAMETERS
                recordProperties = m_BSE_SysParms_Rec;                                                                              // record properties
                annotations      = m_BSE_SysParms_Notes;                                                                            // logfile annotations

                // check whether to add program option columns to BSE_SYSTEM_PARAMETERS file and add them if required
                if ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::ALWAYS) ||                                         // always add option columns?                   
                   ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::GRID) &&                                           // add for grids?
                   (!OPTIONS->GridFilename().empty() || OPTIONS->CommandLineGrid()))) {                                             // have grid file or ranges/sets?
                                                                                                                                    // yes - add program options
                    // iterate over the PROGRAM_OPTION_DETAIL map and add each entry to the recordProperties vector
                    // unfortunately no guarantee of order since it is an unordered map
                    for (auto& iter: PROGRAM_OPTION_DETAIL) {                                                                       // for each entry
                        T_ANY_PROPERTY thisProp = iter.first;                                                                       // program option
                        if (std::find(recordProperties.begin(), recordProperties.end(), thisProp) == recordProperties.end()) {      // already exists in recordProperties vector?
                            recordProperties.push_back(thisProp);                                                                   // no - add it
                        }
                    }
                }
                break;

            case LOGFILE::SSE_DETAILED_OUTPUT:                                                                                      // SSE_DETAILED_OUTPUT
                recordProperties = m_SSE_SNE_Rec;                                                                                   // record properties
                annotations      = m_SSE_SNE_Notes;                                                                                 // logfile annotations
                break;

            case LOGFILE::SSE_SUPERNOVAE:                                                                                           // SSE_SUPERNOVAE
                recordProperties = m_SSE_SNE_Rec;                                                                                   // record properties
                annotations      = m_SSE_SNE_Notes;                                                                                 // logfile annotations
                break;

            case LOGFILE::SSE_SWITCH_LOG:                                                                                           // SSE_SWITCH_LOG
                recordProperties = m_SSE_Switch_Rec;                                                                                // record properties
                annotations      = m_SSE_Switch_Notes;                                                                              // logfile annotations
                break;

            case LOGFILE::SSE_SYSTEM_PARAMETERS:                                                                                    // SSE_SYSTEM_PARAMETERS
                recordProperties = m_SSE_SysParms_Rec;                                                                              // record properties
                annotations      = m_SSE_SysParms_Notes;                                                                            // logfile annotations

                // check whether to add program option columns to SSE_SYSTEM_PARAMETERS file and add them if required
                if ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::ALWAYS) ||                                         // always add option columns?                   
                   ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::GRID) &&                                           // add for grids?
                   (!OPTIONS->GridFilename().empty() || OPTIONS->CommandLineGrid()))) {                                             // have grid file or ranges/sets?
                                                                                                                                    // yes - add program options
                    // iterate over the PROGRAM_OPTION_DETAIL map and add each entry to the recordProperties vector
                    // unfortunately no guarantee of order since it is an unordered map
                    for (auto& iter: PROGRAM_OPTION_DETAIL) {                                                                       // for each entry
                        T_ANY_PROPERTY thisProp = iter.first;                                                                       // program option
                        if (std::find(recordProperties.begin(), recordProperties.end(), thisProp) == recordProperties.end()) {      // already exists in recordProperties vector?
                            recordProperties.push_back(thisProp);                                                                   // no - add it
                        }
                    }
                }
                break;

            default:                                                                                                                // unknown logfile
                recordProperties = {};                                                                                              // no record properties
                annotations      = {};                                                                                              // no annotations
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
        fmtVector        = {};                                                                                                      // no format vector
        annotations      = {};                                                                                                      // no annotations
    }

    return std::make_tuple(recordProperties, fmtVector, annotations);
}


/*
 * Determine HDF5 datatype from COMPAS datatype
 * 
 * 
 * hid_t Log::GetHDF5DataType(const TYPENAME p_COMPASdatatype, const int p_FieldWidth)
 *
 * @param   [IN]    p_COMPASdatatype            COMPAS datatype
 * @param   [IN]    p_FieldWidth                Field width for strings (i.e. length of string).
 * @param   [IN]    p_StringQualifier           String qualifier for p_COMPASdatatype == TYPENAME::STRING
 *                                              Indicates whether the string is fixed or variable length.  Default is fixed-length.
 * @return                                      HDF5 datatype
 */
hid_t Log::GetHDF5DataType(const TYPENAME p_COMPASdatatype, const int p_FieldWidth, const STRING_QUALIFIER p_StringQualifier) {
    
    hid_t h5DataType = -1;                                                                                          // HDF5 datatype - return value

    switch (p_COMPASdatatype) {                                                                                     // which COMPAS datatype?
        case TYPENAME::SHORTINT    : h5DataType = H5T_NATIVE_SHORT; break;
        case TYPENAME::INT         : h5DataType = H5T_NATIVE_INT; break;
        case TYPENAME::LONGINT     : h5DataType = H5T_NATIVE_LONG; break;
        case TYPENAME::USHORTINT   : h5DataType = H5T_NATIVE_USHORT; break;
        case TYPENAME::UINT        : h5DataType = H5T_NATIVE_UINT; break;
        case TYPENAME::ULONGINT    : h5DataType = H5T_NATIVE_ULONG; break;
        case TYPENAME::FLOAT       : h5DataType = H5T_NATIVE_FLOAT; break;
        case TYPENAME::DOUBLE      : h5DataType = H5T_NATIVE_DOUBLE; break;
        case TYPENAME::LONGDOUBLE  : h5DataType = H5T_NATIVE_LDOUBLE; break;
        case TYPENAME::OBJECT_ID   : h5DataType = H5T_NATIVE_ULONG; break;
        case TYPENAME::ERROR       : h5DataType = H5T_NATIVE_INT; break;
        case TYPENAME::STELLAR_TYPE: h5DataType = H5T_NATIVE_INT; break;
        case TYPENAME::MT_CASE     : h5DataType = H5T_NATIVE_INT; break;
        case TYPENAME::MT_TRACKING : h5DataType = H5T_NATIVE_INT; break;
        case TYPENAME::SN_EVENT    : h5DataType = H5T_NATIVE_INT; break;
        case TYPENAME::SN_STATE    : h5DataType = H5T_NATIVE_INT; break;
        case TYPENAME::STRING: {
            hid_t h5DType = H5Tcopy(H5T_C_S1);                                                                      // HDF5 c-string datatype
            size_t size = p_StringQualifier == STRING_QUALIFIER::FIXED_LENGTH ? p_FieldWidth + 1 : H5T_VARIABLE;    // size is dependent upon string type (fixed or variable length)
            (void)H5Tset_size(h5DType, size);                                                                       // size is field width + 1 (for NULL terminator) if fixed; variable if not limited
            (void)H5Tset_cset(h5DType, H5T_CSET_ASCII);                                                             // ASCII (rather than UTF-8)
            h5DataType = h5DType;
            } break;
        case TYPENAME::BOOL: {
            if (OPTIONS->PrintBoolAsString()) {                                                                     // print bool values as strings "TRUE" or "FALSE"?
                hid_t h5DType = H5Tcopy(H5T_C_S1);                                                                  // yes - HDF5 c-string datatype
                (void)H5Tset_size(h5DType, 6);                                                                      // len("FALSE") + 1 (for NULL terminator)
                (void)H5Tset_cset(h5DType, H5T_CSET_ASCII);                                                         // ASCII (rather than UTF-8)
                h5DataType = h5DType;
            }
            else {                                                                                                  // no - print bool values as 1 or 0
                h5DataType = H5T_NATIVE_UCHAR;
            }
            } break;
        default:                                                                                                    // unknown property type
            Squawk(ERR_MSG(ERROR::UNKNOWN_DATA_TYPE));                                                              // announce error
    }

    return h5DataType;                                                                                              // HDF5 datatype
}


/*
 * Create a dataset subordinate to a group in an HDF5 file
 * 
 * 
 * hid_t Log::CreateHDF5Dataset(const string p_Filename, const hid_t p_GroupId, const string p_DatasetName, const hid_t p_H5DataType, const string p_UnitsStr, const size_t p_HDF5ChunkSize)
 *
 * @param   [IN]    p_Filename                  The filename of the HDF5 file (for error logging)
 * @param   [IN]    p_GroupId                   The group id under which the dataset should be created
 * @param   [IN]    p_DatasetName               The dataset name (this (generally) corresponds to the COMPAS column header)
 * @param   [IN]    p_H5DataType                The HDF5 datatype for the dataset
 * @param   [IN]    p_UnitsStr                  The units string to be associated with the dataset (generally corresponds to COMPAS units string)
 * @param   [IN]    p_HDF5ChunkSize             Chunk size for this dataset
 * @return                                      HDF5 dataset id (-1 indicates failure)
 */
hid_t Log::CreateHDF5Dataset(const string p_Filename, const hid_t p_GroupId, const string p_DatasetName, const hid_t p_H5DataType, const string p_UnitsStr, const size_t p_HDF5ChunkSize) {

    hid_t h5Dset = -1;                                                                                              // datset id - return value

    // create a 1-d HDF5 dataspace
    hsize_t h5Dims[1]    = {0};                                                                                     // initially 0, but...
    hsize_t h5MaxDims[1] = {H5S_UNLIMITED};                                                                         // ... unlimited
    hid_t   h5Dspace     = H5Screate_simple(1, h5Dims, h5MaxDims);                                                  // create the dataspace
    hid_t   h5CPlist     = H5Pcreate(H5P_DATASET_CREATE);                                                           // create the dataset creation property list
    (void)H5Pset_alloc_time(h5CPlist, H5D_ALLOC_TIME_INCR);                                                         // allocate space on disk incrementally
    (void)H5Pset_layout(h5CPlist, H5D_CHUNKED);                                                                     // must be chunked when using unlimited dimensions

    hsize_t h5ChunkDims[1] = {p_HDF5ChunkSize};                                                                     // chunk size - affects performance
    herr_t h5Result = H5Pset_chunk(h5CPlist, 1, h5ChunkDims);                                                       // set chunk size
    if (h5Result < 0) {                                                                                             // ok?
        Squawk("ERROR: Unable to set chunk size for HDF5 container file " + p_Filename);                            // no - announce error
    }
    else {                                                                                                          // yes - chunk size set ok
        // create HDF5 dataset
        string h5DsetName = p_DatasetName;                                                                          // dataset name 
        h5DsetName        = utils::trim(h5DsetName);                                                                // remove leading and trailing blanks
        h5Dset            = H5Dcreate(p_GroupId,                                                                    // create the dataset in group p_GroupId
                                      h5DsetName.c_str(),                                                           // dataset name
                                      p_H5DataType,                                                                 // datatype
                                      h5Dspace,                                                                     // dataspace
                                      H5P_DEFAULT,                                                                  // dataset link property list                                                                     
                                      h5CPlist,                                                                     // dataset creation property list
                                      H5P_DEFAULT);                                                                 // dataset access property list
        if (h5Dset < 0) {                                                                                           // dataset created ok?
            Squawk("ERROR: Unable to create HDF5 dataSet " + h5DsetName + " for file " + p_Filename);               // no - announce error
        }
        else {                                                                                                      // yes - dataset created ok
            // create attribute for units
            hid_t h5Dspace = H5Screate(H5S_SCALAR);                                                                 // HDF5 scalar dataspace
            hid_t h5DType  = H5Tcopy(H5T_C_S1);                                                                     // HDF5 c-string datatype

            (void)H5Tset_size(h5DType, p_UnitsStr.length() + 1);                                                    // size is strlen + 1 (for NULL terminator)
            (void)H5Tset_cset(h5DType, H5T_CSET_ASCII);                                                             // ASCII (rather than UTF-8)
            hid_t h5Attr = H5Acreate(h5Dset, "units", h5DType, h5Dspace, H5P_DEFAULT, H5P_DEFAULT);                 // create attribute for units
            if (h5Attr < 0) {                                                                                       // attribute created ok?
                Squawk("ERROR: Unable to create HDF5 attribute " + p_UnitsStr + " for dataSet " + h5DsetName);      // no - announce error
                h5Dset = -1;                                                                                        // fail
            }
            else {                                                                                                  // yes - attribute created ok
                if (H5Awrite(h5Attr, h5DType, (const void *)p_UnitsStr.c_str()) < 0) {                              // write units attributes to file - ok?
                    Squawk("ERROR: Unable to write HDF5 attribute " + p_UnitsStr + " for dataSet " + h5DsetName);   // no - announce error
                    h5Dset = -1;                                                                                    // fail
                }
            }
            (void)H5Aclose(h5Attr);                                                                                 // close attribute 
        }
    }
    (void)H5Sclose(h5CPlist);                                                                                       // close creation property list
    (void)H5Sclose(h5Dspace);                                                                                       // close scalar dataspace

    return h5Dset;                                                                                                  // return dataset id: < 0 = fail
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
 * logged (for example, when looping through a population of binary stars, p_FileSuffix would indicate
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
    LogfileDetailsT      retVal = {-1, "", {}, {}, {}, {}, {}, {}, {}, {}};                                                                     // default return value

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
                    fileDetails.annotations      = m_BSE_BE_Binaries_Notes;
                    break;

                case LOGFILE::BSE_COMMON_ENVELOPES:                                                                                             // BSE_COMMON_ENVELOPES
                    fileDetails.filename         = OPTIONS->LogfileCommonEnvelopes();
                    fileDetails.recordProperties = m_BSE_CEE_Rec;
                    fileDetails.annotations      = m_BSE_CEE_Notes;
                    break;

                case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS:                                                                                       // BSE_DOUBLE_COMPACT_OBJECTS
                    fileDetails.filename         = OPTIONS->LogfileDoubleCompactObjects();
                    fileDetails.recordProperties = m_BSE_DCO_Rec;
                    fileDetails.annotations      = m_BSE_DCO_Notes;
                    break;
               
                case LOGFILE::BSE_PULSAR_EVOLUTION:                                                                                             // BSE_PULSAR_EVOLUTION
                    fileDetails.filename         = OPTIONS->LogfilePulsarEvolution();
                    fileDetails.recordProperties = m_BSE_Pulsars_Rec;
                    fileDetails.annotations      = m_BSE_Pulsars_Notes;
                    break;

                case LOGFILE::BSE_RLOF_PARAMETERS:                                                                                              // BSE_RLOF_PARAMETERS
                    fileDetails.filename         = OPTIONS->LogfileRLOFParameters();
                    fileDetails.recordProperties = m_BSE_RLOF_Rec;
                    fileDetails.annotations      = m_BSE_RLOF_Notes;
                    break;

                case LOGFILE::BSE_SUPERNOVAE:                                                                                                   // BSE_SUPERNOVAE
                    fileDetails.filename         = OPTIONS->LogfileSupernovae();
                    fileDetails.recordProperties = m_BSE_SNE_Rec;
                    fileDetails.annotations      = m_BSE_SNE_Notes;
                    break;

                case LOGFILE::BSE_SWITCH_LOG:                                                                                                   // BSE_SWITCH_LOG
                    fileDetails.filename         = OPTIONS->LogfileSwitchLog();
                    fileDetails.recordProperties = m_BSE_Switch_Rec;
                    fileDetails.annotations      = m_BSE_Switch_Notes;
                    break;

                case LOGFILE::BSE_SYSTEM_PARAMETERS:                                                                                            // BSE_SYSTEM_PARAMETERS
                    fileDetails.filename         = OPTIONS->LogfileSystemParameters();
                    fileDetails.recordProperties = m_BSE_SysParms_Rec;
                    fileDetails.annotations      = m_BSE_SysParms_Notes;

                    // check whether to add program option columns to BSE_SYSTEM_PARAMETERS file and add them if required
                    if ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::ALWAYS) ||                                                 // always add option columns?                   
                       ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::GRID) &&                                                   // add for grids?
                       (!OPTIONS->GridFilename().empty() || OPTIONS->CommandLineGrid()))) {                                                     // have grid file or ranges/sets?
                                                                                                                                                // yes - add program options
                        // iterate over the PROGRAM_OPTION_DETAIL map and add each entry to the recordProperties vector
                        // unfortunately no guarantee of order since it is an unordered map
                        for (auto& iter: PROGRAM_OPTION_DETAIL) {                                                                               // for each entry
                            T_ANY_PROPERTY thisProp = iter.first;                                                                               // program option
                            if (std::find(fileDetails.recordProperties.begin(), fileDetails.recordProperties.end(), thisProp) == fileDetails.recordProperties.end()) { // already exists in recordProperties vector?
                                fileDetails.recordProperties.push_back(thisProp);                                                               // no - add it
                            }
                        }
                    }
                    break;

                case LOGFILE::SSE_SUPERNOVAE:                                                                                                   // SSE_SUPERNOVAE
                    fileDetails.filename         = OPTIONS->LogfileSupernovae();
                    fileDetails.recordProperties = m_SSE_SNE_Rec;
                    fileDetails.annotations      = m_SSE_SNE_Notes;
                    break;

                case LOGFILE::SSE_SWITCH_LOG:                                                                                                   // SSE_SWITCH_LOG
                    fileDetails.filename         = OPTIONS->LogfileSwitchLog();
                    fileDetails.recordProperties = m_SSE_Switch_Rec;
                    fileDetails.annotations      = m_SSE_Switch_Notes;
                    break;

                case LOGFILE::SSE_SYSTEM_PARAMETERS:                                                                                            // SSE_SYSTEM_PARAMETERS
                    fileDetails.filename         = OPTIONS->LogfileSystemParameters();
                    fileDetails.recordProperties = m_SSE_SysParms_Rec;
                    fileDetails.annotations      = m_SSE_SysParms_Notes;

                    // check whether to add program option columns to SSE_SYSTEM_PARAMETERS file and add them if required
                    if ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::ALWAYS) ||                                                 // always add option columns?                   
                       ((OPTIONS->AddOptionsToSysParms() == ADD_OPTIONS_TO_SYSPARMS::GRID) &&                                                   // add for grids?
                       (!OPTIONS->GridFilename().empty() || OPTIONS->CommandLineGrid()))) {                                                     // have grid file or ranges/sets?
                                                                                                                                                // yes - add program options
                        // iterate over the PROGRAM_OPTION_DETAIL map and add each entry to the recordProperties vector
                        // unfortunately no guarantee of order since it is an unordered map
                        for (auto& iter: PROGRAM_OPTION_DETAIL) {                                                                               // for each entry
                            T_ANY_PROPERTY thisProp = iter.first;                                                                               // program option
                            if (std::find(fileDetails.recordProperties.begin(), fileDetails.recordProperties.end(), thisProp) == fileDetails.recordProperties.end()) { // already exists in recordProperties vector?
                                fileDetails.recordProperties.push_back(thisProp);                                                               // no - add it
                            }
                        }
                    }
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
                                fileDetails.annotations      = m_SSE_Detailed_Notes;
                                break;

                            case LOGFILE::BSE_DETAILED_OUTPUT:                                                                                  // BSE_DETAILED_OUTPUT
                                fileDetails.filename         = DETAILED_OUTPUT_DIRECTORY_NAME + "/" + OPTIONS->LogfileDetailedOutput();         // logfile filename with directory
                                fileDetails.recordProperties = m_BSE_Detailed_Rec;                                                              // record properties
                                fileDetails.annotations      = m_BSE_Detailed_Notes;
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
                //    - we should have an HDF5 container (opened in Log::Start())
                //    - all logfiles except detailed output files are included in a single HDF5 file
                //    - detailed output files are created individually as HDF5 files inside their containing directory
                if (m_LogfileType == LOGFILETYPE::HDF5) {                                                                                       // logging to HDF5 files?
                                                                                                                                                // yes
                    string fileExt = "." + LOGFILETYPEFileExt.at(OPTIONS->LogfileType());                                                       // file extension for HDF5 files

                    if (p_Logfile == LOGFILE::SSE_DETAILED_OUTPUT || p_Logfile == LOGFILE::BSE_DETAILED_OUTPUT) {                               // yes - detailed output file (SSE or BSE)?
                        if (m_HDF5DetailedId < 0) {                                                                                             // have HDF5 detailed file?
                                                                                                                                                // no - create it
                            string h5Filename = m_LogBasePath + "/" + m_LogContainerName + "/" + fileDetails.filename + fileExt;                // full filename with path, container, and extension ("/" works on Uni*x and Windows)

                            // check if file already exists - if it does, add a version number before creating new file
                            // no append for detailed output files, so no need to open existing files for appending
            
                            int version = 0;                                                                                                    // logfile version number if required - start at 1
                            while (utils::FileExists(h5Filename)) {                                                                             // file already exists?
                                h5Filename = m_LogBasePath + "/" + m_LogContainerName + "/" + fileDetails.filename + "_" + std::to_string(++version) + fileExt; // yes - add a version number and generate new filename
                            }

                            m_HDF5DetailedId = H5Fcreate(h5Filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);                           // create HDF5 detailed file
                            if (m_HDF5DetailedId < 0) {                                                                                         // created ok?                        
                                Squawk("ERROR: Unable to create HDF5 detailed file with file name " + h5Filename);                              // no - announce error
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

                        for (auto iter = fileDetails.recordProperties.begin(); iter != fileDetails.recordProperties.end();) {                   // for each property to be included in the log record

                            T_ANY_PROPERTY property = *iter;                                                                                    // this record proerty

                            string headerStr = "";
                            string unitsStr  = "";
                            string typeStr   = "";
                            string fmtStr    = "";

                            bool push = true;
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

                                    // check if we're adding PROGRAM_OPTION::NOTES - handled differently

                                    if (programOption == PROGRAM_OPTION::NOTES) {                                                               // PROGRAM_OPTION::NOTES?
                                                                                                                                                // yes
                                        // processing PROGRAM_OPTION::NOTES
                                        //
                                        // program option NOTES is special - the property is actually a vector of strings (OPTIONS->Notes()),
                                        // each of which is printed to a separate column, so we need to add as many entries to fileDetails
                                        // as there are notes columns to add to this file.  The notes required to be written to this file
                                        // are record in fileDetails.annotations (already set).

                                        if (fileDetails.annotations.size() == 0) {                                                              // annotations present?
                                            iter = fileDetails.recordProperties.erase(iter);                                                    // no - remove NOTES from record properties
                                        }
                                        else {                                                                                                  // have annotations
                                            for (size_t idx = 0; idx < fileDetails.annotations.size(); idx ++) {                                // for each user-specified annotation
                                                if (fileDetails.annotations[idx]) {                                                             // include it?
                                                                                                                                                // yes
                                                    details = ProgramOptionDetails(programOption, idx);                                         // property details
                                                    std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details);               // format the headers

                                                    fileDetails.propertyTypes.push_back(std::get<0>(details));                                  // append property typename
                                                    fileDetails.stringTypes.push_back(STRING_QUALIFIER::VARIABLE_LENGTH);                       // append string type - annotations are variable length
                                                    fileDetails.hdrStrings.push_back(headerStr);                                                // append header string for field
                                                    fileDetails.unitsStrings.push_back(unitsStr);                                               // append units string for field
                                                    fileDetails.typeStrings.push_back(typeStr);                                                 // append type string for field
                                                    fileDetails.fmtStrings.push_back(fmtStr);                                                   // append format string for field
                                                }
                                            }
                                            iter++;                                                                                             // next property
                                        }
                                        push = false;                                                                                           // don't need to add property details later - already done
                                    }
                                    else {                                                                                                      // not PROGRAM_OPTION::NOTES
                                        details = ProgramOptionDetails(programOption);                                                          // property details
                                        std::tie(headerStr, unitsStr, typeStr, fmtStr) = FormatFieldHeaders(details);                           // format the headers
                                    }
                                    } break;

                                default:                                                                                                        // unknown property type
                                    push = false;                                                                                               // that's not ok - don't add property details
                                    Squawk(ERR_MSG(ERROR::UNKNOWN_PROPERTY_TYPE));                                                              // show warning
                            }

                            if (push) {                                                                                                         // record details?
                                                                                                                                                // yes
                                fileDetails.propertyTypes.push_back(std::get<0>(details));                                                      // append property typename
                                fileDetails.stringTypes.push_back(STRING_QUALIFIER::FIXED_LENGTH);                                              // append string type - default is fixed length
                                fileDetails.hdrStrings.push_back(headerStr);                                                                    // append header string for field
                                fileDetails.unitsStrings.push_back(unitsStr);                                                                   // append units string for field
                                fileDetails.typeStrings.push_back(typeStr);                                                                     // append type string for field
                                fileDetails.fmtStrings.push_back(fmtStr);                                                                       // append format string for field

                                iter++;                                                                                                         // next property
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
                        //
                        // These are hard-coded here rather than in the *_PROPERTY_DETAIL maps in
                        // constants.h so that they will always be present in the switch file -
                        // this way users can't add or remove them at runtime via the logfile-definitions
                        // option.

                        if (p_Logfile == LOGFILE::BSE_SWITCH_LOG) {                                                                             // BSE Switch Log
                            fileDetails.propertyTypes.push_back(TYPENAME::INT);                                                                 // append property typename
                            fileDetails.hdrStrings.push_back("STAR_SWITCHING");                                                                 // append header string for field
                            fileDetails.unitsStrings.push_back("-");                                                                            // append units string for field
                            fileDetails.typeStrings.push_back("INT");                                                                           // append type string for field
                            fileDetails.fmtStrings.push_back("4.1");                                                                            // append format string for field (size accomodates header string)
                        }

                        if (p_Logfile == LOGFILE::BSE_SWITCH_LOG || p_Logfile == LOGFILE::SSE_SWITCH_LOG) {                                     // BSE Switch Log or SSE Switch Log
                            fileDetails.propertyTypes.push_back(TYPENAME::STELLAR_TYPE);                                                                 // append property typename
                            fileDetails.propertyTypes.push_back(TYPENAME::STELLAR_TYPE);                                                                 // append property typename

                            fileDetails.hdrStrings.push_back("SWITCHING_FROM");                                                                 // append header string for field
                            fileDetails.hdrStrings.push_back("SWITCHING_TO");                                                                   // append header string for field

                            fileDetails.unitsStrings.push_back("-");                                                                            // append units string for field
                            fileDetails.unitsStrings.push_back("-");                                                                            // append units string for field

                            fileDetails.typeStrings.push_back("INT");                                                                           // append type string for field
                            fileDetails.typeStrings.push_back("INT");                                                                           // append type string for field

                            fileDetails.fmtStrings.push_back("4.1");                                                                            // append fromat string for field (size accomodates header string)
                            fileDetails.fmtStrings.push_back("4.1");                                                                            // append format string for field (size accomodates header string)
                        }
                    }

                    // record new open file details
                    m_OpenStandardLogFileIds.insert({p_Logfile, fileDetails});                                                                  // record the new file details and format strings

                    // initialise files:
                    //    - create datasets for HDF5 files
                    //    - write header/units/types strings for CSV/TSV/TXT files
                    if (OPTIONS->LogfileType() == LOGFILETYPE::HDF5) {                                                                          // logging to HDF5 files?
                        for (size_t idx = 0; idx < fileDetails.hdrStrings.size(); idx++) {                                                      // for each property
                            
                            size_t chunkSize = OPTIONS->nObjectsToEvolve() < HDF5_MINIMUM_CHUNK_SIZE || 
                                               p_Logfile == LOGFILE::SSE_DETAILED_OUTPUT             || 
                                               p_Logfile == LOGFILE::BSE_DETAILED_OUTPUT ? HDF5_MINIMUM_CHUNK_SIZE : OPTIONS->HDF5ChunkSize();  // chunk size

                            size_t IOBufSize = OPTIONS->HDF5BufferSize() * chunkSize;                                                           // IO buffer size
                
                            m_Logfiles[fileDetails.id].h5File.chunkSize = chunkSize;                                                            // record chunk size for file
                            m_Logfiles[fileDetails.id].h5File.IOBufSize = IOBufSize;                                                            // record IO buf size for file

                            m_Logfiles[fileDetails.id].h5File.dataSets.push_back({-1, -1, TYPENAME::NONE, {}});                                 // create new dataset

                            // set datatypes
                            int fieldWidth = fileDetails.fmtStrings[idx].empty() ? 0 : (int)std::stod(fileDetails.fmtStrings[idx]);
                            hid_t h5DataType = GetHDF5DataType(fileDetails.propertyTypes[idx], fieldWidth, fileDetails.stringTypes[idx]);       // get HDF5 data type from COMPAS data type
                            if (h5DataType < 0) {                                                                                               // ok?
                                ok = false;                                                                                                     // no - fail
                            }
                            else {                                                                                                              // yes - ok
                                m_Logfiles[fileDetails.id].h5File.dataSets[idx].dataType   = fileDetails.propertyTypes[idx];                    // record COMPAS data type
                                m_Logfiles[fileDetails.id].h5File.dataSets[idx].stringType = fileDetails.stringTypes[idx];                      // record COMPAS string type
                                m_Logfiles[fileDetails.id].h5File.dataSets[idx].h5DataType = h5DataType;                                        // record HDF5 data type

                                // create HDF5 dataset
                                hid_t h5Dset = CreateHDF5Dataset(fileDetails.filename, 
                                                                 m_Logfiles[fileDetails.id].h5File.groupId, 
                                                                 utils::trim(fileDetails.hdrStrings[idx]), 
                                                                 h5DataType, 
                                                                 utils::trim(fileDetails.unitsStrings[idx]),
                                                                 chunkSize);
                                if (h5Dset < 0) {                                                                                               // created ok?
                                    ok = false;                                                                                                 // no - fail
                                }
                                
                                m_Logfiles[fileDetails.id].h5File.dataSets[idx].dataSetId = h5Dset;                                             // record dataset id
                            }
                        }
                    }
                    else {                                                                                                                      // no - FS file
                        string fullHdrsStr  = "";                                                                                               // initialise full headers string
                        string fullUnitsStr = "";                                                                                               // initialise full units string
                        string fullTypesStr = "";                                                                                               // initialise full types string
                        for (size_t idx = 0; idx < fileDetails.hdrStrings.size(); idx++) {                                                      // for each property
                            fullHdrsStr  += fileDetails.hdrStrings[idx] + delimiter;                                                            // append field header string to full header string
                            fullUnitsStr += fileDetails.unitsStrings[idx] + delimiter;                                                          // append field units string to full units string
                            fullTypesStr += fileDetails.typeStrings[idx] + delimiter;                                                           // append field type string to full type string
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
    else {                                                                                                                                      // logfile already exists
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
 * relevant if both p_AddProps and p_SubtractProps are non-empty.  The order is:
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
 *                               const ANY_PROPERTY_VECTOR p_SubtractProps,
 *                               const std::vector<bool>   p_AddNotes,
 *                               const std::vector<bool>   p_SubtractNotes)
 *
 * 
 * @param   [IN]    p_Logfile                   the logfile for which the record specifier should be updated
 * @param   [IN]    p_UseDefaultProps           indicates whether the default properties of the given logfile should be
 *                                              be used as the base set of properties.  If p_UseDefaultProps is true,
 *                                              the base set of properties is set to the current set for the file indicated
 *                                              by p_LogFile (initially the default set from constants.h).
 *                                              If p_UseDefaultProps is false, the base set of of properties is set
 *                                              empty
 * @param   [IN]    p_AddProps                  vector containing the properties to be added to the given logfile properties
 * @param   [IN]    p_SubtractProps             vector containing the properties to be subtracted from the given logfile properties
 * @param   [IN]    p_AddNotes                  vector of booleans indicating which notes are to be added, if PROGRAM_OPTION::NOTES is in p_AddProps
 * @param   [IN]    p_SubtractNotes             vector of booleans indicating which notes are to be subtracted, if PROGRAM_OPTION::NOTES is in p_SubtractProps
 */
void Log::UpdateLogfileRecordSpecs(const LOGFILE             p_Logfile,
                                   bool                      p_UseDefaultProps,
                                   const ANY_PROPERTY_VECTOR p_AddProps,
                                   const ANY_PROPERTY_VECTOR p_SubtractProps,
                                   const std::vector<bool>   p_AddNotes,
                                   const std::vector<bool>   p_SubtractNotes) {

    ANY_PROPERTY_VECTOR baseProps = {};                                                                                 // base props for the given logfile
    std::vector<bool>   baseNotes;                                                                                      // base annotations for the given logfile

    ANY_PROPERTY_VECTOR newProps = {};                                                                                  // new props for the given logfile
    std::vector<bool>   newNotes;                                                                                       // new annotations for the given logfile

    // initialise baseProps, baseNotes, and newNotes
    // if p_UseDefaultProps is TRUE, baseProps is initialised to the default record specifier for the file,
    // otherwise it is initialised to an empty vector.
    // baseNotes and newNotes are both initialised to the current annottions vector for the file - these
    // are only used if the record specifier includes PROGRAM_OPTION::NOTES, and initialising to the current
    // annotations vector for the file makes updating the vector correctly easier.

    switch (p_Logfile) {
        case LOGFILE::BSE_BE_BINARIES:
            if (p_UseDefaultProps) baseProps = BSE_BE_BINARIES_REC;
            baseNotes = m_BSE_BE_Binaries_Notes;
            break;
        case LOGFILE::BSE_COMMON_ENVELOPES:
            if (p_UseDefaultProps) baseProps = BSE_COMMON_ENVELOPES_REC;
            baseNotes = m_BSE_CEE_Notes;
            break;
        case LOGFILE::BSE_DETAILED_OUTPUT:
            if (p_UseDefaultProps) baseProps = BSE_DETAILED_OUTPUT_REC;
            baseNotes = m_BSE_Detailed_Notes;
            break;
        case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS:
            if (p_UseDefaultProps) baseProps = BSE_DOUBLE_COMPACT_OBJECTS_REC;
            baseNotes = m_BSE_DCO_Notes;
            break;
        case LOGFILE::BSE_PULSAR_EVOLUTION:
            if (p_UseDefaultProps) baseProps = BSE_PULSAR_EVOLUTION_REC;
            baseNotes = m_BSE_Pulsars_Notes;
            break;
        case LOGFILE::BSE_RLOF_PARAMETERS:
            if (p_UseDefaultProps) baseProps = BSE_RLOF_PARAMETERS_REC;
            baseNotes = m_BSE_RLOF_Notes;
            break;
        case LOGFILE::BSE_SUPERNOVAE:
            if (p_UseDefaultProps) baseProps = BSE_SUPERNOVAE_REC;
            baseNotes = m_BSE_SNE_Notes;
            break;
        case LOGFILE::BSE_SWITCH_LOG:
            if (p_UseDefaultProps) baseProps = BSE_SWITCH_LOG_REC;
            baseNotes = m_BSE_Switch_Notes;
            break;
        case LOGFILE::BSE_SYSTEM_PARAMETERS:
            if (p_UseDefaultProps) baseProps = BSE_SYSTEM_PARAMETERS_REC;
            baseNotes = m_BSE_SysParms_Notes;
            break;
        case LOGFILE::SSE_DETAILED_OUTPUT:
            if (p_UseDefaultProps) baseProps = SSE_DETAILED_OUTPUT_REC;
            baseNotes = m_SSE_Detailed_Notes;
            break;
        case LOGFILE::SSE_SUPERNOVAE:
            if (p_UseDefaultProps) baseProps = SSE_SUPERNOVAE_REC;
            baseNotes = m_SSE_Detailed_Notes;
            break;
        case LOGFILE::SSE_SWITCH_LOG:
            if (p_UseDefaultProps) baseProps = SSE_SWITCH_LOG_REC;
            baseNotes = m_SSE_Switch_Notes;
            break;
        case LOGFILE::SSE_SYSTEM_PARAMETERS:
            if (p_UseDefaultProps) baseProps = SSE_SYSTEM_PARAMETERS_REC;
            baseNotes = m_SSE_SysParms_Notes;
            break;
        default: break;                                                                                                 // avoids compiler warning
    }
    newNotes = baseNotes;

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
                        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION    : {                                                                                                                      // PROGRAM_OPTION
                        
                            PROGRAM_OPTION thisBaseProperty     = boost::get<PROGRAM_OPTION>(baseProperty);             // base property
                            PROGRAM_OPTION thisSubtractProperty = boost::get<PROGRAM_OPTION>(subtractProperty);         // property to be subtracted
                            
                            isSubtract = thisBaseProperty == thisSubtractProperty;                                      // should subtract (nominally)?

                            if (isSubtract && thisBaseProperty == PROGRAM_OPTION::NOTES) {                              // subtracting PROGRAM_OPTION::NOTES?
                                                                                                                        // yes - update the annotations vector for this logfile
                                size_t notesCount = 0;                                                                  // count of notes user wants
                                for (size_t idx = 0; idx < std::min(newNotes.size(), p_SubtractNotes.size()); idx++) {  // min of sizes for safety - should be same
                                    newNotes[idx] = p_SubtractNotes[idx] ? false : baseNotes[idx];                      // subtract notes index requested
                                    if (newNotes[idx]) notesCount++;
                                }
                                isSubtract = notesCount == 0;                                                           //  only remove PROGRAM_OPTION::NOTES if the user wants none of the notes
                            }}
                            break;
                        default: break;                                                                                 // avoids compiler warning
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
            ANY_PROPERTY_TYPE addPropertyType = boost::apply_visitor(VariantPropertyType(), addProperty);               // add property type

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
                        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION    : {                                                                                                               // PROGRAM_OPTION                        
                        
                            PROGRAM_OPTION thisNewProperty = boost::get<PROGRAM_OPTION>(newProperty);                   // new property
                            PROGRAM_OPTION thisAddProperty  = boost::get<PROGRAM_OPTION>(addProperty);                  // property to be added
                            
                            isAlready = thisNewProperty == thisAddProperty;                                             // should add (nominally)?

                            if (isAlready && thisNewProperty == PROGRAM_OPTION::NOTES) {                                // adding PROGRAM_OPTION::NOTES?
                                                                                                                        // yes - update the annotations vector for this logfile
                                size_t notesCount = 0;                                                                  // count of notes user wants
                                for (size_t idx = 0; idx < std::min(newNotes.size(), p_AddNotes.size()); idx++) {       // min of sizes for safety - should be same
                                    newNotes[idx] = p_AddNotes[idx] ? true : newNotes[idx];                             // add notes index requested
                                    if (newNotes[idx]) notesCount++;
                                }
                                isAlready = notesCount > 0;                                                             //  only add PROGRAM_OPTION::NOTES if the user wants at least one of the notes
                            }}
                            break;
                        default: break;                                                                                 // avoids compiler warning
                    }
                    if (isAlready) break;
                }
            }
            if (!isAlready) newProps.push_back(addProperty);                                                            // add property to newProps if necessary
        }
    }

    // replace existing props and annotations vector for given logfile
    switch (p_Logfile) {
        case LOGFILE::BSE_BE_BINARIES           : m_BSE_BE_Binaries_Rec = newProps; m_BSE_BE_Binaries_Notes = newNotes; break;
        case LOGFILE::BSE_COMMON_ENVELOPES      : m_BSE_CEE_Rec         = newProps; m_BSE_CEE_Notes         = newNotes; break;
        case LOGFILE::BSE_DETAILED_OUTPUT       : m_BSE_Detailed_Rec    = newProps; m_BSE_Detailed_Notes    = newNotes; break;
        case LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS: m_BSE_DCO_Rec         = newProps; m_BSE_DCO_Notes         = newNotes; break;
        case LOGFILE::BSE_PULSAR_EVOLUTION      : m_BSE_Pulsars_Rec     = newProps; m_BSE_Pulsars_Notes     = newNotes; break;
        case LOGFILE::BSE_RLOF_PARAMETERS       : m_BSE_RLOF_Rec        = newProps; m_BSE_RLOF_Notes        = newNotes; break;
        case LOGFILE::BSE_SUPERNOVAE            : m_BSE_SNE_Rec         = newProps; m_BSE_SNE_Notes         = newNotes; break;
        case LOGFILE::BSE_SWITCH_LOG            : m_BSE_Switch_Rec      = newProps; m_BSE_Switch_Notes      = newNotes; break;
        case LOGFILE::BSE_SYSTEM_PARAMETERS     : m_BSE_SysParms_Rec    = newProps; m_BSE_SysParms_Notes    = newNotes; break;
        case LOGFILE::SSE_DETAILED_OUTPUT       : m_SSE_Detailed_Rec    = newProps; m_SSE_Detailed_Notes    = newNotes; break;
        case LOGFILE::SSE_SUPERNOVAE            : m_SSE_SNE_Rec         = newProps; m_SSE_SNE_Notes         = newNotes; break;
        case LOGFILE::SSE_SWITCH_LOG            : m_SSE_Switch_Rec      = newProps; m_SSE_Switch_Notes      = newNotes; break;
        case LOGFILE::SSE_SYSTEM_PARAMETERS     : m_SSE_SysParms_Rec    = newProps; m_SSE_SysParms_Notes    = newNotes; break;
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
 *       n    means integer number
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
 * <prop_spec>  ::= <prop_type> "::" <prop_name> [ <prop_index> ] <prop_delim>
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
 * <prop_index> ::= "[" n "]"
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
 * @return                                      boolean indicating if logfile record specifications were updated successfully:
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
    std::vector<bool>   addNotes = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);                                       // annotations (notes) user wants added to the base properties
    std::vector<bool>   subtractNotes = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);                                  // annotations (notes) user wants subtracted from the base properties

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
                        addNotes      = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);                                  // start with no annotations to be added       
                        subtractNotes = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);                                  // start with no annotations to be subtracted

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
                        UpdateLogfileRecordSpecs(currentLogfile, useDefaultProps, addProps, subtractProps, addNotes, subtractNotes); // update current logfile record specifications

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
                        UpdateLogfileRecordSpecs(currentLogfile, useDefaultProps, addProps, subtractProps, addNotes, subtractNotes); // update current logfile record specifications

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
                        // check for indexed property name - only allowed for PROGRAM_OPTION::NOTES
                        // check that property type and name are a valid, known property type and name
                        // check that subscript for PROGRAM_OPTION::NOTES, if present, is valid

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

                                        // check for indexed PROGRAM_OPTION::NOTES

                                        int notesIdx = -1;                                                                      // index value for PROGRAM_OPTION::NOTES
                                        if (propNameStr.length() > 5 && utils::Equals(propNameStr.substr(0, 5), "notes")) {     // PROGRAM_OPTION::NOTES?
                                            if (propNameStr[5] == '[' && propNameStr[propNameStr.length() - 1] != ']') {        // wanted index, but failed?
                                                error = ERROR::EXPECTED_POSITIVE_INTEGER;                                       // set error - expected an integer index > 0
                                                errorPos += 22;                                                                 // caret position for error
                                            }
                                            else if (propNameStr[5] == '[' && propNameStr[propNameStr.length() - 1] == ']') {   // possibly - indexed?
                                                size_t idxLen = propNameStr.length() - 7;                                       // possibly...
                                                if (idxLen > 0) {                                                               // length of index value > 0?    
                                                    // indexed - check for valid index
                                                    try {
                                                        size_t lastChar;                                                        // for conversion
                                                        notesIdx = std::stoi(propNameStr.substr(6, idxLen), &lastChar);         // try conversion
                                                        if (lastChar != idxLen) {                                               // valid INT only if propNameStr completely consumed
                                                            error = ERROR::EXPECTED_POSITIVE_INTEGER;                           // not valid - set error - expected an integer index > 0
                                                            errorPos += 22;                                                     // caret position for error
                                                            notesIdx = -1;                                                      // don't want to use invalid index value
                                                        }
                                                        else {                                                                  // valid integer - valid index?
                                                            // the index for PROGRAM_OPTION::NOTES is 1-based (i.e. the first
                                                            // note is index 1 (not 0))
                                                            // to be a valid index for PROGRAM_OPTION::NOTES, the index value
                                                            // must be from 1..number of notes headers (inclusive)

                                                            if (notesIdx < 1 or notesIdx > (int)OPTIONS->NotesHdrs().size()) {  // index in valid range?
                                                                error = ERROR::OUT_OF_BOUNDS;                                   // no - set error
                                                                errorPos += 22;                                                 // caret position for error
                                                            }
                                                            else {                                                              // yes - valid index
                                                                propNameStr = propNameStr.substr(0, 5);                         // remove index from property name
                                                            }
                                                        }
                                                    }
                                                    catch (const std::out_of_range& e) {                                        // conversion failed
                                                        error = ERROR::EXPECTED_POSITIVE_INTEGER;                               // set error - expected an integer index > 0
                                                        errorPos += 22;                                                         // caret position for error
                                                    }
                                                    catch (const std::invalid_argument& e) {                                    // conversion failed
                                                        error = ERROR::EXPECTED_POSITIVE_INTEGER;                               // set error - expected an integer index > 0
                                                        errorPos += 22;                                                         // caret position for error
                                                    }
                                                }                              
                                            }
                                        }

                                        if (error == ERROR::NONE) {                                                             // ok so far?
                                                                                                                                // yes
                                            bool found;
                                            PROGRAM_OPTION property;                                                            // lookup property name
                                            std::tie(found, property) = utils::GetMapKey(propNameStr, PROGRAM_OPTION_LABEL, PROGRAM_OPTION::RANDOM_SEED);
                                            if (!found) {                                                                       // property name found?
                                                error = ERROR::UNKNOWN_PROGRAM_OPTION;                                          // no - set error
                                            }
                                            else {                                                                              // found known property name
                                                if (addAssign) {                                                                // add property?
                                                    addProps.push_back(property);                                               // yes - add it to 'add' list

                                                    if (property == PROGRAM_OPTION::NOTES) {                                    // PROGRAM_OPTION::NOTES?
                                                        if (notesIdx < 1) {                                                     // yes - indexed?
                                                            for (size_t idx = 0; idx < addNotes.size(); idx++) addNotes[idx] = true; // no - add all notes
                                                        }
                                                        else {                                                                  // yes, indexed
                                                            addNotes[notesIdx - 1] = true;                                      // add specified note only (0 - based)
                                                        }
                                                    }
                                                }
                                                else {                                                                          // no - subtract property
                                                    subtractProps.push_back(property);                                          // add it to 'subtract' list

                                                    if (property == PROGRAM_OPTION::NOTES) {                                    // PROGRAM_OPTION::NOTES?
                                                        if (notesIdx < 1) {                                                     // yes - indexed?
                                                            for (size_t idx = 0; idx < subtractNotes.size(); idx++) subtractNotes[idx] = true; // no - subtract all notes
                                                        }
                                                        else {                                                                  // yes, indexed
                                                            subtractNotes[notesIdx - 1] = true;                                 // subtract specified note only (0 - based)
                                                        }
                                                    }
                                                }
                                            }
                                            expecting = TOKEN_TYPE::COMMA;                                                      // now expecting a comma or close brace
                                        }
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
