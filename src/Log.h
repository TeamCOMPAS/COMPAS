#ifndef __Log_h__
#define __Log_h__

#define LOGGING Log::Instance()

#include <fstream>
#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/variant.hpp>

#include "hdf5.h"

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Options.h"
#include "LogMacros.h"

using std::string;

/*
 * Log Singleton
 *
 * Provides global logging functionality
 *
 * Singletons are sometimes frowned-upon - but Logging functionality is one area
 * for which they are considered acceptable.  Doing it this way means the logging
 * class instance doesn't need to be passed around to all and sundry.
 *
 * JR, May 2019
 */


/*
 * Note: I wrote this class when I refactored the SSE code, without reference to the BSE code.
 * While refactoring the BSE code (later) I realised that the logging I had implemented for SSE
 * wasn't sufficient for the logging requirements of the BSE code.  To provide for the logging
 * needs of the BSE code I added new functionality almost as a wrapper around the original, SSE,
 * logging functionality.  Some of the original SSE logging functionality has almost been
 * rendered redundant by the new BSE code, but I have left it here (almost) in its entirety
 * because it may still be useful.  We'll see.
 *
 * When I wrote the SSE logging functionality I provided debugging functionality with it, as
 * well as a set of macros to make debugging and the issuing of warning messages easier.  I
 * also wrote a set of logging macros to make logging easier.  The debug macros are still
 * useful, and I would encourage their use.  As well as writing the BSE logging code, I also
 * wrote the Errors class, which is my attempt at making error handling easier.  Some of the
 * functionality in the Errors class supersedes the DBG_WARN* macros provided here, but the
 * DBG_WARN* macros are still useful in some circumstances (and in fact are still used in
 * various places in the code).  The LOG* macros are somewhat less useful, but are left here
 * in case the original SSE logging functionality (that which underlies the BSE logging
 * functionality) is used in the future (as mentioned above, it could still be useful in
 * some circumstances).
 *
 * JR, September 2019
 * 
 * 
 * 
 * HDF5 File Support
 * =================
 * 
 * Further note: support for HDF5 was added to the exiting log functionality, which was written
 * assuming separate files for each of the log file (system parameters, DCOs, SNe, etc.).  That
 * existing functionality has been maintained, and HDF5 added as an option available to the user.
 * HDF5 support in the code has been shoe-horned into the existing functionality, so although
 * there is only a single HDF5 output file (except for detailed output files), each of the groups
 * within the HDF5 file is treated in the logging code as a separate file - because that's how it 
 * was originally written, and we just use the existing infrastructure and fit the HDF5 support 
 * into that framework.  That makes the HDF5 logging code a bit awkward at times, but it works...
 *
 *
 * A brief description of HDF5 files and chunking, in the COMPAS context:
 * (With the caveat that all I know about HDF5 I learned so that I could write this code - it is
 * very likely that there are better ways to do things, so if anyone knows a better way, please
 * either change the code or tell me how to improve it and I'll change it.  Most of what follows
 * is just a brain dump of my reading/research over the past week or so, and it could very well
 * be based on my misunderstanding of what I have read - so uf anyone notices something I've
 * misunderstood please let me know so I can improve the code.)
 * 
 *
 * Data in HDF5 files are arranged in groups and datasets.
 *
 *     A COMPAS output file (e.g. BSE_System_Parameters, BSE_RLOF, etc.) maps to an HDF5 group,
 *     where the group name is the name of the COMPAS output file.
 *
 *     A column in a COMPAS output file (e.g. SEED, Mass(1), Radius(2), etc.) maps to an HDF5 dataset,
 *     where the dataset name is the column heading string.
 *
 *     COMPAS column datatype strings are encoded in the dataset meta-details (dataset.dtype).
 *     COMPAS column units strings are attached to HDF5 datasets as attributes.
 *
 * Each dataset in an HDF5 files is broken into "chunks", where a chunk is defined as a number of dataset
 * entries ("chunks" and "chunking" are HDF5 terms).  In COMPAS, all datasets are 1-d arrays (columns), so
 * a chunk is defined as a number of values in the 1-d array (or column).  Chunking can be enabled or not, 
 * but if chunking is not* enabled a dataset cannot be resized - so if chunking is not enabled the size of 
 * the dataset must be known at the time of creation, and the entire datset created in one go.  That doesn't 
 * work for COMPAS - even though we know the number of systems being evolved, we don't know the number of 
 * entries we'll have in each of the output log files (and therefore the HDF5 datasests if we're logging to 
 * HDF5 files).  So, we enable chunking.
 * 
 * Chunking can improve, or degrade, performance depending upon how it is implemented - mostly related to
 * the chunk size chosen.
 * 
 * Datasets are stored inside an HDF5 file as a number of chunks - the chunks are not guaranteed (not even
 * likely) to be contiguous in the file or on the storage media (HDD, SSD etc.).  Chunks are mapped/indexed
 * in the HDF5 file using a B-tree, and the size of the B-tree, and therefore the traversal time, depends 
 * directly upon the number of chunks allocated for a dataset - so the access time for a chunk increases as
 * the number of chunks in the dataset increases.  So many small chunks will degrade performance.
 * 
 * Chunks are the unit of IO for HDF5 files - all IO to HDF5 is performed on the basis of chunks.  This means
 * that whenever dataset values are accessed (read or written (i.e. changed)), if the value is not already in
 * memory, the entire chunck containing the value must be read from, or written to, the storage media - even
 * if the dataset value being accessed is the only value in the chunk.  So few large chunks could cause 
 * empty, "wasted", space in the HDF5 files (at the end of datasets) - but they could also adversely affect
 * performance by causing unecessary IO traffic (although probably not much in the way we access data in COMPAS
 * files).
 * 
 * HDF5 files implement a chunk cache on a per-dataset basis.  The default size of the chunk cache is 1MB, and
 * its maximum size is 32MB.  The purpose of the chunk cache is to reduce storage media IO - even with SSDs,
 * memory access is way faster than storage media access, so the more of the file data that can be kept in
 * memory and maipulated there, the better.  Assuming the datatype of a particular dataset is DOUBLE, and
 * therefore consumes 8 bytes of storage space, at its maximum size the chunk cache for that dataset could hold
 * 4,000,000 values - so a single chuck with 4,000,000 values, two chunks with 2,000,000 values, four with 
 * 1,000,000, and so on.  Caching a single chunk defeats the purpose of the cache, so chunk sizes somewhat less
 * that 4,000,000 would be most appropriate if the chunk cache is to be utilised.  Chunks too big to fit in the
 * cache simply bypass the cache and are read from, or written to, the storage media directly.
 * 
 * However, the chunk cache is really only useful for random access of the dataset.  Most, if not all, of the
 * access in the COMPAS context (including post-creation analyses) is serial - the COMPAS code writes the
 * datasets from top to bottom, and later analyses (generally) read the datasets the same way.  Caching the
 * chunks for serial access just introduces overhead that costs memory (not much, to be sure: up to 32MB per 
 * open dataset), and degrades performace (albeit it a tiny bit).  For that reason I disable the chunk cache
 * in COMPAS - so all IO to/from an HDF5 file in COMPAS is directly to/from the storage media.  (To be clear,
 * post-creation analysis software can disable the cache or not when accessing HDF5 files created by COMPAS - 
 * disabling the cache here does not affect how other software accesses the files post-creation).
 * 
 * So many small chunks is not so good, and neither is just a few very large chunks.  So what's the optimum
 * chunk size?  It depends upon several things, and probably the most important of those are the final size
 * of the dataset and the access pattern.
 * 
 * As mentioned above, we tend to access datasets serially, and generally from to to bottom, so larger chunks 
 * would seem appropriate, but not so large that we generate HDF5 files with lots of unused space.  However, 
 * disk space, even SSD space, is cheap, so trading space against performance is probably a good trade.
 * 
 * Also as mentioned above, we don't know the final size of (most of) the datasets when creating the HDF5 in
 * COMPAS - though we do know the number of systems being generated, which allows us to determine an upper
 * bound for at least some of the datasets (though not for groups such as BSE_RLOF).
 * 
 * One thing we need to keep in mind is that when we create the HDF5 file we write each dataset of a group
 * in the same iteration - this is analagous to writing a single record in (e.g.) a CSV log file (the HDF5
 * group corresponds to the CSV file, and the HDF5 datasets in the group correspond to the columns in the
 * CSV file).  So for each iteration - typically each system evolved (though each timestep for detailed
 * output files) we do as many IOs to the HDF5 file as there are datasets in the group (columns in the file).
 * We are not bound to reading or writing a single chunk at a time - but we are bound to reading or writing
 * an integral multiple of whole chunks at a time.
 * 
 * We want to reduce the number of storage media accesses when writing (or later reading) the HDF5 files, so
 * larger chunk sizes are appropriate, but not so large that we create excessively large HDF5 files that have
 * lots of unused space (bearing in mind the trade-off mentioned above), especially when were evolving just
 * a few systems (rather than millions).
 * 
 * To really optimise IO performance for HDF5 files we'd choose chunk sizes that are close to multiples of
 * storage media block sizes, but I chose not to go down that rathole...
 * 
 * Based on everything written above, and some tests, I've chosen a default chunk size of 100,000 (dataset
 * entries) for all datasets (HDF5_DEFAULT_CHUNK_SIZE in constants.h).  This clearly trades performance against
 * storage space.  For the (current) default logfile record specifications, per-binary logfile space is about
 * 1K bytes, so in the very worst case we will waste some space at the end of an HDF5 output file, but the 
 * performance gain, especially for post-creation analyses, is significant.  Ensuring the number of systems 
 * evolved is an integral multiple of this fixed chunk size will minimise storage space waste.
 * 
 * I have added the program option --hdf5-chunk-size to allow users to specify the chunk size - the option
 * value defaults to HDF5_DEFAULT_CHUNK_SIZE.
 * 
 * I have chosen a minimum chunk size of 1000 (HDF5_MINIMUM_CHUNK_SIZE in constants.h).  If the number of
 * systems being evolved is >= HDF5_MINIMUM_CHUNK_SIZE the chunk size used will be the value of the hdf5-chunk-size
 * program option (either HDF5_DEFAULT_CHUNK_SIZE or a value specified by the user), but if the number of 
 * systems being evolved is < HDF5_MINIMUM_CHUNK_SIZE the chunk size used will be HDF5_MINIMUM_CHUNK_SIZE.
 * This is just so we don't waste too much storage space when running small tests - and if they are that small
 * performance is probably not going to be much of an issue, so no real trade-off against storage space.  
 * 
 * IO to HDF5 files is buffered in COMPAS - we buffer a number of chunks for each open dataset and write the
 * buffer to the file when the buffer fills (or a partial buffer upon file close if the buffer is not full).
 * This IO buffering is not HDF5 or filesystem buffering - this is a COMPAS-internal implementation to improve
 * performance.  We could do the same for the outher logfile types one day, but I suspect HDF5 is mostly
 * what people will use - and since each record written to a logfile of type other than HDF5 includes values
 * for all columns in the file, there are fewer IO operations to logfiles of type other than HDF5 (1 per system
 * being evolved) and so IO to logfiles of type other than HDF5 has a far less significant impact on performance
 * than does IO to HDF5 logfiles (where column (dataset) data is written individually).
 * 
 * The default HDF5 IO buffer size is defined in constants.h - HDF5_DEFAULT_IO_BUFFER_SIZE - and I have set it
 * to just 1 - so by default no buffering happens.  I have added the program option --hdf5-buffer-size to allow 
 * users to specify the buffer size - the option value defaults to HDF5_DEFAULT_IO_BUFFER_SIZE.  The HDF5 IO 
 * buffer size is specified as a number of chunks.  Users should increase the buffer size for better performance
 * if memory space allows it.
 * 
 * Users should bear in mind that the combination of HDF5 chunk size and HDF5 IO buffer size affect performance,
 * storage space, and memory usage - so they may need to experiment to find a balance that suits their needs.
 * 
 * 
 * String values stored in HDF5 files
 * ==================================
 *
 * COMPAS writes string data to HDF5 output files as C-type strings.  Python interprets C-type
 * strings in HDF5 files as byte arrays - regardless of the specified datatype when written (COMPAS
 * writes the strings as ASCII data (as can be seen with h5dump), but Python ignores that).  Note that
 * this affects the values in datasets (and attributes) only, not the dataset names (or group names,
 * attribute names, etc.).
 * 
 * The only real impact of this is that if the byte array is printed by Python, it will be displayed
 * as (e.g.) "b'abcde'" rather than just "abcde".  All operations on the data work as expected - it
 * is just the output that is impacted.  If that's an issue, use .decode('utf-8') to decode the byte
 * array as a Python string variable.  E.g.
 * 
 *     str = h5File[Group][Dataset][0]
 *     str is a byte array and print(str) will display (e.g.) b'abcde'
 * 
 *     but
 * 
 *     str = h5File[Group][Dataset][0].decode('utf-8')
 *     str is a Python string and print(str) will display (e.g.) abcde
 *
 * 
 * JR, January 2021
 */



/*
 * Format a boost::variant value using a format specification passed as a parameter
 *
 * This is defined as a class for use with boost::apply_visitor().
 * It is only ever used by the Log class, hence the reason it is defined here.
 *
 * This function is applied to a boost::variant value via the boost::apply_visitor() function.
 * The function extracts the underlying primitive value stored in the boost::variant value and formats
 * it using the format string constructed here specifically for the underlying primitive type.
 *
 */
class FormatVariantValue: public boost::static_visitor<string> {
public:
    string operator()(const bool               v, const string fmtStr) const {
                                                                           string fmt = OPTIONS->PrintBoolAsString() ? "%5s" : "%1s";
                                                                           string vS  = OPTIONS->PrintBoolAsString() ? (v ? "TRUE " : "FALSE") : (v ? "1" : "0");
                                                                           return utils::vFormat(fmt.c_str(), vS.c_str());
                                                                       }
    string operator()(const int                v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const short int          v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long int           v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned int       v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned short int v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned long int  v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); } // also handles OBJECT_ID (typedef)
    string operator()(const float              v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const double             v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long double        v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const string             v, const string fmtStr) const { string fmt = fmtStr; fmt = "%-" + fmt + "s"; return utils::vFormat(fmt.c_str(), v.c_str()); }
    string operator()(const ERROR              v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const STELLAR_TYPE       v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_CASE            v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_TRACKING        v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_EVENT           v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_STATE           v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
};


/*
 * Format a boost::variant value using a type-based default format specification
 *
 * This is defined as a class for use with boost::apply_visitor().
 * It is only ever used by the Log class, hence the reason it is defined here.
 *
 * This function is applied to a boost::variant value via the boost::apply_visitor() function.
 * The function extracts the underlying primitive value stored in the boost::variant value and formats
 * it using a type-based default format specification.
 *
 */
class FormatVariantValueDefault: public boost::static_visitor<string> {
public:
    string operator()(const bool               v) const {
                                                    string fmt = OPTIONS->PrintBoolAsString() ? "%5s" : "%1s";
                                                    string vS  = OPTIONS->PrintBoolAsString() ? (v ? "TRUE " : "FALSE") : (v ? "1" : "0");
                                                    return utils::vFormat(fmt.c_str(), vS.c_str());
                                                  }
    string operator()(const int                v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const short int          v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long int           v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned int       v) const { string fmt = "%14.1u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned short int v) const { string fmt = "%14.1u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned long int  v) const { string fmt = "%14.1u"; return utils::vFormat(fmt.c_str(), v); } // also handles OBJECT_ID (typedef)
    string operator()(const float              v) const { string fmt = "%16.8e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const double             v) const { string fmt = "%16.8e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long double        v) const { string fmt = "%16.8e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const string             v) const { string fmt = "%-30s";  return utils::vFormat(fmt.c_str(), v.c_str()); }
    string operator()(const ERROR              v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const STELLAR_TYPE       v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_CASE            v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_TRACKING        v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_EVENT           v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_STATE           v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
};


class Log {

private:

    Log() {                                                                         // constructor - initialise variables
        m_Enabled = false;                                                          // logging disabled initially
        m_HDF5ContainerId = -1;                                                     // no HDF5 container file open initially
        m_HDF5DetailedId = -1;                                                      // no HDF5 detailed file open initially
        m_HDF5ChunkSize = HDF5_DEFAULT_CHUNK_SIZE;                                  // initially just the defined constant
        m_HDF5IOBufSize = HDF5_DEFAULT_IO_BUFFER_SIZE * HDF5_DEFAULT_CHUNK_SIZE;    // initially just the defined constant (scaled)
        m_LogBasePath = ".";                                                        // default log file base path
        m_LogContainerName = DEFAULT_OUTPUT_CONTAINER_NAME;                         // default log file container name                        
        m_LogNamePrefix = "";                                                       // default log file name prefix
        m_LogfileType = DEFAULT_LOGFILE_TYPE;                                       // default log file type
        m_LogLevel = 0;                                                             // default log level - log everything
        m_LogClasses = {};                                                          // no default log classes
        m_DbgLevel = 0;                                                             // default debug level - debug everything
        m_DbgClasses = {};                                                          // no default debug classes
        m_DbgToLogfile = false;                                                     // default is not to log debug records to the log file
        m_DbgLogfileId = -1;                                                        // default is not valid
        m_Logfiles.empty();                                                         // default is no log files
        m_OpenStandardLogFileIds = {};                                              // no open COMPAS standard log files

        m_ObjectIdSwitching = -1L;                                                  // object id of Star object swithcing stellar type - default none
        m_TypeSwitchingFrom = STELLAR_TYPE::NONE;                                   // stellar type from which Star object is switching - default NONE
        m_TypeSwitchingTo   = STELLAR_TYPE::NONE;                                   // stellar type to which Star object is switching - default NONE
        m_PrimarySwitching  = false;                                                // Star swithcing is primary star of binary - default false

        m_SSESupernova_DelayedLogRecord    = "";                                    // delayed log record for SSE_Supernova file - initially empty
        m_SSESupernova_LogRecordProperties = {};                                    // SSE Supernova logfile record properties - initially empty
        m_SSESupernova_LogRecordFmtVector  = {};                                    // SSE Supernova logfile format vector - initially empty
    };
    Log(Log const&) = delete;                                                       // copy constructor does nothing, and not exposed publicly
    Log& operator = (Log const&) = delete;                                          // operator = does nothing, and not exposed publicly

    // instance variable
    static Log       *m_Instance;                                                   // pointer to the instance


    // member variables
    bool                 m_Enabled;                                                 // is logging enabled?

    string               m_HDF5ContainerName;                                       // HDF5 container name
    hid_t                m_HDF5ContainerId;                                         // HDF5 container id
    hid_t                m_HDF5DetailedId;                                          // HDF5 detailed output id
    size_t               m_HDF5ChunkSize;                                           // HDF5 chunk size
    size_t               m_HDF5IOBufSize;                                           // HDF5 IO buffer size

    string               m_LogBasePath;                                             // base path for log files
    string               m_LogContainerName;                                        // container (directory) name for log files
    string               m_LogNamePrefix;                                           // prefix for log files

    LOGFILETYPE          m_LogfileType;                                             // logfile type
    int                  m_LogLevel;                                                // log level
    std::vector <string> m_LogClasses;                                              // log classes

    int                  m_DbgLevel;                                                // debug level
    std::vector <string> m_DbgClasses;                                              // debug classes

    bool                 m_DbgToLogfile;                                            // log debug records to log file?
    int                  m_DbgLogfileId;                                            // log file id of file to which debug statements should be written

    bool                 m_ErrToLogfile;                                            // log error records to log file?
    int                  m_ErrLogfileId;                                            // log file id of file to which error statements should be written


    struct logfileAttrT {
        bool        active;                                                         // currently logging?  (ie log file open)
        LOGFILE     logfiletype;                                                    // logfile type
        LOGFILETYPE filetype;                                                       // file type
        string      name;                                                           // name of log file
        bool        timestamp;                                                      // time stamp enabled?
        bool        label;                                                          // record labels enabled?

        std::ofstream file;                                                         // file pointer for CSV, TSV, TXT files

        struct h5AttrT {                                                            // attributes of HDF5 files
            hid_t fileId;                                                           //    - file id
            hid_t groupId;                                                          //    - group id

            struct h5DataSetsT {                                                    // attributes of HDF5 datasets
                hid_t    dataSetId;                                                 //    - HDF5 dataset id
                hid_t    h5DataType;                                                //    - HDF5 datatype
                TYPENAME dataType;                                                  //    - COMPAS data type
                std::vector<COMPAS_VARIABLE_TYPE> buf;                              //    - write buffer - for chunking
            };

            std::vector<h5DataSetsT> dataSets;                                      // details of datasets
        };

        h5AttrT h5File;                                                             // file details for HDF5 files
    };

    std::vector<logfileAttrT> m_Logfiles;                                           // logfiles - in use and not

    COMPASUnorderedMap<LOGFILE, LogfileDetailsT> m_OpenStandardLogFileIds;          // currently open standard logfiles: id, filename, property details, field format strings

    ANY_PROPERTY_VECTOR m_BSE_BE_Binaries_Rec = BSE_BE_BINARIES_REC;                // default specification
    ANY_PROPERTY_VECTOR m_BSE_CEE_Rec         = BSE_COMMON_ENVELOPES_REC;           // default specification
    ANY_PROPERTY_VECTOR m_BSE_DCO_Rec         = BSE_DOUBLE_COMPACT_OBJECTS_REC;     // default specification
    ANY_PROPERTY_VECTOR m_BSE_Detailed_Rec    = BSE_DETAILED_OUTPUT_REC;            // default specification
    ANY_PROPERTY_VECTOR m_BSE_Pulsars_Rec     = BSE_PULSAR_EVOLUTION_REC;           // default specification
    ANY_PROPERTY_VECTOR m_BSE_RLOF_Rec        = BSE_RLOF_PARAMETERS_REC;            // default specification
    ANY_PROPERTY_VECTOR m_BSE_SNE_Rec         = BSE_SUPERNOVAE_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_BSE_Switch_Rec      = BSE_SWITCH_LOG_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_BSE_SysParms_Rec    = BSE_SYSTEM_PARAMETERS_REC;          // default specification

    ANY_PROPERTY_VECTOR m_SSE_Detailed_Rec    = SSE_DETAILED_OUTPUT_REC;            // default specification
    ANY_PROPERTY_VECTOR m_SSE_SNE_Rec         = SSE_SUPERNOVAE_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_SSE_Switch_Rec      = SSE_SWITCH_LOG_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_SSE_SysParms_Rec    = SSE_SYSTEM_PARAMETERS_REC;          // default specification


    // the following block of variables support the BSE Switch Log file
    
    OBJECT_ID    m_ObjectIdSwitching;                                               // the object id of the Star object switching stellar type
    STELLAR_TYPE m_TypeSwitchingFrom;                                               // the stellar type from which the Star object is switching
    STELLAR_TYPE m_TypeSwitchingTo;                                                 // the stellar type to which the Star object is switching
    bool         m_PrimarySwitching;                                                // flag to indicate whether the primary star of the binary is switching


    // the following block of variables support delayed writes to logfiles
    //
    // For now, only delayed writes to the SSE Supernova file is implemented, and only a single record can be delayed.
    // This functionality was introduced specifically to allow queueing a delayed write the to SSE Supernova file - see
    // the discussion in the description of Log::GetStandardLogFileRecordDetails() in Log.cpp.
    // This functionality probably shouldn't be extended to allow queueing/delaying multiple records for later writing
    // (I don't think we need it, it would probably soak up too much memory if over-used, and it might just cause confusion)
    
    string              m_SSESupernova_DelayedLogRecord;                            // log record to be written to SSE Supernova log file in delayed write
    ANY_PROPERTY_VECTOR m_SSESupernova_LogRecordProperties;                         // SSE Supernova logfile record properties
    std::vector<string> m_SSESupernova_LogRecordFmtVector;                          // SSE Supernova logfile format vector
    
  
    // the following block of variables support the run details file

    std::ofstream                                      m_RunDetailsFile;            // run details file
    std::chrono::time_point<std::chrono::system_clock> m_WallStart;                 // for run details file
    std::chrono::time_point<std::chrono::system_clock> m_WallEnd;                   // for run details file
    clock_t                                            m_ClockStart;                // for run details file

    // member functions

    bool IsValidId(const int p_LogfileId)    { return ((p_LogfileId >= 0) && ((unsigned int)p_LogfileId < m_Logfiles.size())); }
    bool IsActiveId(const int p_LogfileId)   { return IsValidId(p_LogfileId) && m_Logfiles[p_LogfileId].active; }

    void ClearEntry(const int p_LogfileId) {
        if (IsValidId(p_LogfileId)) {
            m_Logfiles[p_LogfileId].active          = false;                       // not active
            m_Logfiles[p_LogfileId].logfiletype     = LOGFILE::NONE;
            m_Logfiles[p_LogfileId].filetype        = LOGFILETYPE::NONE;
            m_Logfiles[p_LogfileId].name            = "";
            m_Logfiles[p_LogfileId].timestamp       = false;
            m_Logfiles[p_LogfileId].label           = false;
            m_Logfiles[p_LogfileId].h5File.fileId   = -1;
            m_Logfiles[p_LogfileId].h5File.groupId  = -1;
            m_Logfiles[p_LogfileId].h5File.dataSets = {};
        }
    }

    bool DoIt(const string p_Class, const int p_Level, const std::vector<string> p_EnabledClasses, const int p_EnabledLevel);
    void Say_(const string p_SayStr);
    bool Write_(const int p_LogfileId, const string p_LogStr);
    bool Write_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const bool p_Flush = false);
    bool Flush_(const int p_LogfileId) { return Write_(p_LogfileId, {}, true); }
    bool Put_(const int p_LogfileId, const string p_LogStr, const string p_Label = "");
    bool Put_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues);
    bool Debug_(const string p_DbgStr);
    bool Close_(const int p_LogfileId);

    PROPERTY_DETAILS StellarPropertyDetails(ANY_STAR_PROPERTY p_Property);
    PROPERTY_DETAILS BinaryPropertyDetails(BINARY_PROPERTY p_Property);
    PROPERTY_DETAILS ProgramOptionDetails(PROGRAM_OPTION p_Property);
    STR_STR_STR_STR  FormatFieldHeaders(PROPERTY_DETAILS p_Details, string p_HeaderSuffix = "");
    LogfileDetailsT  StandardLogFileDetails(const LOGFILE p_Logfile, const string p_FileSuffix);


    std::tuple<bool, LOGFILE> GetLogfileDescriptorKey(const string p_Value);
    std::tuple<bool, LOGFILE> GetStandardLogfileKey(const int p_FileId);


    /*
     * Construct a log record to be written a standard logfile
     *
     * This function constructs a log record to be written to one of the standard COMPAS logfiles.
     * The record to be constructed is identified by the logfile to which it is to be written, and the data
     * is assembled on-the-fly - except possibly for the parameter p_SpecifiedStellarType (see below).  The 
     * star from which the data should be gathered is passed as a parameter, as is the logfile for which the
     * record should be constructed.  The logfile record properties and format vector appropriate to the
     * specified logfile are also passed as parameters - this function does not need to retrieve/construct
     * them.
     *
     * This function will, if the value of parameter p_UseSpecifiedValue is true, replace the value of any 
     * property in the record matching the value of the parameter p_SpecifiedProperty with the value specified
     * by the parameter p_SpecifiedPropertyValue.  (See the variant function GetLogStandardRecord() for an
     * explanation of why p_UseSpecifiedValue is required).
     * 
     *
     * template <class T1, typename T2>
     * string GetLogStandardRecord(const LOGFILE             p_LogFile,
     *                             const T1* const           p_Star,
     *                             const ANY_PROPERTY_VECTOR p_RecordProperties,
     *                             const std::vector<string> p_FmtVector,
     *                             const bool                p_UseSpecifiedValue,
     *                             const ANY_STAR_PROPERTY   p_SpecifiedProperty,
     *                             const T2                  p_SpecifiedPropertyValue)
     *
     * @param   [IN]    p_LogFile                   The logfile for which the record should be constructed
     * @param   [IN]    p_Star                      The object from which the values to be used to construct the record should be retrieved
     * @param   [IN]    p_RecordProperties          The logfile record properties pertaining to p_LogFile
     * @param   [IN]    p_FmtVector                 The logfile format vector pertaining to p_LogFile
     * @param   [IN]    p_UseSpecifiedValue         Flag to specify whether the parameters p_SpecifiedProperty and p_SpecifiedPropertyValue
     *                                              should be used to replace values retrieved from p_Star
     * @param   [IN]    p_SpecifiedProperty         The property type of the value to be replaced by p_SpecifiedPropertyValue
     * @param   [IN]    p_SpecifiedPropertyValue    The value of the property to be replaced
     * 
     * @return                                      String formatted as log file record - empty string if an error occurred
     */
    template <class T1, typename T2>
    string GetLogStandardRecord(const LOGFILE             p_LogFile,
                                const T1* const           p_Star,
                                const ANY_PROPERTY_VECTOR p_RecordProperties,
                                const std::vector<string> p_FmtVector,
                                const bool                p_UseSpecifiedValue,
                                const ANY_STAR_PROPERTY   p_SpecifiedProperty,
                                const T2                  p_SpecifiedPropertyValue) {

        // construct log record from current data

        string logRecord = "";                                                              // the formatted logfile record

        // set delimiter based on logfile type
        string delimiter = "";                                                              // default
        switch (OPTIONS->LogfileType()) {
            case LOGFILETYPE::HDF5: delimiter = ""; break;                                  // HDF5
            case LOGFILETYPE::CSV : delimiter = DELIMITERValue.at(DELIMITER::COMMA); break; // CSV
            case LOGFILETYPE::TSV : delimiter = DELIMITERValue.at(DELIMITER::TAB); break;   // TSV
            case LOGFILETYPE::TXT : delimiter = DELIMITERValue.at(DELIMITER::SPACE); break; // TXT
            default               : delimiter = ""; break;                                  // default
        }

        // get and format values for printing
        bool                 ok;                                                            // flag to indicate property value retrieved ok (or not)
        COMPAS_VARIABLE_TYPE value;                                                         // property value
        string               valueStr;                                                      // string for formatted value

        int index = 0;
        for (auto &property : p_RecordProperties) {                                         // for each property to be included in the log record

            // determine if this property is p_SpecifiedStellarProperty - if it is, replace the value with p_SpecifiedStellarType
            // (this is specific to replacing a STELLAR_TYPE property - may need to generalise this to other types later - cross that bridge then)

            ok = true;
            ANY_STAR_PROPERTY thisProperty;
            switch (boost::apply_visitor(VariantPropertyType(), property)) {
                case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : { STAR_PROPERTY      prop = boost::get<STAR_PROPERTY>(property);      thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : { STAR_1_PROPERTY    prop = boost::get<STAR_1_PROPERTY>(property);    thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : { STAR_2_PROPERTY    prop = boost::get<STAR_2_PROPERTY>(property);    thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: { SUPERNOVA_PROPERTY prop = boost::get<SUPERNOVA_PROPERTY>(property); thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: { COMPANION_PROPERTY prop = boost::get<COMPANION_PROPERTY>(property); thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                default: ok = false; // unknown property type - that's not ok...
            }

            // if not ok at this point we'll just take the current value - if there is really a problem
            // that will fail (as it would have normally if we weren't trying to replace a value)
            
            boost::variant<string> fmtStr(p_FmtVector[index++]);                            // format string for this property

            if (ok && p_UseSpecifiedValue && (thisProperty == p_SpecifiedProperty)) {

                // use value passed as parameter

                value    = p_SpecifiedPropertyValue;
                valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);       // format value
            }
            else {

                // use current value

                std::tie(ok, value) = p_Star->PropertyValue(property);                      // get property flag and value
                if (ok) {                                                                   // have valid property value
                    valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);   // format value
                    logRecord += valueStr + delimiter;                                      // add value string to log record - with delimiter
                }
                else {                                                                      // error formatting value 
                    logRecord = "";                                                         // empty record
                    break;                                                                  // stop now
                }
            }
        }

        if (ok) {
            // if we are constructing an SSE Switch file we add two pre-defined
            // columns to the end of the log record.  These are:
            //
            // ( i) the stellar type from which the star is switching
            // (ii) the stellar type to which the star is switching
            //
            // if we are constructing a BSE Switch file record we add three pre-defined
            // columns to the end of the log record.  These are:
            //
            // (  i) the star switching - 1 = primary, 2 = secondary
            // ( ii) the stellar type from which the star is switching
            // (iii) the stellar type to which the star is switching
            //
            // These are hard-coded here rather than in the *_PROPERTY_DETAIL maps in
            // constants.h so that they will always be present in the switch file -
            // this way users can't add or remove them at runtime via the logfile-definitions
            // option.

            if (p_LogFile == LOGFILE::BSE_SWITCH_LOG) {
                int starSwitching = m_PrimarySwitching ? 1 : 2;                             // primary (1) or secondary (2)
                string fmt        = "%14.1d";                                               // format specifier
                logRecord += utils::vFormat(fmt.c_str(), starSwitching) + delimiter;        // star switching
            }

            if (p_LogFile == LOGFILE::BSE_SWITCH_LOG || p_LogFile == LOGFILE::SSE_SWITCH_LOG) {
                string fmt = "%14.1d";                                                      // format specifier
                logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingFrom) + delimiter;  // switching from
                fmt        = "%12.1d";                                                      // format specifier
                logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingTo) + delimiter;    // switching to
            }

            logRecord = logRecord.substr(0, logRecord.size() - 1);                          // remove the last character - extraneous delimiter
        }

        return logRecord;
    }

    /*
     * This variant of GetLogStandardRecord() is here because I can't readily figure out how to 
     * have an optional parameter in a template function if not providing that optional parameter 
     * doesn't give the compiler enough context to deduce the type of the template type.
     * 
     * Clear as mud?  Well, for the declaration of GetLogStandardRecord() above, template type T2
     * can't be deduced if the parameter p_SpecifiedPropertyValue is not supplied - even if it
     * is given a default value in the parameter list.  There must be a way of making it work, but
     * I can't see it - and I don't want to spend any more time trying to figure it out.  Maybe 
     * later.  In the meantime, this method works.
     * 
     * See description of GetLogStandardRecord() for functionality and parameter descriptions.
     * 
     */
    template <class T>
    string GetLogStandardRecord(const LOGFILE             p_LogFile,
                                const T* const            p_Star,
                                const ANY_PROPERTY_VECTOR p_RecordProperties,
                                const std::vector<string> p_FmtVector) {

        return GetLogStandardRecord(p_LogFile, 
                                    p_Star, 
                                    p_RecordProperties, 
                                    p_FmtVector, 
                                    false, 
                                    ANY_STAR_PROPERTY::STELLAR_TYPE, 
                                    STELLAR_TYPE::NONE);
    }


    /*
     * Write a record to a standard logfile
     *
     * This function writes a log record to one of the standard COMPAS logfiles.
     * The record to be written is identified by the logfile to which it is to be written, and the data
     * is assembled on-the-fly.  The star from which the data should be gathered is passed as a parameter,
     * as is the logfile to which the record should be written, as well as a string suffix for the logfile
     * name.
     *
     * If the logfile is not already open, it will be opened by calling this function.
     *
     * The logging classes and level set by program options are honoured by this function.  The log class
     * and level for this record are passed as parameters.  If the log class for this record is not
     * configured to be logged, this function returns without performing any logging action (i.e. the log
     * file status (open, closed) and contents, will be unchanged by this function.  Similarly, if the
     * log level for this record indicates that the record should not be logged, no changes will be made
     * by this function.
     * 
     *
     * template <class T>
     * bool LogStandardRecord(const string   p_LogClass,
     *                        const int      p_LogLevel,
     *                        const LOGFILE  p_LogFile,
     *                        const T* const p_Star,
     *                        const string   p_LogRecord,
     *                        const string   p_FileSuffix = "")
     *
     * @param   [IN]    p_LogClass                  Class to determine if record should be written
     * @param   [IN]    p_LogLevel                  Level to determine if record should be written
     * @param   [IN]    p_LogFile                   The logfile to which the record should be written
     * @param   [IN]    p_Star                      The star object from which the field values should be retrieved
     * @param   [IN]    p_LogRecord                 Formatted record to be written to the logfile
     *                                                  If this parameter is supplied as an empty string the function will construct
     *                                                  a record from current data.
     *                                                  If this parameter is supplied as a non-empty string it will be written to the
     *                                                  file as is (the function will not construct a record from current data)
     *                                              This parameter is currently ignored for HDF5 files - all standard logfiles are
     *                                              comprised of columns of values that can be written as a single, delimited, string
     *                                              for CSV, TSV, and TXT files, but not for HDF5 files.  HDF5 files require the 
     *                                              actual values to be written into datasets of the correct datatype.  So until we
     *                                              introduce a standard logfile that is actually just a single string (single column),
     *                                              p_LogRecord is ignored for HDF5 files.
     * @param   [IN]    p_FileSuffix                String suffix to be added to the logfile name (optional, default = "")
     * @return                                      Boolean status (true = success, false = failure)
     */
    template <class T>
    bool LogStandardRecord(const string   p_LogClass,
                           const int      p_LogLevel,
                           const LOGFILE  p_LogFile,
                           const T* const p_Star,
                           const string   p_LogRecord,
                           const string   p_FileSuffix = "") {

        bool ok = true;

        LogfileDetailsT fileDetails;                                                                                                    // file details

        fileDetails = StandardLogFileDetails(p_LogFile, p_FileSuffix);                                                                  // get record details - open file (if necessary)
        if (fileDetails.id >= 0) {                                                                                                      // file open?
                                                                                                                                        // yes
            // set delimiter based on logfile type
            string delimiter = "";                                                                                                      // default
            switch (m_Logfiles[fileDetails.id].filetype) {
                case LOGFILETYPE::HDF5: delimiter = ""; break;                                                                          // HDF5
                case LOGFILETYPE::CSV : delimiter = DELIMITERValue.at(DELIMITER::COMMA); break;                                         // CSV
                case LOGFILETYPE::TSV : delimiter = DELIMITERValue.at(DELIMITER::TAB); break;                                           // TSV
                case LOGFILETYPE::TXT : delimiter = DELIMITERValue.at(DELIMITER::SPACE); break;                                         // TXT
                default               : delimiter = ""; break;                                                                          // default
            }

            std::vector<COMPAS_VARIABLE_TYPE> logRecordValues = {};                                                                     // for HDF5 files: vector of values to be written
            string logRecord = "";                                                                                                      // for CVS, TSV, TXT files: the record to be written to the log file

            if (p_LogRecord.empty() || m_Logfiles[fileDetails.id].filetype == LOGFILETYPE::HDF5) {                                      // logfile record passed in is empty, or ignored for HDF5 files
                                                                                                                                        // construct log record from current data
                ANY_PROPERTY_VECTOR properties = fileDetails.recordProperties;                                                          // vector of properties to be printed

                // get values
                //    - format for printing for CSV, TSV and TXT files
                //    - record for HDF5 files
                COMPAS_VARIABLE_TYPE value;                                                                                             // property value
                string               valueStr;                                                                                          // string for formatted value

                int index = 0;
                for (auto &property : properties) {                                                                                     // for each property to be included in the log record
                    std::tie(ok, value) = p_Star->PropertyValue(property);                                                              // get property flag and value
                    if (ok) {                                                                                                           // have valid property value
                        if (m_Logfiles[fileDetails.id].filetype == LOGFILETYPE::HDF5) {                                                 // yes - HDF5 file?
                            logRecordValues.push_back(value);                                                                           // add value to vector of values
                        }
                        else {                                                                                                          // no - CSV, TSV, or TXT file
                            boost::variant<string> fmtStr(fileDetails.fmtStrings[index++]);                                             // get format string
                            valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);                                       // format value
                            logRecord += valueStr + delimiter;                                                                          // add value string to log record - with delimiter
                        }
                    }
                    else {                                                                                                              // unknown property type - should never happen
                        Squawk(ERR_MSG(ERROR::UNKNOWN_PROPERTY_TYPE) + " while writing to logfile " + fileDetails.filename);            // show warning
                        ok = false;                                                                                                     // fail
                        break;                                                                                                          // stop now
                    }
                }

                if (ok) {
                    // if we are writing to the SSE Switch file we add two pre-defined columns
                    // to the end of the log record.  These are:
                    //
                    // ( i) the stellar type from which the star is switching
                    // (ii) the stellar type to which the star is switching
                    //
                    // if we are writing to the BSE Switch file we add three pre-defined columns
                    // to the end of the log record.  These are:
                    //
                    // (  i) the star switching - 1 = primary, 2 = secondary
                    // ( ii) the stellar type from which the star is switching
                    // (iii) the stellar type to which the star is switching

                    string fmt = "%4.1d";                                                                                               // format - all integers here
                    if (p_LogFile == LOGFILE::BSE_SWITCH_LOG) {
                        int starSwitching = m_PrimarySwitching ? 1 : 2;                                                                 // primary (1) or secondary (2)
                        logRecord += utils::vFormat(fmt.c_str(), starSwitching) + delimiter;                                            // star switching
                    }

                    if (p_LogFile == LOGFILE::BSE_SWITCH_LOG || p_LogFile == LOGFILE::SSE_SWITCH_LOG) {
                        logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingFrom) + delimiter;                                      // switching from
                        logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingTo) + delimiter;                                        // switching to
                    }

                    logRecord = logRecord.substr(0, logRecord.size()-1);                                                                // remove the last character - extraneous delimiter
                }
            }
            else {                                                                                                                      // logfile record passed in is not empty
                logRecord = p_LogRecord;                                                                                                // use logfile record passed in
            }

            if (ok) {                                                                                                                   // if all ok, write the record
                if (m_Logfiles[fileDetails.id].filetype == LOGFILETYPE::HDF5) {                                                         // HDF5 file?
                    ok = Put_(fileDetails.id, logRecordValues);                                                                         // yes - write the record
                }
                else {                                                                                                                  // no - CSV, TSV, or TXT file
                    ok = Put_(fileDetails.id, logRecord);                                                                               // write the record
                }
                if (!ok) Squawk(ERR_MSG(ERROR::FILE_WRITE_ERROR) + " while writing to logfile " + fileDetails.filename);                // show warning if record not written ok
            }
        }

        return ok;
    }

    void PrintLogfileRecordDetails(const ANY_PROPERTY_VECTOR& p_LogfileRecord, const string p_LogfileRecordName);

    void UpdateLogfileRecordSpecs(const LOGFILE             p_Logfile,
                                  bool                      p_UseDefaultProps,
                                  const ANY_PROPERTY_VECTOR p_AddProps,
                                  const ANY_PROPERTY_VECTOR p_SubtractProps);

    bool UpdateAllLogfileRecordSpecs();


public:

    ~Log() {                                                                                                                            // destructor
        Stop();                                                                                                                         // stop logging - flushes and closes all open logfiles
        delete m_Instance;                                                                                                              // delete the instance variable
    }

    // instance function
    static Log* Instance();                                                                                                             // the singleton instance is exposed, not the constructor


    // member functions
    void   Start(const string              p_LogBasePath,
                 const string              p_LogContainerName,
                 const string              p_LogNamePrefix,
                 const int                 p_LogLevel,
                 const std::vector<string> p_LogClasses,
                 const int                 p_DbgLevel,
                 const std::vector<string> p_DbgClasses,
                 const bool                p_DbgToFile,
                 const bool                p_ErrorsToFile,
                 const LOGFILETYPE         p_LogfileType);

    void   Stop(std::tuple<int, int> p_ObjectStats = std::make_tuple(0, 0));

    bool   Enabled() const { return m_Enabled; }

    int    Open(const string p_LogFileName, const bool p_Append, const bool p_TimeStamp, const bool p_Label, const LOGFILE p_StandardLogfile = LOGFILE::NONE);
    bool   Close(const int p_LogfileId);

    bool   Write(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr);
    bool   Put(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr);

    bool   Debug(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr);
    bool   DebugWait(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr);

    bool   Error(const string p_ErrStr);

    void   Squawk(const string squawkStr);

    void   Say(const string p_SayClass, const int p_SayLevel, const string p_SayStr);


    // SetSwitchParameters is called by Star::SwitchTo to set the parameters 
    // to be written to the BSE Switch Log file
    void   SetSwitchParameters(const OBJECT_ID    p_ObjectIdSwitching, 
                               const STELLAR_TYPE p_TypeSwitchingFrom, 
                               const STELLAR_TYPE p_TypeSwitchingTo) {
        m_ObjectIdSwitching = p_ObjectIdSwitching;                          // the object id of the Star object switching stellar type
        m_TypeSwitchingFrom = p_TypeSwitchingFrom;                          // the stellar type from which the Star object is switching
        m_TypeSwitchingTo   = p_TypeSwitchingTo;                            // the stellar type to which the Star object is switching
    }

    OBJECT_ID ObjectIdSwitching() { return m_ObjectIdSwitching; }



    // standard logfile logging functions

    bool CloseStandardFile(const LOGFILE p_LogFile, const bool p_Erase = true);
    bool CloseAllStandardFiles();

    std::tuple<ANY_PROPERTY_VECTOR, std::vector<string>> GetStandardLogFileRecordDetails(const LOGFILE p_Logfile);

    template <class T>
    bool LogBeBinary(const T* const p_Binary, const string p_Rec)                               { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_BE_BINARIES)), 0, LOGFILE::BSE_BE_BINARIES, p_Binary, p_Rec); }

    template <class T>
    bool LogBSEDetailedOutput(const T* const p_Binary, const long int p_Id, const string p_Rec) { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT)), 0, LOGFILE::BSE_DETAILED_OUTPUT, p_Binary, p_Rec, "_" + std::to_string(abs(p_Id))); }

    template <class T>
    bool LogBSEPulsarEvolutionParameters(const T* const p_Binary, const string p_Rec)           { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION)), 0, LOGFILE::BSE_PULSAR_EVOLUTION, p_Binary, p_Rec); }

    template <class T>
    bool LogBSESupernovaDetails(const T* const p_Binary, const string p_Rec)                    { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE)), 0, LOGFILE::BSE_SUPERNOVAE, p_Binary, p_Rec); }
    
    template <class T>
    bool LogBSESwitchLog(const T* const p_Binary, const long int p_Id, const bool p_PrimarySwitching) {
        m_PrimarySwitching = p_PrimarySwitching;        
        return LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG)), 0, LOGFILE::BSE_SWITCH_LOG, p_Binary, "");
    }

    template <class T>
    bool LogBSESystemParameters(const T* const p_Binary, const string p_Rec)                    { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SYSTEM_PARAMETERS)), 0, LOGFILE::BSE_SYSTEM_PARAMETERS, p_Binary, p_Rec); }

    template <class T>
    bool LogCommonEnvelope(const T* const p_Binary, const string p_Rec)                         { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_COMMON_ENVELOPES)), 0, LOGFILE::BSE_COMMON_ENVELOPES, p_Binary, p_Rec); }

    template <class T>
    bool LogDoubleCompactObject(const T* const p_Binary, const string p_Rec)                    { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS)), 0, LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS, p_Binary, p_Rec); }

    template <class T>
    bool LogRLOFParameters(const T* const p_Binary, const string p_Rec)                         { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_RLOF_PARAMETERS)), 0, LOGFILE::BSE_RLOF_PARAMETERS, p_Binary, p_Rec); }

    template <class T>
    bool LogSSEDetailedOutput(const T* const p_Star, const int p_Id, const string p_Rec)        { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_DETAILED_OUTPUT)), 0, LOGFILE::SSE_DETAILED_OUTPUT, p_Star, p_Rec, "_" + std::to_string(abs(p_Id))); }

    template <class T>
    bool LogSSESupernovaDetails(const T* const p_Star, const string p_Rec)                      { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE)), 0, LOGFILE::SSE_SUPERNOVAE, p_Star, p_Rec); }

    template <class T>
    bool LogSSESwitchLog(const T* const p_Star, const int p_Id, const string p_Rec)             { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SWITCH_LOG)), 0, LOGFILE::SSE_SWITCH_LOG, p_Star, p_Rec); }

    template <class T>
    bool LogSSESystemParameters(const T* const p_Star, const string p_Rec)                      { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SYSTEM_PARAMETERS)), 0, LOGFILE::SSE_SYSTEM_PARAMETERS, p_Star, p_Rec); }


    template <class T>
    void StashSSESupernovaDetails(const T* const p_Star, const STELLAR_TYPE p_StellarType) {

        // if we don't already have the SSE Supernova log record properties and format vector, get them
        // this will only need to be done one per run
        if (m_SSESupernova_LogRecordProperties.empty() || m_SSESupernova_LogRecordFmtVector.empty()) {
            std::tie(m_SSESupernova_LogRecordProperties, m_SSESupernova_LogRecordFmtVector) = LOGGING->GetStandardLogFileRecordDetails(LOGFILE::SSE_SUPERNOVAE);
        }

        // get a formatted record with current data
        // this will replace any existing stashed record - no queue here
        m_SSESupernova_DelayedLogRecord = GetLogStandardRecord(LOGFILE::SSE_SUPERNOVAE, 
                                                               p_Star, 
                                                               m_SSESupernova_LogRecordProperties, 
                                                               m_SSESupernova_LogRecordFmtVector, 
                                                               true, 
                                                               (ANY_STAR_PROPERTY)STAR_PROPERTY::STELLAR_TYPE, 
                                                               p_StellarType);
    }

    template <class T>
    bool LogStashedSSESupernovaDetails(const T* const p_Star) { 

        bool result = true;

        // if the stashed SSE Supernova record is non-empty, print it, then clear it - otherwise do nothing
        if (!m_SSESupernova_DelayedLogRecord.empty()) {

            result = LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE)), 0, LOGFILE::SSE_SUPERNOVAE, p_Star, m_SSESupernova_DelayedLogRecord);
            m_SSESupernova_DelayedLogRecord = "";
        }

        return result;
    }
};

#endif // __Log_h__
