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

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Options.h"
#include "LogMacros.h"

using std::string;
using std::get;


/*
 * Log Singleton
 *
 * Provides global logging functionality
 *
 * Singletons are sometimes frowned-upon - but Logging functionality is one area
 * for which they are considered acceptable.  Doing it this way means the logging
 * class instance doesn't need to be passed around to all and sundry.
 *
 * Jeff Riley, May 2019
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
  * Jeff Riley, September 2019
  */



/*
 * Format a boost::variant value
 *
 * This is defined as a class for use with boost::apply_visitor().
 * It is only ever used by the Log class, hence the reason it is define here.
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
    string operator()(const unsigned long int  v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); }     // also handles OBJECT_ID (typedef)
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


class Log {

private:

    Log() {                                                                         // constructor - initialise variables
        m_Enabled = false;                                                          // logging disabled initially
        m_LogBasePath = ".";                                                        // default log file base path
        m_LogContainerName = DEFAULT_OUTPUT_CONTAINER_NAME;                         // default log file container name                        
        m_LogNamePrefix = "";                                                       // default log file name prefix
        m_Delimiter = DELIMITERValue.at(DELIMITER::TAB);                            // default delimiter is TAB
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

    string               m_LogBasePath;                                             // base path for log files
    string               m_LogContainerName;                                        // container (directory) name for log files
    string               m_LogNamePrefix;                                           // prefix for log files

    string               m_Delimiter;                                               // filed delimiter for logging

    int                  m_LogLevel;                                                // log level
    std::vector <string> m_LogClasses;                                              // log classes

    int                  m_DbgLevel;                                                // debug level
    std::vector <string> m_DbgClasses;                                              // debug classes

    bool                 m_DbgToLogfile;                                            // log debug records to log file?
    int                  m_DbgLogfileId;                                            // log file id of file to which debug statements should be written

    bool                 m_ErrToLogfile;                                            // log error records to log file?
    int                  m_ErrLogfileId;                                            // log file id of file to which error statements should be written


    struct logfileAttr {
        bool          active;                                                       // currently logging?  (ie log file open)
        std::ofstream file;                                                         // file pointer
        string        name;                                                         // name of log file
        string        delimiter;                                                    // field delimiter string
        bool          timestamp;                                                    // time stamp enabled?
        bool          label;                                                        // record labels enabled?
    };

    std::vector<logfileAttr> m_Logfiles;                                            // logfiles - in use and not

    COMPASUnorderedMap<LOGFILE, LOGFILE_DETAILS> m_OpenStandardLogFileIds;          // currently open standard logfiles: fileId, property details, field format strings

    ANY_PROPERTY_VECTOR m_SSE_Parms_Rec       = SSE_PARAMETERS_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_SSE_Switch_Rec      = SSE_SWITCH_LOG_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_SSE_SN_Rec          = SSE_SUPERNOVA_REC;                  // default specification
    ANY_PROPERTY_VECTOR m_BSE_SysParms_Rec    = BSE_SYSTEM_PARAMETERS_REC;          // default specification
    ANY_PROPERTY_VECTOR m_BSE_DCO_Rec         = BSE_DOUBLE_COMPACT_OBJECTS_REC;     // default specification
    ANY_PROPERTY_VECTOR m_BSE_SNE_Rec         = BSE_SUPERNOVAE_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_BSE_CEE_Rec         = BSE_COMMON_ENVELOPES_REC;           // default specification
    ANY_PROPERTY_VECTOR m_BSE_RLOF_Rec        = BSE_RLOF_PARAMETERS_REC;            // default specification
    ANY_PROPERTY_VECTOR m_BSE_BE_Binaries_Rec = BSE_BE_BINARIES_REC;                // default specification
    ANY_PROPERTY_VECTOR m_BSE_Pulsars_Rec     = BSE_PULSAR_EVOLUTION_REC;           // default specification
    ANY_PROPERTY_VECTOR m_BSE_Detailed_Rec    = BSE_DETAILED_OUTPUT_REC;            // default specification
    ANY_PROPERTY_VECTOR m_BSE_Switch_Rec      = BSE_SWITCH_LOG_REC;                 // default specification


    // the following block of variables support the BSE Switch Log file
    
    OBJECT_ID    m_ObjectIdSwitching;                                               // the object id of the Star object switching stellar type
    STELLAR_TYPE m_TypeSwitchingFrom;                                               // the stellar type from which the Star object is switching
    STELLAR_TYPE m_TypeSwitchingTo;                                                 // the stellar type to which the Star object is switching
    bool         m_PrimarySwitching;                                                // flag to indicate whether the primary star of the binary is switching


    // the following block of variables support delayed writes to logfiles
    //
    // For now, only delayed writes to the SSE Supernova file is implemented, and only a single record can be delayed.
    // This functionality was introduced specifically to allow queueing a delayed write the to SSE Supernova file - see
    // the discusiion in the description of Log::GetStandardLogFileRecordDetails() in Log.cpp.
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
            m_Logfiles[p_LogfileId].active    = false;                              // not active
            m_Logfiles[p_LogfileId].delimiter = "";                                 // clear delimiter
            m_Logfiles[p_LogfileId].name      = "";                                 // clear name
            m_Logfiles[p_LogfileId].timestamp = false;                              // clear timestamp
            m_Logfiles[p_LogfileId].label     = false;                              // clear label
        }
    }

    bool DoIt(const string p_Class, const int p_Level, const std::vector<string> p_EnabledClasses, const int p_EnabledLevel);
    void Say_(const string p_SayStr);
    bool Write_(const int p_LogfileId, const string p_LogStr);
    bool Put_(const int p_LogfileId, const string p_LogStr, const string p_Label = "");
    bool Debug_(const string p_DbgStr);
    bool Close_(const int p_LogfileId);

    PROPERTY_DETAILS StellarPropertyDetails(ANY_STAR_PROPERTY p_Property);
    PROPERTY_DETAILS BinaryPropertyDetails(BINARY_PROPERTY p_Property);
    PROPERTY_DETAILS ProgramOptionDetails(PROGRAM_OPTION p_Property);
    STR_STR_STR_STR  FormatFieldHeaders(PROPERTY_DETAILS p_Details, string p_HeaderSuffix = "");
    LOGFILE_DETAILS  StandardLogFileDetails(const LOGFILE p_Logfile, const string p_FileSuffix);


    std::tuple<bool, LOGFILE> GetLogfileDescriptorKey(const std::string p_Value);
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
     * @return                                      String formatted as log file record
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

        string delimiter = DELIMITERValue.at(OPTIONS->LogfileDelimiter());                  // delimiter specified by OPTIONS (always for standard log files)

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
                default: ok = false; // unknown property tye - that's not ok...
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
                }
                else valueStr = "ERROR!";                                                   // error formatting value
            }

            logRecord += valueStr + delimiter;                                              // add value string to log record - with delimiter
        }

        // if we are constructing an SSE Switch file we add two pre-defined
        // columns to the end of the log record.  These are:
        //
        // ( i) the stellar type from which the star is switching
        // (ii) the stellar type to which the star is switching

        if (p_LogFile == LOGFILE::SSE_SWITCH_LOG) {
            string fmt = "%4.1d";                                                           // format - all integers here
            logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingFrom) + delimiter;      // switching from
            logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingTo) + delimiter;        // switching to
        }

        // if we are constructing a BSE Switch file record we add three pre-defined
        // columns to the end of the log record.  These are:
        //
        // (  i) the star switching - 1 = primary, 2 = secondary
        // ( ii) the stellar type from which the star is switching
        // (iii) the stellar type to which the star is switching

        if (p_LogFile == LOGFILE::BSE_SWITCH_LOG) {
            string fmt = "%4.1d";                                                           // format - all integers here
            int starSwitching = m_PrimarySwitching ? 1 : 2;                                 // primary (1) or secondary (2)
            logRecord += utils::vFormat(fmt.c_str(), starSwitching) + delimiter;            // star switching
            logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingFrom) + delimiter;      // switching from
            logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingTo) + delimiter;        // switching to
        }

        logRecord = logRecord.substr(0, logRecord.size()-1);                                // remove the last character - extraneous delimiter

        return logRecord;
    }

    /*
     * This variant of GetLogStandardRecord() is here because I can't readily figure out how to 
     * have an optional parameter in a template function if not providing that oprtional parameter 
     * doesn't give the compiler enough context to deduce the type of the template type.
     * 
     * Clear as mud?  Well, for the declaration of GetLogStandardRecord() above, template type T2
     * can't be decduced if the parameter p_SpecifiedPropertyValue is not supplied - even if it
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
     * void LogStandardRecord(const string   p_LogClass,
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
     * @param   [IN]    p_FileSuffix                String suffix to be added to the logfile name (optional, default = "")
     */
    template <class T>
    void LogStandardRecord(const string   p_LogClass,
                           const int      p_LogLevel,
                           const LOGFILE  p_LogFile,
                           const T* const p_Star,
                           const string   p_LogRecord,
                           const string   p_FileSuffix = "") {

        LOGFILE_DETAILS fileDetails;                                                                                                    // file details

        fileDetails = StandardLogFileDetails(p_LogFile, p_FileSuffix);                                                                  // get record details - open file (if necessary)
        if (get<0>(fileDetails) >= 0) {                                                                                                 // file is open

            int    fileId    = get<0>(fileDetails);                                                                                     // log file id

            string logRecord = "";                                                                                                      // the record to be written to the log file

            if (p_LogRecord.empty()) {                                                                                                  // logfile record passed in is empty
                                                                                                                                        // construct log record from current data
                ANY_PROPERTY_VECTOR properties = get<2>(fileDetails);                                                                   // vector of properties to be printed

                // get and format values for printing
                bool                 ok;                                                                                                // flag to indicate property value retrieved ok (or not)
                COMPAS_VARIABLE_TYPE value;                                                                                             // property value
                string               valueStr;                                                                                          // string for formatted value

                int index = 0;
                for (auto &property : properties) {                                                                                     // for each property to be included in the log record

                    std::tie(ok, value) = p_Star->PropertyValue(property);                                                              // get property flag and value
                    if (ok) {                                                                                                           // have valid property value
                        boost::variant<string> fmtStr(get<3>(fileDetails)[index++]);                                                    // format string

                        valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);                                           // format value
                    }
                    else valueStr = "ERROR!";                                                                                           // error formatting value

                    logRecord += valueStr + m_Logfiles[fileId].delimiter;                                                               // add value string to log record - with delimiter
                }

                // if we are writing to the SSE Switch file we add two pre-defined columns
                // to the end of the log record.  These are:
                //
                // ( i) the stellar type from which the star is switching
                // (ii) the stellar type to which the star is switching

                if (p_LogFile == LOGFILE::SSE_SWITCH_LOG) {
                    string fmt = "%4.1d";                                                                                               // format - all integers here
                    logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingFrom) + m_Logfiles[fileId].delimiter;                       // switching from
                    logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingTo) + m_Logfiles[fileId].delimiter;                         // switching to
                }

                // if we are writing to the BSE Switch file we add three pre-defined columns
                // to the end of the log record.  These are:
                //
                // (  i) the star switching - 1 = primary, 2 = secondary
                // ( ii) the stellar type from which the star is switching
                // (iii) the stellar type to which the star is switching

                if (p_LogFile == LOGFILE::BSE_SWITCH_LOG) {
                    string fmt = "%4.1d";                                                                                               // format - all integers here
                    int starSwitching = m_PrimarySwitching ? 1 : 2;                                                                     // primary (1) or secondary (2)
                    logRecord += utils::vFormat(fmt.c_str(), starSwitching) + m_Logfiles[fileId].delimiter;                             // star switching
                    logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingFrom) + m_Logfiles[fileId].delimiter;                       // switching from
                    logRecord += utils::vFormat(fmt.c_str(), m_TypeSwitchingTo) + m_Logfiles[fileId].delimiter;                         // switching to
                }

                logRecord = logRecord.substr(0, logRecord.size()-1);                                                                    // remove the last character - extraneous delimiter
            }
            else {                                                                                                                      // logfile record passed in is not empty
                logRecord = p_LogRecord;                                                                                                // use logfile record passed in
            }

            if (!Put_(fileId, logRecord)) DBG_WARN(ERR_MSG(ERROR::FILE_WRITE_ERROR) + ": Log Record");                                  // write the record - show warning if failure
        }
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
                 const string              p_Delimiter);

    void   Stop(std::tuple<int, int> p_ObjectStats = std::make_tuple(0, 0));

    bool   Enabled() const { return m_Enabled; }

    int    Open(const string p_LogFileName, const bool p_Append, const bool p_TimeStamp, const bool p_Label, const string p_Delimiter = "");
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

    std::tuple<ANY_PROPERTY_VECTOR, std::vector<string>> GetStandardLogFileRecordDetails(const LOGFILE p_Logfile);

    template <class T>
    void LogSSEParameters(const T* const p_Star, const int p_Id, const string p_Rec)            { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_PARAMETERS)), 0, LOGFILE::SSE_PARAMETERS, p_Star, p_Rec, "_" + std::to_string(abs(p_Id))); }
    template <class T>
    void LogSSESwitchLog(const T* const p_Star, const int p_Id, const string p_Rec)             { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SWITCH_LOG)), 0, LOGFILE::SSE_SWITCH_LOG, p_Star, p_Rec, "_" + std::to_string(abs(p_Id))); }
    template <class T>
    void LogSSESupernovaDetails(const T* const p_Star, const string p_Rec)                      { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVA)), 0, LOGFILE::SSE_SUPERNOVA, p_Star, p_Rec); }
    template <class T>
    void LogBSEParameters(const T* const p_Binary, const string p_Rec)                          { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SYSTEM_PARAMETERS)), 0, LOGFILE::BSE_SYSTEM_PARAMETERS, p_Binary, p_Rec); }
    template <class T>
    void LogBSEDoubleCompactObject(const T* const p_Binary, const string p_Rec)                 { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS)), 0, LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS, p_Binary, p_Rec); }
    template <class T>
    void LogBSECommonEnvelope(const T* const p_Binary, const string p_Rec)                      { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_COMMON_ENVELOPES)), 0, LOGFILE::BSE_COMMON_ENVELOPES, p_Binary, p_Rec); }
    template <class T>
    void LogBSERLOFParameters(const T* const p_Binary, const string p_Rec)                      { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_RLOF_PARAMETERS)), 0, LOGFILE::BSE_RLOF_PARAMETERS, p_Binary, p_Rec); }
    template <class T>
    void LogBSEBeBinary(const T* const p_Binary, const string p_Rec)                            { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_BE_BINARIES)), 0, LOGFILE::BSE_BE_BINARIES, p_Binary, p_Rec); }
    template <class T>
    void LogBSEDetailedOutput(const T* const p_Binary, const long int p_Id, const string p_Rec) { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT)), 0, LOGFILE::BSE_DETAILED_OUTPUT, p_Binary, p_Rec, "_" + std::to_string(abs(p_Id))); }
    template <class T>
    void LogBSEPulsarEvolutionParameters(const T* const p_Binary, const string p_Rec)           { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION)), 0, LOGFILE::BSE_PULSAR_EVOLUTION, p_Binary, p_Rec); }
    template <class T>
    void LogBSESupernovaDetails(const T* const p_Binary, const string p_Rec)                    { LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE)), 0, LOGFILE::BSE_SUPERNOVAE, p_Binary, p_Rec); }
    
    template <class T>
    void LogBSESwitchLog(const T* const p_Binary, const long int p_Id, const bool p_PrimarySwitching) {
        m_PrimarySwitching = p_PrimarySwitching;        
        LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG)), 0, LOGFILE::BSE_SWITCH_LOG, p_Binary, "", "_" + std::to_string(abs(p_Id)));
    }

    template <class T>
    void StashSSESupernovaDetails(const T* const p_Star, const STELLAR_TYPE p_StellarType) {

        // if we don't already have the SSE Supernova log record properties and format vector, get them
        // this will only need to be done one per run
        if (m_SSESupernova_LogRecordProperties.empty() || m_SSESupernova_LogRecordFmtVector.empty()) {
            std::tie(m_SSESupernova_LogRecordProperties, m_SSESupernova_LogRecordFmtVector) = LOGGING->GetStandardLogFileRecordDetails(LOGFILE::SSE_SUPERNOVA);
        }

        // get a formatted record with current data
        // this will replace any existing stashed record - no queue here
        m_SSESupernova_DelayedLogRecord = GetLogStandardRecord(LOGFILE::SSE_SUPERNOVA, 
                                                               p_Star, 
                                                               m_SSESupernova_LogRecordProperties, 
                                                               m_SSESupernova_LogRecordFmtVector, 
                                                               true, 
                                                               (ANY_STAR_PROPERTY)STAR_PROPERTY::STELLAR_TYPE, 
                                                               p_StellarType);
    }

    template <class T>
    void LogStashedSSESupernovaDetails(const T* const p_Star) { 

        // if the stashed SSE Supernova record is non-empty, print it, then clear it - otherwise do nothing
        if (!m_SSESupernova_DelayedLogRecord.empty()) {

            LogStandardRecord(get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVA)), 0, LOGFILE::SSE_SUPERNOVA, p_Star, m_SSESupernova_DelayedLogRecord);
            m_SSESupernova_DelayedLogRecord = "";
        }
    }


    bool CloseStandardFile(const LOGFILE p_LogFile, const bool p_Erase = true);
    bool CloseAllStandardFiles();
};

#endif // __Log_h__
