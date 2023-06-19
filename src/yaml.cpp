#include <iostream>
#include <stdarg.h>
#include <fstream>
#include <string>
#include <iomanip>

#include "Options.h"
#include "yaml.h"
#include "utils.h"


/*
 * YAML file creation
 *
 * Code to create a new YAML file based on current COMPAS defaults and any user-input option values.
 * See yaml.h for description of the YAML template.
 * 
 * There is no doubt that some of the code here is not the most efficient way of doing this, but realistically
 * this will be run once or twice every not very often - and it runs in a fraction of a second - so I'm not going
 * to engage in premature optimisation...
 */

namespace yaml {

    /*
     * Read content from YAML template file
     *
     * Read YAML template from file passed in p_YAMLtemplateName.
     * Replace content of default template 'yamlTemplate' declared in yaml.h if read ok.
     * 
     * 
     * int ReadYAMLtemplate(const std::string p_YAMLtemplateName)
     * 
     * @param   [IN]    p_YAMLtemplateName          Filename to be read - should be fully qualified
     * @return                                      Integer result:
     *                                                -1 indicates IO error
     *                                                 0 indicates non-existent file
     *                                                 1 indicates file exists but contains no content
     *                                                 2 indicates contents read ok (and replace the default template declared in yaml.h)
     */
    int ReadYAMLtemplate(const std::string p_YAMLtemplateName) {
        int result = -1;                                                                                        // result - initially error

        if (utils::FileExists(p_YAMLtemplateName)) {                                                            // template file exists?

            std::ifstream yamlTemplateFile(p_YAMLtemplateName);                                                 // yes - open the file
            if (yamlTemplateFile.is_open()) {                                                                   // open ok?
 
                std::vector<std::string> yamlTemplateContent;                                                   // yes - file content

                std::string rec;                                                                                // record read from file
                while (std::getline(yamlTemplateFile, rec)) {                                                   // for all records in file - get record from template file
                    yamlTemplateContent.push_back(rec);                                                         // add record to content vector
                }

                result = 1;                                                                                     // template file read ok - no content
                if (yamlTemplateContent.size() > 0) {                                                           // template file has content?
                    yamlTemplate = yamlTemplateContent;                                                         // yes - replace yamlTemplate content
                    result = 2;                                                                                 // template file read ok - content ok
                }
            }
        }
        else result = 0;                                                                                        // no - YAML template file does not exist

        return result;     
    }


    /*
     * Write YAML content to YAML file
     *
     * Writes YAML file content passed in p_YAMLcontent to file passed in p_YAMLname.
     * Will prompt user for overwrite if file already exists.
     * 
     * 
     * int WriteYAMLfile(const std::string p_YAMLname, const std::vector<std::string> p_YAMLcontent)
     * 
     * @param   [IN]    p_YAMLname                  Filename to be written - should be fully qualified
     * @param   [IN]    p_YAMLcontent               YAML content to be written - vector of strings
     * @return                                      Integer result:
     *                                                -1 indicates IO error
     *                                                 0 indicates file already exists - no overwrite
     *                                                 1 indicates success - file written (possibly overwritten) ok
     */
    int WriteYAMLfile(const std::string p_YAMLname, const std::vector<std::string> p_YAMLcontent) {
        int result = -1;                                                                                        // result - initially error

        bool doWrite = true;                                                                                    // should write (overwrite)?
        if (utils::FileExists(p_YAMLname)) {                                                                    // file exists?
                                                                                                                // yes - ask user whether to overwrite
            std::string prompt = "YAML file '" + p_YAMLname + "' already exists - overwrite (Y/N)?";            // set prompt string for overwrite
            std::string response = "";                                                                          // user response

            do {
                std::cout << prompt << std::flush;                                                              // prompt use for input
                std::string ch = "";                                                                            // user response character
                std::string lastCh = "";                                                                        // last user response character (will always be newline)
                size_t charsInput = 0;                                                                          // number of characters input (including newline)
                do {
                    lastCh = ch;                                                                                // set last response character
                    ch = std::cin.get();                                                                        // get response character (one character at a time)
                    charsInput++;                                                                               // count of characters input
                } while (ch != "\n");                                                                           // while not newline
                ch = charsInput == 2 ? lastCh : "";                                                             // we want just 1 character + newline - character prior to newline          
                if (!ch.empty()) response = utils::ToLower(ch);                                                 // set user response
            } while (response != "y" && response != "n");                                                       // while not valid response
            doWrite = response == "y";                                                                          // set write flag
        }

        result = 0;                                                                                             // file not written (no overwrite) result
        if (doWrite) {                                                                                          // should write?
            result = 1;                                                                                         // yes - set file written result
            std::ofstream yamlFile(p_YAMLname);                                                                 // open the file
            if (yamlFile.is_open()) {                                                                           // open ok?
                for (std::size_t i = 0; i < p_YAMLcontent.size(); i++) {                                        // yes - for each record
                    if (!yamlFile.write(p_YAMLcontent[i].substr(3).c_str(), p_YAMLcontent[i].length() - 3)) {   // remove the marker ("HDR" / "OPT") and write
                        result = false;                                                                         // failed - set result ...
                        break;                                                                                  // ... and stop
                    }
                }
            }
            else result = -1;                                                                                   // no - not open ok - set IO error
        }
        return result;
    }


    /*
     * Construct YAML file content and write to file passed in p_YAMLfilename
     *
     * Constructs YAML file content based on YAML template read from file passed in
     * p_YAMLtemplateName.  If no filename passed in p_YAMLtemplateName, the default
     * template from yaml.h is used.
     * 
     * Constructed content is written to file passed in p_YAMLfilename.
     * 
     * 
     * void MakeYAMLfile(const std::string p_YAMLfilename, const std::string p_YAMLtemplateName)
     * 
     * @param   [IN]    p_YAMLfilename          Filename to be written - should be fully qualified
     * @param   [IN]    p_YAMLtemplateName      Template filename to be read - should be fully qualified
     */

    void MakeYAMLfile(const std::string p_YAMLfilename, const std::string p_YAMLtemplateName) {
    // following macro is for convenience and readability - undefined at end of function
    // sets strings to be written to the new YAML file
    #define SET_STRINGS(preamble, option, value, allowed, dfault, comment) {    \
        yamlPreambleStrings.push_back(preamble);                                \
        yamlOptionStrings.push_back(option);                                    \
        yamlValueStrings.push_back(value);                                      \
        yamlAllowedStrings.push_back(allowed);                                  \
        yamlDefaultStrings.push_back(dfault);                                   \
        yamlCommentStrings.push_back(comment); }


        std::vector<std::string> yamlPreambleStrings;                                                                       // new YAML file content - preamble strings
        std::vector<std::string> yamlOptionStrings;                                                                         // new YAML file content - option name strings
        std::vector<std::string> yamlValueStrings;                                                                          // new YAML file content - option value strings
        std::vector<std::string> yamlAllowedStrings;                                                                        // new YAML file content - allowed values strings
        std::vector<std::string> yamlDefaultStrings;                                                                        // new YAML file content - default values strings
        std::vector<std::string> yamlCommentStrings;                                                                        // new YAML file content - comment strings


        // process the option
        // - format strings as required and populate the yaml vectors (declared above)
        auto ProcessOption = [&] (const std::string              p_OptionStr,                                               // option name string
                                  const std::string              p_ValueStr,                                                // option value string
                                  const std::vector<std::string> p_AllowedStr,                                              // option allowed value strings
                                  const std::string              p_DefaultStr,                                              // option default string
                                  const std::string              p_CommentStr,                                              // comment string
                                  const bool                     p_UserSpecified,                                           // user specified option value?
                                  const TYPENAME                 p_DataType,                                                // short data type
                                  const std::string              p_TypeStr) {                                               // detailed data type

            // format option value/default string
            // fix floating point precision
            // fix bool for python (0/true, 1/false -> True, False)
            // add '' for string options if not already present
            // not all data types need to be formatted - some just pass through
            auto FormatString = [&] (const std::string p_Str,                                                               // string to be formatted
                                     const TYPENAME    p_shortType,                                                         // short data type
                                     const std::string p_detailedType) {                                                    // detailed data type
                std::stringstream ss;
                switch (p_shortType) {
                    case TYPENAME::FLOAT     : {       float v = std::stof(p_Str);          ss << std::fixed << std::setprecision(v < 100 ? 6 : 2) << v;     } break;
                    case TYPENAME::DOUBLE    : {      double v = std::stod(p_Str);          ss << std::fixed << std::setprecision(v < 100 ? 6 : 2) << v;     } break;
                    case TYPENAME::LONGDOUBLE: { long double v = std::stold(p_Str);         ss << std::fixed << std::setprecision(v < 100 ? 6 : 2) << v;     } break;
                    case TYPENAME::BOOL      : { std::string s = p_Str; s = utils::trim(s); ss << (s == "1" || utils::Equals(s, "true") ? "True" : "False"); } break;
                    case TYPENAME::STRING    : {
                            bool vecType = (p_detailedType == "VECTOR<STRING>");                                            // vector of strings? 
                            size_t len = p_Str.length();                                                                    // str length
                            if (len == 0) {                                                                                 // empty str?
                                ss << (vecType ? "{ }" : "''");                                                             // yes - just quotes/braces (for vector)
                            }
                            else {                                                                                          // no - str not empty
                                // add quotes if necessary
                                // the assumption here is that if the string already has opening and closing quotes
                                // then we're good - do nothing, but if it has only an opening quote or a closing quote,
                                // or neither, then we need to add both opening and closing quotes
                                // (same assumption for open and close brace for vector strings - if it has braces we leave it alone)

                                if (!vecType) {                                                                             // VECTOR<STRING> ?
                                    if (p_Str[0] != '\'' || p_Str[len - 1] != '\'') { ss << "'"; ss << p_Str;  ss << "'"; } // no - add quotes if necessary
                                    else ss << p_Str;                                                                       // has quotes - just pass through                                                                         
                                }
                                else {                                                                                      // yes - string vector - parse and format string
                                    std::string str = p_Str;                                                                // string to be formatted
                                    str = utils::trim(str);                                                                 // trim whitespace both ends
                                    if (str[0] != '{' || str[len - 1] != '}') {                                             // has open and close brace?
                                        ss << "{ ";                                                                         // no - add open brace
                                        size_t start = str[0] == '{' ? 1 : 0;                                               // start position for find
                                        size_t end   = 0;                                                                   // found position
                                        bool first   = true;                                                                // for comma placement
                                        while ((end = p_Str.find(' ', start)) != std::string::npos) {                       // find space
                                            std::string v = p_Str.substr(start, end - start);                               // got it - extract the substring
                                            if (v.length() > 0) {                                                           // skip empty substring
                                                if (!first) ss << ", ";                                                     // add comma if not first element
                                                if (v[0] != '\'' || v[len - 1] != '\'') { ss << "'"; ss << v; ss << "'"; }  // add quotes if necessary
                                                else ss << v;                                                               // has quotes - just pass through                                 
                                            }
                                            start = end + 1;                                                                // for next space
                                            first = false;                                                                  // not first element
                                        }
                                        std::string v = p_Str.substr(start, end - start);                                   // extract (possible) last substring
                                        if (v.length() > 0) {                                                               // anything there?
                                            if (!first) ss << ", ";                                                         // yes - add comma if not first element (shouldn't be, but just in case...)
                                            if (v[0] != '\'' || v[len - 1] != '\'') { ss << "'"; ss << v; ss << "'"; }      // add quotes if necessary
                                            else ss << v;                                                                   // has quotes - just pass through                                                                         
                                        }
                                        if (str[0] != '{' || str[len - 1] != '}') ss << " }";                               // add close brace if necessary
                                    }
                                    else ss << p_Str;                                                                       // has braces - just pass it through
                                }
                            }   
                        } break;
                    default: ss << p_Str;                                                                                   // just pass it through
                }
                return ss.str();                                                                                            // formatted string
            };

            std::string valueStr = FormatString(p_ValueStr, p_DataType, p_TypeStr);
            std::string defaultStr = FormatString(p_DefaultStr, p_DataType, p_TypeStr);

            std::string preamble = p_UserSpecified ? "     --" : "#    --";                                                 // preamble - all option records are commented, except options that are user specified

            if (p_AllowedStr.size() == 0) {                                                                                 // multiple option values?
                SET_STRINGS(preamble, p_OptionStr, valueStr, "", defaultStr, p_CommentStr);                                 // no - record for YAML file output - value set
            }
            else {                                                                                                          // yes - multiple option values
                // We list the allowed values in the order they are stored in the COMPASUnorderedMap
                // in constants.h.  We could instead list in alphabetical order, but the assumption
                // is that the order in constants.h was deliberate, so we'll maintain it - if that's
                // not true then we could just sort alphabetically here
                std::string allowedValuesStr = "Options: [";                                                                // records allowed option values
                for (size_t idx = 0; idx < p_AllowedStr.size(); idx++) {                                                    // for each allowed value
                    allowedValuesStr += p_AllowedStr[idx];                                                                  // show allowed value
                    if (idx < p_AllowedStr.size() - 1) allowedValuesStr += ",";                                             // add delimiter if necessary
                }
                allowedValuesStr += "]";                                                                                    // close bracket
                SET_STRINGS(preamble, p_OptionStr, valueStr, allowedValuesStr, defaultStr, p_CommentStr);                   // option record for YAML file output
            }
        };


        // begin MakeYAMLfile() body

        // get COMPAS option details
        std::vector<OptionDetailsT> optionDetails = OPTIONS->CmdLineOptionsDetails();

        if (optionDetails.empty()) {                                                                                        // have COMPAS option details?
            std::cerr << "Unable to access COMPAS option details - no YAML file created.\n";                                // no - announce error
        }
        else {                                                                                                              // yes - have options
                                                                                                                            // process options/create YAML file
            auto wallTime           = std::chrono::system_clock::now();                                                     // get wall time
            std::time_t timeNow     = std::chrono::system_clock::to_time_t(wallTime);                                       // current time and date ...
            std::string currentTime = std::string(std::ctime(&timeNow));                                                    // ... as a string ...
            currentTime             = utils::trim(currentTime);                                                             // ... trimmed of whitespace (including newline etc.)

            // read YAML template file if supplied

            if (p_YAMLtemplateName.length() > 0) {                                                                          // have YAML template filename?
                int readResult = ReadYAMLtemplate(p_YAMLtemplateName);                                                      // yes - read YAML template file
                std::string preamble = "*WARNING* File '" + p_YAMLtemplateName + "' ";                                      // preamble for warnings
                std::string warn     = "";                                                                                  // warning string
                if (readResult < 0) warn += "not read: IO error.  COMPAS default template will be used.";                   // IO error
                else if (readResult == 0) warn += "not read: file does not exist.  COMPAS default template will be used.";  // file not found
                else if (readResult == 1) warn += "is empty.  COMPAS default template will be used.";                       // no content in file
                if (warn.length() > 0) std::cout << preamble << warn << "\n";                                               // announce warning if necessary
            }

            // add default records (file headers)

            SET_STRINGS("##~!!~## COMPAS option values", "", "", "", "", "");                                               // header record for YAML file output

            std::string s = "##~!!~## File Created " + currentTime + " by COMPAS v" + VERSION_STRING;                       // time and version stamp
            if (s[23] == ' ') s[23] = '0';                                                                                  // add leading zero for day if necessary
            SET_STRINGS(s, "", "", "", "", "");                                                                             // time and version stamp record for YAML file output

            SET_STRINGS("##~!!~##", "", "", "", "", "");                                                                                    // Blank header line
            SET_STRINGS("##~!!~## The default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has", "", "", "", "", "");   // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output
            SET_STRINGS("##~!!~## all COMPAS option entries commented so that the COMPAS default value for the", "", "", "", "", "");       // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output
            SET_STRINGS("##~!!~## option is used by default. To use a value other than the COMPAS default value,", "", "", "", "", "");     // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output
            SET_STRINGS("##~!!~## users must uncomment the entry and change the option value to the desired value.", "", "", "", "", "");   // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output

            size_t numTemplateRecords = yamlTemplate.size();                                                                // number of records in the YAML template

            // create YAML content

            if (numTemplateRecords == 0) {                                                                                  // have YAML template?
                                                                                                                            // no...
                // we have no template - this should never happen (we've got a default template in the code,
                // but we need to deal with that being missing just in case) - so we'll just write the COMPAS
                // options alphabetically within the datatype categories (booleanChoices, numericalChoices, 
                // stringChoices, listChoices)

                std::cerr << "*WARNING* No YAML template - will write COMPAS options alphabetically.\n";                    // announce warning

                // this is not pretty - we're going to iterate over the COMPAS options 4 times - once
                // for each of the datatype categories mentioned above (see note about premature optimisation
                // at the top of this file)

                std::vector<std::string> categories { "booleanChoices", "numericalChoices", "stringChoices", "listChoices" }; // category names for headers

                for (size_t idx = 0; idx < categories.size(); idx++) {                                                      // for each category

                    std::string category = categories[idx];                                                                 // which category?
                    SET_STRINGS("", "", "", "", "", "");                                                                    // blank record for YAML file output
                    SET_STRINGS(category + ":", "", "", "", "", "");                                                        // category header for YAML file output

                    for (size_t optionIdx = 0; optionIdx < optionDetails.size(); optionIdx++) {                             // for each COMPAS option

                        if (optionDetails[optionIdx].dataType == TYPENAME::BOOL) {                                          // boolean option?
                            if (category != "booleanChoices") continue;                                                     // skip if this category is not booleanChoices 
                        }
                        else if (optionDetails[optionIdx].dataType == TYPENAME::STRING) {                                   // string, or vector of strings, option?
                            if (optionDetails[optionIdx].typeStr == "VECTOR<STRING>") {                                     // ... vector of strings
                                if (category != "listChoices") continue;                                                    // skip if this category is not listChoices 
                            }
                            else {                                                                                          // ... string
                                if (category != "stringChoices") continue;                                                  // skip if this category is not stringChoices
                            }
                        }
                        else {                                                                                              // assume numerical option
                             if (category != "numericalChoices") continue;                                                  // skip if this category is not numericalChoices
                        }

                        // process the option
                        ProcessOption(optionDetails[optionIdx].optionStr,                                                   // option name string
                                      optionDetails[optionIdx].valueStr,                                                    // option value string
                                      optionDetails[optionIdx].allowedStr,                                                  // option allowed value strings
                                      optionDetails[optionIdx].defaultStr,                                                  // option default string
                                      "",                                                                                   // no comment
                                      utils::Equals(optionDetails[optionIdx].sourceStr, "user_supplied"),                   // user specified option value?
                                      optionDetails[optionIdx].dataType,                                                    // option (short) data type  
                                      optionDetails[optionIdx].typeStr);                                                    // option (detailed) data type  
                    }
                }
            }
            else {                                                                                                          // yes - have YAML template
                                                                                                                            // use template to format YAML file
                for (size_t rec = 0; rec < numTemplateRecords; rec++) {                                                     // for each template record

                    std::string thisRec = yamlTemplate[rec];                                                                // template record ...
                    thisRec = utils::trim(thisRec);                                                                         // ... trimmed of whitespace

                    // need to allow for commented options in template, especially since the template
                    // might be an existing YAML file.  Since the template may have the hash character
                    // ("#") throughout, I'll refer to the hash at the start of the record (or as the
                    // first non-whitespace character) as the "leading hash" to differentiate it from
                    // other instances of the hash character in a template record
                    bool   optRec = false;                                                                                  // option template record?
                    size_t startPos = 0;                                                                                    // where to start looking for non-leading hash

                    if (thisRec.length() > 3) {                                                                             // record long enough to be option record?
                        if (thisRec.substr(0, 2) == "--") {                                                                 // yes - has leading "--"?
                            optRec = true;                                                                                  // yes - assume option record
                        }
                        else if (thisRec[0] == '#') {                                                                       // no - not leading "--": leading hash?
                            optRec = false;                                                                                 // for now...
                            if (thisRec.length() > 4) {                                                                     // enough characters to be option record?
                                std::string s = thisRec.substr(1);                                                          // yes... strip leading hash
                                s = utils::ltrim(s);                                                                        // trim whitespace from start
                                if (s.length() > 2 && s.substr(0, 2) == "--") {                                             // more than 2 characters, and option indicater?
                                    optRec = true;                                                                          // yes - looks like it might be an option record
                                    startPos = 1;                                                                           // skip over leading hash in later searches
                                }
                            }
                        }
                    }

                    // parse template record
                    if (!optRec) {                                                                                          // option record?
                        if (yamlTemplate[rec].length() < 8 || yamlTemplate[rec].substr(0, 8) != "##~!!~##") {               // no - skip COMPAS generated headers, and ...
                            SET_STRINGS(yamlTemplate[rec], "", "", "", "", "");                                             // ... output the record as is
                        }
                    }
                    else {                                                                                                  // yes - option record

                        // look for a comment on the option record and, if present, preserve it.
                        // I refer to the hash indiaction a comment on the option record as the
                        // "comment hash" to differentiate it from all the other hashes that might
                        // be there...  Because the template might be an existing YAML file, one or
                        // both of the (special) strings "# Default:" and "# Options:" may be present,
                        // and if so we just skip over them looking for any comment on the option record

                        size_t hashStart = startPos;                                                                        // where to start looking for the comment hash
                        size_t p;                                                                                           // position in string
                        if ((p = yamlTemplate[rec].find("# Default:", hashStart)) != std::string::npos) {                   // "# Default:" present?
                            if (p > hashStart) hashStart = p + 1;                                                           // yes - adjust starting position for comment hash search
                        }

                        if ((p = yamlTemplate[rec].find("# Options:", hashStart)) != std::string::npos) {                   // "# Options:" present?
                            if (p > hashStart) hashStart = p + 1;                                                           // yes - adjust starting position for comment hash search
                        }

                        std::string optionStr;                                                                              // string to hold option name and value
                        std::string commentStr;                                                                             // string to hold optional comment
                        size_t hashLoc = yamlTemplate[rec].find_first_of("#", hashStart);                                   // find the comment if there is one
                        if (hashLoc == std::string::npos) {                                                                 // found hash?
                            optionStr = yamlTemplate[rec].substr(startPos);                                                 // no - no comment present - option string is the record (from startPos)
                            commentStr = "";                                                                                // no comment
                        }
                        else {                                                                                              // yes - hash (comment present)
                            optionStr = yamlTemplate[rec].substr(startPos, hashLoc);                                        // option string is the record up to hash
                            commentStr = yamlTemplate[rec].substr(hashLoc);                                                 // comment string (includes "#")
                        }

                        optionStr = utils::trim(optionStr);                                                                 // trim whitespace from both ends of option string

                        // need to allow for option values in template, especially since the template might
                        // be an existing YAML file. Values in the template are not preserved - we only use
                        // the template to preserve header records, blank records, option groupings, and comments.
                        size_t colonLoc = optionStr.find_first_of(":");                                                     // find the first colon if there is one
                        if (colonLoc != std::string::npos) {                                                                // found colon?
                            optionStr = optionStr.substr(0, colonLoc);                                                      // option string is everything prior to the ":""
                        }

                        // Find the value for the option (COMPAS value, not any value present in the template)
                        // - note: could be user supplied
                        size_t pos = optionStr.find_first_not_of("\t -");                                                   // find option indicator
                        std::string templateOptionStr = optionStr.substr(pos);                                              // strip "--" from template option string

                        bool found = false;                                                                                 // matching COMPAS option string found?  Not initially...
                        for (std::size_t optionIdx = 0; optionIdx < optionDetails.size(); optionIdx++) {                    // for each COMPAS option

                            std::string compasOptionStr = optionDetails[optionIdx].optionStr;                               // COMPAS option name string
                            if (utils::Equals(optionDetails[optionIdx].sourceStr, "calculated")) continue;                  // ignore calculated options                                         
                            if (utils::Equals(compasOptionStr, templateOptionStr)) {                                        // match?

                                found = true;                                                                               // yes
                                ProcessOption(templateOptionStr,                                                            // option name string (as written in the YAML template)
                                              optionDetails[optionIdx].valueStr,                                            // option (COMPAS) value string
                                              optionDetails[optionIdx].allowedStr,                                          // option allowed value strings
                                              optionDetails[optionIdx].defaultStr,                                          // option default string
                                              commentStr,                                                                   // comment string (as written in the YAML template)
                                              utils::Equals(optionDetails[optionIdx].sourceStr, "user_supplied"),           // user specified option value?
                                              optionDetails[optionIdx].dataType,                                            // option (short) data type  
                                              optionDetails[optionIdx].typeStr);                                            // option (detailed) data type
                                break;                                                                                      // we're done
                            }
                        }

                        if (!found) {                                                                                       // matching COMPAS option string found?
                            std::cerr << "*WARNING* Option '" 
                                      << templateOptionStr 
                                      << "' in YAML template is not a valid COMPAS option - ignored.\n";                    // no - announce warning
                        }
                    }
                }
            }

            // add any COMPAS options not in the YAML template - alphabetically at end

            bool extra = false;                                                                                             // extra records?  Initially false
            for (std::size_t optionIdx = 0; optionIdx < optionDetails.size(); optionIdx++) {                                // for each COMPAS option

                std::string compasOptionStr = optionDetails[optionIdx].optionStr;                                           // option name string

                if (utils::Equals(optionDetails[optionIdx].sourceStr, "calculated")) continue;                              // ignore calculated options    
                if (utils::Equals(compasOptionStr, "help")) continue;                                                       // ignore 'help' option          
                if (utils::Equals(compasOptionStr, "version")) continue;                                                    // ignore 'version' option          
                if (utils::Equals(compasOptionStr, "create-YAML-file")) continue;                                           // ignore 'create-yaml-file' option          
                if (utils::Equals(compasOptionStr, "YAML-template")) continue;                                              // ignore 'YAML-template' option          

                bool match = false;
                for (size_t yamlIdx = 0; yamlIdx < yamlOptionStrings.size(); yamlIdx++) {                                   // for each YAML option name
                    if (utils::Equals(compasOptionStr, yamlOptionStrings[yamlIdx])) {                                       // match?
                        match = true;                                                                                       // yes - flag it
                        break;                                                                                              // no need to look further
                    }
                }

                if (!match) {                                                                                               // COMPAS option in YAML options?
                    if (!extra) {                                                                                           // no - first extra record?
                        extra = true;                                                                                       // yes - flag it
                        SET_STRINGS("\n\n### Additional COMPAS options not found in YAML template ###", "", "", "", "", "");// header record for YAML file output                     
                    }
                    ProcessOption(compasOptionStr,                                                                          // option name string (as written in the COMPAS OPTIONS class)
                                  optionDetails[optionIdx].valueStr,                                                        // option value string
                                  optionDetails[optionIdx].allowedStr,                                                      // option allowed value strings
                                  optionDetails[optionIdx].defaultStr,                                                      // option default string
                                  "",                                                                                       // no comment
                                  utils::Equals(optionDetails[optionIdx].sourceStr, "user supplied"),                       // user specified option value?
                                  optionDetails[optionIdx].dataType,                                                        // option (short) data type  
                                  optionDetails[optionIdx].typeStr);                                                        // option (detailed) data type
                }
            }

            // partially build new YAML file records
            // for now, just any header or blank records, and preamble, option name, and option value
            // for option records set maximum default str length and maximum allowed str length
            std::vector<std::string> yamlRecords;                                                                           // new YAML file records
            size_t maxOptionStrLen  = 0;                                                                                    // maximum length of option strings (to locate default string)
            size_t maxDefaultStrLen = 0;                                                                                    // maximum length of default strings (to locate allowed string)
            size_t maxAllowedStrLen = 0;                                                                                    // maximum length of allowed strings (to locate comment string)
            for (std::size_t yamlIdx = 0; yamlIdx < yamlPreambleStrings.size(); yamlIdx++) {                                // for each YAML option
                if (yamlOptionStrings[yamlIdx].empty()) {                                                                   // header or blank record?
                    yamlRecords.push_back("HDR" + yamlPreambleStrings[yamlIdx]);                                            // yes - marker and just preamble
                }
                else {                                                                                                      // no - not header or blank record
                    std::string s = "OPT" + yamlPreambleStrings[yamlIdx] + yamlOptionStrings[yamlIdx] + ": " + yamlValueStrings[yamlIdx]; // construct string including value
                    if (s.length() > maxOptionStrLen) maxOptionStrLen = s.length();                                         // get maximum length of YAML option strings
                    yamlRecords.push_back(s);                                                                               // partial record for YAML file output
                }
                if (yamlDefaultStrings[yamlIdx].length() > maxDefaultStrLen) maxDefaultStrLen = yamlDefaultStrings[yamlIdx].length(); // maximum length of default strings
                if (yamlAllowedStrings[yamlIdx].length() > maxAllowedStrLen) maxAllowedStrLen = yamlAllowedStrings[yamlIdx].length(); // maximum length of allowed strings
            }

            // add default values, allowed values, and template comment as necessary

            size_t defaultPos = maxOptionStrLen + 2;                                                                        // position of default string
            size_t allowedPos = defaultPos + maxDefaultStrLen + 13;                                                         // position of allowed strings
            size_t commentPos = allowedPos + maxAllowedStrLen + 4;                                                          // position of comment string
            for (size_t idx = 0; idx < yamlRecords.size(); idx++) {                                                         // for each YAML record
                std::string s = yamlRecords[idx];                                                                           // start with partially constructed record
                    
                if (yamlRecords[idx].substr(0, 3) == "OPT") {                                                               // option record?
                                                                                                                            // yes - header and blank records skipped here
                    size_t len = s.length();                                                                                // current length of record
                    for (size_t pos = len; pos < defaultPos; pos++) s += " ";                                               // pad with spaces up to option default field
                    s += "# Default: " + yamlDefaultStrings[idx];                                                           // yes add default values string

                    if (!yamlAllowedStrings[idx].empty()) {                                                                 // have allowed values?
                        size_t len = s.length();                                                                            // yes - current length of record
                        for (size_t pos = len; pos < allowedPos; pos++) s += " ";                                           // pad with spaces up to allowed strings field
                        s += "# " + yamlAllowedStrings[idx];                                                                // yes add allowed values string
                    }

                    if (!yamlCommentStrings[idx].empty()) {                                                                 // have comment?
                        size_t len = s.length();                                                                            // yes - current length of record
                        for (size_t pos = len; pos < commentPos; pos++) s += " ";                                           // pad with space ups to comment strings field
                        s += "  " + yamlCommentStrings[idx];                                                                // yes add comment string
                    }
                }
                yamlRecords[idx] = s + "\n";                                                                                // set YAML file record
            }

            // write new YAML file

            int writeResult = WriteYAMLfile(p_YAMLfilename, yamlRecords);                                                   // write records to file
            if (writeResult < 0) std::cerr << "*ERROR* File '" << p_YAMLfilename << "' not written: IO error.\n";           // announce IO error
            else if (writeResult == 0) std::cout << "File '" << p_YAMLfilename << "' not written - file already exists.\n"; // announce already exists
            else std::cout << "File '" << p_YAMLfilename << "' written.\n";                                                 // announce success                
        }
    }
    #undef SET_STRINGS
}
