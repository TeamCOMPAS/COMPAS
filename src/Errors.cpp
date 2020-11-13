/*
 * Class Errors
 *
 * Provides global error handling functionality.
 *
 * Brief description (full documentation comes later...):
 *
 * Errors are defined in the error catalog in constants.h (ERROR_CATALOG).  One day we
 * might move the catalog to a file so it can be changed without changing the code, or
 * even have multiple catalogs provided for internationalisation - but not now.
 *
 * Errors defined in the error catalog have a scope and message text.  The scope is used
 * to determine when/if an error should be printed.
 *
 * The current values for scope are:
 *
 *    NEVER                : The error will not be printed.
 *    ALWAYS               : The error will always be printed.
 *    FIRST                : The error will be printed only on the first time it is encountered anywhere in the program
 *    FIRST_IN_OBJECT_TYPE : The error will be printed only on the first time it is encountered anywhere in objects of
 *                           the same type (e.g. Binary Star objects)
 *    FIRST_IN_STELLAR_TYPE: The error will be printed only on the first time it is encountered anywhere in objects of
 *                           the same stellar type (e.g. HeWD Star obejcts)
 *    FIRST_IN_OBJECT_ID   : The error will be printed only on the first time it is encountered anywhere in an object
 *    FIRST_IN_FUNCTION    : The error will be printed only on the first time it is encountered anywhere in a function of
 *                           an object (i.e. will print twice if encountered in the same function name in different objects)
 *
 *
 * The Errors class provides methods to print both warnings and errors - essentially the same thing,
 * but warning messages are prefixed with "WARNING:", whereas error messages are prefixed with "ERROR:".
 *
 * Errors and warnings are printed by using the macros defined in ErrorsMacros.h.  They are:
 *
 * Error macros:
 *
 *    SHOW_ERROR(error_number)                       : prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog)
 *    SHOW_ERROR(error_number, error_string)         : prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string"
 *    SHOW_ERROR_IF(cond, error_number)              : prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog) if "cond" is TRUE
 *    SHOW_ERROR_IF(cond, error_number, error_string): prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string", if "cond" is TRUE
 *
 * Warning macros:
 *
 *    SHOW_WARN(error_number)                        : prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog)
 *    SHOW_WARN(error_number, error_string)          : prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string"
 *    SHOW_WARN_IF(cond, error_number)               : prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog) if "cond" is TRUE
 *    SHOW_WARN_IF(cond, error_number, error_string) : prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string", if "cond" is TRUE
 *
 * Error and warning message always contain:
 *
 *    The object id of the calling object
 *    The object type of the calling object
 *    The stellar type of the calling object (will be "NONE" if the calling object is not a star-type object)
 *    The function name of the calling function
 *
 * Notes:
 *
 * Any object that uses the Error class (i.e. the SHOW_* macros) *must* expose the following functions:
 *
 *    OBJECT_ID           ObjectId() const                            { return m_ObjectId; }
 *    OBJECT_TYPE         ObjectType() const                          { return m_ObjectType; }
 *    STELLAR_TYPE        StellarType() const                         { return m_StellarType; }
 *
 * These functions are called by the SHOW_* macros.
 * If any of the functions are not applicable to the object, then they must return "*::NONE (all objects should implement ObjectId() correctly)
 */


#include "Errors.h"

Errors* Errors::m_Instance = nullptr;


Errors* Errors::Instance() {
    if (!m_Instance) {
        m_Instance = new Errors();
    }
    return m_Instance;
}


/*
 * Prints the error/warning message if appropriate
 *
 * Checks the scope of the error number, and if necessary (i.e. the scope is not NEVER and not ALWAYS),
 * checks whether the error meets the criteria for being printed.  If it does, this function retreives
 * the associated error message from the catalog, formats the error or warning string, and prints the
 * message.
 *
 *
 * bool ShowIt(const std::string  p_Prefix,
 *             const ERROR        p_Error,
 *             const std::string  p_QualifyingStr,
 *             const OBJECT_ID    p_ObjectId,
 *             const OBJECT_TYPE  p_ObjectType,
 *             const STELLAR_TYPE p_StellarType,
 *             const char*        p_FuncName)
 *
 * @param   [IN]    p_Prefix                    The prefix string (will be "ERROR: " or "WARNING: "
 * @param   [IN]    p_Error                     The error number (actually id)
 * @param   [IN]    p_QualifyingStr             String to be appended to the error/warning string (can be empty)
 * @param   [IN]    p_ObjectId                  The object id of the calling object
 * @param   [IN]    p_ObjectType                The object type of the calling object
 * @param   [IN]    p_StellarType               The stellar type of the calling object (will be "NONE" if not a star-type object)
 * @param   [IN]    p_FuncName                  The nmae of the calling function
 * @return                                      Boolean indicating if the error/warning was printed: true = yes, false = no
 */
bool Errors::ShowIt(const std::string  p_Prefix,
                    const ERROR        p_Error,
                    const std::string  p_QualifyingStr,
                    const OBJECT_ID    p_ObjectId,
                    const OBJECT_TYPE  p_ObjectType,
                    const STELLAR_TYPE p_StellarType,
                    const char*        p_FuncName) {

    if (p_Error == ERROR::NONE) return false;                                                                                                   // nothing to do if no error

    if (p_Prefix == WARNING_PREFIX && !OPTIONS->EnableWarnings()) return false;                                                                      // do nothing

    bool print = false;                                                                                                                         // default - don't print

    COMPASUnorderedMap<
        ERROR,                                                                                                                                  // error id
        std::tuple<                                                                                                                             // details for error id
            ERROR_SCOPE,                                                                                                                        //    scope
            bool,                                                                                                                               //    flag indicating if already printed
            std::vector<OBJECT_TYPE>,                                                                                                           //    object type
            std::vector<STELLAR_TYPE>,                                                                                                          //    stellar type
            std::vector<OBJECT_ID>,                                                                                                             //    vector of non-stellar (main, utils, etc) object ids
            std::vector<OBJECT_ID>,                                                                                                             //    vector of stellar ids
            std::vector<std::string>,                                                                                                           //    vector of function names for non-stellar object ids
            std::vector<std::string>,                                                                                                           //    vector of function names for stellar ids
            std::string                                                                                                                         //    error text
        >
    > ::iterator iter;

	iter = m_ErrorCatalog.find(p_Error);                                                                                                        // look for error in dynamic catalog
	if (iter == m_ErrorCatalog.end()) {                                                                                                         // found?
        COMPASUnorderedMap<ERROR, std::tuple<ERROR_SCOPE, std::string>>::const_iterator staticIter;
        staticIter = ERROR_CATALOG.find(p_Error);                                                                                               // look for error in static catalog
        if (staticIter != ERROR_CATALOG.end()) {                                                                                                // found
            m_ErrorCatalog[p_Error] = { std::get<0>(staticIter->second), false, {}, {}, {}, {}, {}, {}, std::get<1>(staticIter->second) };      // yes - put new entry in dynamic catalog
        }
    }

	iter = m_ErrorCatalog.find(p_Error);                                                                                                        // look for error in dynamic catalog
	if (iter != m_ErrorCatalog.end()) {                                                                                                         // found?
                                                                                                                                                // yes
        std::string funcName(p_FuncName);                                                                                                       // convert p_FuncName to string type
        std::string text;                                                                                                                       // error text to print

        ERROR_SCOPE scope = std::get<0>(iter->second);                                                                                          // scope of error

        if (scope != ERROR_SCOPE::NEVER) {                                                                                                      // nothing to do if scope is NEVER

            bool                      already             = std::get<1>(iter->second);                                                          // already printed once?
            std::vector<OBJECT_TYPE>  objectTypes         = std::get<2>(iter->second);                                                          // object types from which error already printed
            std::vector<STELLAR_TYPE> stellarTypes        = std::get<3>(iter->second);                                                          // stellar types from which error already printed
            std::vector<OBJECT_ID>    nonStellarObjectIds = std::get<4>(iter->second);                                                          // non-stellar object ids from which error already printed
            std::vector<OBJECT_ID>    stellarObjectIds    = std::get<5>(iter->second);                                                          // stellar object ids from which error already printed
            std::vector<std::string>  nonStellarFuncs     = std::get<6>(iter->second);                                                          // non-stellar functions from which error already printed
            std::vector<std::string>  stellarFuncs        = std::get<7>(iter->second);                                                          // stellar functions from which error already printed
                                      text                = std::get<8>(iter->second);                                                          // error text to print

            switch (scope) {
                case ERROR_SCOPE::ALWAYS:               // scope is ALWAYS - no need to check anything else
                    print = true;                                                                                                               // always print it - no need to update entry
                    break;

                case ERROR_SCOPE::FIRST:                // scope is FIRST - only need to check if already printed
                    print = !already;                                                                                                           // print it only if not already printed

                    if (print) {                                                                                                                // if will print, then...
                        iter->second = std::make_tuple(scope, true, objectTypes, stellarTypes, nonStellarObjectIds, stellarObjectIds, nonStellarFuncs, stellarFuncs, text); // update catalog
                    }
                    break;

                case ERROR_SCOPE::FIRST_IN_OBJECT_TYPE: // scope is FIRST in this type of object - need to check if already printed for this object type
                    if (p_ObjectType != OBJECT_TYPE::NONE) {                                                                                    // valid object type?
                        print = (std::find(objectTypes.begin(), objectTypes.end(), p_ObjectType) == objectTypes.end());                         // print it only if not already printed for this object type

                        if (print) {                                                                                                            // if will print, then...
                            objectTypes.push_back(p_ObjectType);                                                                                // add 'p_ObjectType' to 'objectTypes'
                            iter->second = std::make_tuple(scope, already, objectTypes, stellarTypes, nonStellarObjectIds, stellarObjectIds, nonStellarFuncs, stellarFuncs, text); // update catalog
                        }
                    }
                    break;

                case ERROR_SCOPE::FIRST_IN_STELLAR_TYPE: // scope is FIRST in this type of star - need to check if already printed for this stellar type
                    if (p_StellarType != STELLAR_TYPE::NONE) {                                                                                  // valid stellar type?
                        print = (std::find(stellarTypes.begin(), stellarTypes.end(), p_StellarType) == stellarTypes.end());                     // print it only if not already printed for this stellar type

                        if (print) {                                                                                                            // if will print, then...
                            stellarTypes.push_back(p_StellarType);                                                                              // add 'p_StellarType' to 'stellarTypes'
                            iter->second = std::make_tuple(scope, already, objectTypes, stellarTypes, nonStellarObjectIds, stellarObjectIds, nonStellarFuncs, stellarFuncs, text); // update catalog
                        }
                    }
                    break;

                case ERROR_SCOPE::FIRST_IN_OBJECT_ID:   // scope is FIRST in this object - need to check if already printed for this object
                    if (p_ObjectId >= 0) {                                                                                                      // valid object id (0 = main())?

                        print = (std::find(stellarObjectIds.begin(), stellarObjectIds.end(), p_ObjectId) == stellarObjectIds.end());            // look for stellar objectId first
                        if (print) {                                                                                                            // not found, now look for ...
                            print = (std::find(nonStellarObjectIds.begin(), nonStellarObjectIds.end(), p_ObjectId) == nonStellarObjectIds.end()); // ... non-stellar objectId
                        }

                        if (print) {                                                                                                            // if will print, then...
                            if (p_ObjectType == OBJECT_TYPE::MAIN || p_ObjectType == OBJECT_TYPE::UTILS) {                                      // add 'p_ObjectId' and 'funcName' to relevant vectors
                                nonStellarObjectIds.push_back(p_ObjectId);                                                                      // non-stellar objectId
                                nonStellarFuncs.push_back(funcName);                                                                            // non-stellar funcName (required for ERROR_SCOPE::FIRST_IN_FUNCTION)
                            }
                            else {
                                stellarObjectIds.push_back(p_ObjectId);                                                                         // stellar objectId
                                stellarFuncs.push_back(funcName);                                                                               // stellar funcName (required for ERROR_SCOPE::FIRST_IN_FUNCTION)
                            }
                            iter->second = std::make_tuple(scope, already, objectTypes, stellarTypes, nonStellarObjectIds, stellarObjectIds, nonStellarFuncs, stellarFuncs, text); // update catalog
                        }
                    }
                    break;

                case ERROR_SCOPE::FIRST_IN_FUNCTION:    // scope is FIRST in FUNCTION of this object - need to check if already printed for this function for this object
                    if (p_ObjectId >= 0) {                                                                                                      // valid object id (0 = main())?

                        if (!funcName.empty()) {                                                                                                // blank function name?
                                                                                                                                                // no - proceed
                            bool foundFunc = false;                                                                                             // found func for objectId?

                            std::vector<OBJECT_ID>::iterator thisIter = std::find(stellarObjectIds.begin(), stellarObjectIds.end(), p_ObjectId); // look for stellar objectId first
                            if (thisIter != stellarObjectIds.end()) {                                                                           // stellar objectId?
                                                                                                                                                // yes
                                while (thisIter != stellarObjectIds.end() && !foundFunc) {                                                      // while have objectId and not funcName
                                    int thisIndex = std::distance(stellarObjectIds.begin(), thisIter);                                          // get vector index for objectId located
                                    if (stellarFuncs[thisIndex] != funcName) {                                                              // have funcName?
                                        thisIter = std::find(++thisIter, stellarObjectIds.end(), p_ObjectId);                                   // no - find next entry for objectId
                                    }
                                    else foundFunc = true;                                                                                      // yes - already printed
                                }
                            }
                            else {                                                                                                              // not stellar objectId - look for non-stellar objectId

                                std::vector<OBJECT_ID>::iterator thisIter = std::find(nonStellarObjectIds.begin(), nonStellarObjectIds.end(), p_ObjectId); // look for stellar objectId first
                                while (thisIter != nonStellarObjectIds.end() && !foundFunc) {                                                   // while have objectId and not funcName
                                    int thisIndex = std::distance(nonStellarObjectIds.begin(), thisIter);                                       // get vector index for objectId located
                                    if (nonStellarFuncs[thisIndex] != funcName) {                                                               // have funcName?
                                        thisIter = std::find(++thisIter, nonStellarObjectIds.end(), p_ObjectId);                                // no - find next entry for objectId
                                    }
                                    else foundFunc = true;                                                                                      // yes - already printed
                                }
                            }
                            print = !foundFunc;                                                                                                 // print it if not found


                            if (print) {                                                                                                        // if will print, then...
                                if (p_ObjectType == OBJECT_TYPE::MAIN || p_ObjectType == OBJECT_TYPE::UTILS) {                                  // add 'p_ObjectId' and 'funcNmae' to relevant vectors
                                    nonStellarObjectIds.push_back(p_ObjectId);                                                                  // non-stellar objectId
                                    nonStellarFuncs.push_back(funcName);                                                                        // non-stellar funcName (required for ERROR_SCOPE::FIRST_IN_FUNCTION)
                                }
                                else {
                                    stellarObjectIds.push_back(p_ObjectId);                                                                     // stellar objectId
                                    stellarFuncs.push_back(funcName);                                                                           // stellar funcName (required for ERROR_SCOPE::FIRST_IN_FUNCTION)
                                }
                                iter->second = std::make_tuple(scope, already, objectTypes, stellarTypes, nonStellarObjectIds, stellarObjectIds, nonStellarFuncs, stellarFuncs, text); // update catalog
                            }
                        }
                    }
                    break;

                default:                                                                                                                        // unknown scope
                    print = false;                                                                                                              // don't print it
            }
        }

        if (print) {

            std::string       objectTypeStr = (p_ObjectType == OBJECT_TYPE::BASE_STAR) ? STELLAR_TYPE_LABEL.at(p_StellarType) : OBJECT_TYPE_LABEL.at(p_ObjectType);
            std::string       funcNameStr   = funcName.empty() ? "Unknown!" : funcName;
            std::string       qualifyingStr = p_QualifyingStr.empty() ? "" : "; " + p_QualifyingStr;
            std::stringstream objectIdStr;

            objectIdStr << ", (ObjectId = " << p_ObjectId << ")";

            LOGGING->Error(p_Prefix + " in " + objectTypeStr + "::" + funcNameStr + ": " + text + qualifyingStr + objectIdStr.str());
        }
	}

	return print;
}


/*
 * Removes all objectId entries (and associated funcNames) from the Error Catalog
 * This is so entries for deleted objects don't bloat the error catalog
 *
 *
 * void Clean()
 */
void Errors::Clean() {

    if (m_ErrorCatalog.size() == 0) return;                                                                                                                             // nothing to do

    for (auto catalogIter : m_ErrorCatalog) {                                                                                                                           // for each entry in the error catalog
        std::vector<OBJECT_ID> stellarObjectIds = std::get<5>(catalogIter.second);                                                                                      // stellar object ids
        if (stellarObjectIds.size() == 0) continue;                                                                                                                     // no stellar objectIds for this error

        ERROR_SCOPE               scope               = std::get<0>(catalogIter.second);                                                                                // scope of error
        bool                      already             = std::get<1>(catalogIter.second);                                                                                // already
        std::vector<OBJECT_TYPE>  objectTypes         = std::get<2>(catalogIter.second);                                                                                // object types
        std::vector<STELLAR_TYPE> stellarTypes        = std::get<3>(catalogIter.second);                                                                                // stellar types
        std::vector<OBJECT_ID>    nonStellarObjectIds = std::get<4>(catalogIter.second);                                                                                // object ids
        std::vector<std::string>  nonStellarFuncs     = std::get<6>(catalogIter.second);                                                                                // functions
        std::vector<std::string>  stellarFuncs        = std::get<7>(catalogIter.second);                                                                                // functions
        std::string               text                = std::get<8>(catalogIter.second);                                                                                // error text

        stellarObjectIds.clear();                                                                                                                                       // clear stellar objectIds vector
        stellarFuncs.clear();                                                                                                                                           // clear stellar funcs vector
        catalogIter.second = std::make_tuple(scope, already, objectTypes, stellarTypes, nonStellarObjectIds, stellarObjectIds, nonStellarFuncs, stellarFuncs, text);    // update catalog
    }
}
