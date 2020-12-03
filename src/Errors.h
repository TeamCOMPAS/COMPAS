#ifndef __Errors_h__
#define __Errors_h__

#define ERRORS Errors::Instance()

#include <iostream>
#include <unordered_map>
#include <vector>


#include "constants.h"
#include "typedefs.h"

#include "Log.h"
#include "ErrorsMacros.h"
#include "Options.h"

class Log;


/*
 * Errors Singleton
 *
 * Provides global error handling
 *
 * Singletons and global variables are sometimes frowned-upon, but doing it this
 * way means the objects don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Errors {

private:

    Errors() { };
    Errors(Errors const&) = delete;
    Errors& operator = (Errors const&) = delete;

    static Errors *m_Instance;                                                      // pointer to the instance

    // dynamic error catalog
    // this unordered_map records the errors printed and thedetails (object type, object id, function name etc.)
    // there are 2 objectId vectors (non-stellar and stellar object), and similarly 2 funcName vectors
    // it's a bit of a hack, but splitting them out allows me to clean out the stellar vectors after each
    // star/binary is evolved, so the catalog isn't bloated by details for deleted objects
    COMPASUnorderedMap<
        ERROR,                              // error id
        std::tuple<                         // details for error id
            ERROR_SCOPE,                    //    scope
            bool,                           //    flag indicating if already printed
            std::vector<OBJECT_TYPE>,       //    object type
            std::vector<STELLAR_TYPE>,      //    stellar type
            std::vector<OBJECT_ID>,         //    vector of non-stellar (main, utils, etc) object ids
            std::vector<OBJECT_ID>,         //    vector of stellar ids
            std::vector<std::string>,       //    vector of function names for non-stellar object ids
            std::vector<std::string>,       //    vector of function names for stellar ids
            std::string                     //    error text
        >
    > m_ErrorCatalog = {};


public:

    ~Errors() { delete m_Instance; }

    static Errors* Instance();


    bool ShowIt(const std::string  p_Prefix,
                const ERROR        p_Error,
                const std::string  p_QualifyingStr = "",
                const OBJECT_ID    p_ObjectId      = 0,
                const OBJECT_TYPE  p_ObjectType    = OBJECT_TYPE::NONE,
                const STELLAR_TYPE p_StellarType   = STELLAR_TYPE::NONE,
                const char*        p_FuncName      = "Not provided");

    void Clean();

    size_t CatalogSize() { return m_ErrorCatalog.size(); }
};

#endif // __Errors_h_
