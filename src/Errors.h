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

    Errors() {
        // populate error catalog
        for (auto const &it: ERROR_CATALOG) {                                       // iterate over ERROR_CATALOG from constants.h
            m_ErrorCatalog[it.first] = { std::get<0>(it.second), false, {}, {}, {}, {}, std::get<1>(it.second) };
        }
    };
    Errors(Errors const&) = delete;
    Errors& operator = (Errors const&) = delete;

    static Errors *m_Instance;                                                      // pointer to the instance

    std::unordered_map<ERROR, std::tuple<ERROR_SCOPE, bool, std::vector<OBJECT_TYPE>, std::vector<OBJECT_ID>, std::vector<STELLAR_TYPE>, std::vector<std::string>, std::string>> m_ErrorCatalog = {};


public:

    ~Errors() { delete m_Instance; }

    static Errors* Instance();


    bool ShowIt(const std::string  p_Prefix,
                const ERROR        p_Error,
                const std::string  p_QualifyingStr = "",
                const OBJECT_ID    p_ObjectId      = 0,
                const OBJECT_TYPE  p_ObjectType    = OBJECT_TYPE::NONE,
                const STELLAR_TYPE p_StellarType   = STELLAR_TYPE::NONE,
                const char*        p_FuncName      = __builtin_FUNCTION());

};

#endif // __Errors_h_
