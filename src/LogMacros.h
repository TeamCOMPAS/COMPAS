#ifndef __defines_h__
#define __defines_h__

#include <sstream>

//#define DEBUG // comment this line out, or #undef DEBUG, to build production executable (i.e. no DEBUG code)

#define DEBUG_WARNINGS  // comment this line out, or #undef DEBUG_WARNINGS, to build executable without WARNing statements

#define DEBUG_PERTURB // see discussion around FGB::PerturbLuminosityAndRadius()



// DEBUG, DEBUG_WARN and SAY macros

#define GET_MACRO(_0, _1, _2, _3, _4, _5, NAME, ...)                NAME


// Debug (to stdout, and to file if configured)

#ifdef DEBUG

    #define DBG_0()
    #define DBG_1(dbgStr)                                           { std::stringstream _ss; _ss << dbgStr; Log::Instance()->Debug("", 0, _ss.str()); }
    #define DBG_2(dbgLevel, dbgStr)                                 { std::stringstream _ss; _ss << dbgStr; Log::Instance()->Debug("", dbgLevel, _ss.str()); }
    #define DBG_3(dbgClass, dbgLevel, dbgStr)                       { std::stringstream _ss; _ss << dbgStr; Log::Instance()->Debug(dbgClass, dbgLevel, _ss.str()); }
    #define DBG_4()
    #define DBG_5()
    #define DBG(...)                                                GET_MACRO(_0, ##__VA_ARGS__, DBG_5, DBG_4, DBG_3, DBG_2, DBG_1, DBG_0) (__VA_ARGS__)

    #define DBG_ID_0()                                              { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Debug("", 0, _ss.str()); }
    #define DBG_ID_1(dbgStr)                                        { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << dbgStr; Log::Instance()->Debug("", 0, _ss.str()); }
    #define DBG_ID_2(dbgLevel, dbgStr)                              { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << dbgStr; Log::Instance()->Debug("", dbgLevel, _ss.str()); }
    #define DBG_ID_3(dbgClass, dbgLevel, dbgStr)                    { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << dbgStr; Log::Instance()->Debug(dbgClass, dbgLevel, _ss.str()); }
    #define DBG_ID_4()
    #define DBG_ID_5()
    #define DBG_ID(...)                                             GET_MACRO(_0, ##__VA_ARGS__, DBG_ID_5, DBG_ID_4, DBG_ID_3, DBG_ID_2, DBG_ID_1, DBG_ID_0)(__VA_ARGS__)

    #define DBG_IF_0()
    #define DBG_IF_1()
    #define DBG_IF_2(cond, dbgStr)                                  { if (cond) { std::stringstream _ss; _ss << dbgStr; Log::Instance()->Debug("", 0, _ss.str()); }}
    #define DBG_IF_3(cond, dbgLevel, dbgStr)                        { if (cond) { std::stringstream _ss; _ss << dbgStr; Log::Instance()->Debug("", dbgLevel, _ss.str()); }}
    #define DBG_IF_4(cond, dbgClass, dbgLevel, dbgStr)              { if (cond) { std::stringstream _ss; _ss << dbgStr; Log::Instance()->Debug(dbgClass, dbgLevel, _ss.str()); }}
    #define DBG_IF_5()
    #define DBG_IF(...)                                             GET_MACRO(_0, ##__VA_ARGS__, DBG_IF_5, DBG_IF_4, DBG_IF_3, DBG_IF_2, DBG_IF_1, DBG_IF_0)(__VA_ARGS__)

    #define DBG_ID_IF_0()
    #define DBG_ID_IF_1(cond)                                       { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Debug("", 0, _ss.str()); }}
    #define DBG_ID_IF_2(cond, dbgStr)                               { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << dbgStr; Log::Instance()->Debug("", 0, _ss.str()); }}
    #define DBG_ID_IF_3(cond, dbgLevel, dbgStr)                     { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << dbgStr; Log::Instance()->Debug("", dbgLevel, _ss.str()); }}
    #define DBG_ID_IF_4(cond, dbgClass, dbgLevel, dbgStr)           { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << dbgStr; Log::Instance()->Debug(dbgClass, dbgLevel, _ss.str()); }}
    #define DBG_ID_IF_5()
    #define DBG_ID_IF(...)                                          GET_MACRO(_0, ##__VA_ARGS__, DBG_ID_IF_5, DBG_ID_IF_4, DBG_ID_IF_3, DBG_ID_IF_2, DBG_ID_IF_1, DBG_ID_IF_0)(__VA_ARGS__)

    #define DBG_WAIT_0()
    #define DBG_WAIT_1(dbgStr)                                      { std::stringstream _ss; _ss << dbgStr; Log::Instance()->DebugWait("", 0, _ss.str()); }
    #define DBG_WAIT_2(dbgLevel, dbgStr)                            { std::stringstream _ss; _ss << dbgStr; Log::Instance()->DebugWait("", dbgLevel, _ss.str()); }
    #define DBG_WAIT_3(dbgClass, dbgLevel, dbgStr)                  { std::stringstream _ss; _ss << dbgStr; Log::Instance()->DebugWait(dbgClass, dbgLevel, _ss.str()); }
    #define DBG_WAIT_4()
    #define DBG_WAIT_5()
    #define DBG_WAIT(...)                                           GET_MACRO(_0, ##__VA_ARGS__, DBG_WAIT_5, DGB_WAIT_4, DBG_WAIT_3, DBG_WAIT_2, DBG_WAIT_1, DBG_WAIT_0)(__VA_ARGS__)

    #define DBG_WAIT_IF_0()
    #define DBG_WAIT_IF_1()
    #define DBG_WAIT_IF_2(cond, dbgStr)                             { if (cond) { std::stringstream _ss; _ss << dbgStr; Log::Instance()->DebugWait("", 0, _ss.str()); }}
    #define DBG_WAIT_IF_3(cond, dbgLevel, dbgStr)                   { if (cond) { std::stringstream _ss; _ss << dbgStr; Log::Instance()->DebugWait("", dbgLevel, _ss.str()); }}
    #define DBG_WAIT_IF_4(cond, dbgClass, dbgLevel, dbgStr)         { if (cond) { std::stringstream _ss; _ss << dbgStr; Log::Instance()->DebugWait(dbgClass, dbgLevel, _ss.str()); }}
    #define DBG_WAIT_IF_5()
    #define DBG_WAIT_IF(...)                                        GET_MACRO(_0, ##__VA_ARGS__, DBG_WAIT_IF_5, DBG_WAIT_IF_4, DBG_WAIT_IF_3, DBG_WAIT_IF_2, DBG_WAIT_IF_1, DBG_WAIT_IF_0)(__VA_ARGS__)

#else
    #define DBG(...)
    #define DBG_ID(...)
    #define DBG_IF(...)
    #define DBG_ID_IF(...)
    #define DBG_WAIT(...)
    #define DBG_WAIT_IF(...)
#endif


// Warning (to stdout)

#ifdef DEBUG_WARNINGS

    #define DBG_WARN_0()
    #define DBG_WARN_1(warnStr)                                     { std::stringstream _ss; _ss << warnStr; Log::Instance()->Say("", 0, _ss.str()); }
    #define DBG_WARN_2(warnLevel, warnStr)                          { std::stringstream _ss; _ss << warnStr; Log::Instance()->Say("", warnLevel, _ss.str()); }
    #define DBG_WARN_3(warnClass, warnLevel, warnStr)               { std::stringstream _ss; _ss << warnStr; Log::Instance()->Say(warnClass, warnLevel, _ss.str()); }
    #define DBG_WARN_4()
    #define DBG_WARN_5()
    #define DBG_WARN(...)                                           GET_MACRO(_0, ##__VA_ARGS__, DBG_WARN_5, DBG_WARN_4, DBG_WARN_3, DBG_WARN_2, DBG_WARN_1, DBG_WARN_0)(__VA_ARGS__)

    #define DBG_WARN_ID_0()                                         { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";              Log::Instance()->Say("", 0, _ss.str()); }
    #define DBG_WARN_ID_1(warnStr)                                  { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << warnStr; Log::Instance()->Say("", 0, _ss.str()); }
    #define DBG_WARN_ID_2(warnLevel, warnStr)                       { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << warnStr; Log::Instance()->Say("", warnLevel, _ss.str()); }
    #define DBG_WARN_ID_3(warnClass, warnLevel, warnStr)            { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << warnStr; Log::Instance()->Say(warnClass, warnLevel, _ss.str()); }
    #define DBG_WARN_ID_4()
    #define DBG_WARN_ID_5()
    #define DBG_WARN_ID(...)                                        GET_MACRO(_0, ##__VA_ARGS__, DBG_WARN_ID_5, DBG_WARN_ID_4, DBG_WARN_ID_3, DBG_WARN_ID_2, DBG_WARN_ID_1, DBG_WARN_ID_0)(__VA_ARGS__)

    #define DBG_WARN_IF_0()
    #define DBG_WARN_IF_1()
    #define DBG_WARN_IF_2(cond, warnStr)                            { if (cond) { std::stringstream _ss; _ss << warnStr; Log::Instance()->Say("", 0, _ss.str()); }}
    #define DBG_WARN_IF_3(cond, warnLevel, warnStr)                 { if (cond) { std::stringstream _ss; _ss << warnStr; Log::Instance()->Say("", warnLevel, _ss.str()); }}
    #define DBG_WARN_IF_4(cond, warnClass, warnLevel, warnStr)      { if (cond) { std::stringstream _ss; _ss << warnStr; Log::Instance()->Say(warnClass, warnLevel, _ss.str()); }}
    #define DBG_WARN_IF_5()
    #define DBG_WARN_IF(...)                                        GET_MACRO(_0, ##__VA_ARGS__, DBG_WARN_IF_5, DBG_WARN_IF_4, DBG_WARN_IF_3, DBG_WARN_IF_2, DBG_WARN_IF_1, DBG_WARN_IF_0)(__VA_ARGS__)

    #define DBG_WARN_ID_IF_0()
    #define DBG_WARN_ID_IF_1(cond)                                  { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";              Log::Instance()->Say("", 0, _ss.str()); }}
    #define DBG_WARN_ID_IF_2(cond, warnStr)                         { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << warnStr; Log::Instance()->Say("", 0, _ss.str()); }}
    #define DBG_WARN_ID_IF_3(cond, warnLevel, warnStr)              { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << warnStr; Log::Instance()->Say("", warnLevel, _ss.str()); }}
    #define DBG_WARN_ID_IF_4(cond, warnClass, warnLevel, warnStr)   { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << warnStr; Log::Instance()->Say(warnClass, warnLevel, _ss.str()); }}
    #define DBG_WARN_ID_IF_5()
    #define DBG_WARN_ID_IF(...)                                     GET_MACRO(_0, ##__VA_ARGS__, DBG_WARN_ID_IF_5, DBG_WARN_ID_IF_4, DBG_WARN_ID_IF_3, DBG_WARN_ID_IF_2, DBG_WARN_ID_IF_1, DBG_WARN_ID_IF_0)(__VA_ARGS__)

#else
    #define DBG_WARN(...)
    #define DBG_WARN_ID(...)
    #define DBG_WARN_IF(...)
    #define DBG_WARN_ID_IF(...)
#endif


// Messaging (to stdout)

#define SAY_0()
#define SAY_1(sayStr)                                               { std::stringstream _ss; _ss << sayStr; Log::Instance()->Say("", 0, _ss.str()); }
#define SAY_2(sayLevel, sayStr)                                     { std::stringstream _ss; _ss << sayStr; Log::Instance()->Say("", sayLevel, _ss.str()); }
#define SAY_3(sayClass, sayLevel, sayStr)                           { std::stringstream _ss; _ss << sayStr; Log::Instance()->Say(sayClass, sayLevel, _ss.str()); }
#define SAY_4()
#define SAY_5()
#define SAY(...)                                                    GET_MACRO(_0, ##__VA_ARGS__, SAY_5, SAY_4, SAY_3, SAY_2, SAY_1, SAY_0)(__VA_ARGS__)

#define SAY_ID_0()                                                  { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Say("", 0, _ss.str()); }
#define SAY_ID_1(sayStr)                                            { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << sayStr; Log::Instance()->Say("", 0, _ss.str()); }
#define SAY_ID_2(sayLevel, sayStr)                                  { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << sayStr; Log::Instance()->Say("", sayLevel, _ss.str()); }
#define SAY_ID_3(sayClass, sayLevel, sayStr)                        { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << sayStr; Log::Instance()->Say(sayClass, sayLevel, _ss.str()); }
#define SAY_ID_4()
#define SAY_ID_5()
#define SAY_ID(...)                                                 GET_MACRO(_0, ##__VA_ARGS__, SAY_ID_5, SAY_ID_4, SAY_ID_3, SAY_ID_2, SAY_ID_1, SAY_ID_0)(__VA_ARGS__)

#define SAY_IF_0()
#define SAY_IF_1()
#define SAY_IF_2(cond, sayStr)                                      { if (cond) { std::stringstream _ss; _ss << sayStr; Log::Instance()->Say("", 0, _ss.str()); }}
#define SAY_IF_3(cond, sayLevel, sayStr)                            { if (cond) { std::stringstream _ss; _ss << sayStr; Log::Instance()->Say("", sayLevel, _ss.str()); }}
#define SAY_IF_4(cond, sayClass, sayLevel, sayStr)                  { if (cond) { std::stringstream _ss; _ss << sayStr; Log::Instance()->Say(sayClass, sayLevel, _ss.str()); }}
#define SAY_IF_5()
#define SAY_IF(...)                                                 GET_MACRO(_0, ##__VA_ARGS__, SAY_IF_5, SAY_IF_4, SAY_IF_3, SAY_IF_2, SAY_IF_1, SAY_IF_0)(__VA_ARGS__)

#define SAY_ID_IF_0()
#define SAY_ID_IF_1(cond)                                           { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Say("", 0, _ss.str()); }}
#define SAY_ID_IF_2(cond, sayStr)                                   { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << sayStr; Log::Instance()->Say("", 0, _ss.str()); }}
#define SAY_ID_IF_3(cond, sayLevel, sayStr)                         { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << sayStr; Log::Instance()->Say("", sayLevel, _ss.str()); }}
#define SAY_ID_IF_4(cond, sayClass, sayLevel, sayStr)               { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << sayStr; Log::Instance()->Say(sayClass, sayLevel, _ss.str()); }}
#define SAY_ID_IF_5()
#define SAY_ID_IF(...)                                              GET_MACRO(_0, ##__VA_ARGS__, SAY_ID_IF_5, SAY_ID_IF_4, SAY_ID_IF_3, SAY_ID_IF_2, SAY_ID_IF_1, SAY_ID_IF_0)(__VA_ARGS__)


// Logging (to log files with LOG(), and stdout with LOGV())

#define LOG_0()
#define LOG_1()
#define LOG_2(logfileId, logStr)                                    { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, "", 0, _ss.str()); }
#define LOG_3(logfileId, logLevel, logStr)                          { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, "", logLevel, _ss.str()); }
#define LOG_4(logfileId, logClass, logLevel, logStr)                { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); }
#define LOG_5()
#define LOG(...)                                                    GET_MACRO(_0, ##__VA_ARGS__, LOG_5, LOG_4, LOG_3, LOG_2, LOG_1, LOG_0)(__VA_ARGS__)

#define LOG_ID_0()
#define LOG_ID_1(logfileId)                                         { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Put(logfileId, "", 0, _ss.str()); }
#define LOG_ID_2(logfileId, logStr)                                 { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", 0, _ss.str()); }
#define LOG_ID_3(logfileId, logLevel, logStr)                       { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", logLevel, _ss.str()); }
#define LOG_ID_4(logfileId, logClass, logLevel, logStr)             { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); }
#define LOG_ID_5()
#define LOG_ID(...)                                                 GET_MACRO(_0, ##__VA_ARGS__, LOG_ID_5, LOG_ID_4, LOG_ID_3, LOG_ID_2, LOG_ID_1, LOG_ID_0)(__VA_ARGS__)

#define LOG_IF_0()
#define LOG_IF_1()
#define LOG_IF_2()
#define LOG_IF_3(logfileId, cond, logStr)                           { if (cond) { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, "", 0, _ss.str()); }}
#define LOG_IF_4(logfileId, cond, logLevel, logStr)                 { if (cond) { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, "", logLevel, _ss.str()); }}
#define LOG_IF_5(logfileId, cond, logClass, logLevel, logStr)       { if (cond) { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); }}
#define LOG_IF(...)                                                 GET_MACRO(_0, ##__VA_ARGS__, LOG_IF_5, LOG_IF_4, LOG_IF_3, LOG_IF_2, LOG_IF_1, LOG_IF_0)(__VA_ARGS__)

#define LOG_ID_IF_0()
#define LOG_ID_IF_1()
#define LOG_ID_IF_2(logfileId, cond)                                { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Put(logfileId, "", 0, _ss.str()); }}
#define LOG_ID_IF_3(logfileId, cond, logStr)                        { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", 0, _ss.str()); }}
#define LOG_ID_IF_4(logfileId, cond, logLevel, logStr)              { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", logLevel, _ss.str()); }}
#define LOG_ID_IF_5(logfileId, cond, logClass, logLevel, logStr)    { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); }}
#define LOG_ID_IF(...)                                              GET_MACRO(_0, ##__VA_ARGS__, LOG_ID_IF_5, LOG_ID_IF_4, LOG_ID_IF_3, LOG_ID_IF_2, LOG_ID_IF_1, LOG_ID_IF_0)(__VA_ARGS__)

#define LOGV_0()
#define LOGV_1()
#define LOGV_2(logfileId, logStr)                                   { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, "", 0, _ss.str());              Log::Instance()->Say("", 0, _ss.str()); }
#define LOGV_3(logfileId, logLevel, logStr)                         { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, "", logLevel, _ss.str());       Log::Instance()->Say("", logLevel, _ss.str()); }
#define LOGV_4(logfileId, logClass, logLevel, logStr)               { std::stringstream _ss; _ss << logStr; Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); Log::Instance()->Say(logClass, logLevel, _ss.str()); }
#define LOGV_5()
#define LOGV(...)                                                   GET_MACRO(_0, ##__VA_ARGS__, LOGV_5, LOGV_4, LOGV_3, LOGV_2, LOGV_1, LOGV_0)(__VA_ARGS__)

#define LOGV_ID_0()
#define LOGV_ID_1(logfileId)                                        { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Put(logfileId, "", 0, _ss.str());              Log::Instance()->Say("", 0, _ss.str()); }
#define LOGV_ID_2(logfileId, logStr)                                { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", 0, _ss.str());              Log::Instance()->Say("", 0, _ss.str()); }
#define LOGV_ID_3(logfileId, logLevel, logStr)                      { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", logLevel, _ss.str());       Log::Instance()->Say("", logLevel, _ss.str()); }
#define LOGV_ID_4(logfileId, logClass, logLevel, logStr)            { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); Log::Instance()->Say(logClass, logLevel, _ss.str()); }
#define LOGV_ID_5()
#define LOGV_ID(...)                                                GET_MACRO(_0, ##__VA_ARGS__, LOGV_ID_5, LOGV_ID_4, LOGV_ID_3, LOGV_ID_2, LOGV_ID_1, LOGV_ID_0)(__VA_ARGS__)

#define LOGV_IF_0()
#define LOGV_IF_1()
#define LOGV_IF_2()
#define LOGV_IF_3(logfileId, cond, logStr)                          { if (cond) { std::stringstream _ss; _ss << logStr;  Log::Instance()->Put(logfileId, "", 0, _ss.str());              Log::Instance()->Say("", 0, _ss.str()); }}
#define LOGV_IF_4(logfileId, cond, logLevel, logStr)                { if (cond) { std::stringstream _ss; _ss << logStr;  Log::Instance()->Put(logfileId, "", logLevel, _ss.str());       Log::Instance()->Say("", logLevel, _ss.str()); }}
#define LOGV_IF_5(logfileId, cond, logClass, logLevel, logStr)      { if (cond) { std::stringstream _ss; _ss << logStr;  Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); Log::Instance()->Say(logClass, logLevel, _ss.str()); }}
#define LOGV_IF(...)                                                GET_MACRO(_0, ##__VA_ARGS__, LOGV_IF_5, LOGV_IF_4, LOGV_IF_3, LOGV_IF_2, LOGV_IF_1, LOGV_IF_0)(__VA_ARGS__)

#define LOGV_ID_IF_0()
#define LOGV_ID_IF_1()
#define LOGV_ID_IF_2(logfileId, cond)                               { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'";             Log::Instance()->Put(logfileId, "", 0, _ss.str());              Log::Instance()->Say("", 0, _ss.str()); }}
#define LOGV_ID_IF_3(logfileId, cond, logStr)                       { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", 0, _ss.str());              Log::Instance()->Say("", 0, _ss.str()); }}
#define LOGV_ID_IF_4(logfileId, cond, logLevel, logStr)             { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, "", logLevel, _ss.str());       Log::Instance()->Say("", logLevel, _ss.str()); }}
#define LOGV_ID_IF_5(logfileId, cond, logClass, logLevel, logStr)   { if (cond) { std::stringstream _ss; _ss << "IN FUNCTION " << "'" << __PRETTY_FUNCTION__ << "'\n" << logStr; Log::Instance()->Put(logfileId, logClass, logLevel, _ss.str()); Log::Instance()->Say(logClass, logLevel, _ss.str()); }}
#define LOGV_ID_IF(...)                                             GET_MACRO(_0, ##__VA_ARGS__, LOGV_ID_IF_5, LOGV_ID_IF_4, LOGV_ID_IF_3, LOGV_ID_IF_2, LOGV_ID_IF_1, LOGV_ID_IF_0)(__VA_ARGS__)


#endif // __defines_h__
