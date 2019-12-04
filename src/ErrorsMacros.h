
#undef WARNINGS          // define/undef this to turn warnings on/off       JR: todo: make this a program option

// ERRORS macros
#define GET_MACRO(_0, _1, _2, _3, _4, _5, NAME, ...)                NAME

#define SHOW_ERROR_0()
#define SHOW_ERROR_1(error)                         { ERRORS->ShowIt("ERROR: ", error, "", ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_ERROR_2(error, qualifyingStr)          { ERRORS->ShowIt("ERROR: ", error, qualifyingStr, ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_ERROR_3()
#define SHOW_ERROR_4()
#define SHOW_ERROR_5()
#define SHOW_ERROR(...)                             GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_5, SHOW_ERROR_4, SHOW_ERROR_3, SHOW_ERROR_2, SHOW_ERROR_1, SHOW_ERROR_0) (__VA_ARGS__)

#define SHOW_ERROR_IF_0()
#define SHOW_ERROR_IF_1()
#define SHOW_ERROR_IF_2(cond, error)                { if (cond) { ERRORS->ShowIt("ERROR: ", error, "", ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }}
#define SHOW_ERROR_IF_3(cond, error, qualifyingStr) { if (cond) { ERRORS->ShowIt("ERROR: ", error, qualifyingStr, ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }}
#define SHOW_ERROR_IF_4()
#define SHOW_ERROR_IF_5()
#define SHOW_ERROR_IF(...)                          GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_IF_5, SHOW_ERROR_IF_4, SHOW_ERROR_IF_3, SHOW_ERROR_IF_2, SHOW_ERROR_IF_1, SHOW_ERROR_IF_0) (__VA_ARGS__)


#ifdef WARNINGS

#define SHOW_WARN_0()
#define SHOW_WARN_1(error)                          { ERRORS->ShowIt("WARNING: ", error, "", ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_WARN_2(error, qualifyingStr)           { ERRORS->ShowIt("WARNING: ", error, qualifyingStr, ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_WARN_3()
#define SHOW_WARN_4()
#define SHOW_WARN_5()
#define SHOW_WARN(...)                              GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_5, SHOW_WARN_4, SHOW_WARN_3, SHOW_WARN_2, SHOW_WARN_1, SHOW_WARN_0) (__VA_ARGS__)

#define SHOW_WARN_IF_0()
#define SHOW_WARN_IF_1()
#define SHOW_WARN_IF_2(cond, error)                 { if (cond) { ERRORS->ShowIt("WARNING: ", error, "", ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }}
#define SHOW_WARN_IF_3(cond, error, qualifyingStr)  { if (cond) { ERRORS->ShowIt("WARNING: ", error, qualifyingStr, ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }}
#define SHOW_WARN_IF_4()
#define SHOW_WARN_IF_5()
#define SHOW_WARN_IF(...)                           GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_IF_5, SHOW_WARN_IF_4, SHOW_WARN_IF_3, SHOW_WARN_IF_2, SHOW_WARN_IF_1, SHOW_WARN_IF_0) (__VA_ARGS__)

#else
    #define SHOW_WARN(...)
    #define SHOW_WARN_IF(...)
#endif


// macros for static functions
#define SHOW_ERROR_STATIC_0()
#define SHOW_ERROR_STATIC_1()
#define SHOW_ERROR_STATIC_2()
#define SHOW_ERROR_STATIC_3(error, objectType, stellarType)                 { ERRORS->ShowIt("ERROR: ", error, "", -1l, objectType, stellarType, __PRETTY_FUNCTION__); }
#define SHOW_ERROR_STATIC_4(error, qualifyingStr, objectType, stellarType)  { ERRORS->ShowIt("ERROR: ", error, qualifyingStr, -1l, objectType, stellarType, __PRETTY_FUNCTION__); }
#define SHOW_ERROR_STATIC_5()
#define SHOW_ERROR_STATIC(...)                      GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_STATIC_5, SHOW_ERROR_STATIC_4, SHOW_ERROR_STATIC_3, SHOW_ERROR_STATIC_2, SHOW_ERROR_STATIC_1, SHOW_ERROR_STATIC_0) (__VA_ARGS__)

#define SHOW_ERROR_IF_STATIC_0()
#define SHOW_ERROR_IF_STATIC_1()
#define SHOW_ERROR_IF_STATIC_2()
#define SHOW_ERROR_IF_STATIC_3()
#define SHOW_ERROR_IF_STATIC_4(cond, error, objectType, stellarType)                { if (cond) { ERRORS->ShowIt("ERROR: ", error, "", -1l, objectType, stellarType, __PRETTY_FUNCTION__); }}
#define SHOW_ERROR_IF_STATIC_5(cond, error, qualifyingStr, objectType, stellarType) { if (cond) { ERRORS->ShowIt("ERROR: ", error, qualifyingStr, -1l, objectType, stellarType, __PRETTY_FUNCTION__); }}
#define SHOW_ERROR_IF_STATIC(...)                   GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_IF_STATIC_5, SHOW_ERROR_IF_STATIC_4, SHOW_ERROR_IF_STATIC_3, SHOW_ERROR_IF_STATIC_2, SHOW_ERROR_IF_STATIC_1, SHOW_ERROR_IF_STATIC_0) (__VA_ARGS__)


#ifdef WARNINGS

#define SHOW_WARN_STATIC_0()
#define SHOW_WARN_STATIC_1()
#define SHOW_WARN_STATIC_2()
#define SHOW_WARN_STATIC_3(error, objectType, stellarType)                  { ERRORS->ShowIt("WARNING: ", error, "", -1l, objectType, stellarType, __PRETTY_FUNCTION__); }
#define SHOW_WARN_STATIC_4(error, qualifyingStr, objectType, stellarType)   { ERRORS->ShowIt("WARNING: ", error, qualifyingStr, -1l, objectType, stellarType, __PRETTY_FUNCTION__); }
#define SHOW_WARN_STATIC_5()
#define SHOW_WARN_STATIC(...)                       GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_STATIC_5, SHOW_WARN_STATIC_4, SHOW_WARN_STATIC_3, SHOW_WARN_STATIC_2, SHOW_WARN_STATIC_1, SHOW_WARN_STATIC_0) (__VA_ARGS__)

#define SHOW_WARN_IF_STATIC_0()
#define SHOW_WARN_IF_STATIC_1()
#define SHOW_WARN_IF_STATIC_2(cond, error, objectType, stellarType)                 { if (cond) { ERRORS->ShowIt("WARNING: ", error, "", -1l, objectType, stellarType, __PRETTY_FUNCTION__); }}
#define SHOW_WARN_IF_STATIC_3(cond, error, qualifyingStr, objectType, stellarType)  { if (cond) { ERRORS->ShowIt("WARNING: ", error, qualifyingStr, -1l, objectType, stellarType, __PRETTY_FUNCTION__); }}
#define SHOW_WARN_IF_STATIC_4()
#define SHOW_WARN_IF_STATIC_5()
#define SHOW_WARN_IF_STATIC(...)                    GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_IF_STATIC_5, SHOW_WARN_IF_STATIC_4, SHOW_WARN_IF_STATIC_3, SHOW_WARN_IF_STATIC_2, SHOW_WARN_IF_STATIC_1, SHOW_WARN_IF_STATIC_0) (__VA_ARGS__)

#else
    #define SHOW_WARN_STATIC(...)
    #define SHOW_WARN_IF_STATIC(...)
#endif
