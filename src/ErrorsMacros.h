#define ERROR_PREFIX "ERROR: "
#define WARNING_PREFIX "WARNING: "

// ERRORS macros
#define GET_MACRO(_0, _1, _2, _3, _4, _5, NAME, ...) NAME

#define SHOW_ERROR_0()
#define SHOW_ERROR_1(error)                          { ERRORS->ShowIt(ERROR_PREFIX, error, "", ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_ERROR_2(error, qualifyingStr)           { ERRORS->ShowIt(ERROR_PREFIX, error, qualifyingStr, ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_ERROR_3()
#define SHOW_ERROR_4()
#define SHOW_ERROR_5()
#define SHOW_ERROR(...)                              GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_5, SHOW_ERROR_4, SHOW_ERROR_3, SHOW_ERROR_2, SHOW_ERROR_1, SHOW_ERROR_0) (__VA_ARGS__)

#define SHOW_ERROR_IF_0()
#define SHOW_ERROR_IF_1()
#define SHOW_ERROR_IF_2(cond, error)                 { if (cond) { SHOW_ERROR(error); }}
#define SHOW_ERROR_IF_3(cond, error, qualifyingStr)  { if (cond) { SHOW_ERROR(error, qualifyingStr); }}
#define SHOW_ERROR_IF_4()
#define SHOW_ERROR_IF_5()
#define SHOW_ERROR_IF(...)                           GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_IF_5, SHOW_ERROR_IF_4, SHOW_ERROR_IF_3, SHOW_ERROR_IF_2, SHOW_ERROR_IF_1, SHOW_ERROR_IF_0) (__VA_ARGS__)


#define THROW_ERROR_0()
#define THROW_ERROR_1(error)                         { SHOW_ERROR(error); throw static_cast<int>(error); }
#define THROW_ERROR_2(error, qualifyingStr)          { SHOW_ERROR(error, qualifyingStr); throw static_cast<int>(error); }
#define THROW_ERROR_3()
#define THROW_ERROR_4()
#define THROW_ERROR_5()
#define THROW_ERROR(...)                             GET_MACRO(_0, ##__VA_ARGS__, THROW_ERROR_5, THROW_ERROR_4, THROW_ERROR_3, THROW_ERROR_2, THROW_ERROR_1, THROW_ERROR_0) (__VA_ARGS__)

#define THROW_ERROR_IF_0()
#define THROW_ERROR_IF_1()
#define THROW_ERROR_IF_2(cond, error)                { if (cond) { SHOW_ERROR(error); throw static_cast<int>(error); }}
#define THROW_ERROR_IF_3(cond, error, qualifyingStr) { if (cond) { SHOW_ERROR(error, qualifyingStr); throw static_cast<int>(error); }}
#define THROW_ERROR_IF_4()
#define THROW_ERROR_IF_5()
#define THROW_ERROR_IF(...)                          GET_MACRO(_0, ##__VA_ARGS__, THROW_ERROR_IF_5, THROW_ERROR_IF_4, THROW_ERROR_IF_3, THROW_ERROR_IF_2, THROW_ERROR_IF_1, THROW_ERROR_IF_0) (__VA_ARGS__)


#define SHOW_WARN_0()
#define SHOW_WARN_1(error)                           { ERRORS->ShowIt(WARNING_PREFIX, error, "", ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_WARN_2(error, qualifyingStr)            { ERRORS->ShowIt(WARNING_PREFIX, error, qualifyingStr, ObjectId(), ObjectType(), StellarType(), __PRETTY_FUNCTION__); }
#define SHOW_WARN_3()
#define SHOW_WARN_4()
#define SHOW_WARN_5()
#define SHOW_WARN(...)                               GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_5, SHOW_WARN_4, SHOW_WARN_3, SHOW_WARN_2, SHOW_WARN_1, SHOW_WARN_0) (__VA_ARGS__)

#define SHOW_WARN_IF_0()
#define SHOW_WARN_IF_1()
#define SHOW_WARN_IF_2(cond, error)                  { if (cond) { SHOW_WARN(error); }}
#define SHOW_WARN_IF_3(cond, error, qualifyingStr)   { if (cond) { SHOW_WARN(error, qualifyingStr); }}
#define SHOW_WARN_IF_4()
#define SHOW_WARN_IF_5()
#define SHOW_WARN_IF(...)                            GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_IF_5, SHOW_WARN_IF_4, SHOW_WARN_IF_3, SHOW_WARN_IF_2, SHOW_WARN_IF_1, SHOW_WARN_IF_0) (__VA_ARGS__)


// macros for static functions
#define SHOW_ERROR_STATIC_0()
#define SHOW_ERROR_STATIC_1(error)                                                   { ERRORS->ShowIt(ERROR_PREFIX, error, "", -1l, OBJECT_TYPE::NONE, STELLAR_TYPE::NONE, __PRETTY_FUNCTION__); }
#define SHOW_ERROR_STATIC_2(error, qualifyingStr)                                    { ERRORS->ShowIt(ERROR_PREFIX, error, qualifyingStr, -1l, OBJECT_TYPE::NONE, STELLAR_TYPE::NONE, __PRETTY_FUNCTION__); }
#define SHOW_ERROR_STATIC_3()
#define SHOW_ERROR_STATIC_4()
#define SHOW_ERROR_STATIC_5()
#define SHOW_ERROR_STATIC(...)                                                       GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_STATIC_5, SHOW_ERROR_STATIC_4, SHOW_ERROR_STATIC_3, SHOW_ERROR_STATIC_2, SHOW_ERROR_STATIC_1, SHOW_ERROR_STATIC_0) (__VA_ARGS__)

#define SHOW_ERROR_IF_STATIC_0()
#define SHOW_ERROR_IF_STATIC_1()
#define SHOW_ERROR_IF_STATIC_2(cond, error)                                          { if (cond) { SHOW_ERROR_STATIC(error); }}
#define SHOW_ERROR_IF_STATIC_3(cond, error, qualifyingStr)                           { if (cond) { SHOW_ERROR_STATIC(error, qualifyingStr); }}
#define SHOW_ERROR_IF_STATIC_4()
#define SHOW_ERROR_IF_STATIC_5()
#define SHOW_ERROR_IF_STATIC(...)                                                    GET_MACRO(_0, ##__VA_ARGS__, SHOW_ERROR_IF_STATIC_5, SHOW_ERROR_IF_STATIC_4, SHOW_ERROR_IF_STATIC_3, SHOW_ERROR_IF_STATIC_2, SHOW_ERROR_IF_STATIC_1, SHOW_ERROR_IF_STATIC_0) (__VA_ARGS__)


#define THROW_ERROR_STATIC_0()
#define THROW_ERROR_STATIC_1(error)                                                  { SHOW_ERROR_STATIC(error); throw static_cast<int>(error); }
#define THROW_ERROR_STATIC_2(error, qualifyingStr)                                   { SHOW_ERROR_STATIC(error, qualifyingStr); throw static_cast<int>(error); }
#define THROW_ERROR_STATIC_3()
#define THROW_ERROR_STATIC_4()
#define THROW_ERROR_STATIC_5()
#define THROW_ERROR_STATIC(...)                                                      GET_MACRO(_0, ##__VA_ARGS__, THROW_ERROR_STATIC_5, THROW_ERROR_STATIC_4, THROW_ERROR_STATIC_3, THROW_ERROR_STATIC_2, THROW_ERROR_STATIC_1, THROW_ERROR_STATIC_0) (__VA_ARGS__)

#define THROW_ERROR_IF_STATIC_0()
#define THROW_ERROR_IF_STATIC_1()
#define THROW_ERROR_IF_STATIC_2(cond, error)                                         { if (cond) { THROW_ERROR_STATIC(error); }}
#define THROW_ERROR_IF_STATIC_3(cond, error, qualifyingStr)                          { if (cond) { THROW_ERROR_STATIC(error, qualifyingStr); }}
#define THROW_ERROR_IF_STATIC_4()
#define THROW_ERROR_IF_STATIC_5()
#define THROW_ERROR_IF_STATIC(...)                                                   GET_MACRO(_0, ##__VA_ARGS__, THROW_ERROR_IF_STATIC_5, THROW_ERROR_IF_STATIC_4, THROW_ERROR_IF_STATIC_3, THROW_ERROR_IF_STATIC_2, THROW_ERROR_IF_STATIC_1, THROW_ERROR_IF_STATIC_0) (__VA_ARGS__)


#define SHOW_WARN_STATIC_0()
#define SHOW_WARN_STATIC_1(error)                                                    { ERRORS->ShowIt(WARNING_PREFIX, error, "", -1l, OBJECT_TYPE::NONE, STELLAR_TYPE::NONE, __PRETTY_FUNCTION__); }
#define SHOW_WARN_STATIC_2(error, qualifyingStr)                                     { ERRORS->ShowIt(WARNING_PREFIX, error, qualifyingStr, -1l, OBJECT_TYPE::NONE, STELLAR_TYPE::NONE, __PRETTY_FUNCTION__); }
#define SHOW_WARN_STATIC_3()
#define SHOW_WARN_STATIC_4()
#define SHOW_WARN_STATIC_5()
#define SHOW_WARN_STATIC(...)                                                        GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_STATIC_5, SHOW_WARN_STATIC_4, SHOW_WARN_STATIC_3, SHOW_WARN_STATIC_2, SHOW_WARN_STATIC_1, SHOW_WARN_STATIC_0) (__VA_ARGS__)

#define SHOW_WARN_IF_STATIC_0()
#define SHOW_WARN_IF_STATIC_1()
#define SHOW_WARN_IF_STATIC_2(cond, error)                                           { if (cond) { SHOW_WARN_STATIC(error); }}
#define SHOW_WARN_IF_STATIC_3(cond, error, qualifyingStr)                            { if (cond) { SHOW_WARN_STATIC(error, qualifyingStr); }}
#define SHOW_WARN_IF_STATIC_4()
#define SHOW_WARN_IF_STATIC_5()
#define SHOW_WARN_IF_STATIC(...)                                                     GET_MACRO(_0, ##__VA_ARGS__, SHOW_WARN_IF_STATIC_5, SHOW_WARN_IF_STATIC_4, SHOW_WARN_IF_STATIC_3, SHOW_WARN_IF_STATIC_2, SHOW_WARN_IF_STATIC_1, SHOW_WARN_IF_STATIC_0) (__VA_ARGS__)
