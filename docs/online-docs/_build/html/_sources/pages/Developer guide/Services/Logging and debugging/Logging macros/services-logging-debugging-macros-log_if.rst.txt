LOG_IF() macro
==============

::

    LOG_IF(id, cond, ...)

Writes a log record to the log file specified by ``id`` if the condition given by ``cond`` is met. Usage::

    LOG_IF(id, cond, string)                // writes 'string' to the log file if 'cond' is TRUE.
    
    LOG_IF(id, cond, level, string)         // writes 'string' to the log file if 'cond' is TRUE and 
                                            // 'level' <= 'p_LogLevel' in Log::Start().
    
    LOG_IF(id, cond, class, level, string)  // writes 'string' to the log file if 'cond' is TRUE, 
                                            // 'class' is in 'p_LogClasses' in Log::Start(), and 
                                            // 'level' <= 'p_LogLevel' in Log::Start().

Default class is ””; default level is 0.

Examples::

    LOG_IF(MyLogfileId, a > 1.0, ”This is a log record”);
    LOG_IF(SSEfileId, (b == c && a > x), ”The value of x is ” << x << ” km”);
    LOG_IF(CHeBfileId, flag, 2, ”Log string”);
    LOG_IF(SSEfileId, (x >= y), ”CHeB”, 4, ”This is a CHeB only log record”);
