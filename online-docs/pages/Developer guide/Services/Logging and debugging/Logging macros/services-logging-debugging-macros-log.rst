LOG() macro
===========

::

    LOG(id, ...)

Writes a log record to the log file specified by ``id``. Usage::

    LOG(id, string)                // writes 'string' to the log file

    LOG(id, level, string)         // writes 'string' to the log file if 'level' <= 'p_LogLevel' in
                                   // Log::Start()

    LOG(id, class, level, string)  // writes 'string to the log file if 'class' is in 'p_LogClasses' in
                                   // Log::Start(), and 'level' <= 'p_LogLevel' in Log::Start()

Default class is ””; default level is 0.

Examples::

    LOG(SSEfileId, ”This is a log record”);
    LOG(OutputFile2Id, ”The value of x is ” << x << ” km”);
    LOG(MyLogfileId, 2, ”Log string”);
    LOG(SSEfileId, ”CHeB”, 4, ”This is a CHeB only log record”); 