LOG_ID() macro
==============

::

    LOG_ID(id, ...)

Writes a log record prepended with calling function name to the log file specified by ``id``. Usage::

    LOG_ID(id)                        // writes the name of calling function to the log file

    LOG_ID(id, string)                // writes 'string' prepended with name of calling function to the
                                      // log file

    LOG_ID(id, level, string)         // writes 'string' prepended with name of calling function to the 
                                      // log file if 'level' <= 'p_LogLevel' in Log::Start()

    LOG_ID(id, class, level, string)  // writes 'string' prepended with name of calling function to the 
                                      // log file if 'class' is in 'p_LogClasses' in Log::Start(), and 
                                      // 'level' <= 'p_LogLevel' in Log::Start().

Default class is ””; default level is 0.

Examples::

    LOG_ID(Outf1Id);
    LOG_ID(Outf2Id, ”This is a log record”);
    LOG_ID(MyLogfileId, ”The value of x is ” << x << ” km”);
    LOG_ID(OutputFile2Id, 2, ”Log string”);
    LOG_ID(CHeBfileId, ”CHeB”, 4, ”This is a CHeB only log record”);
