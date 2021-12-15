Log::Write(...)
===============

::

    BOOL Log::Write(
        INT    p_LogfileId,  // the identifier of the log file to be written.
        STRING p_LogClass,   // the log class of the record to be written. An empty string (””) satisfies
                             // all checks against enabled classes.
        INT    p_LogLevel,   // the log level of the record to be written. A value of 0 satisfies all 
                             // checks against enabled levels.
        STRING p_LogStr      // the string to be written to the log file.
    )

Writes an unformatted record to the specified log file. If the Log service is enabled, the specified log file is active, and the log 
record ``class`` and ``level`` passed are enabled, the string is written to the file. See :doc:`./services-base-level-logging` for details
regarding log record ``class`` and ``level``.

Returns a boolean indicating whether the record was written successfully. If an error occurred the log file will be disabled.
