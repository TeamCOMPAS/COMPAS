Log::Put(...)
=============

::

    BOOL Log::Put(
        INT    p_LogfileId,  // the identifier of the log file to be written.
        STRING p_LogClass,   // the log class of the record to be written. An empty string (””) satisfies
                             // all checks against enabled classes.
        INT    p_LogLevel,   // the log level of the record to be written. A value of 0 satisfies all 
                             // checks against enabled levels.
        STRING p_LogStr      // the string to be written to the log file.
    )

Writes a minimally formatted record to the specified log file. If the Log service is enabled, the specified log file is active, and the log
record ``class`` and ``level`` passed are enabled, the string is written to the file. See :doc:`./services-base-level-logging` for details
regarding log record ``class`` and ``level``.

If labels are enabled for the log file, a label will be prepended to the record. The label text will be the ``p_LogClass`` parameter.

If timestamps are enabled for the log file, a formatted timestamp is prepended to the record. The timestamp format is **yyyymmdd hh:mm:ss**.

Returns a boolean indicating whether the record was written successfully. If an error occurred the log file will be disabled.
