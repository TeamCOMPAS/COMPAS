Log::Debug(...)
===============

::

  BOOL Log::Debug(
      STRING p_DbgClass,  // the log class of the record to be written. An empty string (””) satisfies
                          // all checks against enabled classes.
      INT    p_DbgLevel,  // the log level of the record to be written. A value of 0 satisfies all 
                          // checks against enabled levels.
      STRING p_DbgStr     // the string to be written to stdout (and optionally to file).
  )

Writes ``p_DbgStr`` to stdout and, if logging is active and so configured (via program option ``--debug-to-file``), writes ``p_DbgStr`` 
to the debug log file.

Returns a boolean indicating whether the record was written successfully. If an error occurred writing to the debug log file, 
the log file will be disabled.
