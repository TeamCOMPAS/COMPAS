Log::Error(p_ErrStr)
====================

::

    BOOL Log::Error(
        STRING p_ErrStr  // the string to be written to stdout (and optionally to file).
    )

Writes ``p_ErrStr`` to stdout and, if logging is active and so configured (via program option ``--errors-to-file``), writes ``p_ErrStr`` 
to the error log file.

Returns a boolean indicating whether the record was written successfully. If an error occurred writing to the error log file, 
the log file will be disabled.
