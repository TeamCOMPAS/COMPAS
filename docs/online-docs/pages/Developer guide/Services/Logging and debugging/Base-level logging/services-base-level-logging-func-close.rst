Log::Close(p_LogfileId)
=======================

::

    BOOL Log::Close(
        INT p_LogfileId  // the identifier of the log file to be closed (as returned by Log::Open()).
    )

Closes the log file specified by the ``p_LogfileId`` parameter. If the log file specified by the logFileId parameter is open, it is flushed
to disk and closed.

Returns a boolean indicating whether the file was closed successfully.
