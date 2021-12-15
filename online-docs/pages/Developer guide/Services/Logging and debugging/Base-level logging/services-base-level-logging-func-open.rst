Log::Open(...)
==============

::

    INT Log::Open(
        STRING  p_LogFileName,     // the name of the log file to be created and opened.
                                   // This should be the filename only – the path, prefix and
                                   // extensions are added by the logging service. If the file 
                                   // already exists, the logging service will append a version 
                                   // number to the name if necessary (see p_Append parameter below).
        BOOL    p_Append,          // flag indicating whether the file should be opened in append mode
                                   // (i.e. existing data is preserved) and new records written to the
                                   // file appended, or whether a new file should be opened (with 
                                   // version number if necessary).
        BOOL    p_TimeStamp,       // flag indicating whether timestamps should be written with each 
                                   // log file record.
        BOOL    p_Label,           // flag indicating whether a label should be written with each log 
                                   // record. This is useful when different types of logging data are 
                                   // being written to the same log file file.
        LOGFILE p_StandardLogfile  // (optional) the standard log file type for this log file.
                                   // If not provided LOGFILE::NONE is used.
    )

Opens a log file. If the ``p_Append`` parameter is TRUE and a file named ``p_LogFileName`` exists, the existing file will be opened and
the existing contents retained, otherwise a new file will be created and opened (not a Standard Log File – see 
:doc:`../services-extended-logging`).

The log file container directory is created at the path specified by the ``p_LogBasePath`` parameter passed to the 
:doc:`./services-base-level-logging-func-start` function. New log files are created in the log file container directory. ``BSE`` detailed 
log files are created in the ``Detailed_Output`` directory, which is created in the log file container directory (if required).

The filename is prefixed by the ``p_LogNamePrefix`` parameter passed to the :doc:`./services-base-level-logging-func-start` function.

If a file with the name as given by the ``p_LogFileName`` parameter already exists, and the ``p_Append`` parameter is false, a version
number will be appended to the filename before the extension (this functionality is largely redundant since the implementation of the 
log file container directory).

The integer log file identifier is returned to the caller. A value of −1 indicates the log file was not opened successfully. Multiple 
log files can be open simultaneously – referenced by the identifier returned.
