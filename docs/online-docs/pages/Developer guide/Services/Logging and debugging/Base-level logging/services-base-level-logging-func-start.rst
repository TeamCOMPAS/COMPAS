Log::Start(...)
===============

::

    VOID Log::Start(**
        STRING      p_LogBasePath,       // the name of the top-level directory in which log files will
                                         // be created.
        STRING      p_LogContainerName,  // the name of the directory to be created at p_LogBasePath to
                                         // hold all log files.
        STRING      p_LogNamePrefix,     // prefix for log file names (can be blank).
        INT         p_LogLevel,          // logging level (see below).
        STRING[]    p_LogClasses,        // enabled logging classes (see below).
        INT         p_DbgLevel,          // debug level (see below).
        STRING[]    p_DbgClasses,        // enabled debug classes (see below).
        BOOL        p_DbgToFile,         // flag indicating whether debug statements should also be 
                                         // written to the log file.
        BOOL        p_ErrorsToFile,      // flag indicating whether error messages should also be 
                                         // written to the log file.
        LOGFILETYPE p_LogfileType        // the log file type (see below).
    )

Initialises and starts the logging and debugging service. Logging parameters are set per the program options specified (using
default values if no options are specified by the user). The log file container directory is created. If a directory with the
name as given by the containerName parameter already exists, a version number will be appended to the directory name. The 
``Run_Details`` file is created within the log file container directory. Any input files specified by the user, such as a grid
file and/or a log file definitions file, are copied to the log file container directory.  Log files to which debug statements
and error messages should be written will be created and opened if required.

The filename to which debug records are written when parameter ``p_ErrorsToFile`` is TRUE is declared in ``constants.h`` – see
enum class ``LOGFILE`` and associated descriptor map ``LOGFILE_DESCRIPTOR``. Currently the name is ”Debug_Log”.

Log file types are defined in the enum class ``LOGFILETYPE`` in ``constants.h``. The log file type and file extension are defined
by the ``p_LogfileType`` parameter:

    .. list-table::
       :widths: 25 75 
       :header-rows: 0
       :class: aligned-text

       * - LOGFILETYPE::TXT
         - will result in a plain text file, delimited by the space character (' '), with a file extension of ”.txt”
       * - LOGFILETYPE::TSV
         - will result in a plain text file, delimited by the tab character ('\\t'), with a file extension of ”.tsv”
       * - LOGFILETYPE::CSV
         - will result in a plain text file, delimited by the comma character (','), with a file extension of ”.csv”
       * - LOGFILETYPE::HDF5
         - will result in an ``HDF5``\ [#f1]_ file, with a file extension of ”.h5”

This function does not return a value.


.. rubric:: Footnotes

.. [#f1] https://www.hdfgroup.org/
