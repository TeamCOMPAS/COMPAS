Extended logging
================

The extended logging service supports standard log files for both Single Star Evolution ``(SSE``) and Binary Star 
Evolution (``BSE``).

The standard log files defined are:

For ``SSE``:

- SSE_System_Parameters log file
- SSE_Supernovae log file
- SSE_Detailed_Output log file
- SSE_Switchlog log file

For ``BSE``:

- BSE_System_Parameters log file
- BSE_Double_Compact Objects log file
- BSE_Common_Envelopes log file
- BSE_Supernovae log file
- BSE_Pulsar_Evolution log file
- BSE_RLOF_Parameters log file
- BSE_Detailed_Output log file
- BSE_Switchlog log file

The Logging service maintains information about each of the standard log files, and will handle creating, opening, writing and
closing the files. For each execution of the COMPAS program, one (and only one) of each of the log files listed above that
pertain to the mode of evolution (``--mode`` option, ``SSE`` or ``BSE``) will be created, except for the ``Detailed_Output`` 
log files, in which case there will be one log file created for each system (single star or binary star) evolved.

The Logging service provides the following public member functions specifically for managing standard log files:

For ``SSE`` log files::

    BOOL LogSSESystemParameters(CONST T* CONST p_Star, CONST string p_Rec)
    BOOL LogSSESupernovaDetails(CONST T* CONST p_Star, CONST string p_Rec)
    BOOL LogSSEDetailedOutput(CONST T* CONST p_Star, CONST int p_Id, CONST string p_Rec)
    BOOL LogSSESwitchLog(CONST T* CONST p_Star, CONST string p_Rec)

Each ``SSE`` function is passed a pointer to the single star for which details are to be logged (``p_Star``), and a string to be 
written to the log file (``p_Rec``). If ``p_Rec`` is an empty string, the function constructs the log record from the current 
attributes of the star and the default record specifier for the log file (see property vectors in ``constants.h``,  e.g. 
``SSE_DETAILED_OUTPUT_REC``). ``LogSSEDetailedOutput()`` is also passed an integer identifier (typically the loop index of the
star) that is appended to the log file name (``p_Id``).

For ``BSE`` log files::

    BOOL LogBSESystemParameters(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogDoubleCompactObject(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogCommonEnvelope(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogBSESupernovaDetails(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogBSEPulsarEvolutionParameters(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogRLOFParameters(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogBSEDetailedOutput(CONST T* CONST p_Binary, CONST long int p_Id, CONST string p_Rec)
    BOOL LogBSESwitchLog(CONST T* CONST p_Binary, CONST bool p_PrimarySwitching)

Each ``BSE`` function is passed a pointer to the binary star for which details are to be logged (``p_Binary``), and a string to 
be  written to the log file (``p_Rec``). If ``p_Rec`` is an empty string, the function constructs the log record from the current 
attributes of the binary and the default record specifier for the log file (see property vectors in ``constants.h``,  e.g. 
``BSE_DETAILED_OUTPUT_REC``). ``LogBSEDetailedOutput()`` is also passed an integer identifier (typically the loop index of the
binary) that is appended to the log file name (``p_Id``).

Each of the functions listed above will, if necessary, create and open the appropriate log file. Internally the Log service opens
(creates first if necessary) once at first use, and keeps the files open for the life of the program.

The Log service provides a further two functions to manage standard log files::

    BOOL CloseStandardFile()
    BOOL CloseAllStandardFiles()

``CloseStandardFile()`` flushes and closes a standard log file. The function returns a boolean indicating whether the log file was
closed successfully.

``CloseAllStandardFiles()`` flushes and closes all currently open standard log files. The function returns a boolean indicating
whether all standard log files were closed successfully.

Standard log file names are supplied via program options (e.g. ``--logfile-system-parameters``), with default values declared in
``constants.h``.

The extended logging service always sets the log record ``class`` to the name of the standard log file being written to, and the
log record ``level`` to 0.  See :doc:`./Base-level logging/services-base-level-logging` for details regarding log record ``class``
and ``level``.

