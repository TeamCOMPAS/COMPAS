Error handling service
======================

An error handling service is provided encapsulated in a singleton object (an instantiation of the ``Errors`` class).

The Errors service provides global error handling functionality.

Errors are defined in the error catalog in ``ErrorCatalog.h`` (see ``ERROR_CATALOG``). |br|

Errors defined in the error catalog have a scope and message text. The scope is used to determine if and when an error should
be printed.

The current values for scope are:

    .. list-table::
       :widths: 30 70 
       :header-rows: 0
       :class: aligned-text

       * - **NEVER**
         - the error will not be printed.
       * - **ALWAYS**
         -  the error will always be printed.
       * - **FIRST**
         -  the error will be printed only on the first time it is encountered anywhere in the program in the current execution of COMPAS.
       * - **FIRST_IN_OBJECT_TYPE**
         - the error will be printed only on the first time it is encountered anywhere in objects of the same type (e.g. Binary Star objects) in the current execution of COMPAS.
       * - **FIRST_IN_STELLAR_TYPE**
         - the error will be printed only on the first time it is encountered anywhere in objects of the same stellar type (e.g. HeWD Star objects) in the current execution of COMPAS.
       * - **FIRST_IN_OBJECT_ID**
         - the error will be printed only on the first time it is encountered anywhere in an object instance in the current execution of COMPAS (useful for debugging).
       * - **FIRST_IN_FUNCTION**
         - the error will be printed only on the first time it is encountered anywhere in the same function of an object instance in the current execution of COMPAS (i.e. will print more than once if encountered in the same function name in different objects; useful for debugging).

The Errors service provides methods to print both warnings and errors – essentially the same thing (for printing), but warning messages are prepended with 
':boldtext:`WARNING:` ', whereas error messages are prepended with ':boldtext:`ERROR:` '.

Errors and warnings are printed by using the macros defined in ``ErrorsMacros.h``. They are:


Error macros
------------

::
    
    SHOW_ERROR(error_number)                         // prints the error message associated with error
                                                     // number (from the error catalog) prepended by
                                                     // 'ERROR: '

    SHOW_ERROR(error_number, error_string)           // prints the error message associated with error
                                                     // number (from the error catalog) prepended by
                                                     // 'ERROR: ', and appends 'error_string'

    SHOW_ERROR_IF(cond, error_number)                // if 'cond' is TRUE, prints the error message
                                                     // associated with error number (from the error
                                                     // catalog) prepended by 'ERROR: '

    SHOW_ERROR_IF(cond, error_number, error_string)  // if 'cond' is TRUE, prints the error message
                                                     // associated with error number (from the error
                                                     // catalog) prepended by 'ERROR: ', and appends 
                                                     // 'error_string'


Warning macros
--------------

The ``WARNING`` macros function in the same way as the ``ERROR`` macros, with the exception that instead of prepending the
message with ':boldtext:`ERROR:` ', the ``WARNING`` macros prepend the message with ':boldtext:`WARNING:` '.

The ``WARNING`` macros are:

::

    SHOW_WARN(error_number)
    SHOW_WARN(error_number, error_string)
    SHOW_WARN_IF(cond, error_number)
    SHOW_WARN_IF(cond, error_number, error_string)


Static macros
-------------

An additional set of macros is provided to be used in static functions and other functions that are not contained within an instantiated 
object (e.g. main()).  

The static ``ERROR`` macros are:

::

   SHOW_ERROR_STATIC(error_number)                       : prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog)
   SHOW_ERROR_STATIC(error_number, error_string)         : prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string"
   SHOW_ERROR_IF_STATIC(cond, error_number)              : if "cond" is TRUE, prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog)
   SHOW_ERROR_IF_STATIC(cond, error_number, error_string): if "cond" is TRUE, prints "ERROR: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string"

The static ``WARNING`` macros are:

::

   SHOW_WARN_STATIC(error_number)                        : prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog)
   SHOW_WARN_STATIC(error_number, error_string)          : prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string"
   SHOW_WARN_IF_STATIC(cond, error_number)               : if "cond" is TRUE, prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog)
   SHOW_WARN_IF_STATIC(cond, error_number, error_string) : if "cond" is TRUE, prints "WARNING: " followed by the error message associated with "error_number" (from the error catalog), and appends "error_string"


Macro output
------------

Error and warning messages printed via the ``SHOW_ERROR(_IF)`` and ``SHOW_WARN(_IF)`` macros will always contain:

::

    The object id of the calling object.
    The object type of the calling object.
    The stellar type of the calling object (will be ”NONE” if the calling object is not a star-type object).
    The function name of the calling function.

Any object that uses the the ``SHOW_ERROR(_IF)`` and ``SHOW_WARN(_IF)`` macros must expose the following functions:

::

    OBJECT_ID ObjectId()** const { return m ObjectId; }
    OBJECT_TYPE ObjectType()** const { return m ObjectType; }
    STELLAR_TYPE StellarType()** const { return m StellarType; }

These functions are called by the ``ERROR`` and ``WARNING`` macros. If any of the functions are not applicable to the object, 
then they must return ':italictext:`type`::NONE' (all objects should implement ObjectId() correctly).

Error and warning messages displayed using the ``SHOW_ERROR(_IF)_STATIC`` and ``SHOW_WARN(_IF)_STATIC`` macros will always contain:

::

   The function name of the calling function

but will not contain:

::

   The object id of the calling object (not available in static functions)
   The object type of the calling object (doesn't add enough information on its own)
   The stellar type of the calling object (not available in static functions)
 

Handling errors at runtime
--------------------------

Early versions of COMPAS (versions prior to v03.00.00) did not have a coherent, robust error-handling strategy. In those versions,
errors were typically displayed as either errors or warnings (depending on the severity) as they occurred, and evolution of the star
(SSE mode) or binary (BSE mode) continued - users were expected to check errors or warnings displayed and use results with appropriate
caution. This was not ideal.

In COMPAS version 03.00.00 the error handling philosophy was changed, and more coherent and robust error-handling code was implemented.
The new error-handling philosophy is to stop evolution of a star or binary if an error occurs, and record in the (SSE/BSE) system
parameters file the fact that an error occurred, and an error number identifying the error that occurred. This way users can check the
system paramers file at the completion of a run for the disposition of a star or binary and, if the evolution of that star or binary was
stopped because an error occurred, the actual error that occurred. Possible dispositions (for both stars and binaries) are given in the
``EVOLUTION_STATUS`` enum class and associate label map in ``typedefs.h``.

The error-handling code implemented in v03.00.00 allows developers to terminate evolution of a star or binary if they determine that a
condition encountered is sufficiently severe that allowing the evolution of the star or binary to continue would produce inconsistent or
untrustable results. In those cases, the developers should terminate the evolution of the star or binary via the use of the ``THROW_ERROR*``
macros (defined in ``ErrorsMacros.h`` and described below).

Developers should use the ``SHOW_WARN*`` macros (defined in ``ErrorsMacros.h`` and described above) to alert users to conditions they want
to bring to the attention of users, but are not sufficiently severe to warrant termination of the evolution of the star or binary.

The ``SHOW_ERROR*`` macros (defined in ``ErrorsMacros.h`` and described above) should be used sparingly - generally only in catch blocks
for errors thrown, or in the (very few) sections of the code not covered by catch blocks.

The class member variable ``m_Error`` (in the ``BaseStar`` class for SSE; ``BaseBinaryStar`` for BSE) should not be set explicitly
throughout the code - it is set in the appropriate error exception catch blocks.  ``m_Error`` is the error value written to the log files.
Note that it is possible that if users choose to add the ``STAR_PROPERTY::ERROR`` (SSE) or ``BINARY_PROPERTY::ERROR`` (BSE) to log files
via the ``logfile-definitions`` file, the value of the error logged to those files may be 0 (``ERROR::NONE``) even for stars or binaries
that ere eventually terminated due to an error - the error value is only set when the error occurs (and is thrown), so some records in some
log files may already have been written prior to the error being identified and evolution terminated.


.. _developer-guide-fp-errors:

Floating-point errors in C++
----------------------------

In C++ implementations that implement the IEEE floating-point standard, in ordinary operation, the division of a finite non-zero floating-point
value by 0[.0] is well-defined and results in ``+infinity`` if the value is greater than zero, ``-infinity`` if the value is less than zero, and
``NaN`` if the value is equal to 0.0, and in each case program execution continues uninterrupted. Integer division by 0 is undefined and results
in a floating-point exception and the process is halted.

The GNU C++ implementation allows us to trap the following floating-point errors:

::

    DIVBYZERO : division by zero, or some other asymptotically infinite result (from finite arguments).
    INEXACT   : a value cannot be represented with exact accuracy (e.g. 0.1, 1.0/3.0, and sqrt(2.0)).
    INVALID   : at least one of the arguments to a floating-point library function is a value for which the function is not defined (e.g. sqrt(-1.0))
    OVERFLOW  : the result of an operation is too large in magnitude to be represented as a value of the return type.
    UNDERFLOW : the result is too small in magnitude to be represented as a value of the return type.

When an enabled floating-point trap is encountered, a ``SIGFPE`` signal is raised. If we don't have a signal handler installed for ``SIGFPE`` the program
is terminated with a floating-point exception. If we do have a signal handler installed for ``SIGFPE``, that signal handler is invoked. Ordinarily, once
the ``SIGFPE`` signal handler is invoked, there is no going back - after doing whetever we need to do to manage the signal, the only valid operations are
to exit the program, or to longjmp to a specific location in the code. Fortunately the GNU C++ designers have given us another option: if we compile with
the ``-fnon-call-exceptions`` compiler flag we can raise an exception safely from the ``SIGFPE`` signal handler, because the throw is just a non-local
transfer of control (just like a longjmp), and then we can just catch the exception raised.


Floating-point errors in COMPAS
-------------------------------

Instrumentation has been implemented in COMPAS v03.00.00 that traps ``DIVBYZERO``, ``INVALID``, ``OVERFLOW``, and ``UNDERFLOW``. Trapping ``INEXACT`` would
mean we'd trap on just about every floating-point calculation (it is really just informational - we know there are many values we can't represent exactly
in base-2).

We have 3 modes for the floating-point error instrumentation:

::

    Instrumentation not active  : This mode is enabled with the option '--fp-error-mode OFF'  (This is the default mode).
    
                                  This mode is just the default behaviour of the C++ compiler, as described above. In this mode,
                                  the program execution will not be interrupted in the event of a floating-point error, but the
                                  error reported in the system parameters file will be set to indicate if a floating-point error
                                  occurred during evolution (and in this mode we can, and do, differentiate between 'DIVBYZERO',
                                  'INVALID', 'OVERFLOW', and 'UNDERFLOW'). Note that an integer divide-by-zero will cause the
                                  execution of the program to halt (and, rather obtusely, will report "Floating point exception").

    Floating-point traps enabled: This mode is enabled with the option '--fp-error-mode ON'.
    
                                  In this mode, floating-point traps 'DIVBYZERO', 'INVALID', 'OVERFLOW', and 'UNDERFLOW' are enabled.
                                  When a floating-point operation traps, a 'SIGFPE' signal is raised and the 'SIGFPE' signal handler
                                  is invoked, and the signal handler raises a 'runtime_error' exception, with a 'what' argument of
                                  "FPE" (and in this mode we cannot, and do not, differentiate between 'DIVBYZERO', 'INVALID', 
                                  'OVERFLOW', and 'UNDERFLOW'). The exception raised will cause the execution of the program to halt
                                  if it is not caught and managed. We catch 'runtime_error' exceptions in 'Star::Evolve()' for SSE
                                  mode, in 'BaseBinaryStar::Evolve()' for BSE mode, and in 'main()' for errors that might occur
                                  outside the evolution of stars or binaries.

    Debug mode                  : This mode is enabled with the option '--fp-error-mode DEBUG'.
    
                                  In this mode, floating-point traps 'DIVBYZERO', 'INVALID', 'OVERFLOW', and 'UNDERFLOW' are enabled.
                                  When a floating-point operation traps, a 'SIGFPE' signal is raised and the 'SIGFPE' signal handler
                                  is called, but instead of raising a 'runtime_error' exception, the signal handler prints the stack
                                  trace that led to the error and halts execution of the program.  In this way, the user can determine
                                  where (to the function level - we do not determine line numbers) the floating-point error occurred.

                                  The construction of the stack trace in debug mode happens inside the signal handler, and the functions
                                  used to do that are generally not signal safe - but we call 'std::exit()' anyway, so that should not
                                  be a problem.
   

Throwing errors
---------------

Another additional set of macros is provided, for both static and non-static functions, that will, after displaying an error (as described above),
throw an exception and cause the ordinary program flow to be interrupted.  These are:

::

   THROW_ERROR(error_number)                              : displays the error (as described above), then throws exception
   THROW_ERROR(error_number, error_string)                : displays the error (as described above), then throws exception
   THROW_ERROR_IF(cond, error_number)                     : if "cond" is TRUE, displays the error (as described above), then throws exception
   THROW_ERROR_IF(cond, error_number, error_string)       : if "cond" is TRUE, displays the error (as described above), then throws exception

   THROW_ERROR_STATIC(error_number)                       : displays the error (as described above), then throws exception
   THROW_ERROR_STATIC(error_number, error_string)         : displays the error (as described above), then throws exception
   THROW_ERROR_IF_STATIC(cond, error_number)              : if "cond" is TRUE, displays the error (as described above), then throws exception
   THROW_ERROR_IF_STATIC(cond, error_number, error_string): if "cond" is TRUE, displays the error (as described above), then throws exception
 
In each case, the exception thrown by the ``THROW*`` macros is the 'error_number' argument cast as an integer, so it can be caught by using 
``catch (int e)`` and inspecting "e".


Writing errors to file
----------------------

When the user sets the ``--errors-to-file`` program option to ``TRUE`` (``FALSE`` by default), errors and, if the user also sets the 
``--enable-warnings`` program option to ``TRUE`` (``FALSE`` by default), warnings will be written to a log file in addition to being
displayed on stdout/stderr as they occur. In this case, the ``Log::Start()`` parameter ``p_ErrorsToFile`` will be ``TRUE``, instructing
the ``LOGGING`` service to write errors and warnings to a log file.

The filename to which error records are written when  is declared in ``LogTypedefs.h`` – see the enum class ``LOGFILE`` and associated
descriptor map ``LOGFILE_DESCRIPTOR``. Currently the name is ``Error_Log``.
