Base-level logging
==================

The Log class member variables are private, and public functions have been created for logging and debugging functionality
required by the code.

The Log service can be accessed by referring to the Log::Instance() object. For example, to check if the logging service is
enabled, call the ``Log::Enabled()`` function::

    Log::Instance()→Enabled();

Since that could become unwieldy, there is a convenience macro to access the Log service. The macro just defines “LOGGING”
as “Log::Instance()”, so calling the ``Log::Enabled()`` function can be written as::

    LOGGING→Enabled();

The Log service must be initialised and started before logging and debugging functionality can be used. Initialise and start 
logging by calling the ``Log::Start()`` function::

    LOGGING→Start(...)

Refer to the description of the ``Log::Start()`` function below for parameter definitions.

The Log service should be stopped before exiting the program – this ensures all open log files are flushed to disk and closed
properly. Stop logging by calling the ``Log::Stop()`` function::

    LOGGING→Stop(...)

Refer to the description of the ``Log::Stop()`` function below for parameter definitions.


Log & debug record filtering
----------------------------

The Log service provides a set of functions and macros to manage log files, and to write log and debug records to the log files,
``stdout``, and ``stderr``. Base-level logging allows developers to tag log and debug records with a string ``class``, and an 
integer ``level``. The Log service will filter log and debug records by ``class`` and ``level``, and only write those records 
that meet the ``class`` and ``level`` filters specified by the users via the ``--log-level``,  ``--log-class``, ``--debug-level``,
and ``--debug-classes`` program options.


Log service functions:
----------------------

The Log service provides the following public member functions:


.. toctree::
   :maxdepth: 1

   ./services-base-level-logging-func-start
   ./services-base-level-logging-func-stop
   ./services-base-level-logging-func-open
   ./services-base-level-logging-func-close
   ./services-base-level-logging-func-write
   ./services-base-level-logging-func-put
   ./services-base-level-logging-func-debug
   ./services-base-level-logging-func-debugWait
   ./services-base-level-logging-func-error
   ./services-base-level-logging-func-squawk
   ./services-base-level-logging-func-say
   ./services-base-level-logging-func-enabled

