Debugging macros
================

A set of macros similar to the :doc:`logging macros <./Logging macros/services-logging-debugging-macros-logging>` is also
provided for debugging purposes. These macros are analogous to their logging counterparts.

The debugging macros write directly to stdout rather than the log file, but their output can also be written to
the log file if desired (see the ``p_DbgToFile`` parameter of ``Log::Start()``, and the ``--debug-to-file`` program
option).

A major difference between the logging macros and the debugging macros is that the debugging macros can be defined
away. The debugging macro definitions are enclosed in an ``#ifdef`` enclosure, and are only present in the source code
if ``#DEBUG`` is defined. This means that if ``#DEBUG`` is not defined (``#undef``), all debugging statements using
the debugging macros will be removed from the source code by the preprocessor before the source is compiled. Un-defining
``#DEBUG`` not only prevents bloat of unused code in the executable, it improves performance. Many of the functions in
the code are called hundreds of thousands, if not millions, of times as the stellar evolution proceeds. Even if the
debugging classes and debugging level are set so that no debug statement is displayed, just checking the debugging level
every time a function is called increases the run time of the program. The suggested use is to enable the debugging
macros (``#define DEBUG``) while developing new code, and disable them (``#undef DEBUG``) to produce a production version
of the executable.

The debugging macros provided are::

    DBG(...)              // analogous to LOG()
    DBG_ID(...)           // analogous to LOG_ID()
    DBG_IF(cond, ...)     // analogous to LOG_IF()
    DBG_ID_IF(cond, ...)  // analogous to LOG_ID_IF()


Two further debugging macros are provided::

    DBG_WAIT(...)
    DBG_WAIT_IF(cond, ...)

The ``DBG_WAIT`` macros function in the same way as their non-wait counterparts (``DBG(...)`` and ``DBG_IF(cond, ...)``),
with the added functionality that they will pause execution of the program and wait for user input before proceeding.

A set of macros for printing warning message is also provided. These are the ``DBG_WARN`` macros::

    DBG_WARN(...)              // analogous to LOG()
    DBG_WARN_ID(...)           // analogous to LOG_ID()
    DBG_WARN_IF(cond, ...)     // analogous to LOG_IF()
    DBG_WARN_ID_IF(cond, ...)  // analogous to LOG_ID_IF()

The ``DBG_WARN`` macros write to stdout via the ``SAY`` macro, so honour the logging classes and level, and are not written
to the debug or errors files.

Note that the ``id`` parameter of the ``LOG`` macros (to specify the logfileId) is not required for the ``DBG`` macros (the
filename to which debug records are written is declared in ``constants.h`` – see the enum class ``LOGFILE`` and associated 
descriptor map ``LOGFILE_DESCRIPTOR``).
