Log::Stop(p_ObjectStats)
========================

::

    VOID Log::Stop(
        TUPLE<INT, INT> p_ObjectStats  // number of stars or binaries requested, count created.
    )

Stops the logging and debugging service. All open log files are flushed to disk and closed (including any Standard Log Files
open - see description of Standard Log Files in :doc:`../services-extended-logging`. The ``Run_Details`` file is populated and closed.

This function does not return a value.
