Logging macros
==============

The logging macros provided by the :doc:`logging service <../services-logging-debugging>` are:


.. toctree::
   :maxdepth: 1

   ./services-logging-debugging-macros-log
   ./services-logging-debugging-macros-log_id
   ./services-logging-debugging-macros-log_if
   ./services-logging-debugging-macros-log_id_if

|br|
The logging macros described above are also provided in a verbose variant. The verbose macros function the same way as their non-verbose
counterparts, with the added functionality that the log records written to the log file will also be written to stdout. The verbose 
logging macros are::

    LOGV(id, ...)
    LOGV_ID(id, ...)
    LOGV_IF(id, cond, ...)
    LOGV_ID_IF(id, cond, ...)


|br|
A further four macros are provided that allow writing directly to stdout rather than a log file. These are::

    SAY( ...)
    SAY_ID( ...)
    SAY_IF(cond, ...)
    SAY_ID IF(cond, ...)

The ``SAY`` macros function the same way as their ``LOG`` counterparts, but write directly to stdout instead of a log file. The ``SAY`` 
macros honour the logging classes and level.
