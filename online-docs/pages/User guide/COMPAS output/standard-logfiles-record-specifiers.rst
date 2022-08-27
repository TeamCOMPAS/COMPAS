Standard log file record specifiers
===================================

Each standard log file has an associated log file record specifier that defines what data are to be written to the log files. Each record 
specifier is a list of known properties that are to be written as the log record for the log file associated with the record specifier. 
Default record specifiers for each of the standard log files are shown in :doc:`standard-logfiles-default-record-specifications`. 
The standard log file record specifiers can be defined by the user at run-time (see :doc:`standard-logfiles-record-specification`).

When specifying known properties, the property name must be prefixed with the property type. The current list of valid property types 
available for use is:

    - STAR_PROPERTY
    - STAR_1_PROPERTY
    - STAR_2_PROPERTY
    - SUPERNOVA_PROPERTY
    - COMPANION_PROPERTY
    - BINARY_PROPERTY
    - PROGRAM_OPTION

The stellar property types (all types except BINARY_PROPERTY AND PROGRAM_OPTION) must be paired with properties from the stellar property list, 
the binary property type BINARY_PROPERTY with properties from the binary property list, and the program option type PROGRAM_OPTION with properties 
from the program option property list.

All standard log files, except the switch log files, have a record type property (column). The record type property is an integer that can be used
to idenitify and filter records within a standard log file. The record type property was introduced primarily to support different types of records
in the detailed output files (BSE and SSE), and is currently unly used in those files, but could be useful for other log files. All record types for
all standard log files are printed to the log files by default, but records of specific record types can be enabled and disabled using the
appropriate program option (see e.g. ``--logfile-detailed-output-record-types`` in :doc:`../Program options/program-options-list-defaults` for the 
detailed output files).

