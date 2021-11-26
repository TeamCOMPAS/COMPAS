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
