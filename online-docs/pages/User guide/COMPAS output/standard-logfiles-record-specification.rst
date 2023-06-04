Standard log file record specification
======================================

The standard log file record specifiers can be changed at run-time by supplying a log file definitions file via the ``--logfile-definitions`` 
program option. This allows users to change what appears in any of the standard COMPAS log files without the need to change the source code
or rebuild the COMPAS executable.

The syntax of the definitions file is fairly simple. The definitions file is expected to contain zero or more log file record specifications, as 
explained below.

Complete lists of available properties selectable for inclusion (or exclusion) from COMPAS log file are available:


.. toctree::
   :maxdepth: 1

   standard-logfiles-record-specification-stellar
   standard-logfiles-record-specification-binary
   standard-logfiles-record-specification-options


For the following specification::

     ::=   means `expands to` or `is defined as`
    { x }  means (possible) repetition: x may appear zero or more times
    [ x ]  means x is optional: x may appear, or not
    <name> is a term (expression)
    "abc"  means literal string "abc"
      n    means integer number
      |    means `OR`
      #    indicates the start of a comment


Log file definitions file specification
---------------------------------------

::

    <def_file>   ::= {<rec_spec>}

    <rec_spec>   ::= <rec_name> <op> "{" { [ <props_list> ] } "}" <spec_delim>

    <rec_name>   ::= "SSE_SYSPARMS_REC"    |   # SSE only
                     "SSE_DETAILED_REC"    |   # SSE only
                     "SSE_SNE_REC"         |   # SSE only
                     "SSE_SWITCH_REC"      |   # SSE only
                     "BSE_SYSPARMS_REC"    |   # BSE only
                     "BSE_SWITCH_REC"      |   # BSE only
                     "BSE_DCO_REC"         |   # BSE only
                     "BSE_SNE_REC"         |   # BSE only
                     "BSE_CEE_REC"         |   # BSE only
                     "BSE_PULSARS_REC"     |   # BSE only
                     "BSE_RLOF_REC"        |   # BSE only
                     "BSE_DETAILED_REC"    |   # BSE only
   
    <op>         ::= "=" | "+=" | "-="

    <props_list> ::= <prop_spec> [ <props_delim> <props_list> ]

    <prop_spec>  ::= <prop_type> "::" <prop_name> [ <prop_index> ] <prop_delim>

    <spec_delim> ::= "::" | "EOL"

    <prop_delim> ::= "," | <spec_delim>

    <prop_type>  ::= "STAR_PROPERTY"       |   # SSE only
                     "STAR1_PROPERTY"      |   # BSE only
                     "STAR2_PROPERTY"      |   # BSE only
                     "SUPERNOVA_PROPERTY"  |   # BSE only
                     "COMPANION_PROPERTY"  |   # BSE only
                     "BINARY_PROPERTY"     |   # BSE only
                     "PROGRAM_OPTION"      |   # SSE_or BSE

    <prop_index> ::= "[" n "]"

    <prop_name>  ::= valid property name for specified property type (see definitions in constants.h)


The log file definitions may contain comments. Comments are denoted by the hash/pound character ("#"). The hash/pound character and any 
text following it on the line in which the hash character appears is ignored by the parser. The hash/pound character can appear anywhere 
on a line - if it is the first character then the entire line is a comment and ignored by the parser, or it can follow valid symbols on 
a line, in which case the symbols before the hash/pound character are parsed and interpreted by the parser.

A log file specification record is initially set to its default value (see :doc:`standard-logfiles-default-record-specifications`). 
The log file definitions file informs COMPAS as to the modifications to the default values the user wants. This means that the log file 
definitions log file is not mandatory, and if the log file definitions file is not present, or contains no valid record specifiers, the log 
file record definitions will remain at their default values.

The assignment operator given in a record specification (``<op>`` in the specification above) can be one of "=", "+=", and "-=".  The meanings of 
these are:

    "=" means the record specifier should be assigned the list of properties specified in the braced-list following the "=" operator. The value of the record specifier prior to the assignment is discarded, and the new value set as described.
    
    "+=" means the list of properties specified in the braced-list following the "+=" operator should be appended to the existing value of the record specifier. Note that the new properties are appended to the existing list, so will appear at the end of the list (properties are printed in the order they appear in the list).
    
    "-=" means the list of properties specified in the braced-list following the "-=" operator should be subtracted from the existing value of the record specifier.

The `prop-index` qualifier may only be used for vector program options.


Example definitions file entries
--------------------------------

Some example log file definitions file entries are::

    SSE_SYSPARMS_REC = { STAR_PROPERTY::RANDOM_SEED, STAR_PROPERTY::RADIUS,
                         STAR_PROPERTY::MASS,
                         STAR_PROPERTY::LUMINOSITY
                       }

    BSE_PULSARS_REC += { STAR_1_PROPERTY::LUMINOSITY, STAR_2_PROPERTY::CORE_MASS, 
                         BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL, COMPANION PROPERTY::RADIUS
                       }

    BSE_PULSARS_REC -= { SUPERNOVA_PROPERTY::TEMPERATURE }

    BSE_PULSARS_REC += { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
                         BINARY_PROPERTY::ORBITAL_VELOCITY
                       }

    BSE_SYSPARMS_REC -= { PROGRAM_OPTION::NOTES }
    BSE_SYSPARMS_REC += { PROGRAM_OPTION::NOTES[1], PROGRAM_OPTION::NOTES[3] }
    
A full example log file record specifications file is shown in :doc:`standard-logfiles-example-definitions-file`.

The record specifications in the definitions file are processed individually in the sequence they appear in the file, and are cumulative: 
for record specifications pertaining to the same record name, the output of earlier specifications is input to later specifications.

For each record specification:

- Properties requested to be added to an existing record specification that already exist in that record specification are ignored. Properties will not appear in a record specification twice.
- Properties requested to be subtracted from an existing record specification that do not exist in that record specification are ignored.

Note that neither of those circumstances will cause a parse error for the definitions file – in both cases the user’s intent is satisfied.
