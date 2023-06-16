Log file annotations
====================

COMPAS provides functionality that allows users to annotate standard log files.

The original motivation for log file annotation functionality was to enable users to describe the contents of custom grid files,
but annotations can be used for any reason.  With the ability to annotate the log files, users can indicate the origin of various
input data - e.g. the user could indicate what IMF was used to draw initial mass values, or what distribution was used to draw 
mass ratio (q), etc.  

Annotations are written to log files as columns of data (aka datasets in ``HDF5`` files). Annotations are specified via program
options, and so can be specified on the command line as well as in grid files. Because annotations are written to log files
as columns of data, functionality is provided to allow the user to specify header strings for the columns, as well as the
actual data for the columns (the annotations).

Two program options are provided:

``--notes-hdrs``

Allows users to specify header strings for annotation columns.  Can only be specified on the command line: not valid in grid
files, and will be ignored, and elicit a warning, if included.

The ``notes-hdrs`` program option is a vector program option, and provides for the specification of one or more annotation
header strings.  Usage is::

    --notes-hdrs headerstr1 headerstr2 headerstr3 ... headerstrN

Header strings are separated by a space character. There is no limit to the number of header strings specified, and the number
specified defines the maximum number of annotations allowed.


``--notes``

Allows users to specify annotations.  Can be specified on the command line and in a grid file.

The ``notes`` program option is a vector program option, and provides for the specification of one or more annotation strings.
Usage is::

    --notes annotation1 annotation2 annotation3 ... annotationN

Annotation strings are separated by a space character.  The number of annotation strings is limited to the number of annotation
header strings specified (via the ``--notes-hdrs`` program option).  If more annotation strings are specified than header strings,
the excess annotation strings will be ignored (and a warning displayed).

When using this notation, all annotation strings must be provided: there is no mechanism to allow a default annotation using the
fallback method for program options (to the command-line value, then to the COMPAS default) - leaving an annotation string blank
would be ambiguous (as to which annotation string had been left blank), and specifying an empty string ("") as an annotation
string would be ambiguous (as to whether the user wanted the annotation string to default, or just be an empty string).

Neither header strings, nor annotation strings, may begin with the dash character ('-'), because the shell parser will parse them
as option names before passing them through to COMPAS.

COMPAS imposes no limit to the length of an individual annotation header string or annotation string (number of characters), but
there may be practical limits imposed by the underlying system. 


Using vector option shorthand notation
--------------------------------------

Because the notation described above could become awkward, and to allow for default annotations, a shorthand notation for vector
program options is provided (see :doc:`../Program options/program-options-vector-options` for details).

Usage using the shorthand notation is::

    --notes-hdrs [headerstr1,headerstr2,headerstr3,...,headerStrN]

    --notes [annotation1,annotation2,annotation3,...,annotationN]

Because the parameters are bounded by the brackets, and delimited by commas (and so are now positional), users can omit specific
annotations::

    --notes [,,annotation3,,annotation5]

In the example above, annotations 1, 2, 4, and those beyond annotation 5 have been omitted. Annotations 1, 2 & 4 will default - if
they are specified in this manner on a grid line they will default to the correspodning annotation specified on the command line; if
they are specified in this manner on the command line they will default to the COMPAS default annotation (the empty string).  In the
example above, if the number of annotations expected, as defined by the number of annotation headers specified via the ``--notes-hdrs``
program option is more than 5, then annotations beyond annotation 5 (the last annotation actually specified by the user) will default
in the same manner as described above.

Spaces in annotation header strings and annotation strings need to be enclosed in quotes, or the shell parser will parse them as 
separate arguments before passing them through to COMPAS.  If the log file type is specified as TXT, then any spaces in annotation
header strings and annotation strings need to be enclosed in quotes to avoid the shell parser parsing them as separate arguments, but
also need to have enclosing quotes propagated to the logfile, or the spaces will be interpreted as delimiters in the log file - in this
case, the user will need to add enclosing escaped quote characters ('\"') before adding the enclosing quotes.  e.g.::

    --notes-hdrs [headerstr1,"\"header str 2\"",headerstr3,...,headerStrN]

The shorthand notation is expanded to the notation described above (the COMPAS code just fills in the omitted annotation strings with
the required defaults), so the caveat mentioned above that neither header strings, nor annotation strings, may begin with the dash 
character ('-') applies to the shorthand notation.

Shorthand notation is optional: users may choose to use the notation described above rather than shorthand notation, but in that case
all annotation strings must be specified (no omissions, no defaults).


Including annotations in log files
----------------------------------

Annotations can be included in the record specifiers for any of the standard log files in the same way that other program options
are included: they can be included in the default log file record specifiers in ``constants.h`` (compile-time configuration: see
:doc:`../../Developer guide/constants-dot-h`), and/or they can be added to, or removed from, the default log file record specifiers
via the use of a log file definitions file (run-time configuration: see :doc:`./standard-logfiles-record-specification`).


Compile-time configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

The property specifier ``PROGRAM_OPTION::NOTES`` is available to be included in any of the default log file record specifiers in
``constants.h``.  If the property specifier ``PROGRAM_OPTION::NOTES`` is included in a default log file record specifier:

    - a column will be included in the log file for each of the annotation header strings specified by the user (via the ``notes-hdrs``
      program option), with each column having the respective header string specified by the user
    - each record in the log file will include an annotation string for each annotation column, the contents of which is determined by
      what the user specified using the ``--notes`` program option.

Adding the property specifier ``PROGRAM_OPTION::NOTES`` to the default record specifier in ``constants.h`` for a log file causes
*all* annotation columns (as specified via the ``notes-hdrs`` program option) to be included in the log file: including annotations
in a log file via compile-time configuration is an all-or-nothing proposition.


Run-time configuration
~~~~~~~~~~~~~~~~~~~~~~

COMPAS provides functionality to allow users to change which properties are to be written to the standard log files at run-time:
see :doc:`./standard-logfiles-record-specification`. The property specifier ``PROGRAM_OPTION::NOTES`` can be added to, or removed from,
any of the log file record specifiers by the use of this functionality.

Furthermore, because at run-time the number of annotation columns is known (information not known at compile-time), the 
``PROGRAM_OPTION::NOTES`` property specifier can (optionally) be indexed to allow the specification of a particular annotation column.
Thus, ``PROGRAM_OPTION::NOTES`` (with no index) indicates *all* annotations columns, whereas ``PROGRAM_OPTION::NOTES[2]`` indicates 
annotation column 2 (1-based: the first annotation column is indicated by ``PROGRAM_OPTION::NOTES[1]``). By using the optional index, 
users can add specific annotations columns to, or remove them from, any of the log files.

For example, this log file definitions file entry::

    SSE_SYSPARMS_REC -= { PROGRAM_OPTION::NOTES }

removes all annotations columns from the SSE System Parameters log file.

This entry::

    BSE_PULSARS_REC += { PROGRAM_OPTION::NOTES[1] }

adds annotations column 1 to the BSE Pulsar Evolution log file.

These entries::

    BSE_SYSPARMS_REC -= { PROGRAM_OPTION::NOTES }
    BSE_SYSPARMS_REC += { PROGRAM_OPTION::NOTES[1], PROGRAM_OPTION::NOTES[3] }

    BSE_SNE_REC += { PROGRAM_OPTION::NOTES[2], PROGRAM_OPTION::NOTES[3],  PROGRAM_OPTION::NOTES[7] }

    BSE_RLOF_REC -= { PROGRAM_OPTION::NOTES[5] }

    BSE_CEE_REC += { PROGRAM_OPTION::NOTES }

will:

    - remove all annotations columns from the BSE System Parameters log file
    - add annotations columns 1 and 3 to the BSE System Parameters log file
    - add annotations columns 2, 3, and 7 to the BSE Supernovae log file
    - remove annotations column 5 from the BSE RLOF log file
    - add all annotation columns to the BSE Common envelopes files

Specifying an index value less than 1 or greater than the number of annotation headers specified will result in an error.

See :doc:`./standard-logfiles-record-specification` for more details.

|br|
The property ``PROGRAM_OPTION::NOTES`` is included in the default record specifier in ``constants.h`` for both the SSE System Parameters
log file, and the BSE System Parameters log file.

Note that whichever configuration method is used to include annotations in log files, if no annotation headers are specified via the 
``notes-hdrs`` program option, no annotations will be included in any log file.


Log file space considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Annotations are strings of (effectively) unlimited length, so if several annotations are recorded in a log file for each system evolved,
the storage space required for the log file will grow accordingly. Sometimes annotations will apply to the population rather than to each
individual system, and in those cases it may not be necessary to record the same annotation string for possibly millions of systems in the
log file(s). In such situations, it is likely that the annotation headers and strings will be specified on the command line rather than
for each system in a grid file, and so will be recorded in the ``Run_Details`` file - so removing the ``PROGRAM_OPTION::NOTES`` property
from the log file record specifications will prevent them from being repeated needlessly in the log files, and they can be retrieved as
required from the ``Run_Details`` file.

