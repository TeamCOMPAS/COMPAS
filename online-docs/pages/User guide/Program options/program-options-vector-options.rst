Vector program options
======================

Most of the program options available in COMPAS allow users to specify a `single` value for the option (e.g. ``--initial-mass-1 7.5``,
or ``--metallicity 0.0142``, etc.), while others allow users to specify multiple values (e.g. ``--log-classes class1 class2 class3``,
``--notes "this is note 1" note2 a_third_note``, etc.). The latter are `vector` options, because they allow users to specify a vector
of values.

Currently, the only vector program options are:

- --debug-classes
- --log-classes
- --notes-hdrs
- --notes


The notation for vector program options provides for the specification of one or more values. e.g.::

    --log-classes class1 class2 class3 ... classN

Option values are separated by a space character. Unless otherwise specified in the documentation for the program option, there is no
limit to the number of values specified.

When using this notation, all options required must be provided: there is no mechanism to allow a default value using the fallback method
for program options (to the command-line value, then to the COMPAS default) - leaving a value blank would be ambiguous (as to which e.g.
`log-class` had been left blank), and specifying an empty string ("") for a value would be ambiguous (as to whether the user wanted the
option value to default, or just be an empty string).

Option values (in general, but also specifically for vector options) may not begin with the dash character ('-'), because the shell parser
will parse them as option names before passing them through to COMPAS.

COMPAS imposes no limit to the length (number of characters) of an individual option values that are specified as strings, but there may 
be practical limits imposed by the underlying system.


Shorthand notation
------------------

Because the notation described above could become awkward, and to allow for default values for vector options, a shorthand notation for
vector options is provided. Usage using the shorthand notation is::

    --debug-classes [class1,headerclass2str2,class3,...,classN]

    --notes [annotation1,annotation2,annotation3,...,annotationN]

Because the parameters are bounded by the brackets, and delimited by commas (and so are now positional), users can omit specific values::

    --notes [,,annotation3,,annotation5]

In the example above, annotations 1, 2, 4, and those beyond annotation 5 have been omitted. Annotations 1, 2 & 4 will default - if they are
specified in this manner on a grid line they will default to the correspodning annotation specified on the command line; if they are specified
in this manner on the command line they will default to the COMPAS default annotation.

Spaces in option values (in general, but also specifically for vector options) strings need to be enclosed in quotes, or the shell parser will
parse them as separate arguments before passing them through to COMPAS.  If the log file type is specified as TXT, then any spaces in option
values need to be enclosed in quotes to avoid the shell parser parsing them as separate arguments, but also need to have enclosing quotes 
propagated to the logfile, or the spaces will be interpreted as delimiters in the log file.  e.g. - in the following example, enclosing escaped
quote characters ('\"') are added before adding the enclosing quotes::

    --notes-hdrs [headerstr1,"\"header str 2\"",headerstr3,...,headerStrN]


The shorthand notation is expanded to the notation described above (the COMPAS code just fills in the omitted values with the required defaults),
so the caveat mentioned above that option values may not begin with the dash character ('-') applies to the shorthand notation.

Shorthand notation is optional: users may choose to use the notation described above rather than shorthand notation, but in that case all option
values must be specified (no omissions, no defaults).

