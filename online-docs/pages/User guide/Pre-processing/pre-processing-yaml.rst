COMPAS YAML file
================

The COMPAS YAML file is an options configuration file for use with :doc:`runSubmit.py<./pre-processing-runSubmit>`.

`YAML files <https://en.wikipedia.org/wiki/YAML>`__ are text files that use a minimal syntax, with Python-style
indentation to indicate nesting, and a colon-centered syntax for expressing key-value pairs.

The COMPAS YAML file contains entries, including **key-value** pairs, for all COMPAS options, where the **key** is the
COMPAS option name (string), and the **value** is the required option value. A COMPAS option entry in the YAML file
also indicates (as comments) the COMPAS default value for the option and, where applicable, the allowed values for the
option.

The default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has all COMPAS option entries commented so
that the COMPAS default value for the option is used by default. To use a value other than the COMPAS default value,
users must uncomment the entry and change the option value to the desired value.

Users can either use the default YAML file provided with the distribution (``compasConfigDefault.yaml``), or create
a custom YAML file.


Creating a custom YAML file
---------------------------

There are two ways to create a custom YAML file for COMPAS:

* modify the default YAML file ``compasConfigDefault.yaml``
* create, and optionally modify, a new YAML file via COMPAS


Creating a YAML file via COMPAS
-------------------------------

Use the COMPAS option ``--create-YAML-file`` to create a new YAML file via COMPAS.  The ``--create-YAML-file``
option requires an argument specifying the name of the YAML file to be created - the file created will be with
the name as provided, including any file extension.

For example, the following commands will create YAML files named ``myYAML.yaml`` and ``newYAML`` respectively:

::

    ./COMPAS --create-YAML-file myYAML.yaml
    ./COMPAS --create-YAML-file newYAML

When COMPAS creates a YAML file, all option values in the resultant YAML file will be the COMPAS default values,
except for options that the user specified on the command line when the YAML file was created.  Furthermore, all
lines in the YAML file specifying a COMPAS option will be commented if the option value is the COMPAS default
value, and not commented if the option value was specified on the command line by the user. This allows users to
create project-specific YAML files if they wish.
 
For example, the following command will create a YAML file with all option values set to the COMPAS default, and
with all lines specifying options commented:

::

    ./COMPAS --create-YAML-file myYAML.yaml

The following command, however,  will create a YAML file with all option values set to the COMPAS default, except
for the ``--logfile-type`` option, which will be set to  ``'csv'``, and all lines specifying options will be
commented except for the line specifying the ``--logfile-type`` option, which will not be commented:

::

    ./COMPAS --create-YAML-file myYAML.yaml --logfile-type csv

If the file specified by the ``--create-YAML-file`` option already exists, the user will be asked if they wish the
existing file to be overwritten.



YAML file format
~~~~~~~~~~~~~~~~

The format of the YAML file created by COMPAS is determined by a template - either the COMPAS default template
(defined in the header file ``yaml.h``), or provided by the user via the ``--YAML-template`` option. An existing
COMPAS YAML file can be used as a template.

The following command will create a YAML file formatted according to the template contained within the file
``myYAMLtemplate.yaml``:

::

    ./COMPAS --create-YAML-file myYAML.yaml --YAML-template myYAMLtemplate.yaml

If the ``--YAML-template`` option is not specified, or if it is and the file specified does not exist or is not
readable, the default COMPAS template will be used.


YAML template rules
~~~~~~~~~~~~~~~~~~~

COMPAS uses the following rules (listed in no particular order) when it creates a YAML file:


    - The following two records will be automatically written to the start of YAML file:

          - ##~!!~## COMPAS option values
          - ##~!!~## Created at ddd MMM DD HH:MM:SS YYYY by COMPAS vxx.yy.zz
    - Lines in the template beginning with *"##~!!~##"* will not be preserved (these are assumed to be COMPAS generated headers, and will be rewritten by COMPAS).
    - Leading *'#'* characters on option definition lines in the template will not be preserved (but they may be rewritten by COMPAS).
    - Option comments in the template must be preceded by *"# "* or they will not be preserved.
    - Strings in the template beginning with *"# Default: "* and up to (but not including) the next *'#'* (or end of line if no *'#'*) will not be preserved.
    - Strings in the template beginning with *"# Options: "* and up to (but not including) the next *'#'* (or end of line if no *'#'*) will not be preserved.
    - Blank lines in the template will be preserved.
    - Option values in the template will not be preserved (but they may be rewritten by COMPAS).
    - Option values written by COMPAS will be the option default values unless COMPAS was run with command-line options set - if the user executed COMPAS and specified options on the command line, the user-specified values will be written to the YAML file, and those option records in the YAML file will not be commented.
    - Options present in the template that are not valid COMPAS options will be ignored and not written to the YAML file.
    - Any COMPAS options that are not present in the template will be written in alphabetical order at the end of the YAML file.

In the following example template:

::

    0001     ##~!!~## COMPAS option values
    0002     ##~!!~## File Created Tue Feb 14 13:09:06 2023 by COMPAS v02.34.06
    0003
    0004     # first comment
    0005
    0006     booleanChoices:
    0007         ### BINARY PROPERTIES
    0008     #    --allow-touching-at-birth          # Default: False                                        # second comment
    0009
    0010         ### STELLAR PROPERTIES
    0011         --mass-loss-prescription: 'HURLEY'  # Default: 'VINK'  # Options: ['VINK','HURLEY','NONE']    third comment

- Lines 0001 and 0002 will not be preserved (but will be replaced by new COMPAS headers).
- The blank line at line 0003 will be preserved.
- The comment *"first comment"* (on line 0004) will be preserved.
- The blank line at line 0005 will be preserved.
- The header *"booleanChoices:"* on line 0006 will be preserved.
- The header *"### BINARY PROPERTIES"* on line 0007 will be preserved.
- The leading *'#'* on line 0008 will not be preserved (but may be rewritten by COMPAS if the option is set to default).
- The string beginning with *"# Default: "* and extending to the next *'#'* on line 0008 will not be preserved (but will be replaced by COMPAS).
- The comment *"second comment"* on line 0008 will be preserved.
- The blank line at line 0009 will be preserved.
- The header *"### STELLAR PROPERTIES"* on line 0010 will be preserved.
- The string beginning with *"# Default: "* and extending to the next *'#'* on line 0011 will not be preserved (but will be replaced by COMPAS).
- The string beginning with *"# Options: "* and extending to the next *'#'* (or, in this case because there is no subsequent *'#'*, the end of the line) on line 0011 will not be preserved (but will be replaced by COMPAS).
- The comment *"third comment"* on line 0011 will not be preserved - there is no *"# "* prefix, so it will be subsumed by the *"# Options: "* string (which extends from *"# Options: "* to the end of the line).
