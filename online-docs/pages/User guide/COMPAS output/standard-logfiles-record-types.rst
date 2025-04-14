Standard log file record types
==============================

All standard COMPAS logfiles, except the switch log files, have a record type property (column).  The record type property
is an unsigned integer that can take any value in the range :math:`0..4294967295`.

The record type is used to differentiate records within a standard log file. The functionality was introduced primarily to
support different types of records in the detailed output files (BSE and SSE), but has been extended to all other log files
(except the switch log files).

The use case for the detailed output log files is to differentiate between records written to the file when the star or binary
and constituent stars are known to be in self-consistent states (that is, the attributes of the star or binary and constituent
stars have been correctly and completely updated), and records written to the file, perhaps mid-timestep, when the star or
binary and/or constituent stars may not be self-consistent.  Since the record type property can take any value in the range 
:math:`0..4294967295` there is scope to identify many different events or situations, in any of the standard log files.  We may,
for example, use different record types to indicate that a detailed output record was written immediately prior to, or immediately
following, a particular event or calculation.  Or we may want to indicate that a supernovae record was written to the supernovae
file prior to the SN event, and another record following the SN event.  Or for the RLOF file we may want to differentiate between
records written pre-MT and post_MT - the possibilities are (almost) limitless.

Each standard log file has its own set of record types - select a file below to see the available record types for that file.

.. toctree::
   :maxdepth: 1
   
   standard-logfiles-record-types-sse-system-parameters
   standard-logfiles-record-types-sse-supernovae
   standard-logfiles-record-types-sse-detailed-output
   standard-logfiles-record-types-sse-pulsars
   standard-logfiles-record-types-sse
   standard-logfiles-record-types-sse-switchlog

   standard-logfiles-record-types-bse-system-parameters
   standard-logfiles-record-types-bse-supernovae
   standard-logfiles-record-types-bse-rlof
   standard-logfiles-record-types-bse-pulsar-evolution
   standard-logfiles-record-types-bse-double-compact-objects
   standard-logfiles-record-types-bse-common-envelopes
   standard-logfiles-record-types-bse-detailed-output
   standard-logfiles-record-types-bse-switchlog

Since the record type property is an unsigned integer, filtering the output files by record type is very simple.  Even so, users may
want to disable the logging of some record types - perhaps to limit the size of the log files produced.  For this reason, program options
are provided to enable/disable the logging of different record types for log files that have a record type property.  These options are:

    - --logfile-system-parameters-record-types
    - --logfile-supernovae-record-types
    - --logfile-detailed-output-record-types
    - --logfile-common-envelopes-record-types
    - --logfile-rlof-parameters-record-types
    - --logfile-double-compact-objects-record-types
    - --logfile-pulsar-evolution-record-types

For these options, the record types to be enabled are specified as a bitmap, with each bit corresponding to a record type. To construct
the bitmap, for each record type to be enabled, raise 2 to the power of (record type - 1), then sum the results - the sum is the bitmap,
and the integer value to be entered for this option.

Example:

To enable record types :math:`1`, :math:`4`, and :math:`9`, the option value should be

.. math::

   2^{(1 - 1)} + 2^{(4 - 1)} + 2^{(9 - 1)} = 2^0 + 2^3 + 2^8 = 1 + 8 + 256 = 265

:math:`265`` as a binary number is written as 0100001001, with the 1st, 4th, and 9th bits enabled (counting 1-based from the least-significant
bit being the right-most), corresponding to the record types :math:`1`, :math:`4`, and :math:`9` being enabled, and all other record types disabled.

A value of -1 for the bitmap is shorthand for all bits enabled - all record types enabled.
