Standard log file format
========================

COMPAS can produce log files in several formats:

    - Hierarchical Data Format version 5 (``HDF5``)\ [#f1]_
    - Comma Separated Values (``CSV``)
    - Tab Separated Values (``TSV``)
    - Plain text: space separated values (``TXT``)
    
The log file type is set using the ``--logfile-type`` program option.

Standard ``CSV``, ``TSV``, and ``TXT`` log files are human-readable files, and formatted in a similar fashion. Each standard
``CSV``, ``TSV``, and ``TXT`` log file consists of three header records followed by data records. Header records and data records
are delimiter separated fields, and the fields as specified by the log file record specifier.

The header records for all standard ``CSV``, ``TSV``, and ``TXT`` log files are:

    Header record 1: Column Data Type Names |br|
    Header record 2: Column Units (where applicable) |br|
    Header record 3: Column Headings

Column Data Type Names are taken from the set **{ BOOL, INT, FLOAT, STRING }**, where

    .. list-table::
       :widths: 12 88 
       :header-rows: 0
       :class: aligned-text

       * - **BOOL**
         - indicates the data value will be a boolean value.
       * - 
         - Boolean data values will be recorded in the log file in either numerical format (1 or 0, where 1 = TRUE and 0 = FALSE), or string format ("TRUE" or "FALSE"), depending upon the value of the ``--print-bool-as-string`` program option.
       * - **INT**
         - indicates the data value will be an integer number.
       * - **FLOAT**
         - indicates the data value will be a floating-point number.
       * - **STRING**
         - indicates the data value will be a text string.

Column Units is a string indicating the units of the corresponding data values (e.g. "Msol\*AU\ :sup:`2`\ \*yr\ :sup:`-1`\ ",
"Msol", "AU", etc.). The Column Units value may be blank where units are not applicable, or may be one of:

    .. list-table::
       :widths: 12 88 
       :header-rows: 0
       :class: aligned-text

       * - **Count**
         - indicates the data value is the total of a counted entity.
       * - **State**
         - indicates the data value describes a state (e.g. "Unbound" state is "TRUE" or "FALSE").
       * - **Event**
         - the data value describes an event status (e.g. "Simultaneous_RLOF" is "TRUE").

Column Headings are string labels that describe the corresponding data values. The heading strings for stellar properties of
constituent stars of a binary will have appropriate identifiers appended. That is, heading strings for:

    .. list-table::
       :widths: 38 62 
       :header-rows: 0
       :class: aligned-text

       * - STAR_1 PROPERTY::properties
         - will have ":boldtext:`(1)`" appended
       * - STAR_2 PROPERTY::properties
         - will have ":boldtext:`(2)`" appended
       * - SUPERNOVA_PROPERTY::properties
         - will have ":boldtext:`(SN)`" appended: any column with a header with a suffix of ":boldtext:`(SN)`" represents an attribute of the star undergoing a supernova event, either just before the supernova [e.g., `Mass_Total@CO(SN)`] or just after the supernovae [e.g., `Mass(SN)`].
       * - COMPANION_PROPERTY::properties
         - will have ":boldtext:`(CP)`" appended: any column with a header with a suffix of ":boldtext:`(CP)`" represents an attribute of the the companion after the supernova event.

``HDF5`` files are not human-readable. The ``HDF5`` file format supports large, complex, heterogeneous data, enabling the data to be stored
in a structured way in a single file. When the ``HDF5`` format is specified for COMPAS log files, a single ``HDF5`` file is produced for 
non-detailed output log files, containing all non-detailed output log files described above. Detailed output files are created, as for
other logfile types, as individual files (in this case, ``HDF5`` files), in the ’Detailed_Output’ container directory.

Each file described above is created as a `group` within the ``HDF5`` file, with the name of the group set to the name of the file
(e.g. "BSE_System_Parameters"). Each column in the files described above is created as a `dataset` within its corresponding group in the
``HDF5`` file, with the name of the datset set to the column header as described above (e.g. "Mass(1)"). Each dataset in an ``HDF5`` file
is typed, and the dataset data types are set to the column data types as described above. The column units described above are attached to
their corresponding datasets in the ``HDF5`` file as `attributes`.


.. rubric:: Footnotes

.. [#f1] https://www.hdfgroup.org/
