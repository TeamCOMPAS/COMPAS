SSE switchlog
=============

Default record definition for the SSE SwitchLog log file::

    const ANY_PROPERTY_VECTOR SSE_SWITCH_LOG_REC = {
        STAR_PROPERTY::RANDOM_SEED,
        STAR_PROPERTY::TIME
    };


The default record specification can be modified at runtime via a logfile record specifications file (program option ``--logfile-definitions``).
See :doc:`standard-logfiles-record-specification` for details.

Note that the SSE SwitchLog file has the following columns automatically appended to each record:

    - The stellar type from which the star is switching.
    - The stellar type to which the star is switching.

|br|
**SWITCHING_FROM**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` ``BaseStar::m_StellarType``
   * - Description:
     - The stellar type of the star immediately prior to the switch.
   * - Header String:
     - "SWITCHING_FROM"

**SWITCHING_TO**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - Not applicable
   * - Description:
     - The stellar type to which the star will switch (i.e. the stellar type immediately following the switch).
   * - Header String:
     - "SWITCHING_TO"

These columns will always be automatically appended to each SSE Switch Log record: they cannot be removed via the log file record 
specifications file.
