BSE switchlog
=============

Default record definition for the BSE SwitchLog log file::

    const ANY_PROPERTY_VECTOR BSE_SWITCH_LOG_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        BINARY_PROPERTY::TIME
    };


The default record specification can be modified at runtime via a logfile record specifications file (program option ``--logfile-definitions``).
See :doc:`standard-logfiles-record-specification` for details.

Note that the BSE SwitchLog file has the following columns automatically appended to each record:

    - The constituent star switching stellar type: 1 = Primary, 2 = Secondary.
    - The stellar type from which the star is switching.
    - The stellar type to which the star is switching.

|br|
**STAR_SWITCHING**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` ``BaseBinaryStar::m_Star1/m_Star2``
   * - Description:
     - The constituent star switching stellar type, where 1 = Primary, and 2 = Secondary.
   * - Header String:
     - "STAR_SWITCHING"

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
     - The stellar type of the constituent star immediately prior to the switch.
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
     - The stellar type to which the constituent star will switch (i.e. the stellar type immediately following the switch).
   * - Header String:
     - "SWITCHING_TO"

These columns will always be automatically appended to each BSE Switch Log record: they cannot be removed via the log file record 
specifications file.
