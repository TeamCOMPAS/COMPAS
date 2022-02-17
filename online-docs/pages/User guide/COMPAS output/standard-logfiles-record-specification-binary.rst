Binary properties
=================

When specifying known properties in a log file record specification record, the property name must be prefixed with 
the property type. Currently there is a single binary property type available for use: BINARY_PROPERTY.

For example, to specify the property ``SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE`` for a binary star being evolved in ``BSE``, use::

    BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE


.. _binary-props-top:

Following is an alphabetical list of binary properties available for inclusion in log file record specifiers.

**Jump to**
:ref:`A <binary-props-A>` :ref:`B <binary-props-B>` :ref:`C <binary-props-C>` :ref:`D <binary-props-D>`
:ref:`E <binary-props-E>` :ref:`F <binary-props-F>` :ref:`G <binary-props-G>` :ref:`H <binary-props-H>`
:ref:`I <binary-props-I>` :ref:`J <binary-props-J>` :ref:`K <binary-props-K>` :ref:`L <binary-props-L>`
:ref:`M <binary-props-M>` :ref:`N <binary-props-N>` :ref:`O <binary-props-O>` :ref:`P <binary-props-P>`
:ref:`Q <binary-props-Q>` :ref:`R <binary-props-R>` :ref:`S <binary-props-S>` :ref:`T <binary-props-T>`
:ref:`U <binary-props-U>` :ref:`V <binary-props-V>` :ref:`W <binary-props-W>` :ref:`X <binary-props-X>`
:ref:`Y <binary-props-Y>` :ref:`Z <binary-props-Z>`


Following is the list of binary properties available for inclusion in log file record specifiers.
Binary Properties


.. _binary-props-A:

.. _binary-props-B:

.. _binary-props-C:

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CIRCULARIZATION_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CircularizationTimescale
   * - Description:
     - Tidal circularisation timescale (Myr)
   * - Header String:
     - Tau_Circ

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_AT_LEAST_ONCE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_CEDetails.CEEcount
   * - Description:
     - Flag to indicate if there has been at least one common envelope event.
   * - Header String:
     - CEE    

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_EVENT_COUNT**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.CEEcount
   * - Description:
     - The number of common envelope events.
   * - Header String:
     - CE_Event_Counter

.. _binary-props-D:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DIMENSIONLESS_KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_uK
   * - Description:
     - Dimensionless kick magnitude supplied by user (see option --fix-dimensionless-kick-magnitude).
   * - Header String:
     - Kick_Magnitude(uK)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DOUBLE_CORE_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.doubleCoreCE
   * - Description:
     - Flag to indicate double-core common envelope.
   * - Header String:
     - Double_Core_CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_Dt
   * - Description:
     - Current timestep (Myr).
   * - Header String:
     - dT

.. _binary-props-E:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_Eccentricity
   * - Description:
     - Orbital eccentricity.
   * - Header String:
     - Eccentricity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_AT_DCO_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_EccentricityAtDCOFormation
   * - Description:
     - Orbital eccentricity at DCO formation.
   * - Header String:
     - Eccentricity@\ DCO

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_INITIAL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_EccentricityInitial
   * - Description:
     - Supplied by user via grid file or sampled from distribution (see ``--eccentricity-distribution`` option) (default).
   * - Header String:
     - Eccentricity@\ ZAMS

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.eccentricity
   * - Description:
     - Eccentricity immediately following common envelope event.
   * - Header String:
     - Eccentricity>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_EccentricityPreSN
   * - Description:
     - Eccentricity of the binary immediately prior to supernova event.
   * - Header String:
     - Eccentricity<SN

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.eccentricity
   * - Description:
     - Eccentricity at the onset of RLOF leading to the CE.
   * - Header String:
     - Eccentricity<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ERROR**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_Error
   * - Description:
     - Error number (if error condition exists, else 0).
   * - Header String:
     - Error

.. _binary-props-F:

.. _binary-props-G:

.. _binary-props-H:

.. _binary-props-I:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ID**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_ObjectId
   * - Description:
     - Unique object identifier for ``C++`` object – used in debugging to identify objects.
   * - Header String:
     - ID

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IMMEDIATE_RLOF_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.immediateRLOFPostCEE
   * - Description:
     - Flag to indicate if either star overflows its Roche lobe immediately following common envelope event.
   * - Header String:
     - Immediate_RLOF>CE

.. _binary-props-J:

.. _binary-props-K:

.. _binary-props-L:

.. _binary-props-M:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.mass
   * - Description:
     - Mass of the primary star immediately following common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(1)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.mass
   * - Description:
     - Mass of the primary star immediately prior to common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(1)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.mass
   * - Description:
     - Mass of the secondary star immediately following common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(2)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.mass
   * - Description:
     - Mass of the secondary star immediately prior to common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(2)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_ENV_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_MassEnv1
   * - Description:
     - Envelope mass of the primary star (\ :math:`M_\odot`).
   * - Header String:
     - Mass_Env(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_ENV_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_MassEnv2
   * - Description:
     - Envelope mass of the secondary star (\ :math:`M_\odot`).
   * - Header String:
     - Mass_Env(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASSES_EQUILIBRATED**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_MassesEquilibrated
   * - Description:
     - Flag to indicate whether chemically homogeneous stars had masses equilibrated and orbit circularised due to Roche lobe overflow during evolution.
   * - Header String:
     - Equilibrated

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASSES_EQUILIBRATED_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_MassesEquilibratedAtBirth
   * - Description:
     - Flag to indicate whether stars had masses equilibrated and orbit circularised at birth due to Roche lobe overflow.
   * - Header String:
     - Equilibrated_At_Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_TRANSFER_TRACKER_HISTORY**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_MassTransferTrackerHistory
   * - Description:
     - Indicator of mass transfer history for the binary. Will be printed as one of:

        .. list-table::
           :widths: 35 5
           :header-rows: 0
           :class: aligned-text
           
           * - NO MASS TRANSFER
             - = 0
           * - STABLE FROM 1 TO 2
             - = 1
           * - STABLE FROM 2 TO 1
             - = 2
           * - CE FROM 1 TO 2
             - = 3
           * - CE FROM 2 TO 1
             - = 4
           * - CE DOUBLE CORE
             - = 5
           * - CE BOTH MS
             - = 6
           * - CE MS WITH CO
             - = 7

   * - Header String:
     - MT_History

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MERGES_IN_HUBBLE_TIME**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_MergesInHubbleTime
   * - Description:
     - Flag to indicate if the binary compact remnants merge within a Hubble time.
   * - Header String:
     - Merges_Hubble_Time

.. _binary-props-N:

.. _binary-props-O:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OPTIMISTIC_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.optimisticCE
   * - Description:
     - Flag that returns TRUE if we have a Hertzsprung-gap star, and we allow it to survive the CE.
   * - Header String:
     - Optimistic_CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_VELOCITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_OrbitalVelocity
   * - Description:
     - Orbital velocity (\ :math:`km s^{-1}`).
   * - Header String:
     - Orbital_Velocity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_VELOCITY_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_OrbitalVelocityPreSN
   * - Description:
     - Orbital velocity immediately prior to supernova event (\ :math:`km s^{-1}`).
   * - Header String:
     - Orbital_Velocity<SN

.. _binary-props-P:

.. _binary-props-Q:

.. _binary-props-R:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.radius
   * - Description:
     - Radius of the primary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - Radius(1)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.radius
   * - Description:
     - Radius of the primary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - Radius(1)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.radius
   * - Description:
     - Radius of the secondary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - Radius(2)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.radius
   * - Description:
     - Radius of the secondary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - Radius(2)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED**
     -
   * - Data type:
     - UNSIGNED LONG
   * - COMPAS variable:
     - BaseBinaryStar::m_RandomSeed
   * - Description:
     - Seed for random number generator for this binary star. Optionally supplied by user via program option ``--random-seed``; default generated from system time.
   * - Header String:
     - SEED

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→isCE
   * - Description:
     - Flag to indicate if the RLOF leads to a common-envelope event 
   * - Header String:
     - CEE>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→eccentricity
   * - Description:
     - Eccentricity immediately after RLOF.
   * - Header String:
     - Eccentricity>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_EVENT_COUNTER**
     -
   * - Data type:
     - UNSIGNED INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→eventCounter
   * - Description:
     - The number of times the binary has overflowed its Roche lobe up to and including this episode
   * - Header String:
     - MT_Event_Counter

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→semiMajorAxis
   * - Description:
     - Semi-major Axis(\ :math:`R_\odot`) immediately after RLOF.
   * - Header String:
     - SemiMajorAxis>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→mass1
   * - Description:
     - Mass (\ :math:`M_\odot`) of the primary immediately after RLOF.
   * - Header String:
     - Mass(1)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→mass2
   * - Description:
     - Mass (\ :math:`M_\odot`) of the secondary immediately after RLOF.
   * - Header String:
     - Mass(2)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→radius1
   * - Description:
     - Radius (\ :math:`R_\odot`) of the primary immediately after RLOF.
   * - Header String:
     - Radius(1)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→radius2
   * - Description:
     - Radius (\ :math:`R_\odot`) of the secondary immediately after RLOF.
   * - Header String:
     - Radius(2)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→isRLOF1
   * - Description:
     - Flag to indicate whether the primary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(1)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→isRLOF2
   * - Description:
     - Flag to indicate whether the secondary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(2)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(1)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(1)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR1_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(2)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(2)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR2_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→time
   * - Description:
     - Time since ZAMS (Myr) immediately after RLOF.
   * - Header String:
     - Time>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→eccentricity
   * - Description:
     - Eccentricity at the onset of RLOF.
   * - Header String:
     - Eccentricity<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→semiMajorAxis
   * - Description:
     - Semi-major Axis (\ :math:`R_\odot`) at the onset of RLOF.
   * - Header String:
     - SemiMajorAxis<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→mass1
   * - Description:
     - Mass (\ :math:`M_\odot`) of the primary at the onset of RLOF.
   * - Header String:
     - Mass(1)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→mass2
   * - Description:
     - Mass (\ :math:`M_\odot`) of the secondary at the onset of RLOF.
   * - Header String:
     - Mass(2)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→radius1
   * - Description:
     - Radius (\ :math:`R_\odot`) of the primary at the onset of RLOF.
   * - Header String:
     - Radius(1)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→radius2
   * - Description:
     - Radius (\ :math:`R_\odot`) of the secondary at the onset of RLOF.
   * - Header String:
     - Radius(2)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→isRLOF1
   * - Description:
     - Flag to indicate whether the primary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(1)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→isRLOF2
   * - Description:
     - Flag to indicate whether the secondary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(2)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(1)<MT

`Note that this property has the same header string as RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(1)<MT

`Note that this property has the same header string as RLOF_PRE_MT_STAR1_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType2
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(2)<MT

`Note that this property has the same header string as RLOF_PRE_MTvSTAR2_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType2
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(2)<MT

`Note that this property has the same header string as RLOF_PRE_MT_STAR2_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→time
   * - Description:
     - Time since ZAMS (Myr) at the onset of RLOF.
   * - Header String:
     - Time<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_RocheLobeRadius
   * - Description:
     - Roche radius of the primary star (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(1)|a

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_RocheLobeRadius
   * - Description:
     - Roche radius of the secondary star (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(2)|a

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.rocheLobe1to2
   * - Description:
     - Roche radius of the primary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(1)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.rocheLobe2to1
   * - Description:
     - Roche radius of the secondary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(2)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.rocheLobe1to2
   * - Description:
     - Roche radius of the primary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(1)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.rocheLobe2to1
   * - Description:
     - Roche radius of the secondary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(2)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STAR_TO_ROCHE_LOBE_RADIUS_RATIO_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Calculated using BinaryConstituentStar::m_StarToRocheLobeRadiusRatio
   * - Description:
     - Ratio of the primary star’s stellar radius to Roche radius (R/RL), evaluated at periapsis.
   * - Header String:
     - Radius(1)|RL

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STAR_TO_ROCHE_LOBE_RADIUS_RATIO_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Calculated using BinaryConstituentStar::m_StarToRocheLobeRadiusRatio
   * - Description:
     - Ratio of the secondary star’s stellar radius to Roche radius (R/RL), evaluated at periapsis.
   * - Header String:
     - Radius(2)|RL

.. _binary-props-S:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_AT_DCO_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxisAtDCOFormation
   * - Description:
     - Semi-major axis at DCO formation (AU).
   * - Header String:
     - SemiMajorAxis@\ DCO

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_INITIAL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxisInitial
   * - Description:
     - Semi-major axis at ZAMS (AU).
   * - Header String:
     - SemiMajorAxis@\ ZAMS

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.semiMajorAxis
   * - Description:
     - Semi-major axis immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxisPreSN
   * - Description:
     - Semi-major axis immediately prior to supernova event (AU).
   * - Header String:
     - SemiMajorAxis<SN

`Note that this property has the same header string as SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_SemiMajorAxisPreSN
   * - Description:
     - Semi-major axis immediately prior to supernova event (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis<SN

`Note that this property has the same header string as SEMI_MAJOR_AXIS_PRE_SUPERNOVA. It is expected that one or the other is printed in any file, but 
not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.semiMajorAxis
   * - Description:
     - Semi-major axis at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxis
   * - Description:
     - Semi-major axis at ZAMS (AU).
   * - Header String:
     - SemiMajorAxis@\ ZAMS

`Note that this property has the same header string as SEMI_MAJOR_AXIS_RSOL. It is expected that one or the other is printed in any file, but not both. 
If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_RSOL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_SemiMajorAxis
   * - Description:
     - Semi-major axis at ZAMS (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis@\ ZAMS

`Note that this property has the same header string as SEMI_MAJOR_AXIS. It is expected that one or the other is printed in any file, but not both. If both 
are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SIMULTANEOUS_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.simultaneousRLOF
   * - Description:
     - Flag to indicate that both stars are undergoing RLOF.
   * - Header String:
     - Simultaneous_RLOF

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STABLE_RLOF_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.stableRLOFPostCEE
   * - Description:
     - Flag to indicate stable mass transfer after common envelope event.
   * - Header String:
     - Stable_RLOF>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_MERGER**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_StellarMerger
   * - Description:
     - Flag to indicate the stars merged (were touching) during evolution.
   * - Header String:
     - Merger

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_MERGER_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_StellarMergerAtBirth
   * - Description:
     - Flag to indicate the stars merged (were touching) at birth.
   * - Header String:
     - Merger_At_Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star immediately following common envelope event.
   * - Header String:
     - Stellar_Type(1)>CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star at the onset of RLOF leading to the common-envelope episode.
   * - Header String:
     - Stellar_Type(1)<CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star immediately following common envelope event.
   * - Header String:
     - Stellar_Type(2)>CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF leading to the common-envelope episode.
   * - Header String:
     - Stellar_Type(2)<CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the primary star immediately following common envelope event. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(1)>CE

`Note that this property has the same header string as STELLAR_TYPE_1_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the primary star at the onset of RLOF leading to the common-envelope episode. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(1)<CE

`Note that this property has the same header string as STELLAR_TYPE_1_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, but not 
both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the secondary star immediately following common envelope event. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(2)>CE

`Note that this property has the same header string as STELLAR_TYPE_2_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, but not 
both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF leading to the common-envelope episode. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(2)<CE

`Note that this property has the same header string as STELLAR_TYPE_2_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, but not 
both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_STATE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_SupernovaState
   * - Description:
     - Indicates which star(s) went supernova. Will be printed as one of:

        .. list-table::
           :widths: 45 5
           :header-rows: 0
           :class: aligned-text

           * - No supernova
             - = 0
           * - Star 1 is the supernova
             - = 1
           * - Star 2 is the supernova
             - = 2
           * - Both stars are supernovae
             - = 3

   * - Header String:
     - Supernova_State

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SYNCHRONIZATION_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SynchronizationTimescale
   * - Description:
     - Tidal synchronisation timescale (Myr).
   * - Header String:
     - Tau_Sync

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SYSTEMIC_SPEED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SystemicVelocity
   * - Description:
     - Post-supernova systemic (centre-of-mass) velocity (\ :math:`km s^{-1}`).
   * - Header String:
     - Systemic_Velocity

.. _binary-props-T:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_Time
   * - Description:
     - Time since ZAMS (Myr).
   * - Header String:
     - Time

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIME_TO_COALESCENCE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_TimeToCoalescence
   * - Description:
     - Time between formation of double compact object and gravitational-wave driven merger (Myr).
   * - Header String:
     - Coalescence_Time

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TOTAL_ANGULAR_MOMENTUM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_TotalAngularMomentum
   * - Description:
     - Total angular momentum calculated using regular conservation of energy (\ :math:`M_\odot  AU^2 y^{r−1}`).
   * - Header String:
     - Ang_Momentum_Total

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TOTAL_ENERGY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_TotalEnergy
   * - Description:
     - Total energy calculated using regular conservation of energy (\ :math:`M_\odot  AU^2 y^{r−1}`).
   * - Header String:
     - Energy_Total

.. _binary-props-U:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **UNBOUND**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_Unbound
   * - Description:
     - Flag to indicate the binary is unbound (or has become unbound after a supernova event).
   * - Header String:
     - Unbound

.. _binary-props-V:

.. _binary-props-W:

.. _binary-props-X:

.. _binary-props-Y:

.. _binary-props-Z:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_LOBE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_ZetaLobe
   * - Description:
     - The logarithmic derivative of Roche lobe radius with respect to donor mass for :math:`q=\frac{Md}{Ma}` at the onset of the RLOF.
   * - Header String:
     - Zeta_Lobe

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_STAR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_ZetaStar
   * - Description:
     - Mass-radius exponent of the star at the onset of the RLOF. Calculated differently based on the value of program option ``--zeta-prescription``
   * - Header String:
     - Zeta_Star

