Stellar properties
==================

When specifying known properties in a log file record specification record, the property name must be prefixed with 
the property type. The current list of valid stellar property types available for use is:

    - STAR_PROPERTY for all stars for ``SSE``
    - STAR_1_PROPERTY for the primary star of a binary for ``BSE``
    - STAR_2_PROPERTY for the secondary star of a binary for ``BSE``
    - SUPERNOVA_PROPERTY for the exploding star in a supernova event for ``BSE``
    - COMPANION_PROPERTY for the companion star in a supernova event for ``BSE``

For example, to specify the property ``TEMPERATURE`` for an individual star being evolved for ``SSE``, use::

    STAR_PROPERTY::TEMPERATURE


.. _stellar-props-top:

Following is an alphabetical list of stellar properties available for inclusion in log file record specifiers.

**Jump to**
:ref:`A <stellar-props-A>` :ref:`B <stellar-props-B>` :ref:`C <stellar-props-C>` :ref:`D <stellar-props-D>`
:ref:`E <stellar-props-E>` :ref:`F <stellar-props-F>` :ref:`G <stellar-props-G>` :ref:`H <stellar-props-H>`
:ref:`I <stellar-props-I>` :ref:`J <stellar-props-J>` :ref:`K <stellar-props-K>` :ref:`L <stellar-props-L>`
:ref:`M <stellar-props-M>` :ref:`N <stellar-props-N>` :ref:`O <stellar-props-O>` :ref:`P <stellar-props-P>`
:ref:`Q <stellar-props-Q>` :ref:`R <stellar-props-R>` :ref:`S <stellar-props-S>` :ref:`T <stellar-props-T>`
:ref:`U <stellar-props-U>` :ref:`V <stellar-props-V>` :ref:`W <stellar-props-W>` :ref:`X <stellar-props-X>`
:ref:`Y <stellar-props-Y>` :ref:`Z <stellar-props-Z>` |_| |_| |_| |_| :ref:`Supernova events/states <supernova-events-states>`

.. _stellar-props-A:

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **AGE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Age
   * - Description:
     - Effective age (changes with mass loss/gain) (Myr)
   * - Header Strings:
     - Age, Age(1), Age(2), Age(SN), Age(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ANGULAR_MOMENTUM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_AngularMomentum
   * - Description:
     - Angular momentum (\ :math:`M_\odot AU^2 yr^{−1}`)
   * - Header Strings:
     - Ang_Momentum, Ang_Momentum(1), Ang_Momentum(2), Ang_Momentum(SN), Ang_Momentum(CP)

.. _stellar-props-B:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.bindingEnergy
   * - Description:
     - Absolute value of the envelope binding energy at the onset of unstable RLOF (erg). Used for calculating post-CE separation.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     -  Binding_Energy@CE(1), Binding_Energy@CE(2), Binding_Energy@CE(SN), Binding_Energy@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_FIXED**
     -
   * - Data type:
     -  DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.fixed
   * - Description:
     - Absolute value of the envelope binding energy calculated using a fixed lambda parameter (erg). Calculated using lambda = m_Lambdas.fixed.
   * - Header Strings:
     - BE_Fixed, BE_Fixed(1), BE_Fixed(2), BE_Fixed(SN), BE_Fixed(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_KRUCKOW**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.kruckow
   * - Description:
     - Absolute value of the envelope binding energy calculated using the fit by :cite:`Vigna-Gomez2018` to :cite:`Kruckow2016` (erg). Calculated using alpha = OPTIONS→CommonEnvelopeSlopeKruckow().
   * - Header Strings:
     - BE_Kruckow, BE_Kruckow(1), BE_Kruckow(2), BE_Kruckow(SN), BE_Kruckow(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_LOVERIDGE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.loveridge
   * - Description:
     - Absolute value of the envelope binding energy calculated as per :cite:`Loveridge2011` (erg). Calculated using lambda = m_Lambdas.loveridge.
   * - Header Strings:
     - BE_Loveridge, BE_Loveridge(1), BE_Loveridge(2), BE_Loveridge(SN), BE_Loveridge(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_LOVERIDGE_WINDS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.loveridgeWinds
   * - Description:
     - Absolute value of the envelope binding energy calculated as per :cite:`Webbink1984` & :cite:`Loveridge2011` including winds (erg). Calculated using lambda = m_Lambdas.loveridgeWinds.
   * - Header Strings:
     - BE_Loveridge_Winds, BE_Loveridge_Winds(1), BE_Loveridge_Winds(2), BE_Loveridge_Winds(SN), BE_Loveridge_Winds(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_NANJING**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.nanjing
   * - Description:
     - Absolute value of the envelope binding energy calculated as per :doc:`Xu & Li (2010) <../../references>` (erg). Calculated using lambda = m_Lambdas.nanjing.
   * - Header Strings:
     - BE_Nanjing, BE_Nanjing(1), BE_Nanjing(2), BE_Nanjing(SN), BE_Nanjing(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.bindingEnergy
   * - Description:
     - Absolute value of the binding energy immediately after CE (erg).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Binding_Energy>CE(1), Binding_Energy>CE(2), Binding_Energy>CE(SN), Binding_Energy>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.bindingEnergy
   * - Description:
     - Absolute value of the binding energy at the onset of unstable RLOF leading to the CE (erg). 
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Binding_Energy<CE(1), Binding_Energy<CE(2), Binding_Energy<CE(SN), Binding_Energy<CE(CP)

.. _stellar-props-C:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable: BaseStar::m_CHE
   * - Description:
     - Flag to indicate whether the star evolved as a ``CH`` star for its entire MS lifetime.
   * -
     - TRUE indicates star evolved as ``CH`` star for entire MS lifetime.
   * -
     - FALSE indicates star spun down and switched from ``CH`` to a ``MS_gt_07`` star.
   * - Header Strings:
     - CH_on_MS, CH_on_MS(1), CH_on_MS(2), CH_on_MS(SN), CH_on_MS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CO_CORE_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_COCoreMass
   * - Description:
     - Carbon-Oxygen core mass (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_CO_Core, Mass_CO_Core(1), Mass_CO_Core(2), Mass_CO_Core(SN), Mass_CO_Core(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CO_CORE_MASS_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.COCoreMass
   * - Description:
     - Carbon-Oxygen core mass at the onset of unstable RLOF leading to the CE (\ :math:`M\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Mass_CO_Core@CE(1), Mass_CO_Core@CE(2), Mass_CO_Core@CE(SN), Mass_CO_Core@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.COCoreMassAtCOFormation
   * - Description:
     - Carbon-Oxygen core mass immediately prior to a supernova (\ :math:`M\odot`).
   * - Header Strings:
     - Mass CO_Core@\ CO, Mass_CO_Core@CO(1), Mass_CO_Core@CO(2), Mass_CO_Core@CO(SN), Mass_CO_Core@CO(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CORE_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_CoreMass
   * - Description:
     - Core mass (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_Core, Mass_Core(1), Mass_Core(2), Mass_Core(SN), Mass_Core(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CORE_MASS_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.CoreMass
   * - Description:
     - Core mass at the onset of unstable RLOF leading to the CE (\ :math:`M\odot`)
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Mass_Core@CE(1), Mass_Core@CE(2), Mass_Core@CE(SN), Mass_Core@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CORE_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.CoreMassAtCOFormation
   * - Description:
     - Core mass immediately prior to a supernova (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_Core@\ CO, Mass_Core@CO(1), Mass_Core@CO(2), Mass_Core@CO(SN), Mass_Core@CO(CP)

.. _stellar-props-D:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DRAWN_KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.drawnKickMagnitude
   * - Description:
     - Magnitude of natal kick without accounting for fallback (\ :math:`km s^{−1}`). Supplied by user in grid file or drawn from distribution (default). This value is used to calculate the actual kick magnitude.
   * - Header Strings:
     - Drawn_Kick_Magnitude, Drawn_Kick_Magnitude(1), Drawn_Kick_Magnitude(2), Drawn_Kick_Magnitude(SN), Drawn_Kick_Magnitude(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DOMINANT_MASS_LOSS_RATE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_DMLR
   * - Description:
     - Current dominant mass loss rate printed as one of

        .. list-table::
           :widths: 85 15
           :header-rows: 0
           :class: aligned-text

           * - None 
             - = 0
           * - Nieuwenhuijzen and de Jager 
             - = 1
           * - Kudritzki and Reimers 
             - = 2
           * - Vassiliadis and Wood 
             - = 3
           * - Wolf-Rayet-like (Hamann, Koesterke and de Koter) 
             - = 4
           * - Vink 
             - = 5
           * - Luminous Blue Variable 
             - = 6

   * - Header Strings:
     - Dominant_Mass_Loss_Rate, Dominant_Mass_Loss_Rate(1), Dominant_Mass_Loss_Rate(2), Dominant_Mass_Loss_Rate(SN), Dominant_Mass_Loss_Rate(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Dt
   * - Description: 
     - Current timestep (Myr).
   * - Header Strings: 
     - dT, dT(1), dT(2), dT(SN), dT(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DYNAMICAL_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_DynamicalTimescale
   * - Description:
     - Dynamical time (Myr).
   * - Header Strings:
     - Tau_Dynamical, Tau_Dynamical(1), Tau_Dynamical(2), Tau_Dynamical(SN), Tau_Dynamical(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.dynamicalTimescale
   * - Description:
     - Dynamical time immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Dynamical>CE(1), Tau_Dynamical>CE(2), Tau_Dynamical>CE(SN), Tau_Dynamical>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.dynamicalTimescale
   * - Description:
     - Dynamical timescale immediately prior to common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Dynamical<CE(1), Tau_Dynamical<CE(2), Tau_Dynamical<CE(SN), Tau_Dynamical<CE(CP)

.. _stellar-props-E:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRIC_ANOMALY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.eccentricAnomaly
   * - Description:
     - Eccentric anomaly calculated using Kepler’s equation.
   * - Header Strings:
     - Eccentric_Anomaly, Eccentric_Anomaly(1), Eccentric_Anomaly(2), Eccentric_Anomaly(SN), Eccentric_Anomaly(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ENV_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_EnvMass
   * - Description:
     - Envelope mass calculated using :cite:`Hurley2000` (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_Env, Mass_Env(1), Mass_Env(2), Mass_Env(SN), Mass_Env(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ERROR**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseStar::m_Error
   * - Description:
     - Error number (if error condition exists, else 0).
   * - Header Strings:
     - Error, Error(1), Error(2), Error(SN), Error(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_AIC**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star underwent an accretion-induced collapse at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_AIC, Experienced_AIC(1), Experienced_AIC(2), Experienced_AIC(SN), Experienced_AIC(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_CCSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as a core-collapse supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_CCSN, Experienced_CCSN(1), Experienced_CCSN(2), Experienced_CCSN(SN), Experienced_CCSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_ECSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as an electron-capture supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_ECSN, Experienced_ECSN(1), Experienced_ECSN(2), Experienced_ECSN(SN), Experienced_ECSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_PISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as an pair-instability supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_PISN, Experienced_PISN(1), Experienced_PISN(2), Experienced_PISN(SN), Experienced_PISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_PPISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as a pulsational pair-instability supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_PPISN, Experienced_PPISN(1), Experienced_PPISN(2), Experienced_PPISN(SN), Experienced_PPISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BinaryConstituentStar::m_RLOFDetails.experiencedRLOF
   * - Description:
     - Flag to indicate whether the star has overflowed its Roche Lobe at any time prior to the current timestep.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Experienced_RLOF(1), Experienced_RLOF(2), Experienced_RLOF(SN), Experienced_RLOF(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_SN_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - The type of supernova event experienced by the star prior to the current timestep. Printed as one of

        .. list-table::
           :widths: 10 5
           :header-rows: 0
           :class: aligned-text

           * - NONE
             - = 0
           * - CCSN
             - = 1
           * - ECSN
             - = 2
           * - PISN
             - = 4
           * - PPISN
             - = 8
           * - USSN
             - = 16
           * - AIC
             - = 32

   * -
     - (see :ref:`Supernova events/states <supernova-events-states>` for explanation).

   * - Header Strings:
     - Experienced_SN_Type, Experienced_SN_Type(1), Experienced_SN_Type(2), Experienced_SN_Type(SN), Experienced_SN_Type(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_USSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as an ultra-stripped supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_USSN, Experienced_USSN(1), Experienced_USSN(2), Experienced_USSN(SN), Experienced_USSN(CP)

.. _stellar-props-F:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **FALLBACK_FRACTION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.fallbackFraction
   * - Description:
     - Fallback fraction during a supernova.
   * - Header Strings:
     - Fallback_Fraction, Fallback_Fraction(1), Fallback_Fraction(2), Fallback_Fraction(SN), Fallback_Fraction(CP)

.. _stellar-props-G:

.. _stellar-props-H:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **HE_CORE_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_HeCoreMass
   * - Description:
     - Helium core mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_He_Core, Mass_He_Core(1), Mass_He_Core(2), Mass_He_Core(SN), Mass_He_Core(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **HE_CORE_MASS_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.HeCoreMass
   * - Description:
     - Helium core mass at the onset of unstable RLOF leading to the CE (\ :math:`M_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Mass_He_Core@CE(1), Mass_He_Core@CE(2), Mass_He_Core@CE(SN), Mass_He_Core@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.HeCoreMassAtCOFormation
   * - Description:
     - Helium core mass immediately prior to a supernova (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_He_Core@\ CO, Mass_He_Core@CO(1), Mass_He_Core@CO(2), Mass_He_Core@CO(SN), Mass_He_Core@CO(CP)

.. _stellar-props-I:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ID**
     -
   * - Data type:
     - UNSIGNED LONG INT
   * - COMPAS variable:
     - BaseStar::m_ObjectId
   * - Description:
     - Unique object identifier for ``C++`` object – used in debugging to identify objects.
   * - Header Strings:
     - ID, ID(1), ID(2), ID(SN), ID(CP)

`Note that this property has the same header string as BINARY_PROPERTY::ID & BINARY_PROPERTY::RLOF_CURRENT_ID. It is expected that one or 
the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_StellarType
   * - Description:
     - Stellar type at zero age main-sequence (per :cite:`Hurley2000`).
   * - Header Strings:
     - Stellar_Type@\ ZAMS, Stellar_Type@ZAMS(1), Stellar_Type@ZAMS(2), Stellar_Type@ZAMS(SN), Stellar_Type@ZAMS(CP)

`Note that this property has the same header string as INITIAL_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any 
file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseStar::m_StellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) at zero age main-sequence. e.g. "First Giant Branch", "Core Helium Burning", "Helium White Dwarf", etc.
   * - Header Strings:
     - Stellar_Type@\ ZAMS, Stellar_Type@ZAMS(1), Stellar_Type@ZAMS(2), Stellar_Type@ZAMS(SN), Stellar_Type@ZAMS(CP)

`Note that this property has the same header string as INITIAL_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_AIC**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently undergoing accretion-induced collapse
   * - Header Strings:
     - AIC, AIC(1), AIC(2), AIC(SN), AIC(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_CCSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently a core-collapse supernova.
   * - Header Strings:
     - CCSN, CCSN(1), CCSN(2), CCSN(SN), CCSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_ECSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently an electron-capture supernova.
   * - Header Strings:
     - ECSN, ECSN(1), ECSN(2), ECSN(SN), ECSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_HYDROGEN_POOR**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.isHydrogenPoor
   * - Description:
     - Flag to indicate if the star is hydrogen poor.
   * - Header Strings:
     - Is_Hydrogen_Poor, Is_Hydrogen_Poor(1), Is_Hydrogen_Poor(2), Is_Hydrogen_Poor(SN), Is_Hydrogen_Poor(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_PISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently a pair-instability supernova.
   * - Header Strings:
     - PISN, PISN(1), PISN(2), PISN(SN), PISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_PPISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently a pulsational pair-instability supernova.
   * - Header Strings:
     - PPISN, PPISN(1), PPISN(2), PPISN(SN), PPISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_RLOFDetails.isRLOF
   * - Description:
     - Flag to indicate whether the star is currently undergoing Roche Lobe overflow.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - RLOF(1), RLOF(2), RLOF(SN), RLOF(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_USSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently an ultra-stripped supernova.
   * - Header Strings:
     - USSN, USSN(1), USSN(2), USSN(SN), USSN(CP)

.. _stellar-props-J:

.. _stellar-props-K:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.kickMagnitude
   * - Description:
     - Magnitude of natal kick received during a supernova (\ :math:`km s^{−1}`). Calculated using the drawn kick magnitude.
   * - Header Strings:
     - Applied_Kick_Magnitude, Applied_Kick_Magnitude(1), Applied_Kick_Magnitude(2), Applied_Kick_Magnitude(SN), Applied_Kick_Magnitude(CP)

.. _stellar-props-L:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.lambda
   * - Description:
     - Common-envelope lambda parameter calculated at the unstable RLOF leading to the CE.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Lambda@CE(1), Lambda@CE(2), Lambda@CE(SN), Lambda@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_DEWI**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.dewi
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Dewi2000` using the fit from Appendix A of :doc:`Claeys et al. (2014) <../../references>`.
   * - Header Strings:
     - Dewi, Dewi(1), Dewi(2), Dewi(SN), Dewi(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_FIXED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.fixed
   * - Description:
     - Universal common envelope lambda parameter specified by the user (program option ``--common-envelope-lambda``).
   * - Header Strings:
     - Lambda_Fixed, Lambda_Fixed(1), Lambda_Fixed(2), Lambda_Fixed(SN), Lambda_Fixed(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckow
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set by program option ``--common-envelope-slope-Kruckow``. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow, Kruckow(1), Kruckow(2), Kruckow(SN), Kruckow(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW_BOTTOM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckowBottom
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set to −1. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow_Bottom, Kruckow_Bottom(1), Kruckow_Bottom(2), Kruckow_Bottom(SN), Kruckow_Bottom(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW_MIDDLE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckowMiddle
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set to :math:`-\frac{4}{5}`. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow_Middle, Kruckow_Middle(1), Kruckow_Middle(2), Kruckow_Middle(SN), Kruckow_Middle(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW_TOP**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckowTop
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set to :math:`-\frac{2}{3}`. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow_Top, Kruckow_Top(1), Kruckow_Top(2), Kruckow_Top(SN), Kruckow_Top(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_LOVERIDGE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.loveridge
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Webbink1984` & :cite:`Loveridge2011`.
   * - Header Strings:
     - Loveridge, Loveridge(1), Loveridge(2), Loveridge(SN), Loveridge(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_LOVERIDGE_WINDS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.loveridgeWinds
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Webbink1984` & :cite:`Loveridge2011` including winds.
   * - Header Strings:
     - Loveridge_Winds, Loveridge_Winds(1), Loveridge_Winds(2), Loveridge_Winds(SN), Loveridge_Winds(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_NANJING**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.nanjing
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :doc:`Xu & Li (2010) <../../references>`.
   * - Header Strings:
     - Lambda_Nanjing, Lambda_Nanjing(1), Lambda_Nanjing(2), Lambda_Nanjing(SN), Lambda_Nanjing(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LBV_PHASE_FLAG**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseStar::m_LBVphaseFlag
   * - Description:
     - Flag to indicate if the star ever entered the luminous blue variable phase.
   * - Header Strings:
     - LBV_Phase_Flag, LBV_Phase_Flag(1), LBV_Phase_Flag(2), LBV_Phase_Flag(SN), LBV_Phase_Flag(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LUMINOSITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Luminosity
   * - Description:
     - Luminosity (\ :math:`L_\odot`).
   * - Header Strings:
     - Luminosity, Luminosity(1), Luminosity(2), Luminosity(SN), Luminosity(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LUMINOSITY_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.luminosity
   * - Description:
     - Luminosity immediately following common envelope event (\ :math:`L_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Luminosity>CE(1), Luminosity>CE(2), Luminosity>CE(SN), Luminosity>CE(CP)


.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LUMINOSITY_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.luminosity
   * - Description:
     - Luminosity at the onset of unstable RLOF leading to the CE (\ :math:`L_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Luminosity<CE(1), Luminosity<CE(2), Luminosity<CE(SN), Luminosity<CE(CP)

.. _stellar-props-M:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Mass
   * - Description:
     - Mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass, Mass(1), Mass(2), Mass(SN), Mass(CP)

`Note that this property has the same header string as RLOF_CURRENT_STAR1_MASS & RLOF_CURRENT_STAR2_MASS. It is expected that one or 
the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_0**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Mass0
   * - Description:
     - Effective initial mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_0, Mass_0(1), Mass_0(2), Mass_0(SN), Mass_0(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_LOSS_DIFF**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_MassLossDiff
   * - Description:
     - The amount of mass lost due to winds (\ :math:`M_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - dmWinds(1), dmWinds(2), dmWinds(SN), dmWinds(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_TRANSFER_DIFF**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_MassTransferDiff
   * - Description:
     - The amount of mass accreted or donated during a mass transfer episode (\ :math:`M_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - dmMT(1), dmMT(2), dmMT(SN), dmMT(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MDOT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Mdot
   * - Description:
     - Mass loss rate (\ :math:`M_\odot yr^{−1}`).
   * - Header Strings:
     - Mdot, Mdot(1), Mdot(2), Mdot(SN), Mdot(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MEAN_ANOMALY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.meanAnomaly
   * - Description:
     - Mean anomaly of supernova kick. Supplied by user in ``grid`` file, default = random number drawn from [0..2π).
   * -
     - See https://en.wikipedia.org/wiki/Mean_anomaly for explanation.
   * - Header Strings:
     - SN Kick Mean Anomaly, SN_Kick_Mean_Anomaly(1), SN_Kick_Mean_Anomaly(2), SN_Kick_Mean_Anomaly(SN), SN_Kick_Mean_Anomaly(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Metallicity
   * - Description:
     - ZAMS Metallicity.
   * - Header Strings:
     - Metallicity@\ ZAMS, Metallicity@ZAMS(1), Metallicity@ZAMS(2), Metallicity@ZAMS(SN), Metallicity@ZAMS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_DONOR_HIST**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - BaseStar::m_MassTransferDonorHistory
   * - Description:
     - A list of all of the stellar types from which the current star was a Mass Transfer donor. This can be readily converted into the different cases of Mass Transfer, depending on the working definition. The output string is formatted as #-#-#... where each # represents a stellar type, in chronological order. E.g, 2-8 means the star was an MT donor as a ``HG`` (stellar type 2) star, and later as a ``HeHG`` (stellar type 8) star.
   * - Header Strings:
     - MT_Donor_Hist, MT_Donor_Hist(1), MT_Donor_Hist(2), MT_Donor_Hist(SN), MT_Donor_Hist(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MZAMS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_MZAMS
   * - Description:
     - ZAMS Mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass@\ ZAMS, Mass@ZAMS(1), Mass@ZAMS(2), Mass@ZAMS(SN), Mass@ZAMS(CP)

.. _stellar-props-N:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NUCLEAR_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_NuclearTimescale
   * - Description:
     - Nuclear timescale (Myr).
   * - Header Strings:
     - Tau_Nuclear, Tau_Nuclear(1), Tau_Nuclear(2), Tau_Nuclear(SN), Tau_Nuclear(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NUCLEAR_TIMESCALE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.nuclearTimescale
   * - Description:
     - Nuclear timescale immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau Nuclear>CE(1), Tau Nuclear>CE(2), Tau Nuclear>CE(SN), Tau Nuclear>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.nuclearTimescale
   * - Description:
     - Nuclear timescale at the onset of unstable RLOF leading to the CE (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Nuclear<CE(1), Tau_Nuclear<CE(2), Tau_Nuclear<CE(SN), Tau_Nuclear<CE(CP)

.. _stellar-props-O:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OMEGA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Omega
   * - Description:
     - Angular frequency (\ :math:`yr^{−1}`).
   * - Header Strings:
     - Omega, Omega(1), Omega(2), Omega(SN), Omega(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OMEGA_BREAK**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_OmegaBreak
   * - Description:
     - Break-up angular frequency (\ :math:`yr^{−1}`).
   * - Header Strings:
     - Omega_Break, Omega_Break(1), Omega_Break(2), Omega_Break(SN), Omega_Break(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OMEGA_ZAMS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_OmegaZAMS
   * - Description:
     - Angular frequency at ZAMS (\ :math:`yr^{−1}`).
   * - Header Strings:
     - Omega@\ ZAMS, Omega@ZAMS(1), Omega@ZAMS(2), Omega@ZAMS(SN), Omega@ZAMS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_ENERGY_POST_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_PostSNeOrbitalEnergy
   * - Description:
     - Absolute value of orbital energy immediately following supernova event (\ :math:`M_\odot AU^2 yr^{−2}`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Orbital_Energy>SN(1), Orbital_Energy>SN(2), Orbital_Energy>SN(SN), Orbital_Energy>SN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_ENERGY_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_PreSNeOrbitalEnergy
   * - Description:
     - Orbital energy immediately prior to supernova event (\ :math:`M_\odot AU^2 yr^{−2}`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Orbital_Energy<SN(1), Orbital_Energy<SN(2), Orbital_Energy<SN(SN), Orbital_Energy<SN(CP)

.. _stellar-props-P:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.magneticField
   * - Description:
     - Pulsar magnetic field strength (G).
   * - Header Strings:
     - Pulsar_Mag_Field, Pulsar_Mag_Field(1), Pulsar_Mag_Field(2), Pulsar_Mag_Field(SN), Pulsar_Mag_Field(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_SPIN_DOWN_RATE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.spinDownRate
   * - Description:
     - Pulsar spin-down rate.
   * - Header Strings:
     - Pulsar_Spin_Down, Pulsar_Spin_Down(1), Pulsar_Spin_Down(2), Pulsar_Spin_Down(SN), Pulsar_Spin_Down(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_SPIN_FREQUENCY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.spinFrequency
   * - Description:
     - Pulsar spin angular frequency (\ :math:`rads s^{−1}`).
   * - Header Strings:
     - Pulsar_Spin_Freq, Pulsar_Spin_Freq(1), Pulsar_Spin_Freq(2), Pulsar_Spin_Freq(SN), Pulsar_Spin_Freq(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_SPIN_PERIOD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.spinPeriod
   * - Description:
     - Pulsar spin period (ms).
   * - Header Strings:
     - Pulsar_Spin_Period, Pulsar_Spin_Period(1), Pulsar_Spin_Period(2), Pulsar_Spin_Period(SN), Pulsar_Spin_Period(CP)

.. _stellar-props-Q:

.. _stellar-props-R:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIAL_EXPANSION_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_RadialExpansionTimescale
   * - Description:
     - e-folding time of stellar radius (Myr).
   * - Header Strings:
     - Tau_Radial, Tau_Radial(1), Tau_Radial(2), Tau_Radial(SN), Tau_Radial(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.radialExpansionTimescale
   * - Description:
     - e-folding time of stellar radius immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Radial<CE(1), Tau_Radial<CE(2), Tau_Radial<CE(SN), Tau_Radial<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.radialExpansionTimescale
   * - Description:
     - e-folding time of stellar radius at the onset of unstable RLOF leading to the CE (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Radial<CE(1), Tau_Radial<CE(2), Tau_Radial<CE(SN), Tau_Radial<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Radius
   * - Description:
     - Radius (\ :math:`R_\odot`).
   * - Header Strings:
     - Radius, Radius(1), Radius(2), Radius(SN), Radius(CP)

`Note that this property has the same header string as RLOF_CURRENT_STAR1_RADIUS & RLOF_CURRENT_STAR2_RADIUS. It is expected that one or 
the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_RandomSeed
   * - Description:
     - Seed for random number generator for this star.
   * - Header Strings:
     - SEED, SEED(1), SEED(2), SEED(SN), SEED(CP)

`Note that this property has the same header string as BINARY_PROPERTY::RANDOM_SEED & BINARY_PROPERTY::RLOF_CURRENT_RANDOM_SEED. It is 
expected that one or the other is printed in any file, but not both. If both are printed then the file will contain two columns with the 
same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RECYCLED_NEUTRON_STAR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the object was a recycled neutron star at any time prior to the current timestep (was a neutron star accreting mass).
   * - Header Strings:
     - Recycled_NS, Recycled_NS(1), Recycled_NS(2), Recycled_NS(SN), Recycled_NS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_ONTO_NS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star transferred mass to a neutron star at any time prior to the current timestep.
   * - Header Strings:
     - RLOF->NS, RLOF->NS(1), RLOF->NS(2), RLOF->NS(SN), RLOF->NS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RUNAWAY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star was unbound by a supernova event at any time prior to the current timestep. (i.e Unbound after supernova event and not a WD, NS, BH or MR).
   * - Header Strings:
     - Runaway, Runaway(1), Runaway(2), Runaway(SN), Runaway(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RZAMS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_RZAMS
   * - Description:
     - ZAMS Radius (\ :math:`R_\odot`).
   * - Header Strings:
     - R@\ ZAMS, R@ZAMS(1), R@ZAMS(2), R@ZAMS(SN), R@ZAMS(CP)

.. _stellar-props-S:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SN_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - The type of supernova event currently being experienced by the star. Printed as one of

        .. list-table::
           :widths: 20 10
           :header-rows: 0
           :class: aligned-text
    
           * - NONE
             - = 0
           * - CCSN
             - = 1
           * - ECSN
             - = 2
           * - PISN
             - = 4
           * - PPISN
             - = 8
           * - USSN
             - = 16
           * - AIC
             - = 32
   * -
     - (see :ref:`Supernova events/states <supernova-events-states>` for explanation).
   * - Header Strings:
     - SN_Type, SN_Type(1), SN_Type(2), SN_Type(SN), SN_Type(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_StellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`).
   * - Header Strings:
     - Stellar_Type, Stellar_Type(1), Stellar_Type(2), Stellar_Type(SN), Stellar_Type(CP)

`Note that this property has the same header string as STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseStar::m_StellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`). e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header Strings:
     - Stellar_Type, Stellar_Type(1), Stellar_Type(2), Stellar_Type(SN), Stellar_Type(CP)

`Note that this property has the same header string as STELLAR_TYPE. It is expected that one or the other is printed in any file, but 
not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_PREV**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_StellarTypePrev
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) at previous timestep.
   * - Header Strings:
     - Stellar_Type_Prev, Stellar_Type_Prev(1), Stellar_Type_Prev(2), Stellar_Type_Prev(SN), Stellar_Type_Prev(CP)

`Note that this property has the same header string as STELLAR_TYPE_PREV_NAME. It is expected that one or the other is printed in any 
file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_PREV_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseStar::m_StellarTypePrev
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) at previous timestep. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header Strings:
     - Stellar_Type_Prev, Stellar_Type_Prev(1), Stellar_Type_Prev(2), Stellar_Type_Prev(SN), Stellar_Type_Prev(CP)

`Note that this property has the same header string as STELLAR_TYPE_PREV. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_KICK_MAGNITUDE_MAGNITUDE_RANDOM_NUMBER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.kickMagnitudeRandom
   * - Description:
     - Random number for drawing the supernova kick magnitude (if required). Supplied by user in grid file, default = random number drawn from [0..1).
   * - Header Strings:
     - SN_Kick_Magnitude_Random_Number, SN_Kick_Magnitude_Random_Number(1), SN_Kick_Magnitude_Random_Number(2), SN_Kick_Magnitude_Random_Number(SN), SN_Kick_Magnitude_Random_Number(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_PHI**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.phi
   * - Description:
     - Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad). Supplied by user in grid file, default = random number drawn from [0..2π).
   * - Header Strings:
     - SN_Kick_Phi, SN_Kick_Phi(1), SN_Kick_Phi(2), SN_Kick_Phi(SN), SN_Kick_Phi(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_THETA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.theta
   * - Description:
     - Angle between the orbital plane and the ’z’ axis of supernovae vector (rad). Supplied by user in grid file, default = drawn from distribution specified by program option ``--kick direction``.
   * - Header Strings:
     - SN_Kick_Theta, SN_Kick_Theta(1), SN_Kick_Theta(2), SN_Kick_Theta(SN), SN_Kick_Theta(CP)

.. _stellar-props-T:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TEMPERATURE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Temperature
   * - Description:
     - Effective temperature (K).
   * - Header Strings:
     - Teff, Teff(1), Teff(2), Teff(SN), Teff(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TEMPERATURE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.temperature
   * - Description:
     - Effective temperature immediately following common envelope event (K).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Teff>CE(1), Teff>CE(2), Teff>CE(SN), Teff>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TEMPERATURE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.temperature
   * - Description:
     - Effective temperature at the unstable RLOF leading to the CE (K).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Teff<CE(1), Teff<CE(2), Teff<CE(SN), Teff<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **THERMAL_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_ThermalTimescale
   * - Description: 
     - Thermal timescale (Myr).
   * - Header Strings:
     - Tau_Thermal, Tau_Thermal(1), Tau_Thermal(2), Tau_Thermal(SN), Tau_Thermal(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **THERMAL_TIMESCALE_POST_COMMON ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.thermalTimescale
   * - Description:
     - Thermal timescale immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Thermal>CE(1), Tau_Thermal>CE(2), Tau_Thermal>CE(SN), Tau_Thermal>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.thermalTimescale
   * - Description:
     - Thermal timescale at the onset of the unstable RLOF leading to the CE (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Thermal<CE(1), Tau_Thermal<CE(2), Tau_Thermal<CE(SN), Tau_Thermal<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Time
   * - Description:
     - Time since ZAMS (Myr).
   * - Header Strings:
     - Time, Time(1), Time(2), Time(SN), Time(CP)

`Note that this property has the same header string as BINARY_PROPERTY::TIME & BINARY_PROPERTY::RLOF_CURRENT_TIME. It is expected that one 
or the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIMESCALE_MS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Timescales[tMS]
   * - Description:
     - Main Sequence timescale (Myr).
   * - Header Strings: tMS, tMS(1), tMS(2), tMS(SN), tMS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.totalMassAtCOFormation
   * - Description:
     - Total mass of the star at the beginning of a supernova event (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_Total@\ CO, Mass_Total@CO(1), Mass_Total@CO(2), Mass_Total@CO(SN), Mass_Total@CO(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TRUE_ANOMALY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.trueAnomaly
   * - Description:
     - True anomaly calculated using Kepler’s equation (rad).
   * -
     - See https://en.wikipedia.org/wiki/True anomaly for explanation.
   * - Header Strings:
     - True_Anomaly(psi), True_Anomaly(psi)(1), True_Anomaly(psi)(2), True_Anomaly(psi)(SN), True_Anomaly(psi)(CP)

.. _stellar-props-U:

.. _stellar-props-V:

.. _stellar-props-W:

.. _stellar-props-X:

.. _stellar-props-Y:

.. _stellar-props-Z:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_HURLEY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.hurley
   * - Description:
     - Adiabatic exponent calculated per :cite:`Hurley2000` using core mass.
   * - Header Strings:
     - Zeta_Hurley, Zeta_Hurley(1), Zeta_Hurley(2), Zeta_Hurley(SN), Zeta_Hurley(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_HURLEY_HE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.hurleyHe
   * - Description:
     - Adiabatic exponent calculated per :cite:`Hurley2000` using He core mass.
   * - Header Strings:
     - Zeta_Hurley_He, Zeta_Hurley_He(1), Zeta_Hurley_He(2), Zeta_Hurley_He(SN), Zeta_Hurley_He(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_SOBERMAN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.soberman
   * - Description:
     - Adiabatic exponent calculated per :cite:`Soberman1997` using core mass.
   * - Header Strings:
     - Zeta_Soberman, Zeta_Soberman(1), Zeta_Soberman(2), Zeta_Soberman(SN), Zeta_Soberman(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_SOBERMAN_HE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.sobermanHe
   * - Description:
     - Adiabatic exponent calculated per :cite:`Soberman1997` using He core mass.
   * - Header Strings:
     - Zeta_Soberman_He, Zeta_Soberman_He(1), Zeta_Soberman_He(2), Zeta_Soberman_He(SN), Zeta_Soberman_He(CP)

.. _supernova-events-states:

:ref:`Back to Top <stellar-props-top>`

Supernova events/states
-----------------------

Supernova events/states, both current ("is") and past ("experienced"), are stored within COMPAS as bitmaps. That means different values 
can be ORed or ANDed into the bit map, so that various events or states can be set concurrently.

The values shown below for the ``SN_EVENT`` type are powers of 2 so that they can be used in a bitmap and manipulated with bit-wise logical 
operators. Any of the individual supernova event/state types that make up the ``SN_EVENT`` type can be set independently of any other event/state.

`constants.h` defines an enum class for ``SN_EVENT``, and an associated label map, ``SN_EVENT_LABEL``, to provide labels for the events.  These 
are shown below::

    enum class SN_EVENT: int {
        NONE         = 0,
        CCSN         = 1,
        ECSN         = 2,
        PISN         = 4,
        PPISN        = 8,
        USSN         = 16,
        AIC          = 32,
    };


    const COMPASUnorderedMap<SN_EVENT, std::string> SN_EVENT_LABEL = {
        { SN EVENT::NONE,         "No Supernova" },
        { SN EVENT::CCSN,         "Core Collapse Supernova" },
        { SN EVENT::ECSN,         "Electron Capture Supernova" },
        { SN EVENT::PISN,         "Pair Instability Supernova" },
        { SN EVENT::PPISN,        "Pulsational Pair Instability Supernova" },
        { SN EVENT::USSN,         "Ultra Stripped Supernova" },
        { SN EVENT::AIC,          "Accretion-Induced Collapse" },
    };

A convenience function (shown below) is provided in ``utils.cpp`` to interpret the bit map.

::

    /*
    * Returns a single SN type based on the SN EVENT parameter passed
    *
    * Returns (in priority order):
    *
    * SN EVENT::CCSN iff CCSN bit is set and USSN bit is not set
    * SN EVENT::ECSN iff ECSN bit is set
    * SN EVENT::PISN iff PISN bit is set
    * SN EVENT::PPISN iff PPISN bit is set
    * SN EVENT::USSN iff USSN bit is set
    * SN EVENT::AIC iff AIC bit is set
    * SN EVENT::NONE otherwise
    *
    *
    * @param [IN] p SNEvent SN EVENT mask to check for SN event type
    * @return SN EVENT
    /
    SN EVENT SNEventType(const SN EVENT p SNEvent) {
        if ((p SNEvent & (SN EVENT::CCSN | SN EVENT::USSN)) == SN EVENT::CCSN ) return SN EVENT::CCSN;
        if ((p SNEvent & SN EVENT::ECSN )                   == SN EVENT::ECSN ) return SN EVENT::ECSN;
        if ((p SNEvent & SN EVENT::PISN )                   == SN EVENT::PISN ) return SN EVENT::PISN;
        if ((p SNEvent & SN EVENT::PPISN)                   == SN EVENT::PPISN) return SN EVENT::PPISN;
        if ((p SNEvent & SN EVENT::USSN )                   == SN EVENT::USSN ) return SN EVENT::USSN;
        if ((p SNEvent & SN EVENT::AIC )                    == SN EVENT::AIC )  return SN EVENT::AIC;

        return SN EVENT::NONE;
    }

