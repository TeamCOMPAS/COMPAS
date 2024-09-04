Program option properties
=========================

When specifying known properties in a log file record specification record, the property name must be prefixed with 
the property type. Currently there is a single binary property type available for use: PROGRAM_OPTION.

For example, to specify the program option property ``RANDOM_SEED``, use::

    PROGRAM OPTION::RANDOM_SEED


.. _spec-options-props-top:

Following is an alphabetical list of program option properties available for inclusion in log file record specifiers.

**Jump to**
:ref:`A <spec-options-props-A>` :ref:`B <spec-options-props-B>` :ref:`C <spec-options-props-C>` :ref:`D <spec-options-props-D>`
:ref:`E <spec-options-props-E>` :ref:`F <spec-options-props-F>` :ref:`G <spec-options-props-G>` :ref:`H <spec-options-props-H>`
:ref:`I <spec-options-props-I>` :ref:`J <spec-options-props-J>` :ref:`K <spec-options-props-K>` :ref:`L <spec-options-props-L>`
:ref:`M <spec-options-props-M>` :ref:`N <spec-options-props-N>` :ref:`O <spec-options-props-O>` :ref:`P <spec-options-props-P>`
:ref:`Q <spec-options-props-Q>` :ref:`R <spec-options-props-R>` :ref:`S <spec-options-props-S>` :ref:`T <spec-options-props-T>`
:ref:`U <spec-options-props-U>` :ref:`V <spec-options-props-V>` :ref:`W <spec-options-props-W>` :ref:`X <spec-options-props-X>`
:ref:`Y <spec-options-props-Y>` :ref:`Z <spec-options-props-Z>`

.. _spec-options-props-A:

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ADD_OPTIONS_TO_SYSPARMS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_AddOptionsToSysParms
   * - Description:
     - Value of program option ``--add-options-to-sysparms``
   * - Header String:
     - PO_Add_Options_To_SysParms

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AllowMainSequenceStarToSurviveCommonEnvelope
   * - Description:
     - Value of program option ``--common-envelope-allow-main-sequence-survive``
   * - Header String:
     - PO_Allow_MS_To_Survive_CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ALLOW_RLOF_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AllowRLOFAtBirth
   * - Description:
     - Value of program option ``--allow-rlof-at-birth``
   * - Header String:
     - PO_Allow_RLOF@\ Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ALLOW_TOUCHING_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AllowTouchingAtBirth
   * - Description:
     - Value of program option ``--allow-touching-at-birth``
   * - Header String:
     - PO_Allow_Touching@\ Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ANG_MOM_CONSERVATION_DURING_CIRCULARISATION**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AngularMomentumConservationDuringCircularisation
   * - Description:
     - Value of program option ``--angular-momentum-conservation-during-circularisation``
   * - Header String:
     - PO_Conserve_AngMom@\ Circ

.. _spec-options-props-B:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BLACK_HOLE_KICKS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_BlackHoleKicks
   * - Description:
     - Value of program option ``--black-hole-kicks`` |br|
   * - Header String:
     - PO_BH_Kicks

DEPRECATION NOTICE: property ``BLACK_HOLE_KICKS`` has been deprecated and will soon be removed. Please use ``BLACK_HOLE_KICKS_MODE`` in future.

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BLACK_HOLE_KICKS_MODE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_BlackHoleKicks
   * - Description:
     - Value of program option ``--black-hole-kicks-mode`` |br|
   * - Header String:
     - PO_BH_Kicks

.. _spec-options-props-C:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CASE_BB_STABILITY_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CaseBBStabilityPrescription
   * - Description:
     - Value of program option ``--case-BB-stability-prescription``
   * - Header String:
     - PO_BB_Mass_xFer_Stblty_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CHE_MODE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CheMode
   * - Description:
     - Value of program option ``--chemically-homogeneous-evolution``
   * - Header String:
     - PO_CHE_Mode

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CIRCULARISE_BINARY_DURING_MT**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_CirculariseBinaryDuringMassTransfer
   * - Description:
     - Value of program option ``--circularise-binary-during-mass-transfer``
   * - Header String:
     - PO_Circularise@\ MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_ALPHA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeAlpha
   * - Description:
     - Value of program option ``--common-envelope-alpha``
   * - Header String:
     - PO_CE_Alpha

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_ALPHA_THERMAL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeAlphaThermal
   * - Description:
     - Value of program option ``--common-envelope-alpha-thermal``
   * - Header String:
     - PO_CE_Alpha_Thermal

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_LAMBDA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeLambda
   * - Description:
     - Value of program option ``--common-envelope-lambda``
   * - Header String:
     - PO_CE_Lambda

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_LAMBDA_MULTIPLIER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeLambdaMultiplier
   * - Description:
     - Value of program option ``--common-envelope-lambda-multiplier``
   * - Header String:
     - PO_CE_Lambda_Multiplier

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_LAMBDA_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CommonEnvelopeLambdaPrescription
   * - Description:
     - Value of program option ``--common-envelope-lambda-prescription``
   * - Header String:
     - PO_CE_Lambda_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionConstant
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-constant``
   * - Header String:
     - PO_CE_Mass_Accr_Constant

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionMax
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-max``
   * - Header String:
     - PO_CE_Mass_Accr_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionMin
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-min``
   * - Header String:
     - PO_CE_Mass_Accr_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionPrescription
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-prescription``
   * - Header String:
     - PO_CE_Mass_Accr_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeRecombinationEnergyDensity
   * - Description:
     - Value of program option ``--common-envelope-recombination-energy-density``
   * - Header String:
     - PO_CE_Recomb_Enrgy_Dnsty

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_SLOPE_KRUCKOW**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeSlopeKruckow
   * - Description:
     - Value of program option ``--common-envelope-slope-kruckow``
   * - Header String:
     - PO_CE_Slope_Kruckow

.. _spec-options-props-D:

.. _spec-options-props-E:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_Eccentricity
   * - Description:
     - Value of program option ``--eccentricity``
   * - Header String:
     - PO_Eccentricity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_EccentricityDistribution
   * - Description:
     - Value of program option ``--eccentricity-distribution``
   * - Header String:
     - PO_Eccentricity_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_EccentricityDistributionMax
   * - Description:
     - Value of program option ``--eccentricity-max``
   * - Header String:
     - PO_Eccentricity_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_EccentricityDistributionMin
   * - Description:
     - Value of program option ``--eccentricity-min``
   * - Header String:
     - PO_Eccentricity_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EDDINGTON_ACCRETION_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_EddingtonAccretionFactor
   * - Description:
     - Value of program option ``--eddington-accretion-factor``
   * - Header String:
     - PO_Eddington_Accr_Factor

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ENVELOPE_STATE_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_EnvelopeStatePrescription
   * - Description:
     - Value of program option ``--envelope-state-prescription``
   * - Header String:
     - PO_Envelope_State_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ENABLE_ROTATIONALLY_ENHANCED_MASS_LOSS**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_EnableRotationallyEnhancedMassLoss
   * - Description:
     - Value of program option ``--enable-rotationally-enhanced-mass-loss``
   * - Header String:
     - PO_Enable_Rotationally_Enhanced_Mass_Loss

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ENHANCE_CHE_LIFETIMES_LUMINOSITIES**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_EnhanceCHELifetimesLuminosities
   * - Description:
     - Value of program option ``--enhance-CHE-lifetimes-luminosities``
   * - Header String:
     - PO_Enhance_CHE_lifetimes_luminosities

.. _spec-options-props-F:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **FRYER_SUPERNOVA_ENGINE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_FryerSupernovaEngine
   * - Description:
     - Value of program option ``--fryer-supernova-engine``
   * - Header String:
     - PO_Fryer_SN_Engine

.. _spec-options-props-G:

.. _spec-options-props-H:

.. _spec-options-props-I:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMass
   * - Description:
     - Value of program option ``--initial-mass``
   * - Header String:
     - PO_Initial_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMass1
   * - Description:
     - Value of program option ``--initial-mass-1``
   * - Header String:
     - PO_Initial_Mass(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMass2
   * - Description:
     - Value of program option ``--initial-mass-2``
   * - Header String:
     - PO_Initial_Mass(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_InitialMassFunction
   * - Description:
     - Value of program option ``--initial-mass-function``
   * - Header String:
     - PO_Initial_Mass_Function

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMassFunctionMax
   * - Description:
     - Value of program option ``--initial-mass-max``
   * - Header String:
     - PO_Initial_Mass_Func_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMassFunctionMin
   * - Description:
     - Value of program option ``--initial-mass-min``
   * - Header String:
     - PO_Initial_Mass_Func_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION_POWER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMassFunctionPower
   * - Description:
     - Value of program option ``--initial-mass-power``
   * - Header String:
     - PO_Initial_Mass_Func_Power

.. _spec-options-props-J:

.. _spec-options-props-K:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_DIRECTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_KickDirectionDistribution
   * - Description:
     - Value of program option ``--kick-direction``
   * - Header String:
     - PO_Kick_Direction_Dstrbtn

DEPRECATION NOTICE: property ``KICK_DIRECTION`` has been deprecated and will soon be removed. Please use ``KICK_DIRECTION_DISTRIBUTION`` in future.

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_DIRECTION_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_KickDirectionDistribution
   * - Description:
     - Value of program option ``--kick-direction-distribution``
   * - Header String:
     - PO_Kick_Direction_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_DIRECTION_POWER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickDirectionPower
   * - Description:
     - Value of program option ``--kick-direction-power``
   * - Header String:
     - PO_Kick_Direction_Power

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_SCALING_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickScalingFactor
   * - Description:
     - Value of program option ``--kick-scaling-factor``
   * - Header String:
     - PO_Kick_Scaling_Factor

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitude
   * - Description:
     - Value of program option ``--kick-magnitude``
   * - Header String:
     - PO_Kick_Magnitude

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitude1
   * - Description:
     - Value of program option ``--kick-magnitude-1``
   * - Header String:
     - PO_Kick_Magnitude(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitude2
   * - Description:
     - Value of program option ``--kick-magnitude-2``
   * - Header String:
     - PO_Kick_Magnitude(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_KickMagnitudeDistribution
   * - Description:
     - Value of program option ``--kick-magnitude-distribution``
   * - Header String:
     - PO_Kick_Magnitude_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeDistributionMaximum
   * - Description:
     - Value of program option ``--kick-magnitude-max``
   * - Header String:
     - PO_Kick_Magnitude_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaCCSN BH
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-CCSN-BH``
   * - Header String:
     - PO_Sigma_Kick_CCSN_BH

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaCCSN NS
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-CCSN-NS``
   * - Header String:
     - PO_Sigma_Kick_CCSN_NS

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaForECSN
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-ECSN``
   * - Header String:
     - PO_Sigma_Kick_ECSN

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaForUSSN
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-USSN``
   * - Header String:
     - PO_Sigma_Kick_USSN

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MEAN_ANOMALY_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMeanAnomaly1
   * - Description:
     - Value of program option ``--kick-mean-anomaly-1``
   * - Header String:
     - PO_Kick_Mean_Anomaly(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MEAN_ANOMALY_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMeanAnomaly2
   * - Description:
     - Value of program option ``--kick-mean-anomaly-2``
   * - Header String:
     - PO_Kick_Mean_Anomaly(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_RANDOM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeRandom
   * - Description:
     - Value of program option ``--kick-magnitude-random``
   * - Header String:
     - PO_Kick_Magnitude_Random

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_RANDOM_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeRandom1
   * - Description:
     - Value of program option ``--kick-magnitude-random-1``
   * - Header String:
     - PO_Kick_Magnitude_Random(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_RANDOM_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeRandom2
   * - Description:
     - Value of program option ``--kick-magnitude-random-2``
   * - Header String:
     - PO_Kick_Magnitude_Random(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_PHI_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickPhi1
   * - Description:
     - Value of program option ``--kick-phi-1``
   * - Header String:
     - PO_Kick_Mean_Phi(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_PHI_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickPhi2
   * - Description:
     - Value of program option ``--kick-phi-2``
   * - Header String:
     - PO_Kick_Mean_Phi(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_THETA_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickTheta1
   * - Description:
     - Value of program option ``--kick-theta-1``
   * - Header String:
     - PO_Kick_Mean_Theta(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_THETA_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickTheta2
   * - Description:
     - Value of program option ``--kick-theta-2``
   * - Header String:
     - PO_Kick_Mean_Theta(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LBV_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_LuminousBlueVariableFactor
   * - Description:
     - Value of program option ``--luminous-blue-variable-multiplier``
   * - Header String:
     - PO_LBV_Factor

.. _spec-options-props-L:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LBV_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_LuminousBlueVariablePrescription
   * - Description:
     - Value of program option ``--luminous-blue-variable-prescription``
   * - Header String:
     - PO_LBV_Mass_Loss_Prscrptn

DEPRECATION NOTICE: property ``LBV_PRESCRIPTION`` has been deprecated and will soon be removed. Please use ``LBV_MASS_LOSS_PRESCRIPTION`` in future.

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LBV_MASS_LOSS_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_LuminousBlueVariablePrescription
   * - Description:
     - Value of program option ``--LBV-mass-loss-prescription``
   * - Header String:
     - PO_LBV_Mass_Loss_Prscrptn

:ref:`Back to Top <spec-options-props-top>`

.. _spec-options-props-M:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_LOSS_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassLossPrescription
   * - Description:
     - Value of program option ``--mass-loss-prescription``
   * - Header String:
     - PO_Mass_Loss_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassRatio
   * - Description:
     - Value of program option ``-``-mass-ratio``
   * - Header String:
     - PO_Mass_Ratio

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassRatioDistribution
   * - Description:
     - Value of program option ``--mass-ratio-distribution``
   * - Header String:
     - PO_Mass_Ratio_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassRatioDistributionMax
   * - Description:
     - Value of program option ``--mass-ratio-max``
   * - Header String:
     - PO_Mass_Ratio_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassRatioDistributionMin
   * - Description:
     - Value of program option ``--mass-ratio-min``
   * - Header String:
     - PO_Mass_Ratio_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_EVOLUTION_TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MaxEvolutionTime
   * - Description:
     - Value of program option ``--maximum-evolution-time``
   * - Header String:
     - PO_Max_Evolution_Time

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_DONOR_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MaximumMassDonorNandezIvanova
   * - Description:
     - Value of program option ``--maximum-mass-donor-nandez-ivanova``
   * - Header String:
     - PO_Max_Donor_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_NEUTRON_STAR_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MaximumNeutronStarMass
   * - Description:
     - Value of program option ``--maximum-neutron-star-mass``
   * - Header String:
     - PO_Max_NS_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_TIMESTEPS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MaxNumberOfTimestepIterations
   * - Description:
     - Value of program option ``--maximum-number-timestep-iterations``
   * - Header String:
     - PO_Max_Timesteps

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MCBUR1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_mCBUR1
   * - Description:
     - Value of program option ``--mcbur1``
   * - Header String:
     - PO_MCBUR1

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_Metallicity
   * - Description:
     - Value of program option ``--metallicity``
   * - Header String:
     - PO_Metallicity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MetallicityDistribution
   * - Description:
     - Value of program option ``--metallicity-distribution``
   * - Header String:
     - PO_Metallicity_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MetallicityDistributionMax
   * - Description:
     - Value of program option ``--metallicity-max``
   * - Header String:
     - PO_Metallicity_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MetallicityDistributionMin
   * - Description:
     - Value of program option ``--metallicity-min``
   * - Header String:
     - PO_Metallicity_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MINIMUM_MASS_SECONDARY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MinimumMassSecondary
   * - Description:
     - Value of program option ``--minimum-secondary-mass``
   * - Header String:
     - PO_Min_Secondary_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_ACCRETION_EFFICIENCY_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferAccretionEfficiencyPrescription
   * - Description:
     - Value of program option ``--mass-transfer-accretion-efficiency-prescription``
   * - Header String:
     - PO_MT_Acc_Efficiency_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_ANG_MOM_LOSS_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferAngularMomentumLossPrescription
   * - Description:
     - Value of program option ``--mass-transfer-angular-momentum-loss-prescription``
   * - Header String:
     - PO_MT_AngMom_Loss_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_FRACTION_ACCRETED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassTransferFractionAccreted
   * - Description:
     - Value of program option ``--mass-transfer-fa``
   * - Header String:
     - PO_MT_Fraction_Accreted

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_JLOSS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassTransferJloss
   * - Description:
     - Value of program option ``--mass-transfer-jloss``
   * - Header String:
     - PO_MT_JLoss

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_THERMAL_LIMIT_C**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassTransferCParameter
   * - Description:
     - Value of program option ``--mass-transfer-thermal-limit-C``
   * - Header String:
     - PO_MT_Thermal_Limit_C

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_REJUVENATION_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferRejuvenationPrescription
   * - Description:
     - Value of program option ``--mass-transfer-rejuvenation-prescription``
   * - Header String:
     - PO_MT_Rejuvenation_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_THERMALLY_LIMITED_VARIATION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferThermallyLimitedVariation
   * - Description:
     - Value of program option ``--mass-transfer-thermal-limit-accretor``
   * - Header String:
     - PO_MT_Thermally_Lmtd_Variation

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MULLER_MANDEL_KICK_MULTIPLIER_BH**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MullerMandelKickBH
   * - Description:
     - Value of program option ``--muller-mandel-kick-multiplier-BH``
   * - Header String:
     - PO_MM_Kick_Multiplier_BH

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MULLER_MANDEL_KICK_MULTIPLIER_NS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MullerMandelKickNS
   * - Description:
     - Value of program option ``--muller-mandel-kick-multiplier-NS``
   * - Header String:
     - PO_MM_Kick_Multiplier_NS

.. _spec-options-props-N:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NEUTRINO_MASS_LOSS_ASSUMPTION_BH**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_NeutrinoMassLossAssumptionBH
   * - Description:
     - Value of program option ``--neutrino-mass-loss-BH-formation``
   * - Header String:
     - PO_Neutrino_Mass_Loss_Assmptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NEUTRINO_MASS_LOSS_VALUE_BH**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_NeutrinoMassLossValueBH
   * - Description:
     - Value of program option ``--neutrino-mass-loss-BH-formation-value``
   * - Header String:
     - PO_Neutrino_Mass_Loss_Value

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NOTES**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - Options::m_Notes
   * - Description:
     - Value of program option ``--notes``
   * - Header String:
     - as specified by program option ``--Notes-Hdrs``

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NS_EOS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_NeutronStarEquationOfState
   * - Description:
     - Value of program option ``--neutron-star-equation-of-state``
   * - Header String:
     - PO_NS_EOS

.. _spec-options-props-O:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_OrbitalPeriod
   * - Description:
     - Value of program option ``--orbital-period``
   * - Header String:
     - PO_Orbital_Period

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_OrbitalPeriodDistribution
   * - Description:
     - Value of program option ``--orbital-period-distribution``
   * - Header String:
     - PO_Orbital_Period_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_OrbitalPeriodDistributionMax
   * - Description:
     - Value of program option ``--orbital-period-max``
   * - Header String:
     - PO_Orbital_Period_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_OrbitalPeriodDistributionMin
   * - Description:
     - Value of program option ``--orbital-period-min``
   * - Header String:
     - PO_Orbital_Period_Min

.. _spec-options-props-P:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PISN_LOWER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PairInstabilityLowerLimit
   * - Description:
     - Value of program option ``--PISN-lower-limit``
   * - Header String:
     - PO_PISN_Lower_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PISN_UPPER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PairInstabilityUpperLimit
   * - Description:
     - Value of program option ``--PISN-upper-limit``
   * - Header String:
     - PO_PISN_Upper_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PPI_CO_CORE_SHIFT_HENDRIKS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsationalPairInstabilityCOCoreShiftHendriks
   * - Description:
     - Value of program option ``--PPI-CO-Core-Shift-Hendriks``
   * - Header String:
     - PO_PPI_CO_CORE_SHIFT_HENDRIKS

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PPI_LOWER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsationalPairInstabilityLowerLimit
   * - Description:
     - Value of program option ``--PPI-lower-limit``
   * - Header String:
     - PO_PPI_Lower_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PPI_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_PulsationalPairInstabilityPrescription
   * - Description:
     - Value of program option ``--pulsational-pair-instability-prescription``
   * - Header String:
     - PO_PPI_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PPI_UPPER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsationalPairInstabilityUpperLimit
   * - Description:
     - Value of program option ``--PPI-upper-limit``
   * - Header String:
     - PO_PPI_Upper_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_PulsarBirthMagneticFieldDistribution
   * - Description:
     - Value of program option ``--pulsar-birth-magnetic-field-distribution``
   * - Header String:
     - PO_Pulsar_Mag_Field_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthMagneticFieldDistributionMax
   * - Description:
     - Value of program option ``--pulsar-birth-magnetic-field-distribution-max``
   * - Header String:
     - PO_Pulsar_Mag_Field_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthMagneticFieldDistributionMin
   * - Description:
     - Value of program option ``--pulsar-birth-magnetic-field-distribution-min``
   * - Header String:
     - PO_Pulsar_Mag_Field_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_PulsarBirthSpinPeriodDistribution
   * - Description:
     - Value of program option ``--pulsar-birth-spin-period-distribution``
   * - Header String:
     - PO_Pulsar_Spin_Period_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_BIRTH_SPIN_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthSpinPeriodDistributionMax
   * - Description:
     - Value of program option ``--pulsar-birth-spin-period-distribution-max``
   * - Header String:
     - PO_Pulsar_Spin_Period_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_BIRTH_SPIN_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthSpinPeriodDistributionMin
   * - Description:
     - Value of program option ``--pulsar-birth-spin-period-distribution-min``
   * - Header String:
     - PO_Pulsar_Spin_Period_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options:m_PulsarMagneticFieldDecayMassscale
   * - Description:
     - Value of program option ``--pulsar-magnetic-field-decay-massscale``
   * - Header String:
     - PO_Pulsar_Mag_Field_Decay_mScale

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options:m_PulsarMagneticFieldDecayTimescale
   * - Description:
     - Value of program option ``--pulsar-magnetic-field-decay-timescale``
   * - Header String:
     - PO_Pulsar_Mag_Field_Decay_tScale

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MINIMUM_MAGNETIC_FIELD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarLog10MinimumMagneticField
   * - Description:
     - Value of program option ``--pulsar-minimum-magnetic-field``
   * - Header String:
     - PO_Pulsar_Minimum_Mag_Field

.. _spec-options-props-Q:

.. _spec-options-props-R:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED**
     -
   * - Data type:
     - UNSIGNED LONG INT
   * - COMPAS variable:
     - Options::m_RandomSeed
   * - Description:
     - Value of program option ``--random-seed``
   * - Header String:
     - PO_SEED(OPTION)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED_CMDLINE**
     -
   * - Data type:
     - UNSIGNED LONG INT
   * - COMPAS variable:
     - Options::m_FixedRandomSeed
   * - Description:
     - Value of program option ``--random-seed`` (specified on the commandline)
   * - Header String:
     - PO_SEED(CMDLINE)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **REMNANT_MASS_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_RemnantMassPrescription
   * - Description:
     - Value of program option ``--remnant-mass-prescription``
   * - Header String:
     - PO_Remnant_Mass_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROTATIONAL_VELOCITY_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_RotationalVelocityDistribution
   * - Description:
     - Value of program option ``--rotational-velocity-distribution``
   * - Header String:
     - PO_Rotational_Velocity_Dstrbtn

.. _spec-options-props-S:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SCALE_CHE_MASS_LOSS_SURF_HE_ABUNDANCE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_ScaleCHEMassLossWithSurfaceHeliumAbundance
   * - Description:
     - Value of program option ``--scale-CHE-mass-loss-with-surface-helium-abundance``
   * - Header String:
     - PO_Scale_CHE_Mass_Loss_Surf_He_Abundance

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SCALE_TERMINAL_WIND_VEL_METALLICITY_POWER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_ScaleTerminalWindVelocityWithMetallicityPower
   * - Description:
     - Value of program option ``--scale-terminal-wind-velocity-with-metallicity-power``
   * - Header String:
     - PO_Scale_Terminal_Wind_Vel_Metallicity_Power
   
.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_SemiMajorAxis
   * - Description:
     - Value of program option ``--semi-major-axis``
   * - Header String:
     - PO_Semi-Major_Axis

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_SemiMajorAxisDistribution
   * - Description:
     - Value of program option ``--semi-major-axis-distribution``
   * - Header String:
     - PO_Semi-Major_Axis_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_SemiMajorAxisDistributionMax
   * - Description:
     - Value of program option ``--semi-major-axis-max``
   * - Header String:
     - PO_Semi-Major_Axis_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_SemiMajorAxisDistributionMin
   * - Description:
     - Value of program option ``--semi-major-axis-min``
   * - Header String:
     - PO_Semi-Major_Axis_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_ZETA_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_StellarZetaPrescription
   * - Description:
     - Value of program option ``--stellar-zeta-prescription``
   * - Header String:
     - PO_Stellar_Zeta_Prscrptn

.. _spec-options-props-T:

.. _spec-options-props-U:

.. _spec-options-props-V:

.. _spec-options-props-W:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **WR_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_WolfRayetFactor
   * - Description:
     - Value of program option ``--wolf-rayet-multiplier``
   * - Header String:
     - PO_WR_Factor

.. _spec-options-props-X:

.. _spec-options-props-Y:

.. _spec-options-props-Z:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_ADIABATIC_ARBITRARY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_ZetaAdiabaticArbitrary
   * - Description:
     - Value of program option ``--zeta-adiabatic-arbitrary``
   * - Header String:
     - PO_Zeta_Adiabatic_Arbitrary

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_MS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_ZetaMainSequence
   * - Description:
     - Value of program option ``--zeta-main-sequence``
   * - Header String:
     - PO_Zeta_Main_Sequence_Giant

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_RADIATIVE_ENVELOPE_GIANT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_ZetaRadiativeEnvelopeGiant
   * - Description:
     - Value of program option ``--zeta-radiative-envelope-giant``
   * - Header String:
     - PO_Zeta_Radiative_Envelope_Giant

