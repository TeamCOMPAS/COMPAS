Standard log files
==================

COMPAS defines several standard log files that may be produced depending upon the simulation mode (Single Star Evolution (SSE), 
or Binary Star Evolution (BSE), see the ``--mode`` program option)), and the value of various program options.

The standard log files are:

    .. list-table::
       :widths: 22 78 
       :header-rows: 0
       :class: aligned-text

       * - System Parameters
         - Records summary information for all stars, or binary stars, during evolution.
       * -
         -
       * - Supernovae
         - Records summary information for all stars that experience a SN event during evolution.
       * -
         -
       * - Detailed Output
         - Records detailed information for a star, or binary star, during evolution.
       * -
         - Enable with program option ``--detailed-output``.
       * -
         -
       * - SwitchLog
         - Records detailed information for all stars, or binary stars, at the time of each stellar type switch during evolution.
       * - 
         - Enable with program option ``--switch-log``.
       * -
         -
       * - Double Compact Objects
         - Records summary information for all binary systems that form DCOs during BSE.
       * -
         -
       * - Common Envelopes
         - Records summary information for all binary systems that experience CEEs during BSE.
       * -
         -
       * - Pulsar Evolution
         - Records detailed Pulsar evolution information during BSE.
       * -
         - Enable with program option ``--evolve-pulsars``.
       * -
         -
       * - RLOF
         - Records detailed information RLOF events during BSE.
       * - 
         - Enable with program option ``--rlof-printing``.

