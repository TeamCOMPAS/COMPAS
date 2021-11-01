Rand::Random(p_Lower, p_Upper)
==============================

::

    DOUBLE Rand::Random(const DOUBLE p_Lower, const DOUBLE p_Upper)

Returns a random floating point number uniformly distributed in the range [``p_Lower``, ``p_Upper``), where ``p_Lower`` 
:math:`\small \leq` ``p_Upper``.

(``p_Lower`` and ``p_Upper`` will be swapped if ``p_Lower`` :math:`\small \gt` ``p_Upper`` as passed)