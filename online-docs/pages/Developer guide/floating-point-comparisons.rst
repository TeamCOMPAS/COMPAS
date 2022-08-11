Floating-point comparisons
==========================

Floating-point comparisons are inherently problematic. Testing floating-point numbers for equality, or even inequality, is fraught with
problems due to the internal representation of floating-point numbers: floatingpoint numbers are stored with a fixed number of binary 
digits, which limits their precision and accuracy. The problems with floating-point comparisons are even more evident if one or both of 
the numbers being compared are the results of (perhaps several) floating-point operations (rather than comparing constants).

To avoid the problems associated with floating-point comparisons it is (almost always) better to do any such comparisons with a tolerance
rather than an absolute comparison. To this end, a floating-point comparison function has been provided, and (almost all of) the 
floating-point comparisons in the code have been changed to use that function. The function uses both an absolute tolerance and a relative 
tolerance, which are both declared in constants.h. Whether the function uses a tolerance or not can be changed by ``#define``-ing or 
``#undef``-ing the ``COMPARE_WITH_TOLERANCE`` flag in ``constants.h`` (so the change is a compile-time change, not run-time).


The compare function is defined in utils.h and is implemented as follows::

    static int Compare(const double p_X, const double p_Y) { |br|
    #ifdef COMPARE WITH TOLERANCE
    
        return (fabs(p X – p Y) <= max( FLOAT_TOLERANCE_ABSOLUTE,
                                        FLOAT_TOLERANCE_RELATIVE * 
                                        max( fabs(p_X), fabs(p Y)))) ? 0 
                                                                     : (p_X < p_Y ? –1 : 1);
    #else

        return (p_X == p_Y) ? 0 : (p_X < p_Y ? –1 : 1);

    #endif



If ``COMPARE_WITH_TOLERANCE`` is defined, ``p_X`` and ``p_Y`` are compared with tolerance values, whereas if ``COMPARE_WITH_TOLERANCE`` is
not defined the comparison is an absolute comparison.

The function returns an integer indicating the result of the comparison:
    .. list-table::
       :widths: 8 92 
       :header-rows: 0
       :class: aligned-text

       * - –1 
         - indicates that ``p_X`` is considered to be less than ``p_Y``
       * - |_| |_| 0
         - indicates ``p_X`` and ``p_Y`` are considered to be equal
       * - +1
         - indicates that ``p_X`` is considered to be greater than ``p_Y``

The comparison is done using both an absolute tolerance and a relative tolerance. The tolerances can be defined to be the same number, or
different numbers. If the relative tolerance is defined as 0.0, the comparison is done using the absolute tolerance only, and if the 
absolute tolerance is defined as 0.0 the comparison is done with the relative tolerance only.

Absolute tolerances are generally more effective when the numbers being compared are small – so using an absolute tolerance of (say) 
0.0000005 is generally effective when comparing single-digit numbers (or so), but is less effective when comparing numbers in the thousands
or millions. For comparisons of larger numbers a relative tolerance is generally more effective (the actual tolerance is wider because the 
relative tolerance is multiplied by the larger absolute value of the numbers being compared).

The tolerances used for the comparison are defined in ``constants.h`` as ``FLOAT_TOLERANCE_ABSOLUTE`` and ``FLOAT_TOLERANCE_RELATIVE``.

There is a little overhead in the comparisons even when the tolerance comparison is disabled, but it shouldn’t be prohibitive.
