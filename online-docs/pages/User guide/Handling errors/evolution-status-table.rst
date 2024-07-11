COMPAS evolution status
=======================

Following is a list of COMPAS evolution status numbers, corresponding symbolic name, and meaning:

1. DONE |br|
   Simulation completed without error
#. ERROR |br|
   Evolution stopped because an error occurred
#. TIMES_UP |br|
   Evolution stopped because the allowed time was exhausted
#. STEPS_UP |br|
   Evolution stopped because the allowed number of timesteps was exhausted
#. NO_TIMESTEPS_READ |br|
   Evolution stopped because there was an error reading the user-provided timesteps
#. TIMESTEPS_EXHAUSTED |br|
   The user-provided timesteps were exhausted before evolution was complete - default timesteps used
#. TIMESTEPS_NOT_CONSUMED |br|
   Evolution completed before all user-provided timesteps were consumed
#. SSE_ERROR |br|
   Evolution stopped because an error occurred while evolving one of the constituent stars of a binary star
#. BINARY_ERROR |br|
   Evolution stopped because an error occurred while evolving a binary star
#. DCO_MERGER_TIME |br|
   Evolution stopped because the the evolution time exceeded DCO merger (formation + coalescence) time
#. STARS_TOUCHING |br|
   Evolution stopped because the constituent stars of a binary star are touching
#. STELLAR_MERGER |br|
   Evolution stopped because the constituent stars of a binary star merged
#. STELLAR_MERGER_AT_BIRTH |br|
   Evolution stopped because the constituent stars of a binary star merged at birth
#. DCO |br|
   Evolution stopped because a double compact object formed
#. WD_WD |br|
   Evolution stopped because a double white dwarf formed
#. MASSLESS_REMNANT |br|
   Evolution stopped because a massless remnant formed
#. UNBOUND |br|
   Evolution stopped because the binary star became unbound
#. NOT_STARTED |br|
   Simulation not started - not complete
#. STARTED |br|
   Simulation started - not complete
