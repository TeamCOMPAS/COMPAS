SSE detailed output file record types
=====================================

Following is a list of the SSE Detailed Output file record type numbers and corresponding symbolic names, and their meaning:

1. INITIAL_STATE |BR|
   Record describes the initial state of the star
#. PRE_MASS_LOSS |BR|
   Record was logged after timestep taken, but before mass loss resolution
#. POST_MASS_LOSS |BR|
   Record was logged after after mass loss resolution
#. TIMESTEP_COMPLETED |BR|
   Record was logged immediately following the completion of the timestep (after all changes to the star)
#. FINAL_STATE |BR|
   Record describes the final state of the star
