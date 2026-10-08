BSE pulsar evolution file record types
======================================

Following is a list of the BSE Pulsar Evolution file record type numbers and corresponding symbolic names, and their meaning:

1. DEFAULT |BR|
   Default BSE_PULSAR_EVOLUTION file record type
#. POST_SN |BR|
   Record was logged immediately following a supernova event
#. POST_BINARY_TIMESTEP |BR|
   Record was logged immediately following binary timestep (i.e. the evolution of the binary system for a single timestep)

Each record describes the state at the end of a timestep, and ``DT`` is the length of that timestep. By default only
the timestep records are written (``--logfile-pulsar-evolution-record-types`` = 4). In the timestep in which a neutron
star forms, the timestep was taken by its progenitor: the record holds the pulsar's birth values and the neutron star's
``Age`` is 0. To weight records by the time each pulsar has spent in that state (e.g. for a constant star formation rate
snapshot), use :math:`\min(\mathrm{DT}, \mathrm{Age})` rather than ``DT``.
