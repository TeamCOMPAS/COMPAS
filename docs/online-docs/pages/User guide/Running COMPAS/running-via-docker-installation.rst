Installing the COMPAS Docker image
----------------------------------

COMPAS ``Docker`` images for all releases of COMPAS are hosted on ``dockerHub``, and tagged\ [#f1]_ with the COMPAS version number.

The latest COMPAS ``Docker`` compiled version of COMPAS (dev branch) can be retrieved by executing either of the following commands::

    docker pull teamcompas/compas
    docker pull teamcompas/compas:latest

Other versions can be retrieved by specifying the tag that corresponds to the COMPAS version required. For example, to retrieve the
image for COMPAS version 2.12.0, type::

    docker pull teamcompas/compas:2.12.0


To see all available versions, go to the `TeamCOMPAS docker hub page <https://hub.docker.com/u/teamcompas>`_.


.. rubric:: Footnotes

.. [#f1] https://docs.docker.com/engine/reference/commandline/tag/
