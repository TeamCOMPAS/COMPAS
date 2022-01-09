COMPAS build methods
====================

We provide two methods to build the COMPAS executable: a method to build COMPAS locally, and a method to create a ``Docker`` image.

Note that each method has a separate, and different, makefile: if functionality is added, modified, or removed from either of the
provided makefiles (``Makefile`` locally, ``Makefile.docker`` for the docker image), the corresponding changes should also be made
to the other makefile.


.. toctree::
   :maxdepth: 1

   ./COMPAS-local-build
   ./docker-developer

   