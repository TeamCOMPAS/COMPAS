Running COMPAS via Docker
=========================

Docker functionality has been added to COMPAS to reduce time and effort required to set up the COMPAS deployment environment.

Instead of having to install and configure several libraries and tools (e.g. ``python``/``pip``, ``numpy``, ``g++``, ``boost``, ``hdf5``) 
which can vary considerably beween operating systems and existing toolchains, users can instead opt to install ``Docker`` and run COMPAS 
with a single command.

This also gives users the ability to run COMPAS on cloud solutions like ``AWS EC2``\ [#f1]_ or ``Google Compute Engine``\ [#f2]_ where hundreds 
of cores can be provisioned without having to manually configure the environment.

``Docker`` works by creating an isolated and standalone environment known as a `container`\ [#f3]_. Containers can be created or destroyed 
without affecting the host machine or other containers (containers can still interact with each other and the host machine through mounted 
directories/files or exposed ports).

Containers are instances of images. An image is a pre-defined setup/environment that is instantiated when started as a container (containers 
are to images what objects are to classes in the OO paradigm)\ [#f4]_. 

Containers are (almost) always run as a ``Linux`` environment. A major benefit of this is the ability to run Linux applications in a ``Windows`` 
or ``macOS`` environment without having to jump through hoops or have a diminished experience.

Image definitions can be defined by users (e.g. ``Dockerfiles``); there are also standard images publicly available on ``dockerHub``\ [#f5]_

This following sections assume ``Docker`` has been installed and is running. For Windows and MacOS users, see 
`Docker Desktop <https://www.docker.com/products/docker-desktop>`_.


.. toctree::
   :maxdepth: 1

   running-via-docker-installation
   running-via-docker-running


.. rubric:: Footnotes

.. [#f1] https://aws.amazon.com/ec2/
.. [#f2] https://cloud.google.com/compute
.. [#f3] https://www.docker.com/resources/what-container
.. [#f4] `https://stackoverflow.com/questions/23735149 <https://stackoverflow.com/questions/23735149/what-is-the-difference-between-a-docker-image-and-a-container#:~:text=An%20instance%20of%20an%20image,of%20layers%20as%20you%20describe.&text=You%20can%20see%20all%20your,an%20image%20is%20a%20container>`_
.. [#f5] https://hub.docker.com/

