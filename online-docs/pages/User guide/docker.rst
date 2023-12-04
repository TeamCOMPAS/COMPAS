COMPAS and Docker
=================

Docker has been added to COMPAS to reduce time and effort required to
set up the COMPAS deployment environment.

Instead of having to install and configure several libraries and tools
(e.g. python/pip, numpy, g++, boost) which can vary considerably between
operating systems and existing toolchains, users can instead opt to
install Docker and run COMPAS with a single command.

This also gives users the ability to run COMPAS on cloud solutions like
`AWS EC2 <https://aws.amazon.com/ec2/>`__ or `Google Compute
Engine <https://cloud.google.com/compute>`__ where hundreds of cores can
be provisioned without having to manually configure the environment.

Docker works by creating an isolated and standalone environment known
as a `container <https://www.docker.com/resources/what-container>`__.
Containers can be created or destroyed without affecting the host
machine or other containers\*.

Containers are instances of images. An image is a pre-defined
setup/environment that is instantiated when started as a container
(containers are to images what objects are to classes). More
`here <https://stackoverflow.com/questions/23735149/what-is-the-difference-between-a-docker-image-and-a-container#:~:text=An%20instance%20of%20an%20image,of%20layers%20as%20you%20describe.&text=You%20can%20see%20all%20your,an%20image%20is%20a%20container.>`__
on the relationship between images and container.

Containers are (almost) always run as a Linux environment. A major
benefit of this is the ability to run Linux applications in a Windows or
MacOS environment without having to jump through hoops or have a
diminished experience.

Image definitions can be defined by users (e.g. Dockerfiles); there are
also standard images publicly available on `Docker
Hub <https://hub.docker.com/>`__

All that is required to start using COMPAS with Docker is the "Usage"
section (the "CI/CD" section is also highly recommended).
The other sections are provided for extra info.

\* Containers can still interact with each other and the host machine
through mounted directories/files or exposed ports.



Usage


N.B. This section assumes `Docker <https://www.docker.com/>`__ has
been installed and is running.
For Windows and MacOS users, see
`here <https://www.docker.com/products/docker-desktop>`__.

Installing


The latest compiled version of COMPAS (dev branch) can be retrieved by
running
``docker pull teamcompas/compas``

Other versions can be used by adding a version
`tag <https://docs.docker.com/engine/reference/commandline/tag/>`__.
For example, COMPAS version 2.12.0 would be
``teamcompas/compas:2.12.0``.
To see all available versions, go to the TeamCOMPAS docker hub page
`here <https://hub.docker.com/u/teamcompas>`__.

Running


COMPAS can still be configured via command line arguments passed to the
COMPAS executable or via a ``runSubmit.py`` file.

Running runSubmit.py


To run COMPAS via a ``runSubmit.py`` file, the command is a little
more complex.

::

    docker run                                                  \
    --rm                                                    \
    -it                                                     \
    -v $(pwd)/compas-logs:/app/COMPAS/logs                  \
    -v $(pwd)/runSubmit.py:/app/starts/runSubmit.py   \
    -e COMPAS_EXECUTABLE_PATH=/app/COMPAS/bin/COMPAS        \
    -e COMPAS_LOGS_OUTPUT_DIR_PATH=/app/COMPAS/logs         \
    teamcompas/compas                                       \
    python3 /app/starts/runSubmit.py                     

Breaking down this command:

``docker run``
creates a container

``--rm``
`Clean
up <https://docs.docker.com/engine/reference/run/#clean-up---rm>`__
destroy the container once it finishes running the command

``-it``
short for `-i and
-t <https://docs.docker.com/engine/reference/run/#foreground>`__ -
provides an interactive terminal

``-v <path-on-host>:<path-in-container>``
`Bind mounts <https://docs.docker.com/storage/bind-mounts/>`__
mount ``<path-on-host>`` to ``<path-in-container``>
This time we not only want to get the output from COMPAS on the host
machine, we also want to supply a ``runSubmit.py`` to the container
from the host machine.

NOTE: if you decide to execute using ``runSubmit.py``, you will need 
a ``compasConfigDefault.yaml``  file in the same directory. This file 
can be find in the same directory as the ``runSubmit.py``, and contains
the default COMPAS choices for stellar and binary physics. These choices
can be changed by modifying the options available in the ``compasConfigDefault.yaml`` 
file.

``-e VAR_NAME=value``
`Environment
variables <https://docs.docker.com/engine/reference/run/#env-environment-variables>`__
set the environment variable ``VAR_VAME`` to ``value``

``teamcompas/compas``
the image to run

``python3 /app/starts/runSubmit.py``
the command to run when the container starts

Run the COMPAS executable


To run the COMPAS executable directly (i.e. without ``runSubmit.py``)

::

    docker run                                  \
    --rm                                    \
    -it                                     \
    -v $(pwd)/compas-logs:/app/COMPAS/logs  \
    teamcompas/compas                       \
    bin/COMPAS                              \
    --number-of-binaries=5                  \
    --outputPath=/app/COMPAS/logs

Breaking down this command:

``docker run``
creates a container

``--rm``
`Clean
up <https://docs.docker.com/engine/reference/run/#clean-up---rm>`__
destroy the container once it finishes running the command

``-it``
short for `-i and
-t <https://docs.docker.com/engine/reference/run/#foreground>`__ -
provides an interactive terminal

``-v <path-on-host>:<path-in-container>``
`Bind mounts <https://docs.docker.com/storage/bind-mounts/>`__
mount ``<path-on-host>`` to ``<path-in-container>``
In this instance, make it so
``$(pwd)/compas-logs on my machine is the same as``/app/COMPAS/logs\`
inside the container

``teamcompas/compas``
the image to run

``bin/COMPAS``
the command to run when the container starts

``--number-of-binaries``
anything after the given start command is passed to that command, in
this case, the flag to set the number of binaries

``--outputPath /app/COMPAS/logs``
same as above, anything after the start command is given to that start
command, here it forces logs to go to the directory that is mapped to
the host machine

More info on ``docker run``
`here <https://docs.docker.com/engine/reference/run/>`__

NOTE 1:

Two new environment variables have been added, both of these apply to
``runSubmit.py`` only and are non-breaking changes.

``COMPAS_EXECUTABLE_PATH`` is an addition to the default
``runSubmit.py`` that overrides where ``runSubmit.py`` looks for
the compiled COMPAS.
This override exists purely for ease-of-use from the command line.

``COMPAS_LOGS_OUTPUT_DIR_PATH`` is also an addition to the default
``runSubmit.py`` that overrides where logs are placed.
The override exists because the mounted directory (option ``-v``) is
created before COMPAS runs. COMPAS sees that the directory where it's
supposed to put logs already exists, so it created a different (i.e.
non-mapped) directory to deposit logs in.

NOTE 2:

The ``docker run ...`` examples above both use the ``-it`` options.
If you want to run multiple instances of COMPAS, I would highly
recommend using `detached
mode <https://docs.docker.com/engine/reference/run/#detached--d>`__
(``-d``) instead.
All container output will be hidden.

An example where this would be useful is if you were running 4
instances of COMPAS at once.
You could copy/paste the following into the terminal...

::

    docker run --rm -d -v $(pwd)/compas-logs/run_0:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_01.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    
    docker run --rm -d -v $(pwd)/compas-logs/run_1:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_02.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    
    docker run --rm -d -v $(pwd)/compas-logs/run_2:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_03.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    
    docker run --rm -d -v $(pwd)/compas-logs/run_3:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_04.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py

...which would run 4 separate instances of COMPAS, each with its own
``runSubmit.py`` file and logging directory, and all console output
suppressed.

You may want to check the console output to see how far into the run
COMPAS is.
The command for this is ``docker logs <container_id>``.
You can get the container id by running ``docker ps``.



CI/CD


The latest version of COMPAS (dev branch) is available at
``teamcompas/compas``.
This is provided automatically by CI/CD.

Whenever a push to
`TeamCOMPAS/dev <https://github.com/TeamCOMPAS/COMPAS/tree/dev>`__ a
continuous deployment process automatically
`builds <https://docs.docker.com/engine/reference/commandline/build/>`__
a new image and deploys it to DockerHub with a ``tag`` that corresponds
to the value of ``VERSION_STRING`` in ``constants.h``.

At time of writing, `GitHub
Actions <https://github.com/features/actions>`__ is facilitating the
above process. While this is convenient (because it's free and well
supported) it is quite slow. I have plans to create a
`runner <https://help.github.com/en/actions/getting-started-with-github-actions/core-concepts-for-github-actions#runner>`__
locally with a high core count that can be used to compile COMPAS
quickly, but haven't gotten around to it yet.

You can realistically expect the latest COMPAS docker image to be
available 5 - 10 minutes after pushing/merging.

The Github Actions configuration is in
``/.github/workflows/dockerhub-ci.yml``.

Atlassian has a `good
writeup <https://www.atlassian.com/continuous-delivery/principles/continuous-integration-vs-delivery-vs-deployment>`__
about what CI/CD is.



Bonus Info


Dockerfile


The `Dockerfile <https://docs.docker.com/engine/reference/builder/>`__
defines how the docker image is constructed.

Images are created as a combination of layers.
During the build process each layer is cached and only updated on
subsequent builds if that layer would change.

The Dockerfile for COMPAS is made up of 8 layers.

``FROM ubuntu:18.04``
Use `Ubuntu 18.04 <https://hub.docker.com/_/ubuntu>`__ as a base
(provided by Docker Hub)
`https://docs.docker.com/engine/reference/builder/#from <FROM>`__ docs

``WORKDIR /app/COMPAS``
Effectively ``cd /app/COMPAS`` within the container.
`WORKDIR <https://docs.docker.com/engine/reference/builder/#workdir>`__
docs

``RUN apt-get update && apt-get install -y ...``
Install the required dependencies.
``-y`` so there's no prompt to install any of the packages.
``update`` and ``install`` are in the same layer because now if there
are any updates, it will force all of the dependencies to be
re-installed
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

``RUN pip3 install numpy``
Install numpy.
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

``COPY src/ src/``
Copy ``./src/`` directory from the local machine to ``./src`` in the
container (remembering that ``WORKDIR`` changes the cwd).
`COPY <https://docs.docker.com/engine/reference/builder/#copy>`__ docs

``RUN mkdir obj bin logs``
Create the directories required by COMPAS.
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

``ENV COMPAS_ROOT_DIR /app/COMPAS``
Set the required environment variable(s).
`ENV <https://docs.docker.com/engine/reference/builder/#env>`__ docs

``RUN cd src && make -f Makefile.docker -j $(nproc)``
Make COMPAS using a specific makefile (more below) and as many cores
as possible.
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

Dockerfiles will usually end with a ``CMD`` directive that specifies
what command should run when the container is started.
COMPAS doesn't have a ``CMD`` directive because some users will want
to run the executable directly and some will want to use
``runSubmit.``.
`CMD <https://docs.docker.com/engine/reference/builder/#cmd>`__ docs

Makefile.docker


A separate makefile is required for Docker in this scenario for two
reasons.

#. To separate compiled files from source files
#. To prevent the usage of ``-march=native``

``-march=native`` is a fantastic optimisation for users who compile
and run COMPAS on the same machine, however it causes fatal errors when
running COMPAS on a machine that it was not compiled for.
`Docs <https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html>`__ for
``-march``.

This selects the CPU to generate code for at compilation time by
determining the processor type of the **compiling machine**.

Using -march=native enables all instruction subsets supported by the
local machine (hence the result might not run on different
machines).



