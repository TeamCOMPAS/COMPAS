=================
COMPAS and Docker
=================

`Docker <https://www.docker.com/>`__ has been integrated into COMPAS to simplify the setup of the COMPAS deployment environment.
Instead of manually installing and configuring various libraries and tools (e.g., python/pip, numpy, g++, boost) which can differ across operating systems and toolchains, users can install Docker and run COMPAS with a single command.
This also allows users to run COMPAS on cloud platforms like `AWS EC2 <https://aws.amazon.com/ec2/>`__ or `Google Compute Engine <https://cloud.google.com/compute>`__, where hundreds of cores can be provisioned without manual environment configuration.

-----
Usage
-----

**Note:** This section assumes `Docker <https://www.docker.com/>`__ is installed and running.
For Windows and MacOS users, see `here <https://www.docker.com/products/docker-desktop>`__.

Installing
~~~~~~~~~~

Retrieve the latest compiled version of COMPAS (dev branch) by running:

.. code-block:: bash

    docker pull teamcompas/compas:latest

To use other versions, add a version `tag <https://docs.docker.com/engine/reference/commandline/tag/>`__. For example, for COMPAS version 2.12.0:

.. code-block:: bash

    docker pull teamcompas/compas:2.12.0

To see all available versions, visit the TeamCOMPAS Docker Hub page `here <https://hub.docker.com/u/teamcompas>`__.

Running
~~~~~~~

COMPAS can be configured via command line arguments passed to the COMPAS executable or via a `runSubmit.py` file.

Running `runSubmit.py`
^^^^^^^^^^^^^^^^^^^^^^

To run COMPAS using a `runSubmit.py` file, use the following command:

.. code-block:: bash

    docker run
    --rm
    -it
    -v $(pwd)/compas-logs:/app/COMPAS/logs
    -v $(pwd)/runSubmit.py:/app/starts/runSubmit.py
    -e COMPAS_EXECUTABLE_PATH=/app/COMPAS/bin/COMPAS
    -e COMPAS_LOGS_OUTPUT_DIR_PATH=/app/COMPAS/logs
    teamcompas/compas
    python3 /app/starts/runSubmit.py


Breaking down this command:
    - `docker run`: Creates a container.
    - `--rm`: Cleans up and destroys the container once it finishes running the command.
    - `-it`: Provides an interactive terminal.
    - `-v <path-on-host>:<path-in-container>`: Mounts `<path-on-host>` to `<path-in-container>`.
    - `-e VAR_NAME=value`: Sets the environment variable `VAR_NAME` to `value`.
    - `teamcompas/compas`: The image to run.
    - `python3 /app/starts/runSubmit.py`: The command to run when the container starts.

**Note:** If using `runSubmit.py`, ensure a `compasConfigDefault.yaml` file is in the same directory. This file contains default COMPAS choices for stellar and binary physics, which can be modified.

Running the COMPAS executable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run the COMPAS executable directly:

.. code-block:: bash

    docker run --rm -it
    -v $(pwd)/compas-logs:/app/COMPAS/logs teamcompas/compas bin/COMPAS
    --number-of-binaries=5
    --outputPath=/app/COMPAS/logs


Breaking down this command:
    - `docker run`: Creates a container.
    - `--rm`: Cleans up and destroys the container once it finishes running the command.
    - `-it`: Provides an interactive terminal.
    - `-v <path-on-host>:<path-in-container>`: Mounts `<path-on-host>` to `<path-in-container>`.
    - `teamcompas/compas`: The image to run.
    - `bin/COMPAS`: The command to run when the container starts.
    - `--number-of-binaries`: Sets the number of binaries.
    - `--outputPath /app/COMPAS/logs`: Forces logs to go to the specified directory.

More information on `docker run` can be found `here <https://docs.docker.com/engine/reference/run/>`__.

    - **Note 1:**
        Two new environment variables have been added for `runSubmit.py` only:
            - `COMPAS_EXECUTABLE_PATH`: Overrides where `runSubmit.py` looks for the compiled COMPAS.
            - `COMPAS_LOGS_OUTPUT_DIR_PATH`: Overrides where logs are placed.

    - **Note 2:**
        The `docker run ...` examples use the `-it` options.
        For running multiple instances of COMPAS, use `detached mode <https://docs.docker.com/engine/reference/run/#detached--d>`__ (`-d`).
        This hides all container output.

Example for running 4 instances of COMPAS:

.. code-block:: bash

    docker run --rm -d
    -v $(pwd)/compas-logs/run_0:/app/COMPAS/logs
    -v $(pwd)/runSubmitMMsolar_01.py:/app/starts/runSubmit.py
    teamcompas/compas python3 /app/starts/runSubmit.py

    &  docker run --rm -d
    -v $(pwd)/compas-logs/run_1:/app/COMPAS/logs
    -v $(pwd)/runSubmitMMsolar_02.py:/app/starts/runSubmit.py
    teamcompas/compas python3 /app/starts/runSubmit.py

    &  docker run --rm -d
    -v $(pwd)/compas-logs/run_2:/app/COMPAS/logs
    -v $(pwd)/runSubmitMMsolar_03.py:/app/starts/runSubmit.py
    teamcompas/compas python3 /app/starts/runSubmit.py

    &  docker run --rm -d
    -v $(pwd)/compas-logs/run_3:/app/COMPAS/logs
    -v $(pwd)/runSubmitMMsolar_04.py:/app/starts/runSubmit.py
    teamcompas/compas python3 /app/starts/runSubmit.py


To check the console output, use `docker logs <container_id>`. Get the container id by running `docker ps`.

-----
CI/CD
-----

Whenever a push is made to `TeamCOMPAS/dev <https://github.com/TeamCOMPAS/COMPAS/tree/dev>`__, a continuous deployment process automatically builds a new image and deploys it to DockerHub with a `tag` corresponding to the value of `VERSION_STRING` in `constants.h`.
Expect the latest COMPAS docker image to be available 5-10 minutes after pushing/merging.

The GitHub Actions configuration is in `.github/workflows/dockerhub-ci.yml <https://github.com/TeamCOMPAS/COMPAS/blob/dev/.github/workflows/dockerhub-ci.yml>`__.

----------
Bonus Info
----------

Dockerfile
^^^^^^^^^^

The `Dockerfile <https://docs.docker.com/engine/reference/builder/>`__ defines how the docker image is constructed.

Images are created as a combination of layers. Each layer is cached and only updated on subsequent builds if that layer changes.

The Dockerfile for COMPAS consists of 8 layers:
    - `FROM ubuntu:18.04`: Uses `Ubuntu 18.04 <https://hub.docker.com/_/ubuntu>`__ as a base.
    - `WORKDIR /app/COMPAS`: Changes the working directory to `/app/COMPAS`.
    - `RUN apt-get update && apt-get install -y ...`: Installs required dependencies.
    - `RUN pip3 install numpy`: Installs numpy.
    - `COPY src/ src/`: Copies the `./src/` directory from the local machine to `./src` in the container.
    - `RUN mkdir obj bin logs`: Creates the required directories.
    - `ENV COMPAS_ROOT_DIR /app/COMPAS`: Sets the required environment variable(s).
    - `RUN cd src && make -f Makefile.docker -j $(nproc)`: Compiles COMPAS using a specific makefile and as many cores as possible.

Dockerfiles usually end with a `CMD` directive specifying the command to run when the container starts. COMPAS does not have a `CMD` directive because some users will run the executable directly, while others will use `runSubmit.py`.

Makefile.docker
^^^^^^^^^^^^^^^

A separate makefile is required for Docker to:
    1. Separate compiled files from source files.
    2. Prevent the usage of `-march=native`.

`-march=native` optimizes for the compiling machine's CPU, causing errors when running on a different machine. More information on `-march` can be found `here <https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html>`__.