Running the COMPAS Docker image
===============================

COMPAS can be configured as usual via command line arguments passed to the COMPAS executable or via a ``runSubmit.py`` file in the 
``Docker`` environment.


via runSubmit.py
-------------------

To run COMPAS via a ``runSubmit.py`` file, type::

    docker run                                                  \
        --rm                                                    \
        -it                                                     \
        -v $(pwd)/compas-input:/app/COMPAS/config               \
        -v $(pwd)/compas-logs:/app/COMPAS/logs                  \
        -v $(pwd)/runSubmit.py:/app/starts/runSubmit.py   \
        -e COMPAS_EXECUTABLE_PATH=/app/COMPAS/bin/COMPAS        \
        -e COMPAS_INPUT_DIR_PATH=/app/COMPAS/config             \
        -e COMPAS_LOGS_OUTPUT_DIR_PATH=/app/COMPAS/logs         \
        teamcompas/compas                                       \
        python3 /app/starts/runSubmit.py                     


NOTE: if you decide to execute using ``runSubmit.py``, you will need 
a ``compasConfigDefault.yaml``  file in the same directory. This file 
can be find in the same directory as the ``runSubmit.py``, and contains
the default COMPAS choices for stellar and binary physics. These choices
can be changed by modifying the options availabe in the ``compasConfigDefault.yaml`` 
file.

Breaking down this command:

**docker run** |br|
creates a container.

**--rm** |br|
destroy the container once it finishes running the command\ [#f1]_.

**-it** |br|
short for [-i and -t] - provides an interactive terminal\ [#f2]_.

**-v <path-on-host>:<path-in-container>** |br|
mount ``<path-on-host>`` to ``<path-in-container>``\ [#f3]_. |br|

This time we not only want to read the COMPAS input files (i.e. grid file and/or logfile-definitions file) on the
host from the container, and get the output from COMPAS in the container onto the host machine, we also want to 
supply a ``runSubmit.py`` to the container from the host machine.

**-e VAR_NAME=value** |br|
set the environment variable ``VAR_VAME`` to `value`\ [#f4]_.

**teamcompas/compas** |br|
the image to run.

**python3 /app/starts/runSubmit.py** |br|
the command to run when the container starts.


via the command line
--------------------

To run the COMPAS executable from the command line (i.e. without ``runSubmit.py``), type::

    docker run                                      \
        --rm                                        \
        -it                                         \
        -v $(pwd)/compas-input:/app/COMPAS/config   \
        -v $(pwd)/compas-logs:/app/COMPAS/logs      \
        teamcompas/compas                           \
        bin/COMPAS                                  \
        --number-of-systems=5                       \
        --output-path=/app/COMPAS/logs


Breaking down this command:

**docker run** |br|
creates a container\ [#f5]_.

**--rm** |br|
destroy the container once it finishes running the command\ [#f1]_.

**-it** |br|
short for [-i and -t] - provides an interactive terminal\ [#f2]_.

**-v <path-on-host>:<path-in-container>** |br|
mount ``<path-on-host>`` to ``<path-in-container>``\ [#f3]_. |br|

In this instance, make it so |br|
   `$(pwd)/compas-input` on my machine is the same as `/app/COMPAS/config` inside the container. |br|
   `$(pwd)/compas-logs` on my machine is the same as `/app/COMPAS/logs` inside the container.

**teamcompas/compas** |br|
the image to run.

**bin/COMPAS** |br|
the command to run when the container starts.

**--number-of-systems=5** |br|
the flag to set the number of binaries.

**--output-path=/app/COMPAS/logs** |br|
forces logs to go to the directory that is mapped to the host machine.



Environment variables
---------------------

Three new environment variables are used in ``runSubmit.py``.  These environment variables are used primarily in the ``Docker``
environment, and are non-breaking changes (i.e. benign to other environments).

``COMPAS_EXECUTABLE_PATH`` specifies where ``runSubmit.py`` looks for the COMPAS executable. This override exists purely for 
ease-of-use from the command line.

`COMPAS_LOGS_OUTPUT_DIR_PATH` specifies where COMPAS output log files are created. The override exists because the mounted directory 
(option `-v`) is created before COMPAS runs. COMPAS sees that the directory where it's supposed to put logs already exists, so it 
creates a different (i.e. non-mapped) directory for the output log files.

`COMPAS_INPUT_DIR_PATH` specifies where input files (such as the ``grid`` file, or ``logfile-definitions`` file are located.


Detached mode
-------------

The ``docker run`` examples above use the ``-it`` option.
To run multiple instances of COMPAS, an alternative is to use detached mode (`-d`)\ [#f6]_. In detached mode, containers are run in 
the background of the current shell - they do not receive input or display output.

Typing::

    docker run --rm -d -v $(pwd)/compas-logs/run_0:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_01.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    docker run --rm -d -v $(pwd)/compas-logs/run_1:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_02.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    docker run --rm -d -v $(pwd)/compas-logs/run_2:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_03.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    docker run --rm -d -v $(pwd)/compas-logs/run_3:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_04.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &

runs 4 separate instances of COMPAS, each with its own ``runSubmit.py`` file and logging directory, and all local console output supressed.

To see the console output of detached containers to check progress, first type::

  docker ps

to get the container id of interest, then type::

    docker logs container_id


.. rubric:: Footnotes

.. [#f1] https://docs.docker.com/engine/reference/run/#clean-up---rm
.. [#f2] https://docs.docker.com/engine/reference/run/#foreground
.. [#f3] https://docs.docker.com/storage/bind-mounts/
.. [#f4] https://docs.docker.com/engine/reference/run/#env-environment-variables
.. [#f5] https://docs.docker.com/engine/reference/run/
.. [#f6] https://docs.docker.com/engine/reference/run/#detached--d

   
