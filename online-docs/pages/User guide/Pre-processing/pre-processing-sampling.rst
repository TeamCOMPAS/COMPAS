Sampling in COMPAS
==================



Here are some basic instructions for efficient sampling of the COMPAS
input parameters, using the python sampling package Stroopwafel.
Below that are instructions for how to sample from the correlated 
parameter distributions outlined in Moe & DiStefano 2017. 

Note that the intended Stroopwafel functionality for "Adaptive
Importance Sampling" is not yet implemented, but is currently in
development.

Requirements




If you have not already, you will need to install Stroopwafel. If you
have admin rights, Stroopwafel can be installed on your system with
``pip install stroopwafel``.

Instructions




To use Stroopwafel sampling, copy
``preProcessing/stroopwafelInterface.py`` into your working directory.

Settings
~~~~~~~~

NOTE: This sampling method is currently being updated as part of an upgrade
in our method to parse user-defined options. We plan to address this shortly. 
Please bear with us and contact the COMPAS team if an urgent solution is needed.

1. runSubmit


If you are running COMPAS on default settings, skip this section.

If you have many non-default COMPAS arguments, you may want to set
them in the ``compasConfigDefault.yaml``, that is read and executed by the 
``runSubmit.py`` file in the same directory. For now, the file must
be named this way and placed in the same directory as the ``stroopwafelInterface.py``
file.

A configurable runSubmit file can be found in the ``preProcessing/``
directory.

Set your desired options, then set the ``userunSubmit`` parameter to ``True``
in the ``stroopwafelInterface.py``.

2. Stroopwafel inputs


The lines below ``userunSubmit`` represent stroopwafel inputs.
These are treated as
defaults, but can be overridden by command-line arguments to
stroopwafel.
See ``python3 stroopwafelInterface.py --help``.

``num_systems`` is the total number of binaries you would like to
evolve.
This value overrides the value set in the ``compasConfigDefault.yaml`` file.

``num_cores`` is the number of cores you would like to use to
parallelize your run. More cores means your run will finish sooner, but
may reduce your ability to run other tasks while you wait. On linux
systems, the command ``echo $(nproc)`` will tell you how many (virtual)
CPUs you have available.

``num_per_core`` is the number of systems to run on a core at a given
time. This translates to the number of systems in a single batch file.
This is more relevant for adaptive importance sampling.

``mc_only`` specifies if you would like to do naive MC sampling only.
Currently, this option must be set to True

``run_on_hpc`` specifies if you are running on a High-Performance
Computer (HPC).
If so, see `docs/compasHPC.md <compasHPC.md>`__ for assistance.

``output_folder`` a string specifying the output folder. Relative paths
will be appended onto the current directory path.

``output_filename`` a string specifying the name of the output samples
file.

``debug`` whether to print the COMPAS output/error.

3. Sampling parameters


Sampled parameters will be combined into grid files which COMPAS then
reads in.
Users should choose which parameters they would like to be sampled
over, as well as
the relevant distributions.

See the `COMPAS
Documentation <https://github.com/TeamCOMPAS/COMPAS/blob/Documentation/COMPAS_Documentation.pdf>`__
for details on which sets of
parameters are allowed/required.

See `Stroopwafel
Documentation <https://github.com/lokiysh/stroopwafel>`__ for details on
which distributions are available.

4. Run Stroopwafel


When your satisfied with your settings, simply run with
``python3 stroopwafelInterface.py``. The output will be collected into
batch containers in your output folder.
To postprocess the output, see
`getting\_started.md <getting_started.md>`__



Moe & DiStefano
~~~~~~~~~~~~~~~

To sample from the Moe & DiStefano 2017 distributions, the sampler script
can be found in ``preProcessing/sampleMoeDiStefano.py``. As described in the 
header, only the number of systems and upper and lower mass bounds can be
set (the parameter correlations break if you try to set the other bounds).
These values are set at the bottom of the script.

