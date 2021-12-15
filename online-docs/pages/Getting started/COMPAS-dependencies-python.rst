Installing Python
=================

Python and some selected libraries are required for interfacing with the code, and also for post-processing. We recommend using ``python3``. The 
``matplotlib`` and ``numpy`` libraries should also be installed. The libraries ``scipy``, ``astropy``, and ``pandas`` are also used in some other scripts.

First check if you have ``python3`` installed. If you do, typing the following should give you the version number::

    python3 --version

If you do not have ``python3`` installed, install it by following the instructions below for your OS:

- For macOS, We recommend installing ``Python`` and its libraries using MacPorts. You can follow the instructions on `MacPorts Python installation on Mac <https://astrofrog.github.io/macports-python/>`__.
- For Linux, install `python3` using your package manager (e.g. in Ubuntu, run `sudo apt-get install python3`). We recommend installing the required python libraries using the package installer ``pip``. E.g. To install ``numpy``, run `pip install numpy`; for ``h5py``, run `pip install h5py`.

