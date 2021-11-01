Tutorial: simple COMPAS run
===========================

This tutorial assumes that you have already built the COMPAS executable as described in :doc:`../../Getting started/building-COMPAS`.

For this example you will need the python script ``pythonSubmitDemo.py``, which specifies all the program options (physics assumptions, 
output types) and runs COMPAS in the terminal. Although the primary functionality of COMPAS is to evolve a whole population of binary 
stars rapidly, for now, let's focus on evolving a single stellar system and examining the detailed output.

If you haven't yet defined the ``COMPAS_ROOT_DIR`` environment variable, do that now::

    export COMPAS_ROOT_DIR=path-to-compas

where `path-to-compas` should be replaced with the path to the parent directory of the COMPAS `src` directory. Depending upon your system,
for the ``export`` command to take effect, it may be necessary to either restart your session or execute the following command::

    source ~/.bashrc

To start, change to the ``examples/methods_paper_plots/detailed_evolution/`` directory::

  cd $COMPAS_ROOT_DIR/examples/methods_paper_plots/detailed_evolution/

where you will find the script ``pythonSubmitDemo.py`` for this demo.


.. toctree::
   :maxdepth: 1

   ./example-compas-run-grid
   ./example-compas-run-detailed-output
   
