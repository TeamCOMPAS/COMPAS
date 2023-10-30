Post-processing tools
=====================

The COMPAS suite includes some useful post-processing tools that are located in the `postProcessing` directory.
In here, you will find scripts used to parse COMPAS ``HDF5``\ [#f1]_ output and calculate rates/weights/samples
to place the results of the simulation into a scientific context.


HDF5 tools
----------

.. toctree::
   :maxdepth: 1

   HDF5 basics <hdf5/post-processing-hdf5-info>
   H5 Copying/concatenating <hdf5/post-processing-h5copy>
   H5 View  <hdf5/post-processing-h5view>
   H5 Sample <hdf5/post-processing-h5sample>
   Python Demo <hdf5/WorkingWithHDF5.ipynb>

Post-processing Demos
---------------------

.. toctree::
    :maxdepth: 1

    Basic COMPAS Data Analysis <notebooks/DataAnalysis.ipynb>
    Cosmic Integration <notebooks/CosmicIntegration.ipynb>
    Spin prescriptions <notebooks/spin_prescriptions/calculate_black_hole_spin_pop_synth.ipynb>
    CHE Tutorial <CHE_paper_tutorial/CHE_evolution_demo_ANSWERS.ipynb>


.. rubric:: Footnotes

.. [#f1] https://www.hdfgroup.org/