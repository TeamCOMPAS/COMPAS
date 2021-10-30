Cosmic integration
==================


This section contains a series of tutorials on how to use and get the most out of the COMPAS Cosmic Integration tools. The COMPAS
Cosmic Integration tools are based on work done in :cite:`Neijssel2019`.

The duration from the birth of the binary until the merger as a double compact object (DCO) can range from a few million years 
(lifetime of the stars) to more than 100 Gigayears depending on the evolution of the system.

This means two things:

    #. DCOs merge at different redshifts
    #. Multiple DCOs merging at a specific redshift could have formed at different times

We thus need to know the star formation that went into forming a single system. However, the star formation rate is non-constant over 
the lifetime of the universe. Furthermore, star formation is heavily dependent on the metallicity of the star forming gas, which also 
changes over the lifetime of the universe. Combined, we call this the metallicity-specific star formation rate (MSSFR).

The cosmic-integration pipeline tries to predict the population of DCO mergers assuming a model for the MSSFR. Here we show how to use 
the combination of scripts: we show the purpose of each script in the pipeline, and how to combine them to derive our results. These 
notes show how to call the functions, and how to construct the pipeline: they do not offer any derivations. We assume that you have a 
COMPAS data set in ``HDF5`` file format. If not, see [COMPAS h5 files.](../H5/0_H5Overview.ipynb). 

Note that although these pipelines are tailored for a COMPAS simulation, with some slight adjustments one could use these same 
pipelines for other (mock-)simulations, provided all the information is given (see section 7 of `class COMPAS <./notebooks/cosmicIntegration/1_ClassCOMPAS.ipynb>`__).


.. toctree::
   :maxdepth: 1

   Class COMPAS: setting the data. <./notebooks/cosmicIntegration/1_ClassCOMPAS.ipynb>
   Class MSSFR: choose a model and plot it. <./notebooks/cosmicIntegration/2_MSSFR-or-SFRD-prescriptions.ipynb>
   Selection effects module: choose the sensitivity and estimate detection probability. <./notebooks/cosmicIntegration/3_SelectionEffects.ipynb>
   Rate at a single redshift: combine data, MSSFR, and selection effects. <./notebooks/cosmicIntegration/4_MergersAtSingleRedshift.ipynb>
   Class CosmicIntegrator: rate as a function of redshift. <./notebooks/cosmicIntegration/5_MergersAsFunctionOfRedhisft.ipynb>
   Fast cosmic integration module: do it all in one go! <./notebooks/cosmicIntegration/6_FastCosmicIntegrator.ipynb>


If you make use of any these pipelines, we would appreciate you citing the following papers:

    - For Cosmic integration: :cite:`Neijssel2019` |br|
    - For Selection effects: :cite:`Barrett2018` 


