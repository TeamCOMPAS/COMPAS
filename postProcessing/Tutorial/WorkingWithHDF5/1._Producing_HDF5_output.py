# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Introduction
#
# Here we show an example of how to run a basic COMPAS simulation in order to produce output in an HDF5 format. 
#
# *Note:* It is required to have some COMPAS output to complete the rest of the post-processing notebooks, but if you already have an output file, you can skip this section.

# ## To run COMPAS, you need to set the environment variable COMPAS_ROOT_DIR to the top level COMPAS directory. Ensure that you have compiled COMPAS, and that there is a COMPAS executable.
#
# ## We will store the output in the Tutorial directory for easy access later on



# ## For simplicity, we will run COMPAS on all default settings except for the number of binaries produced. Here, we will run 1k binaries, which is relatively small. Later, we will look at production of double compact objects (DCOs), but these objects are sufficiently rare that in order to study them, we will need to run a simulation of ~1M binaries, but this can take some time depending on your hardware resources.
#
# ## *Note:* This will produce output in the COMPAS_Output directory only the first time it is run. If you run this multiple times, new directories will be created as COMPAS_Output_X. Ensure that you are using the data in the desired directory.

# To run terminal commands in a jupyter notebook, prepend the command with an exclamation mark
# !COMPAS -n 1000 -o {outputDir}

# ## Non-default settings can be seen by running `COMPAS --help`

# !COMPAS --help
