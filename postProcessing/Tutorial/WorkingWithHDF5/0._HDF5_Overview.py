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
# COMPAS simulations produce all output by default in the form of an [HDF5 file](https://www.hdfgroup.org/solutions/hdf5/), which is a compact and memory efficient, but non-human-readable data file format. In order to interact with these output files, we use the python module `h5py`. 
#
# Within each file are a variety of HDF5 Groups, representing a specific common event in binary evolution, e.g Roche-Lobe Overflow or Supernovae. These are described throughout the post-processing jupyter notebooks.

# ## Material
#
# ### [1. Producing HDF5 output:](./1._Producing_HDF5_output.ipynb)
# How to run default COMPAS and produce a simple output file.
#         
#         
# ### [2. Reading HDF5 files:](./2._Reading_HDF5_files.ipynb)
# The basics and syntax of loading an HDF5 file.
#         
#         
# ### [3. Re-writing HDF5 files:](./3._Rewriting_HDF5.ipynb)
# How to rewrite/reduce the HDF5 data.


