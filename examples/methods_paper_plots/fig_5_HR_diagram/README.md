# How to recreate/edit the HR diagram from the COMPAS methods paper

This folder contains everything you need to reproduce (or edit!) the HRD from the methods paper. You'll need to do the following:

1. Run `python create_fig_5_grid.py` to create a grid of stars for the plot (you can change this to a custom range of masses or metallicities if you like)
2. Run `python pythonSubmit.py` to run COMPAS for this grid
3. Open `make_fig_5.ipynb` in Jupyter Lab/notebook and run the whole thing to create the figure (and learn how to change it)

For reference, the changes from the *default* pythonSubmit.py are just:

- use detailed output
- use smaller time steps (to get smoother lines for the plot)