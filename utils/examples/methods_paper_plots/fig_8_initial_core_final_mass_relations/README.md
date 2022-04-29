# How to recreate/edit the initial-core-final mass relation plot from the COMPAS methods paper

This folder contains everything you need to reproduce (or edit!) the initial-core-final mass relation plot from the methods paper. You'll need to do the following:

1. Run `python create_fig_8_grids.py` to create a grid of stars for the plot
2. Run three different python submit files: `python pythonSubmitDefaults.py`,`python pythonSubmitRapid.py` and `python pythonSubmitMandelMueller.py` to run COMPAS for these three grids(note this is probably going to take quite a few minutes!)
3. Open `make_fig_8.ipynb` in Jupyter Lab/notebook and run the whole thing to create the figure (and learn how to change it)

For reference, the changes from the *default* pythonSubmit.py are just:

- change the output file name (to distinguish the 3 grids)
- change the logfiles definitions so that only the variables needed are included
- change grid name
- (for some) change the remnant mass prescription