#make clean
#make -j 4
$COMPAS_ROOT_DIR/src/COMPAS -v
cd $COMPAS_ROOT_DIR/examples/methods_paper_plots/detailed_evolution/
rm -r COMPAS_Output
python3 pythonSubmitDemo.py
python3 detailed_evol_plotter.py