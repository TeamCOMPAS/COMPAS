echo "Run example COMPAS job using the config example_bbh_compas_config.yaml"
compas_run_submit example_bbh_compas_config.yaml > example_bbh.log
cat example_bbh.log
echo "Generating detailed evolution plot"
compas_plot_detailed_evolution "./COMPAS_Output/Detailed_Output/BSE_Detailed_Output_0.h5" --dont-show >> example_bbh.log
echo "Out files:"
ls -l


