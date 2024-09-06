#bin/bash
echo ">>>  GENERATING TEST COMPAS DATA <<<"
# if $COMPAS_EXECUTABLE_PATH not set, set to ../../src/COMPAS
COMPAS_EXECUTABLE_PATH=${COMPAS_EXECUTABLE_PATH:-../../src/COMPAS}
$COMPAS_EXECUTABLE_PATH \
  -n 2 \
  --initial-mass-1 35 \
  --initial-mass-2 31 \
  -a 3.5 \
  --random-seed 0 \
  --metallicity 0.001 \
  --detailed-output \
  > compas_run.log
cat compas_run.log
echo "Generating detailed evolution plot"
compas_plot_detailed_evolution "./COMPAS_Output/Detailed_Output/BSE_Detailed_Output_0.h5" --dont-show >> detailed_evolution.log
echo "Out files:"
ls -l
echo ">>> DONE <<<"
