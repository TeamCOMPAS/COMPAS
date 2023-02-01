#! /bin/sh

export COMPAS_ROOT_DIR=${GITHUB_WORKSPACE}
cd ${GITHUB_WORKSPACE}/misc/examples/methods_paper_plots/detailed_evolution
python3 runSubmitDemo.py
