import os


DATA_PATH = os.path.join(
        os.path.dirname(__file__),
        "../utils/examples/methods_paper_plots/detailed_evolution/COMPAS_Output/COMPAS_Output.h5"
        )

assert os.path.exists(DATA_PATH)  # Check if path exists
