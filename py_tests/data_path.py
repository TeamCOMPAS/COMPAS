import os


DETAILED_EVOLUTION_PATH = os.path.join(
        os.path.dirname(__file__),
        "../misc/examples/methods_paper_plots/detailed_evolution"
        )

def get_data_path():
    DATA_PATH = os.path.join(DETAILED_EVOLUTION_PATH, "COMPAS_Output/COMPAS_Output.h5")

    if not os.path.exists(DATA_PATH):  # Check if path exists
        curr_dir = os.getcwd()
        os.chdir(DETAILED_EVOLUTION_PATH)
        os.system('python runSubmitDemo.py')
        os.system(f"cd {curr_dir}")
        print("Generated COMPAS test data")

    return DATA_PATH
