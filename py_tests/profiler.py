from compas_python_utils.cosmic_integration import FastCosmicIntegration
import time

if __name__ == '__main__':
    t0 = time.time()
    example_compas_output_path = "/home/avaj040/Documents/projects/data/COMPAS_data/jeff_data/h5out_5M.h5"
    (
        detection_rate,
        formation_rate,
        merger_rate,
        redshifts,
        COMPAS,
    ) =  FastCosmicIntegration.find_detection_rate(
        path=example_compas_output_path,
    )
    runtime = time.time() - t0
    print(f"runtime: {runtime:.2f}s")