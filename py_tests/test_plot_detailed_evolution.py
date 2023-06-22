import os.path
import time

from compas_python_utils.detailed_evolution_plotter import plot_detailed_evolution


def test_plotter(example_compas_output_path, capsys, test_archive_dir):
    data_path = example_compas_output_path
    bse_detailed_out_path = os.path.join(
        os.path.dirname(data_path), "Detailed_Output/BSE_Detailed_Output_0.h5"
    )
    t0 = time.time()
    plot_detailed_evolution.run_main_plotter(
        bse_detailed_out_path, outdir=test_archive_dir, show=False
    )
    runtime = time.time() - t0

    sout = capsys.readouterr().out

    assert runtime < 30
    assert os.path.exists(os.path.join(test_archive_dir, "vanDenHeuvelPlot.eps"))
    assert os.path.exists(os.path.join(test_archive_dir, "detailedEvolutionPlot.eps"))

    with open(os.path.join(test_archive_dir, "test_plotter.log"), "w") as f:
        f.write(sout)
