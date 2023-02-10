import os.path
import unittest
import io
from contextlib import redirect_stdout

from pytest_utils import get_compas_output_path, time_func_and_capture_sout
from compas_python_utils.detailed_evolution_plotter.plot_detailed_evolution import run_main_plotter

OUTDIR = os.path.join(os.path.dirname(__file__), "output_test")
class TestPlotDetailedEvolution(unittest.TestCase):

    def setUp(self) -> None:
        if not os.path.exists(OUTDIR):
            os.mkdir(OUTDIR)

    def test_plotter(self):
        data_path = get_compas_output_path()
        bse_detailed_out_path = os.path.join(
            os.path.dirname(data_path),
            "Detailed_Output/BSE_Detailed_Output_0.h5"
        )
        _, runtime, sout = time_func_and_capture_sout(
            run_main_plotter,
            bse_detailed_out_path, outdir=OUTDIR, show=False
        )
        self.assertLess(runtime, 30)
        self.assertTrue(os.path.exists(os.path.join(OUTDIR, "vanDenHeuvelPlot.eps")))
        self.assertTrue(os.path.exists(os.path.join(OUTDIR, "detailedEvolutionPlot.eps")))

        if not os.path.exists(OUTDIR):
            os.mkdir(OUTDIR)
        with open(os.path.join(OUTDIR, "test_plotter.log"), "w") as f:
            f.write(sout)
