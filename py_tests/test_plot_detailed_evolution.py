import os.path

from compas_python_utils.detailed_evolution_plotter.plot_detailed_evolution import run_main_plotter
from contextlib import redirect_stdout

from data_path import get_data_path
import unittest


class TestPlotDetailedEvolution(unittest.TestCase):

    def test_plotter(self):
        data_path = get_data_path()
        bse_detailed_out_path = os.path.join(
            os.path.dirname(data_path),
            "Detailed_Output/BSE_Detailed_Output_0.h5"
        )
        run_main_plotter(bse_detailed_out_path)