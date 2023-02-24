import unittest

from compas_python_utils.preprocessing.runSubmit import runSubmit
from pytest_utils import time_func_and_capture_sout

class TestRunSubmit(unittest.TestCase):
    def test_run_submit(self):
        """ Capture the stdout from runSubmit(execute=False) and check that it is as expected """
        _, _, sout = time_func_and_capture_sout(runSubmit, execute=False)
        # test that sout has some of the following kwargs
        expected_kwargs = [
            "COMPAS",
            "debug-to-file",
            "detailed-output",
            "quiet",
        ]
        for kwarg in expected_kwargs:
            self.assertTrue(
                kwarg in sout,
                f"Expected kwarg {kwarg} not found in sout: {sout}"
            )
