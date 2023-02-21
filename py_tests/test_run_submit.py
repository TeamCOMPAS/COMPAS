import unittest

#import os.path
#from pytest_utils import get_compas_output_path, time_func_and_capture_sout
## BEGIN Delete later
#import sys
#sys.path.append( os.path.expandvars(os.environ['COMPAS_ROOT_DIR']))
#import subprocess
## END Delete later
#from compas_python_utils.preprocessing.runSubmit import main as runSubmitMain
#
#OUTDIR = os.path.join(os.path.dirname(__file__), "output_test/")
OUTDIR = "output_test"
print(OUTDIR)

from compas_python_utils.preprocessing.runSubmit import runSubmit
from pytest_utils import time_func_and_capture_sout

class TestRunSubmit(unittest.TestCase):
    
    def setUp(self) -> None:
        if not os.path.exists(OUTDIR):
            os.mkdir(OUTDIR)


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


#        _, runtime, sout = time_func_and_capture_sout(
#            runSubmitMain,
#            output_directory=OUTDIR
#        )
#
#        #print(sout)
#        print(runtime)
#
#        xx = os.listdir(OUTDIR)
#        print(xx)
#
#
#
#
#        print('huh?')
#        self.assertLess(runtime, 1)
#        self.assertTrue(os.path.exists(os.path.join(OUTDIR, "COMPAS_Output.h5")))
#
#
#if __name__ == "__main__":
#    test = TestRunSubmit()
#    test.test_run_submit()
