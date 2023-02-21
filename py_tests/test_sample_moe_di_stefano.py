import unittest
import shutil, os
from compas_python_utils.preprocessing.sampleMoeDiStefano import main


class TestSampleMoeDiStefano(unittest.TestCase):

    def setUp(self) -> None:
        self.outdir = "out_test"
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self) -> None:
        shutil.rmtree(self.outdir)

    def test_main(self):
        fname = os.path.join(self.outdir, "gridfn.txt")
        # This can take ~1 minute to run
        main(cli_args=['--gridname', fname, '--nSamples', '10'])
        self.assertTrue(os.path.exists(fname))
        # TODO: ask reinhold if there is a way to check the contents
        # TODO: ask reinhold how to allow createParameterDistributionsAndSampler to be run with smaller grid-points