import pytest
import os
from compas_python_utils.preprocessing import sampleMoeDiStefano


@pytest.mark.slow
def test_generate_MoeDiStefano_gridfile(tmp_path):
    fname = os.path.join(tmp_path, "gridfn.txt")
    # This can take ~1 minute to run
    sampleMoeDiStefano.main(cli_args=['--gridname', fname, '--nSamples', '10'])
    assert os.path.exists(fname)
    # TODO: ask reinhold if there is a way to check the contents
    # TODO: ask reinhold how to allow createParameterDistributionsAndSampler to be run with smaller grid-points
