import os

import numpy as np
import pytest
from conftest import get_compas_data
from deepdiff import DeepDiff

from compas_python_utils import h5sample


def test_sample(tmp_path, example_compas_output_path):
    """Test that h5sample can sample a file"""
    np.random.seed(0)
    init_file = example_compas_output_path

    test_kwargs = [
        dict(n=100, replace=True),
        dict(frac=0.5, replace=False),
        dict(frac=1.2, replace=True),
        dict(seed_group="BSE_Double_Compact_Objects", frac=2.0, replace=True),
    ]

    for kwg in test_kwargs:
        # creating new h5 file to copy contents to
        new_file = f"{tmp_path}/sampled_compas_out.h5"
        if os.path.exists(new_file):
            os.remove(new_file)
        h5sample.sample_h5(init_file, new_file, **kwg)
        assert os.path.exists(new_file), f"File {new_file} does not exist"
        init_data = get_compas_data(init_file)
        new_data = get_compas_data(new_file)
        diff = DeepDiff(init_data, new_data)
        assert (
            len(diff) > 0
        ), f"The sampled file is the same as the original when using kwgs: {kwg}"


def test_argparser(capsys):
    with pytest.raises(SystemExit):
        h5sample.parse_args(["-h"])
    output = capsys.readouterr().out
    assert "Sample an COMPAS h5 file" in output
