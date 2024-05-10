import os

import numpy as np
from conftest import get_compas_data
from deepdiff import DeepDiff

from compas_python_utils import h5sample
import pytest


@pytest.mark.parametrize(
"test_kwargs",
[
    dict(n=100, replace=True),
    dict(frac=0.5, replace=False),
    dict(frac=1.2, replace=True),
    dict(seed_group="BSE_Double_Compact_Objects", frac=2.0, replace=True),
],
)
def test_sample(test_kwargs, tmp_path, fake_compas_output):
    """Test that h5sample can sample a file"""
    np.random.seed(0)
    init_file = fake_compas_output

    new_file = f"{tmp_path}/sampled_compas_out.h5"
    if os.path.exists(new_file):
        os.remove(new_file)
    h5sample.sample_h5(init_file, new_file, **test_kwargs)
    assert os.path.exists(new_file), f"File {new_file} does not exist"
    init_data = get_compas_data(init_file)
    new_data = get_compas_data(new_file)
    diff = DeepDiff(init_data, new_data)
    assert (
        len(diff) > 0
    ), f"The sampled file is the same as the original when using kwgs: {test_kwargs}"


def test_argparser(capsys):
    with pytest.raises(SystemExit):
        h5sample.parse_args(["-h"])
    output = capsys.readouterr().out
    assert "Sample an COMPAS h5 file" in output
