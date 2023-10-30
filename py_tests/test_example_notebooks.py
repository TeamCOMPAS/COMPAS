import glob
from tqdm.auto import tqdm
import jupytext
import os
import nbformat
from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor
import pytest

HERE = os.path.dirname(__file__)
NB_DIR = os.path.join(HERE, "../online-docs/pages/User guide/Post-processing/notebooks/")

EXAMPLE_NB = [
    os.path.join(NB_DIR, "CosmicIntegration.py"),
    # add others here
]


@pytest.mark.webtest
def test_run_examples():
    """Test that all examples run without error"""
    for example_py_file in tqdm(EXAMPLE_NB, desc="Running examples"):
        ipynb = __convert_to_ipynb(example_py_file)
        success = __execute_ipynb(ipynb, execute_dir=NB_DIR)
        if not success:
            os.remove(ipynb)
            assert False, f"{example_py_file} failed to run"


def __convert_to_ipynb(py_file):
    """Convert a python file to a notebook"""
    ipynb_fn = py_file.replace(".py", ".ipynb")
    if os.path.exists(ipynb_fn):
        os.remove(ipynb_fn)
    nb = jupytext.read(py_file, fmt="py:light")
    jupytext.write(nb, ipynb_fn, fmt="ipynb")
    return ipynb_fn


def __execute_ipynb(notebook_filename: str, execute_dir: str = None) -> bool:
    """
    :param notebook_filename: path of notebook to process
    :return: bool if notebook-preprocessing successful/unsuccessful
    """
    success = True
    with open(notebook_filename) as f:
        notebook = nbformat.read(f, as_version=4)

    ep = ExecutePreprocessor(timeout=-1)
    try:
        ep.preprocess(notebook, {"metadata": {"path": execute_dir}})
    except CellExecutionError as e:
        print(
            f"Preprocessing {notebook_filename} failed:\n\n {e.traceback}"
        )
        success = False
    finally:
        with open(notebook_filename, mode="wt") as f:
            nbformat.write(notebook, f)
    return success
