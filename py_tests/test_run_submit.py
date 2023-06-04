from typing import List

from yaml.scanner import ScannerError

from compas_python_utils.preprocessing.runSubmit import (
    DEFAULT_CONFIG_FILE, runSubmit
)


def test_run_submit(tmp_path, capsys):
    """Capture the stdout from runSubmit(execute=False) and check that it is as expected"""

    runSubmit(execute=False)
    _check_if_expected_kwgs_in_stdout(capsys, ["COMPAS", "output-path"])

    temp_ini, ini_contents = _make_tmp_ini(tmp_path)
    try:
        runSubmit(execute=False, cli_args=[temp_ini])
    except ScannerError:
        raise Exception(
            f"runSubmit failed to parse the ini file. Ini file:\n{ini_contents}"
        )
    _check_if_expected_kwgs_in_stdout(capsys, ["log-level", "debug-to_file"])


def _make_tmp_ini(tmp_path) -> str:
    """
    Make a temporary ini file for testing
    uncomment all default options to ensure they can be parsed
    """
    with open(DEFAULT_CONFIG_FILE, "r") as f:
        lines = f.read()
        lines = lines.replace("\n#", "\n")
    ini = f"{tmp_path}/test.ini"
    with open(ini, "w") as f:
        f.writelines(lines)
    return ini, lines


def _check_if_expected_kwgs_in_stdout(capsys, expected_kwgs: List[str]):
    """Check that all expected_kwgs are in stdout"""
    stdout = capsys.readouterr().out
    for expected_kwarg in expected_kwgs:
        assert (
            expected_kwarg in stdout,
            f"Expected kwarg {expected_kwarg} not found in sout: {stdout}",
        )
