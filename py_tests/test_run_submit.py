from compas_python_utils.preprocessing.runSubmit import runSubmit


def test_run_submit(capsys):
    """ Capture the stdout from runSubmit(execute=False) and check that it is as expected """
    runSubmit(execute=False)
    # test that sout has some of the following kwargs
    expected_kwargs = [
        "COMPAS",
        "debug-to-file",
        "detailed-output",
        "quiet",
    ]
    sout = capsys.readouterr().out
    for kwarg in expected_kwargs:
        msg = f"Expected kwarg {kwarg} not found in sout: {sout}"
        assert kwarg in sout, msg
