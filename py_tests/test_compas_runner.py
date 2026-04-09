import os

from compas_python_utils.compas_runner import main, resolve_compas_executable, run_compas


def _make_fake_executable(path, contents=None):
    executable_contents = contents or "#!/usr/bin/env bash\nprintf 'fake-compas %s\\n' \"$*\"\n"
    path.write_text(executable_contents)
    os.chmod(path, 0o755)


def test_resolve_compas_executable_from_env(monkeypatch, tmp_path):
    fake_executable = tmp_path / "COMPAS"
    _make_fake_executable(fake_executable)

    monkeypatch.setenv("COMPAS_EXECUTABLE_PATH", str(fake_executable))
    monkeypatch.delenv("COMPAS_BUNDLE_ROOT", raising=False)

    assert resolve_compas_executable() == str(fake_executable)


def test_resolve_compas_executable_from_bundle_root(monkeypatch, tmp_path):
    bundle_root = tmp_path / "COMPAS-linux-x86_64"
    bundle_root.mkdir()
    launcher = bundle_root / "run_compas.sh"
    _make_fake_executable(launcher)

    monkeypatch.delenv("COMPAS_EXECUTABLE_PATH", raising=False)
    monkeypatch.setenv("COMPAS_BUNDLE_ROOT", str(bundle_root))

    assert resolve_compas_executable() == str(launcher)


def test_run_compas_executes_resolved_binary(monkeypatch, tmp_path):
    fake_executable = tmp_path / "COMPAS"
    _make_fake_executable(fake_executable)

    monkeypatch.setenv("COMPAS_EXECUTABLE_PATH", str(fake_executable))
    completed_process = run_compas(["-v"], capture_output=True, text=True)

    assert completed_process.returncode == 0
    assert "fake-compas -v" in completed_process.stdout


def test_main_print_path(monkeypatch, tmp_path, capsys):
    fake_executable = tmp_path / "COMPAS"
    _make_fake_executable(fake_executable)

    monkeypatch.setenv("COMPAS_EXECUTABLE_PATH", str(fake_executable))

    assert main(["--print-path"]) == 0
    assert str(fake_executable) in capsys.readouterr().out
