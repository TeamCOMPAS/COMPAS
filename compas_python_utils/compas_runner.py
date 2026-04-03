import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Optional, Sequence


PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent


def _is_runnable_file(path: Path) -> bool:
    return path.is_file() and (os.access(path, os.X_OK) or path.suffix == ".sh")


def _validate_explicit_path(path: Path, variable_name: str) -> str:
    if not _is_runnable_file(path):
        raise FileNotFoundError(
            f"{variable_name} points to a non-runnable path: {path}"
        )
    return str(path)


def _candidate_paths() -> Iterable[Path]:
    bundle_root = os.environ.get("COMPAS_BUNDLE_ROOT")
    if bundle_root:
        bundle_path = Path(bundle_root)
        yield bundle_path / "run_compas.sh"
        yield bundle_path / "bin" / "COMPAS"

    yield PACKAGE_ROOT / "bundled" / "COMPAS-linux-x86_64" / "run_compas.sh"
    yield PACKAGE_ROOT / "bundled" / "COMPAS-linux-x86_64" / "bin" / "COMPAS"

    compas_root = Path(os.environ.get("COMPAS_ROOT_DIR", REPO_ROOT))
    yield compas_root / "src" / "COMPAS"
    yield compas_root / "bin" / "COMPAS"

    yield REPO_ROOT / "src" / "COMPAS"
    yield REPO_ROOT / "bin" / "COMPAS"


def resolve_compas_executable() -> str:
    explicit_path = os.environ.get("COMPAS_EXECUTABLE_PATH")
    if explicit_path:
        return _validate_explicit_path(Path(explicit_path), "COMPAS_EXECUTABLE_PATH")

    bundle_root = os.environ.get("COMPAS_BUNDLE_ROOT")
    if bundle_root:
        candidates = [Path(bundle_root) / "run_compas.sh", Path(bundle_root) / "bin" / "COMPAS"]
        for candidate in candidates:
            if _is_runnable_file(candidate):
                return str(candidate)
        raise FileNotFoundError(
            "COMPAS_BUNDLE_ROOT is set, but neither run_compas.sh nor bin/COMPAS "
            f"exists under {bundle_root}"
        )

    for candidate in _candidate_paths():
        if _is_runnable_file(candidate):
            return str(candidate)

    raise FileNotFoundError(
        "Unable to locate a COMPAS executable. Set COMPAS_EXECUTABLE_PATH to an "
        "existing executable, or set COMPAS_BUNDLE_ROOT to an extracted bundle directory."
    )


def run_compas(
    compas_args: Optional[Sequence[str]] = None,
    executable: Optional[str] = None,
    check: bool = True,
    **kwargs,
) -> subprocess.CompletedProcess:
    resolved_executable = executable or resolve_compas_executable()
    executable_path = Path(resolved_executable)
    if executable_path.suffix == ".sh":
        command = ["bash", resolved_executable, *(compas_args or [])]
    else:
        command = [resolved_executable, *(compas_args or [])]
    return subprocess.run(command, check=check, **kwargs)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the COMPAS executable from Python.",
    )
    parser.add_argument(
        "--print-path",
        action="store_true",
        help="Print the resolved COMPAS executable path and exit.",
    )
    args, compas_args = parser.parse_known_args(argv)

    executable = resolve_compas_executable()

    if args.print_path:
        print(executable)
        return 0

    completed_process = run_compas(compas_args=compas_args, executable=executable, check=False)
    return completed_process.returncode


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
