
import re
from pathlib import Path

try:
    from importlib.metadata import PackageNotFoundError, version as package_version
except ImportError:  # pragma: no cover
    from importlib_metadata import PackageNotFoundError, version as package_version


def _version_from_changelog() -> str:
    changelog_path = Path(__file__).resolve().parent.parent / "src" / "changelog.h"
    version_match = re.search(
        r'VERSION_STRING = ["\']([^"\']+)["\']',
        changelog_path.read_text(encoding="utf-8"),
    )
    if not version_match:
        return "0.0.0"
    return ".".join(str(int(part)) for part in version_match.group(1).split("."))


__all__ = []

__author__ = "Team COMPAS"
__email__ = "teamcompas@users.noreply.github.com"
__uri__ = "https://github.com/TeamCOMPAS/COMPAS"
__license__ = "MIT"
__description__ = "COMPAS"
__copyright__ = "Copyright 2022 COMPAS developers"
__contributors__ = "https://github.com/TeamCOMPAS/COMPAS/graphs/contributors"

try:
    __version__ = package_version("compas-popsynth")
except PackageNotFoundError:
    __version__ = _version_from_changelog()
