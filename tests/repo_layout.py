"""Find the process folders, scheme folders, and example scripts in the repository."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / "quacs"

# Folder names that hold support files, not schemes.
NON_SCHEME_DIRS = {"data", "docs", "diagrams", "examples", "inputs", "output", "tests"}


def _is_content_dir(path: Path) -> bool:
    return path.is_dir() and not path.name.startswith((".", "_"))


def process_dirs():
    """Return each process folder, for example ``quacs/dust``."""
    return sorted(p for p in PACKAGE_ROOT.iterdir() if _is_content_dir(p))


def scheme_dirs():
    """Return each scheme folder, for example ``quacs/dust/ginoux``."""
    return sorted(
        s
        for p in process_dirs()
        for s in p.iterdir()
        if _is_content_dir(s) and s.name not in NON_SCHEME_DIRS
    )


def example_scripts():
    """Return each ``examples/*.py`` script in a scheme folder."""
    return sorted(
        f for s in scheme_dirs() for f in (s / "examples").glob("*.py") if f.name != "__init__.py"
    )
