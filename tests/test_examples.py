"""Run each example script in a scheme folder and check that it exits with no error.

To exclude a script, put this line in the script, with a reason:

    # quacs: skip-example <reason>
"""

import os
import subprocess
import sys

import pytest
from repo_layout import REPO_ROOT, example_scripts

SKIP_MARKER = "# quacs: skip-example"
TIMEOUT_SECONDS = 600


def _rel(path):
    return str(path.relative_to(REPO_ROOT))


@pytest.mark.parametrize("script", example_scripts(), ids=_rel)
def test_example_runs(script, tmp_path):
    for line in script.read_text().splitlines():
        if line.strip().startswith(SKIP_MARKER):
            pytest.skip(line.strip()[len(SKIP_MARKER) :].strip() or "skip marker")

    env = dict(os.environ, MPLBACKEND="Agg")
    # Run in a temporary folder so that output files do not go into the repository.
    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
    )
    assert result.returncode == 0, (
        f"{_rel(script)} failed with exit code {result.returncode}\n"
        f"--- stdout ---\n{result.stdout[-4000:]}\n--- stderr ---\n{result.stderr[-4000:]}"
    )
