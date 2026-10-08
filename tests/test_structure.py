"""Check that each process and scheme folder follows the layout in CONTRIBUTING.md."""

import re

import pytest
from repo_layout import REPO_ROOT, process_dirs, scheme_dirs

REQUIRED_SECTIONS = [
    "## Target parameterization template",
    "## Classification",
    "## Inputs",
    "## Outputs",
    "## References",
]
STATUS_PATTERN = re.compile(
    r"\|\s*\*\*Status\*\*\s*\|\s*(implemented|in progress|placeholder)\s*\|"
)


def _rel(path):
    return str(path.relative_to(REPO_ROOT))


@pytest.mark.parametrize("process", process_dirs(), ids=_rel)
def test_process_has_readme_that_links_each_scheme(process):
    readme = process / "README.md"
    assert readme.is_file(), f"{_rel(process)} has no README.md"
    text = readme.read_text()
    for scheme in scheme_dirs():
        if scheme.parent == process:
            assert f"({scheme.name}/)" in text, (
                f"{_rel(readme)} does not link to the scheme folder {scheme.name}/"
            )


@pytest.mark.parametrize("scheme", scheme_dirs(), ids=_rel)
def test_scheme_readme_follows_template(scheme):
    readme = scheme / "README.md"
    assert readme.is_file(), f"{_rel(scheme)} has no README.md"
    text = readme.read_text()
    match = STATUS_PATTERN.search(text)
    assert match, f"{_rel(readme)} has no valid '| **Status** | ... |' row"
    missing = [s for s in REQUIRED_SECTIONS if s not in text]
    assert not missing, f"{_rel(readme)} does not have these sections: {missing}"


@pytest.mark.parametrize("scheme", scheme_dirs(), ids=_rel)
def test_scheme_with_code_has_tests_or_examples(scheme):
    """A scheme that is not a placeholder must have code that CI runs."""
    status = STATUS_PATTERN.search((scheme / "README.md").read_text())
    if status is None or status.group(1) == "placeholder":
        pytest.skip("placeholder scheme")
    has_tests = any((scheme / "tests").glob("test_*.py"))
    has_examples = any((scheme / "examples").glob("*.py")) or any(
        (scheme / "examples").glob("*.ipynb")
    )
    assert has_tests or has_examples, f"{_rel(scheme)} has no tests and no examples"
