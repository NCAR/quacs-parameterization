## Description

<!-- Describe the parameterization: what physical process does it represent,
     what model/scheme is implemented, and what it produces. -->

## Scientific references

<!-- List the key papers this parameterization is based on. -->

## Checklist

- [ ] Parameterization lives in its own folder under `quacs/` and is importable as `from quacs.<name> import ...`
- [ ] `quacs/<name>/__init__.py` exports the public API
- [ ] `quacs/<name>/README.md` describes the science, inputs/outputs, and references
- [ ] Tests in `quacs/<name>/tests/` — `pytest` passes from the repo root
- [ ] New dependencies added to the top-level `pyproject.toml`
- [ ] `ruff check .` passes
- [ ] `ruff format --check .` passes
