## Description

<!-- Which target parameterization does this PR add or change?
     Tell what physical or chemical process it represents and what it produces. -->

## Scientific references

<!-- List the key papers and the original code that this parameterization follows. -->

## Checklist

- [ ] The scheme is in `quacs/<process>/<scheme>/` and imports as `from quacs.<process>.<scheme> import ...`
- [ ] `quacs/<process>/<scheme>/README.md` follows `templates/scheme/README.md`, with complete **Inputs** and **Outputs** tables
- [ ] The process `README.md` lists the scheme and its status
- [ ] At least one example script or notebook in `examples/` runs
- [ ] New dependencies are in `pyproject.toml`
- [ ] `pytest` passes from the repository root
- [ ] `ruff check .` and `ruff format --check .` pass
