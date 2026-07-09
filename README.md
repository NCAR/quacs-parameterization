# quacs-parameterization

A collection of scientific parameterizations for QUACS. Each parameterization
lives in its own sub-package under `quacs/` and can be imported, tested, and
extended independently.

## Structure

```
quacs/
└── <parameterization>/
    ├── __init__.py       # public API
    ├── README.md         # science description, inputs/outputs, references
    ├── data/             # any data files required
    ├── examples/         # runnable example scripts
    └── tests/            # pytest tests
```

## Installation

```bash
pip install -e ".[dev]"
```

## Running tests

```bash
pytest
```

## Style

Code is formatted and linted with [ruff](https://docs.astral.sh/ruff/):

```bash
ruff check .
ruff format .
```

## Contributing

See [`.github/PULL_REQUEST_TEMPLATE.md`](.github/PULL_REQUEST_TEMPLATE.md) for
the contribution checklist. Every parameterization must include a `README.md`
and tests that pass under `pytest` from the repo root.

## Parameterizations

| Name | Description | Reference |
|------|-------------|-----------|
| [`drydep`](quacs/drydep/README.md) | GEOS-Chem offline dry deposition scheme for Hg0 | Wesely (1989); Feinberg et al. (2022) |
