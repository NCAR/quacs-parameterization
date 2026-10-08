# Contributing a parameterization

This repository collects the target parameterizations for the QUACS project.
Each parameterization must make its inputs and outputs clear, and CI must run it.

## Layout

Each parameterization is a *scheme* inside a *process* folder:

```
quacs/
└── <process>/                 # for example: dust, drydep, lightning
    ├── README.md              # table of the schemes in this process
    ├── __init__.py
    └── <scheme>/              # for example: ginoux, wesely
        ├── README.md          # copy of templates/scheme/README.md, filled in
        ├── __init__.py        # public API: from quacs.<process>.<scheme> import ...
        ├── <scheme>.py        # the physics, as functions with no side effects
        ├── data/              # small data files that the scheme needs (optional)
        ├── examples/
        │   ├── musica_box_model.py   # runs the scheme with MUSICA
        │   └── <scheme>.ipynb        # plots and input configuration (optional)
        └── tests/
            ├── __init__.py
            └── test_<scheme>.py      # unit tests (optional, but welcome)
```

Use lowercase names with underscores for the folders. Name the scheme folder after the
scheme or its first author, for example `ginoux` or `lnox_decaria`.

## Steps to add a scheme

1. Find the process folder for your scheme. If none exists, make one with a `README.md` and an `__init__.py`.
2. If a placeholder folder exists for your scheme, use it. If not, make a new scheme folder.
3. Copy `templates/scheme/README.md` to the scheme folder and fill in every section.
4. Put the physics in `.py` modules. Put plots and input configuration in notebooks.
5. Add at least one example script or notebook in `examples/`.
6. Add the scheme to the table in the process `README.md`. Set the status.
7. Add new dependencies to `pyproject.toml`.
8. Run the checks below.

## The scheme README

The README is the contract of the scheme. It must have these parts:

- A header table with **Status**, **Target class**, **Contributors**, and **Original source**.
  The status is `implemented`, `in progress`, or `placeholder`.
- The target parameterization template, with the same fields as the target list.
- The classification, with the 5 dimensions from the proposal: complexity, solving strategy,
  aerosol representation, grid, and interdependence.
- An **Inputs** table. Each row gives the name, description, units, shape, and source.
  The source is MPAS-A, external dataset, constant, or other parameterization.
- An **Outputs** table. Each row gives the name, description, units, shape, and destination.
  The destination is MPAS-A state, MICM rate parameter, or diagnostic.
- **References**.

Use the same names in the Inputs and Outputs tables as the arguments and return values in the code.

## Code and notebooks

- Put all physics in `.py` modules so that other code can import it and tests can call it.
- Use notebooks only to configure inputs, call the module, and plot results.
- Import the scheme with its full path, for example `from quacs.dust.ginoux import ...`.
  Do not change `sys.path`.

## What CI checks

- `ruff check .` and `ruff format --check .` pass. This includes notebooks.
- `tests/test_structure.py`: each process and scheme folder has a README in the correct form.
- `tests/test_examples.py`: each `examples/*.py` script runs with no error.
  To skip a script, add the line `# quacs: skip-example <reason>` to the script.
- `pytest --nbmake`: each notebook runs with no error.
- Each `tests/test_*.py` file in a scheme folder passes.

## Run the checks

```bash
pip install -e ".[dev]"
ruff check .
ruff format --check .
pytest
```
