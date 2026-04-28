# drydep/simple

Offline GEOS-Chem dry deposition scheme integrated with a MUSICA MICM box model for Hg0. Computes a dry deposition velocity from meteorology and Olson land-cover parameters, converts it to a first-order loss rate, and time-integrates Hg0 concentration forward with the Rosenbrock solver.

This uses the [GEOS-Chem standalone drypdep scheme](https://github.com/arifein/offline-drydep).

## Usage

**Simple 0-D calculation** (no MUSICA required):
```
python examples/ex_drydep_simple.py
```
Prints the dry deposition velocity for a single atmospheric column with default tropical mixed-forest conditions.

**MUSICA box model** (requires MUSICA installed):
```
python examples/musica_drydep_box_model.py
```
Integrates Hg0 concentration over ~5 e-folding lifetimes and prints a time series of concentration decay.
