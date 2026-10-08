# Dry deposition

Schemes that compute the dry deposition velocity of gases and the loss that it causes.

| Scheme | Description | Status |
|---|---|---|
| [`wesely`](wesely/) | Wesely (1989), from CAM mo_drydep.F90 | implemented |
| [`simple`](simple/) | GEOS-Chem offline dry deposition, Fortran translation | implemented |
| [`lhs`](lhs/) | GEOS-Chem offline dry deposition with Latin hypercube sampling | implemented |

The `data/` folder holds the Olson 2001 land-cover files that `simple` and `lhs` use.

Each scheme folder follows the layout in [CONTRIBUTING.md](../../CONTRIBUTING.md).
