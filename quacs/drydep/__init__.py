"""
GEOS-Chem offline dry deposition scheme for Hg0 (and other trace gases).

Two sub-packages are available:

``quacs.drydep.simple``
    Original FORTRAN-translated implementation.  Single-species, hardcoded
    tropical-forest land cover.  Minimal dependencies (numpy, xarray).

    >>> from quacs.drydep.simple import compute_drydep_rate
    >>> dvel, k = compute_drydep_rate()

``quacs.drydep.lhs``
    Pythonic rewrite with multi-species support, ``LandCoverPatch`` objects
    for flexible land cover, vectorized multi-cell computation, and a Latin
    Hypercube Sampling sensitivity driver.  Requires musica.

    >>> from quacs.drydep.lhs import compute_drydep_rate, olson_land_cover
    >>> from quacs.drydep.lhs.musica_box_model import hg0, so2, o3
    >>> land_cover = [olson_land_cover.tropical_forest()]
    >>> result = compute_drydep_rate([hg0], met={...}, land_cover=land_cover,
    ...                              coefficients=olson_land_cover.coefficients)
"""

from quacs.drydep import lhs, simple
