"""CMAQ LS99 FENGSHA dust emissions on two-dimensional cells.

Soil tables contain all 16 CMAQ rows. Soil type remapping occurs only when
``soil_type_is_premapped`` is false (the CMAQ ``.not. WRF_V4P`` branch).
Land-use columns are percentages of each grid cell in the selected scheme.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field, fields

import numpy as np

# ============================================================================
# Constants (from CMAQ DUST_EMIS.F and SUBST_CONST)
# ============================================================================
GRAV = 9.80622  # gravitational acceleration [m/s2]
KARMAN = 0.4  # von Karman constant
RHOP = 2650.0  # soil particle density [kg/m3]
PI = 3.14159265  # CONST.EXT:48, single-precision CMAQ literal
GPKG = 1.0e3  # AERO_DATA.F:737
DUST_DENS_CMAQ = np.array([2200.0, 2156.55, 2536.92])  # AERO_DATA.F:539
TINY_REAL4 = float(np.float32(1.1754944e-38))  # TINY(0.0) as a 4-byte real

# Wind-blown dust species table dust_spc, AERO_DATA.F:505-532.
# Each row: (CMAQ name, molecular weight [g/mol], fine fraction, coarse fraction [g/g]).
# DUST_EMIS.F:1042-1045 splits each mode's mass with the two fractions. The molecular
# weight feeds only the DESID setup (DUST_EMIS.F:238-239), which this port leaves out.
DUST_SPC = (
    ("SO4", 96.0, 0.02250, 0.02655),
    ("NO3", 62.0, 0.00020, 0.00160),
    ("CL", 35.5, 0.00945, 0.01190),
    ("NH4", 18.0, 0.00005, 0.0),
    ("NA", 23.0, 0.03935, 0.0),
    ("CA", 40.1, 0.07940, 0.0),
    ("MG", 24.3, 0.01900, 0.0),
    ("K", 39.1, 0.03770, 0.0),
    ("POC", 220.0, 0.01075, 0.0),
    ("PNCOM", 220.0, 0.00430, 0.0),
    ("LVPO1", 218.0, 0.01501, 0.0),
    ("LVOO1", 136.0, 2.23e-5, 0.0),
    ("EC", 12.0, 0.0, 0.0),
    ("FE", 55.8, 0.03355, 0.0),
    ("AL", 27.0, 0.05695, 0.0),
    ("SI", 28.1, 0.19425, 0.0),
    ("TI", 47.9, 0.00280, 0.0),
    ("MN", 54.9, 0.00115, 0.0),
    ("H2O", 18.0, 0.00541, 0.00637),
    ("OTHR", 200.0, 0.48319, 0.0),
    ("SOIL", 100.0, 0.0, 0.95358),
    ("MN_HAPS", 54.9, 0.00041, 0.00041),
    ("NI", 58.7, 2.0e-05, 2.0e-05),
    ("CR_III", 52.0, 5.6e-05, 5.6e-05),
    ("AS", 74.92, 5.5e-06, 5.5e-06),
    ("PB", 207.2, 1.6e-05, 1.6e-05),
    ("CD", 112.4, 9.7e-08, 9.7e-08),
    ("PHG", 200.5, 1.4e-08, 1.4e-08),
)
DUST_SPC_NAMES = tuple(row[0] for row in DUST_SPC)
DUST_SPC_MW = tuple(row[1] for row in DUST_SPC)
SPCFAC = tuple((row[2], row[3]) for row in DUST_SPC)

# Shao & Lu (2000) threshold literals from the pinned CMAQ source
A_THRESH = 260.60061
B_THRESH = 1.6540342e-6

# Marticorena & Bergametti (1997) constant for saltation flux
AMEN = 1.0
CFAC = 1000.0 * AMEN / GRAV

# Erodibility potential weights, in (clay, silt, sand) order
# DUST_EMIS.F: eropot(3) = (/ 0.08 clay, 1.00 silt, 0.12 sand /)
EROPOT = np.array([0.08, 1.00, 0.12])

# Dust mass fractions: Kok (2011), Nabat et al. (2012)
FRACMJ = 0.07  # accumulation mode
FRACMK = 0.93  # coarse mode

# Drag partitioning constants: Xi & Sokolik (2015)
MV, SIGV, BETAV = 0.16, 1.45, 202.0
MB, SIGB, BETAB = 0.5, 1.0, 90.0

# Mean mass median particle diameters [m]
# Chatenet et al. (1996), Menut et al. (2013)
DP = np.array(
    [
        690.0e-6,  # Coarse sand
        210.0e-6,  # Fine-medium sand
        125.0e-6,  # Silt
        2.0e-6,
    ]
)  # Clay

# ============================================================================
# Soil texture lookup table (from CMAQ DUST_EMIS.F, lines ~698-716)
# ============================================================================
# 16-category CMAQ soil classification, [coarse sand,
# fine-medium sand, silt, clay], verified line-by-line against the
# `soiltxt` array literal in DUST_EMIS.F. Rows 4-11 were previously wrong
# in this script (replaced with fabricated interpolated values); corrected
# here to match source exactly.
NSLTYP = 16

SOILTXT = np.array(
    [
        [0.46, 0.46, 0.05, 0.03],  #  1 Sand
        [0.41, 0.41, 0.18, 0.00],  #  2 Loamy Sand
        [0.29, 0.29, 0.32, 0.10],  #  3 Sandy Loam
        [0.00, 0.17, 0.70, 0.13],  #  4 Silt Loam
        [0.00, 0.10, 0.85, 0.05],  #  5 Silt
        [0.00, 0.43, 0.39, 0.18],  #  6 Loam
        [0.29, 0.29, 0.15, 0.27],  #  7 Sandy Clay Loam
        [0.00, 0.10, 0.56, 0.34],  #  8 Silty Clay Loam
        [0.00, 0.32, 0.34, 0.34],  #  9 Clay Loam
        [0.00, 0.52, 0.06, 0.42],  # 10 Sandy Clay
        [0.00, 0.06, 0.47, 0.47],  # 11 Silty Clay
        [0.00, 0.22, 0.20, 0.58],  # 12 Clay
        [0.00, 0.00, 0.00, 0.00],  # 13 Organic Material
        [0.00, 0.00, 0.00, 0.00],  # 14 Water
        [0.00, 0.00, 0.00, 0.00],  # 15 Bedrock
        [0.00, 0.00, 0.00, 0.00],  # 16 Other
    ]
)

# Soil moisture upper limits [kg/kg gravimetric] per soil type
# DUST_EMIS.F soilml1(nsltyp), all 16 rows
SOILML1 = np.array(
    [
        0.242,
        0.257,
        0.286,
        0.350,
        0.350,
        0.307,
        0.277,
        0.350,
        0.332,
        0.284,
        0.357,
        0.344,
        0.329,
        0.000,
        0.170,
        0.280,
    ]
)


def _remap_soil_type(j):
    """
    Apply the pre-WRF_V4P CMAQ soil remap, per DUST_EMIS.F 829-833:

        if j > 4:  j = j + 1
        if j > 13: j = 13

    Parameters
    ----------
    j : int or ndarray of int   1-based soil type index (raw WRF/STATSGO)

    Returns
    -------
    int or ndarray of int       1-based CMAQ table row.
    """
    j = np.where(j > 4, j + 1, j)
    j = np.where(j > 13, 13, j)
    return j


# ============================================================================
# Sandblasting efficiency lookup (from CMAQ DUST_EMIS.F, lines 900-943)
# ============================================================================
def _get_flxfac(soil_type_idx, settings=None):
    """Return (flxfac1, flxfac2) for a 1-based CMAQ soil table row."""
    settings = Ls99Settings() if settings is None else settings
    j = soil_type_idx  # 1-based CMAQ table row
    if j in (1, 2):  # sand
        return settings.flxfac1[0], settings.flxfac2[0]
    elif j in (3, 4, 6, 8, 9):  # loam
        return settings.flxfac1[1], settings.flxfac2[1]
    elif j == 7:  # sandy clay loam
        return settings.flxfac1[2], settings.flxfac2[2]
    elif j in (5, 10, 11, 12):  # clay
        return settings.flxfac1[3], settings.flxfac2[3]
    else:  # others (13) - no dust
        return settings.flxfac1[4], settings.flxfac2[4]


# ============================================================================
# Land-use scheme registry (from CMAQ lus_data_module.F)
# ============================================================================
# Each scheme maps n_dlcat desert columns to BELD3 classes. CMAQ appends
# one cropland slot after those columns and sets emap(n_dlcat+1)=4. That
# slot has zero coverage in this routine but remains in qam and diagnostics.
# `categories` and `dmap` describe only the supplied desert columns.
BELD3_AG_PLACEHOLDER = 4

LANDUSE_SCHEMES = {
    "usgs24": {
        "n_dlcat": 4,
        "categories": [
            "Shrubland",
            "Barren or Sparsely Vegetated",
            "Mixed Tundra",
            "Bare Ground Tundra",
        ],
        "dmap": [1, 3, 3, 3],
    },
    "modis": {
        "n_dlcat": 1,
        "categories": ["(no named desert category - placeholder scheme)"],
        "dmap": [None],
    },
    "modis_noah": {
        "n_dlcat": 4,
        "categories": [
            "Closed Shrublands",
            "Open Shrublands",
            "Barren or Sparsely Vegetated",
            "Barren Tundra",
        ],
        "dmap": [1, 1, 3, 3],
    },
    "nlcd_modis": {
        "n_dlcat": 7,
        "categories": [
            "Closed Shrublands",
            "Open Shrublands",
            "Barren Land (Rock-Sand-Clay)",
            "Dwarf Scrub",
            "Shrub-Scrub",
            "Tundra",
            "Barren or Sparsely Vegetated",
        ],
        "dmap": [1, 1, 3, 3, 3, 3, 3],
    },
    "nlcd40": {
        "n_dlcat": 6,
        "categories": [
            "Closed Shrublands",
            "Open Shrublands",
            "Barren Land (Rock-Sand-Clay)",
            "Dwarf Scrub",
            "Shrub-Scrub",
            "Barren or Sparsely Vegetated",
        ],
        "dmap": [1, 1, 3, 3, 3, 3],
    },
}

# 1-based CMAQ LUFRAC indices: lus_data_module.F desert maps and LUS_DEFN.F:100-300.
HOST_LANDUSE = {
    "usgs24": (
        (8, 19, 22, 23),
        (
            (8, 16, 19, 22, 23, 24),
            (2, 3, 4, 5, 6, 7, 9, 10, 17, 20),
            (1,),
            (11, 12, 13, 14, 15, 18, 21),
        ),
        24,
    ),
    "modis": ((0,), ((), (), (), ()), 1),
    "modis_noah": (
        (6, 7, 16, 20),
        ((6, 7, 15, 16, 17, 19, 20), (8, 9, 10, 12, 14, 18), (13,), (1, 2, 3, 4, 5, 11)),
        20,
    ),
    "nlcd_modis": (
        (37, 38, 7, 12, 13, 18, 47),
        (
            (1, 2, 7, 12, 13, 18, 31, 46, 47, 48),
            (14, 15, 16, 17, 19, 20, 37, 38, 40, 41, 43, 45),
            (4, 5, 6, 44),
            (9, 10, 11, 21, 22, 24, 32, 33, 34, 35, 36, 39, 42),
        ),
        50,
    ),
    "nlcd40": (
        (6, 7, 27, 31, 32, 16),
        (
            (21, 22, 27, 31, 32, 15, 16, 17),
            (33, 34, 35, 36, 37, 38, 6, 7, 9, 10, 12, 14),
            (24, 25, 26, 13),
            (28, 29, 30, 39, 40, 1, 2, 3, 4, 5, 8, 11),
        ),
        40,
    ),
}


def _host_landuse(luf, scheme, shape):
    """Pick CMAQ desert fractions and merged diagnostic cover from full LUFRAC."""
    if scheme not in HOST_LANDUSE:
        raise ValueError(f"Unknown land-use scheme: {scheme}")
    indices, groups, ncat = HOST_LANDUSE[scheme]
    luf = np.asarray(luf)
    if luf.shape != (*shape, ncat):
        raise ValueError(f"luf must have {ncat} land-use columns")
    # LUS_DEFN.F:140-147 comments out MODIS merging; its allocated uland is undefined, so use zero.
    desert = np.zeros((*shape, 1)) if scheme == "modis" else luf[..., np.array(indices) - 1] * 100.0
    cover = np.zeros((*shape, 4))
    # LUS_DEFN.F:99,104-280 converts each LUFRAC before adding in source order.
    for group_number, columns in enumerate(groups):
        for column in columns:
            cover[..., group_number] += luf[..., column - 1] * 100.0
    cover = np.clip(cover, 0.0, 100.0)
    total = np.zeros(shape)
    for group_number in range(4):
        total += cover[..., group_number]
    scale = np.ones_like(total)
    np.divide(100.0, total, out=scale, where=total > 100.0)
    cover = cover * scale.astype(np.float32)[..., None]
    return desert, cover


# Physical constants indexed by BELD3 class (1-4). Class 2 has no source
# definition in any scheme's dmap and is retained only as an array slot
# for shape compatibility with the Fortran 4-element LAMBDAB/HB_LAMBDAB/HV
# arrays; no desert land-use column maps to class 2.
# Index 0 unused (BELD3 classes are 1-based); indices 1-4 below map to
# BELD3 classes 1-4 respectively.
# From Darmenova et al. (2009), Xi & Sokolik (2015)
LAMBDAB = {1: 0.03, 2: 0.04, 3: 0.0001, 4: 0.15}
HB_LAMBDAB = {1: 6.0e-4, 2: 8.0e-4, 3: 2.0e-6, 4: 3.0e-3}

# Seasonal vegetation heights [m] by BELD3 class, indexed by month
HV_SEASONAL = {
    "nov_feb": {1: 0.05, 2: 0.05, 3: 0.05, 4: 0.05},
    "mar": {1: 0.15, 2: 0.05, 3: 0.10, 4: 0.05},
    "apr": {1: 0.15, 2: 0.10, 3: 0.10, 4: 0.05},
    "may": {1: 0.12, 2: 0.20, 3: 0.10, 4: 0.10},
    "jun": {1: 0.12, 2: 0.15, 3: 0.10, 4: 0.30},
    "jul": {1: 0.10, 2: 0.12, 3: 0.10, 4: 0.50},
    "aug": {1: 0.10, 2: 0.12, 3: 0.10, 4: 0.50},
    "sep": {1: 0.10, 2: 0.10, 3: 0.10, 4: 0.30},
    "oct": {1: 0.05, 2: 0.08, 3: 0.10, 4: 0.10},
}

FLXFAC1 = (5.886e-05, 5.2974e-05, 9.4176e-05, 2.3544e-05, 0.0)
FLXFAC2 = (1.5215430, 1.0758933, 1.0758933, 0.1964303, 0.3402273)
KVH_OFFSET = 0.24
WMAX_CLAY, WMAX_BASE = 14.0, 17.0
SOIMT_BASE, SOIMT_SAND = 0.511, 0.126
MOIST_CUTOFF, MOIST_SCALE, MOIST_POWER = 0.01, 1.21, 0.68
LAMBDA_V = 0.35
VEG_MAX, VEG_MIN = 0.95, 0.005
Z0_LOW_SCALE, Z0_LOW_POWER, Z0_SPLIT, Z0_HIGH_SCALE, Z0_HIGH_POWER = 0.96, 1.07, 0.2, 0.083, -0.46
FRUF_CAP = 10.0
WIND_HEIGHT = 10.0
RAIN_MAX, SNOW_MAX = 0.0254, 0.001
DGV = (1.3914, 5.2590)


@dataclass(frozen=True)
class Ls99Settings:
    rhop: float = field(default=RHOP, metadata={"source": "DE.F:879", "what": "particle density"})
    a_thresh: float = field(
        default=A_THRESH, metadata={"source": "DE.F:1243", "what": "dry threshold size slope"}
    )
    b_thresh: float = field(
        default=B_THRESH, metadata={"source": "DE.F:1244", "what": "dry threshold cohesion"}
    )
    amen: float = field(default=AMEN, metadata={"source": "DE.F:1241", "what": "saltation scale"})
    cfac: float = field(
        default=CFAC, metadata={"source": "DE.F:1242", "what": "saltation coefficient"}
    )
    dp: tuple = field(
        default_factory=lambda: tuple(DP),
        metadata={"source": "DE.F:720-724", "what": "texture particle diameters"},
    )
    eropot: tuple = field(
        default_factory=lambda: tuple(EROPOT),
        metadata={"source": "DE.F:167", "what": "erodibility texture weights"},
    )
    soiltxt: tuple = field(
        default_factory=lambda: tuple(map(tuple, SOILTXT)),
        metadata={"source": "DE.F:698", "what": "soil texture table"},
    )
    soilml1: tuple = field(
        default_factory=lambda: tuple(SOILML1),
        metadata={"source": "DE.F:677", "what": "soil moisture limit table"},
    )
    flxfac1: tuple = field(
        default=FLXFAC1, metadata={"source": "DE.F:908-942", "what": "sandblasting intercepts"}
    )
    flxfac2: tuple = field(
        default=FLXFAC2, metadata={"source": "DE.F:908-942", "what": "sandblasting slopes"}
    )
    kvh_offset: float = field(
        default=KVH_OFFSET, metadata={"source": "DE.F:991", "what": "sandblasting offset"}
    )
    wmax_clay: float = field(
        default=WMAX_CLAY, metadata={"source": "DE.F:873", "what": "dry moisture clay slope"}
    )
    wmax_base: float = field(
        default=WMAX_BASE, metadata={"source": "DE.F:873", "what": "dry moisture base"}
    )
    soimt_base: float = field(
        default=SOIMT_BASE, metadata={"source": "DE.F:879", "what": "moisture conversion base"}
    )
    soimt_sand: float = field(
        default=SOIMT_SAND, metadata={"source": "DE.F:879", "what": "moisture conversion sand"}
    )
    moist_cutoff: float = field(
        default=MOIST_CUTOFF, metadata={"source": "DE.F:889", "what": "dry moisture cutoff"}
    )
    moist_scale: float = field(
        default=MOIST_SCALE, metadata={"source": "DE.F:892", "what": "wet threshold scale"}
    )
    moist_power: float = field(
        default=MOIST_POWER, metadata={"source": "DE.F:893", "what": "wet threshold exponent"}
    )
    mv: float = field(
        default=MV, metadata={"source": "DE.F:607-614", "what": "vegetation drag factor"}
    )
    sigv: float = field(
        default=SIGV, metadata={"source": "DE.F:607-614", "what": "vegetation shelter factor"}
    )
    betav: float = field(
        default=BETAV, metadata={"source": "DE.F:607-614", "what": "vegetation drag ratio"}
    )
    mb: float = field(default=MB, metadata={"source": "DE.F:607-614", "what": "bare drag factor"})
    sigb: float = field(
        default=SIGB, metadata={"source": "DE.F:607-614", "what": "bare shelter factor"}
    )
    betab: float = field(
        default=BETAB, metadata={"source": "DE.F:607-614", "what": "bare drag ratio"}
    )
    lambda_v: float = field(
        default=LAMBDA_V, metadata={"source": "DE.F:854", "what": "frontal area coefficient"}
    )
    veg_max: float = field(
        default=VEG_MAX, metadata={"source": "DE.F:852", "what": "vegetation upper limit"}
    )
    veg_min: float = field(
        default=VEG_MIN, metadata={"source": "DE.F:852", "what": "vegetation lower limit"}
    )
    lambdab: tuple = field(
        default_factory=lambda: tuple(LAMBDAB.values()),
        metadata={"source": "DE.F:663-675", "what": "bare roughness by class"},
    )
    hb_lambdab: tuple = field(
        default_factory=lambda: tuple(HB_LAMBDAB.values()),
        metadata={"source": "DE.F:663-675", "what": "bare height by class"},
    )
    hv_seasonal: tuple = field(
        default_factory=lambda: tuple(tuple(values.values()) for values in HV_SEASONAL.values()),
        metadata={"source": "DE.F:757-780", "what": "seasonal vegetation heights"},
    )
    z0_low_scale: float = field(
        default=Z0_LOW_SCALE, metadata={"source": "DE.F:966", "what": "low roughness scale"}
    )
    z0_low_power: float = field(
        default=Z0_LOW_POWER, metadata={"source": "DE.F:966", "what": "low roughness exponent"}
    )
    z0_split: float = field(
        default=Z0_SPLIT, metadata={"source": "DE.F:966", "what": "roughness branch"}
    )
    z0_high_scale: float = field(
        default=Z0_HIGH_SCALE, metadata={"source": "DE.F:968", "what": "high roughness scale"}
    )
    z0_high_power: float = field(
        default=Z0_HIGH_POWER, metadata={"source": "DE.F:968", "what": "high roughness exponent"}
    )
    fruf_cap: float = field(
        default=FRUF_CAP, metadata={"source": "DE.F:985", "what": "drag factor cap"}
    )
    wind_height: float = field(
        default=WIND_HEIGHT,
        metadata={
            "source": "DE.F:972",
            "what": "height of the wind input in m; default 10 m because the input is the 10 m wind",
        },
    )
    rain_max: float = field(
        default=RAIN_MAX, metadata={"source": "DE.F:865", "what": "rain cutoff"}
    )
    snow_max: float = field(
        default=SNOW_MAX, metadata={"source": "DE.F:866", "what": "snow cutoff"}
    )
    fracmj: float = field(default=FRACMJ, metadata={"source": "DE.F:84", "what": "fine dust split"})
    fracmk: float = field(
        default=FRACMK, metadata={"source": "DE.F:85", "what": "coarse dust split"}
    )
    dust_dens: tuple = field(
        default_factory=lambda: tuple(DUST_DENS_CMAQ[:2]),
        metadata={"source": "AERO_DATA.F:539", "what": "dust mode densities"},
    )
    spcfac: tuple | None = field(
        default=SPCFAC,
        metadata={
            "source": "AERO_DATA.F:505-532 (dust_spc table; rows named in DUST_SPC_NAMES)",
            "what": "species mode fractions",
        },
    )
    dgv: tuple = field(
        default=DGV, metadata={"source": "DE.F:90-91", "what": "mode geometric diameters"}
    )

    def __post_init__(self):
        for name in (
            "dp",
            "eropot",
            "soiltxt",
            "soilml1",
            "flxfac1",
            "flxfac2",
            "lambdab",
            "hb_lambdab",
            "hv_seasonal",
            "dust_dens",
            "spcfac",
            "dgv",
        ):
            value = getattr(self, name)
            if isinstance(value, np.ndarray):
                value = value.copy()
                value.setflags(write=False)
                object.__setattr__(self, name, value)


def tunable_names(settings_class):
    return tuple(item.name for item in fields(settings_class))


def _get_hv(jday, settings=None):
    """Return the seasonal-vegetation-height dict (by BELD3 class) for a
    given Julian day of year."""
    if 59 < jday <= 90:
        key, index = "mar", 1
    elif 90 < jday <= 120:
        key, index = "apr", 2
    elif 120 < jday <= 151:
        key, index = "may", 3
    elif 151 < jday <= 181:
        key, index = "jun", 4
    elif 181 < jday <= 212:
        key, index = "jul", 5
    elif 212 < jday <= 243:
        key, index = "aug", 6
    elif 243 < jday <= 273:
        key, index = "sep", 7
    elif 273 < jday <= 304:
        key, index = "oct", 8
    else:
        key, index = "nov_feb", 0
    return HV_SEASONAL[key] if settings is None else settings.hv_seasonal[index]


# ============================================================================
# Core parameterization
# ============================================================================
def dust_hflux(dp, soiltxt_cell, fmoit, fruf, ustr, sd_ep, dens, settings=None):
    """
    Horizontal saltation flux, summed over 4 particle sizes.
    Port of the CMAQ dust_hflux function (DUST_EMIS.F lines 1220-1278).
    A comparison with the compiled pinned code shows no difference.

    Parameters
    ----------
    dp       : array [4]    particle diameters [m]
    soiltxt_cell : array [4] soil texture fractions at this cell
    fmoit    : float        soil moisture correction factor on u*t
    fruf     : float        drag partitioning correction factor on u*t
    ustr     : float        friction velocity [m/s]
    sd_ep    : float        soil-dust erodibility potential
    dens     : float        air density [kg/m3]

    Returns
    -------
    hflux    : float        total horizontal dust flux [g/m/s]
    """
    settings = Ls99Settings() if settings is None else settings
    fac = settings.cfac * settings.amen * dens * sd_ep
    hflux = 0.0
    for n in range(4):
        # Shao & Lu (2000) threshold x moisture x roughness
        utstar = np.sqrt(settings.a_thresh * dp[n] + settings.b_thresh / dp[n]) * fmoit * fruf
        if ustr > utstar:
            # White (1979) saltation flux
            utem = (ustr + utstar) * (ustr**2 - utstar**2)
            hflux += fac * utem * soiltxt_cell[n]
    return hflux


class DustInputWarning(UserWarning):
    """Input is outside its physical range; CMAQ continues the calculation."""


def _surface_rate(
    u10m,
    soim_vol,
    dens,
    snocov,
    rn,
    rc,
    vegfrac,
    soil_type,
    luf,
    jday,
    scheme,
    soil_type_is_premapped,
    lwmask,
    layer_height,
    rjacm,
    clay,
    csand,
    fmsand,
    wmax_last,
    sd_ep_last,
    settings=None,
):
    """Return CMAQ's single flux calculation and the fields it diagnoses."""
    settings = Ls99Settings() if settings is None else settings
    if scheme not in LANDUSE_SCHEMES:
        raise ValueError(f"Unknown land-use scheme: {scheme}")
    shape = np.shape(u10m)
    n = LANDUSE_SCHEMES[scheme]["n_dlcat"]
    if np.shape(luf) != (*shape, n):
        raise ValueError(f"luf must have {n} land-use columns")
    for name, value, low, high in (
        ("vegfrac", vegfrac, 0, 1),
        ("soim_vol", soim_vol, 0, 1),
        ("snocov", snocov, 0, 1),
        ("rn", rn, 0, np.inf),
        ("rc", rc, 0, np.inf),
        ("u10m", u10m, 0, np.inf),
        ("dens", dens, 0, np.inf),
        ("layer_height", layer_height, 0, np.inf),
        ("rjacm", rjacm, 0, np.inf),
    ):
        if np.any((np.asarray(value) < low) | (np.asarray(value) > high)):
            warnings.warn(f"{name} outside physical range", DustInputWarning, stacklevel=2)
    qam = np.zeros((*shape, n + 1))
    elus = np.zeros_like(qam)
    ustr = np.zeros_like(qam)
    kvh = np.zeros_like(qam)
    fruf = np.zeros_like(qam)
    soimt = np.zeros(shape)
    wmax = np.array(np.broadcast_to(wmax_last, shape), dtype=float, copy=True)
    fmoit = np.zeros(shape)
    sd_ep = np.array(np.broadcast_to(sd_ep_last, shape), dtype=float, copy=True)
    vf_out = np.zeros(shape)
    hv = _get_hv(jday, settings)
    for c in range(shape[0]):
        for r in range(shape[1]):
            j = int(soil_type[c, r])
            if not soil_type_is_premapped:
                j = int(_remap_soil_type(j))
            if not 1 <= j <= NSLTYP:
                raise ValueError(f"soil type {j} outside CMAQ table")
            if clay is None:
                cs, fm, cl = np.asarray(settings.soiltxt)[j - 1, (0, 1, 3)]
            else:
                cs, fm, cl = csand[c, r], fmsand[c, r], clay[c, r]
            sandf = cs + fm
            siltf = 1.0 - cl - sandf
            texture = np.array([cs, fm, siltf, cl])
            if np.any((texture < 0) | (texture > 1)):
                warnings.warn("soil texture outside physical range", DustInputWarning, stacklevel=2)
            vf = max(min(vegfrac[c, r], settings.veg_max), settings.veg_min)
            vf_out[c, r] = vf
            vegfree = 1.0 - vf
            lambdav = -settings.lambda_v * np.log(vegfree)
            if not (
                lwmask[c, r] > 0
                and rn[c, r] + rc[c, r] <= settings.rain_max
                and snocov[c, r] < settings.snow_max
            ):
                continue
            wmax[c, r] = (settings.wmax_clay * cl + settings.wmax_base) * cl
            soimt[c, r] = (
                soim_vol[c, r]
                * 1000.0
                / (settings.rhop * (settings.soimt_base + settings.soimt_sand * sandf))
            )
            if soimt[c, r] > settings.soilml1[j - 1]:
                continue
            if soimt[c, r] <= settings.moist_cutoff * wmax[c, r]:
                fmoit[c, r] = 1.0
            else:
                fmoit[c, r] = np.sqrt(
                    1.0
                    + settings.moist_scale
                    * (100.0 * soimt[c, r] - wmax[c, r]) ** settings.moist_power
                )
            sd_ep[c, r] = (
                cl * settings.eropot[0] + siltf * settings.eropot[1] + sandf * settings.eropot[2]
            )
            flxfac1, flxfac2 = _get_flxfac(j, settings)
            elus[c, r, :n] = luf[c, r, :] * vegfree
            for m in range(n + 1):
                coverage = elus[c, r, m]
                if coverage < 0 or coverage > 100:
                    raise ValueError(f"elus outside 0 to 100 at {(c, r, m)}: {coverage}")
                if coverage > 0:
                    cls = LANDUSE_SCHEMES[scheme]["dmap"][m] if m < n else BELD3_AG_PLACEHOLDER
                    if cls is None:
                        raise ValueError("positive coverage has no CMAQ land-use class")
                    lam = settings.lambdab[cls - 1] + lambdav
                    height = (settings.hb_lambdab[cls - 1] + hv[cls - 1] * lambdav) / lam
                    z0 = (
                        settings.z0_low_scale * lam**settings.z0_low_power
                        if lam <= settings.z0_split
                        else settings.z0_high_scale * lam**settings.z0_high_power
                    ) * height
                    ustr[c, r, m] = KARMAN * u10m[c, r] / np.log(settings.wind_height / z0)
                    fruf2 = (
                        (1.0 - settings.sigv * settings.mv * lambdav)
                        * (1.0 + settings.betav * settings.mv * lambdav)
                        * (1.0 - settings.sigb * settings.mb * settings.lambdab[cls - 1] / vegfree)
                        * (1.0 + settings.betab * settings.mb * settings.lambdab[cls - 1] / vegfree)
                    )
                    fruf[c, r, m] = np.sqrt(fruf2) if fruf2 > 1.0 else settings.fruf_cap
                    kvh[c, r, m] = flxfac1 * (settings.kvh_offset + flxfac2 * ustr[c, r, m])
                    hflux = dust_hflux(
                        settings.dp,
                        texture,
                        fmoit[c, r],
                        fruf[c, r, m],
                        ustr[c, r, m],
                        sd_ep[c, r],
                        dens[c, r],
                        settings,
                    )
                    vflux = hflux * kvh[c, r, m]
                    qam[c, r, m] += vflux * (rjacm[c, r] / layer_height[c, r]) * (coverage * 0.01)
                if coverage == 0 and qam[c, r, m] != 0:
                    raise ValueError("elus is zero but qam is nonzero")
    dust_em = qam.sum(axis=-1)
    return dict(
        qam=qam,
        dust_em=dust_em,
        elus=elus,
        ustr=ustr,
        kvh=kvh,
        fruf=fruf,
        soimt=soimt,
        wmax=wmax,
        fmoit=fmoit,
        sd_ep=sd_ep,
        vegfrac=vf_out,
        surface_rate=dust_em * layer_height / rjacm,
    )


def fengsha_dust_driver(
    u10m,
    soim_vol,
    dens,
    snocov,
    rn,
    rc,
    vegfrac,
    soil_type,
    luf,
    jdate,
    lwmask,
    layer_height,
    jacobm,
    wmax_last=None,
    sd_ep_last=None,
    scheme="usgs24",
    soil_type_is_premapped=False,
    clay=None,
    csand=None,
    fmsand=None,
    diagnostic_totals=None,
    sync_seconds=1.0,
    mpas=False,
    settings=None,
):
    """Port CMAQ F:821-1147; return named rates and accumulated diagnostics.

    ``spcfac`` is [species, mode]; ``dust_dens`` is [mode] in kg/m3.
    ``jacobm`` is the host Jacobian field. The port forms ``rjacm`` as
    ASX_DATA_MOD.F:875-882 does: ``1 / jacobm``, or 1.0 when ``mpas`` is
    true. With ``mpas`` the field is never read, so ``jacobm`` may be None.
    Without ``mpas`` it is required.
    ``wmax_last`` and ``sd_ep_last`` are the prior CMAQ saved fields.
    None means zero fields, as at the first CMAQ step (DUST_EMIS.F:281-282).
    Pass returned fields back on the next call. ``diagnostic_totals`` is
    the prior sum in rate-seconds.
    """
    settings = Ls99Settings() if settings is None else settings
    spcfac, dust_dens = settings.spcfac, settings.dust_dens
    shape = np.shape(u10m)
    luf, cover = _host_landuse(luf, scheme, shape)
    dens = np.asarray(dens)[:, 0, :]
    if mpas:
        soim_vol = np.asarray(soim_vol)[:, 0, :]
    if np.any(np.asarray(snocov) < 0):
        warnings.warn("snocov outside physical range", DustInputWarning, stacklevel=2)
    snocov = np.maximum(snocov, 0.0)
    jday = int(jdate) % 1000
    layer_height = np.broadcast_to(layer_height, shape)
    lwmask = np.broadcast_to(lwmask, shape)
    if mpas:
        # ASX_DATA_MOD.F:875-876: RJACM = 1.0, JACOBM is not read.
        rjacm = np.ones(shape)
    else:
        # ASX_DATA_MOD.F:878-882: RJACM = 1.0 / JACOBM.
        if jacobm is None:
            raise ValueError("jacobm is required unless mpas is true")
        jacobm = np.broadcast_to(np.asarray(jacobm, dtype=float), shape)
        if np.any(jacobm <= 0):
            warnings.warn("jacobm outside physical range", DustInputWarning, stacklevel=2)
        with np.errstate(divide="ignore"):
            rjacm = 1.0 / jacobm
    if wmax_last is None:
        wmax_last = np.zeros(shape)
    if sd_ep_last is None:
        sd_ep_last = np.zeros(shape)
    if any(x is not None for x in (clay, csand, fmsand)) and not all(
        x is not None for x in (clay, csand, fmsand)
    ):
        raise ValueError("PX texture needs clay, csand, and fmsand")
    if spcfac is None or dust_dens is None:
        raise ValueError("spcfac and dust_dens are required")
    core = _surface_rate(
        u10m,
        soim_vol,
        dens,
        snocov,
        rn,
        rc,
        vegfrac,
        soil_type,
        luf,
        jday,
        scheme,
        soil_type_is_premapped,
        lwmask,
        layer_height,
        rjacm,
        clay,
        csand,
        fmsand,
        wmax_last,
        sd_ep_last,
        settings,
    )
    edust = core["dust_em"][..., None] * np.array([settings.fracmj, settings.fracmk])
    core["mode_mass"] = edust
    spcfac = np.asarray(spcfac)
    if spcfac.ndim != 2 or spcfac.shape[1] != 2:
        raise ValueError("spcfac must have shape [species, 2]")
    core["species_rates"] = edust[..., :, None] * spcfac.T
    dust_dens = np.asarray(dust_dens)
    if dust_dens.shape != (2,):
        raise ValueError("dust_dens must have two mode densities")
    f6dpi = 6.0 / PI
    m3 = edust * f6dpi / (GPKG * dust_dens) + TINY_REAL4
    l2sg = np.log(2.0) * np.log(2.0)
    dgv = np.asarray(settings.dgv)
    factnum = 1.0e18 * np.exp(4.5 * l2sg) / dgv**3
    factm2 = 1.0e6 * np.exp(0.5 * l2sg) / dgv
    factsrf = PI * factm2
    core["number_rates"] = m3 * factnum
    core["surface_area_rates"] = m3 * factsrf
    urban, forest = cover[..., 2], cover[..., 3]
    parts = [
        core["qam"],
        core["dust_em"][..., None],
        core["elus"],
        core["elus"].sum(axis=-1)[..., None],
        core["ustr"],
        core["kvh"],
        core["fruf"],
        *[
            np.asarray(x)[..., None]
            for x in (
                core["fmoit"],
                core["sd_ep"],
                core["wmax"],
                core["vegfrac"],
                np.broadcast_to(urban, shape),
                np.broadcast_to(forest, shape),
            )
        ],
        core["number_rates"],
        core["surface_area_rates"],
    ]
    for mode in range(2):
        for species in range(spcfac.shape[0]):
            if spcfac[species, mode] > 0:
                parts.append(core["species_rates"][..., mode, species, None])
    core["diagnostics"] = np.concatenate(parts, axis=-1)
    previous = 0 if diagnostic_totals is None else diagnostic_totals
    core["diagnostic_totals"] = previous + core["diagnostics"] * sync_seconds
    return core


def dust_emission(
    u10m,
    soim_vol,
    dens,
    snocov,
    rain=None,
    vegfrac=None,
    soil_type=None,
    luf=None,
    jday=None,
    scheme="usgs24",
    soil_type_is_premapped=False,
    *,
    rn=None,
    rc=None,
    lwmask=None,
    layer_height=1.0,
    rjacm=1.0,
    clay=None,
    csand=None,
    fmsand=None,
    settings=None,
):
    """Return two surface mass rates [g/m2/s] from the shared CMAQ calculation."""
    settings = Ls99Settings() if settings is None else settings
    shape = np.shape(u10m)
    if rain is not None:
        if rn is not None or rc is not None:
            raise ValueError("supply rain or rn and rc")
        rn, rc = rain, np.zeros(shape)
    if rn is None or rc is None:
        raise ValueError("both rn and rc are required")
    if lwmask is None:
        lwmask = np.ones(shape)
    if any(x is not None for x in (clay, csand, fmsand)) and not all(
        x is not None for x in (clay, csand, fmsand)
    ):
        raise ValueError("PX texture needs clay, csand, and fmsand")
    core = _surface_rate(
        u10m,
        soim_vol,
        dens,
        snocov,
        rn,
        rc,
        vegfrac,
        soil_type,
        luf,
        jday,
        scheme,
        soil_type_is_premapped,
        np.broadcast_to(lwmask, shape),
        np.broadcast_to(layer_height, shape),
        np.broadcast_to(rjacm, shape),
        clay,
        csand,
        fmsand,
        np.zeros(shape),
        np.zeros(shape),
        settings,
    )
    return core["surface_rate"][..., None] * np.array([settings.fracmj, settings.fracmk])


# ============================================================================
# Standalone validation: box-model demo
# ============================================================================
