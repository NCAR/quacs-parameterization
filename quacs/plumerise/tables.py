"""Default lookup tables used by the workflow."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class VegetationSplit:
    smoldering_fraction: float
    flaming_fraction: float
    prm_veg_group: str
    source_mapping: str
    notes: str = ""


@dataclass(frozen=True)
class HeatFlux:
    heat_flux_kw_m2: float
    source_mapping: str
    notes: str = ""


DEFAULT_VEGETATION_SPLIT_TABLE: dict[str, VegetationSplit] = {
    "forest": VegetationSplit(0.30, 0.70, "forest", "quacs_placeholder"),
    "savanna": VegetationSplit(0.45, 0.55, "savanna", "quacs_placeholder"),
    "grassland": VegetationSplit(0.55, 0.45, "grassland", "quacs_placeholder"),
    "shrubland": VegetationSplit(0.50, 0.50, "shrubland", "quacs_placeholder"),
    "default": VegetationSplit(0.40, 0.60, "default", "fallback"),
}

DEFAULT_HEAT_FLUX_TABLE: dict[str, HeatFlux] = {
    "forest": HeatFlux(80.0, "quacs_placeholder"),
    "savanna": HeatFlux(55.0, "quacs_placeholder"),
    "grassland": HeatFlux(35.0, "quacs_placeholder"),
    "shrubland": HeatFlux(45.0, "quacs_placeholder"),
    "default": HeatFlux(50.0, "fallback"),
}


def normalize_vegetation_class(vegetation_class):
    return vegetation_class.strip().lower()


def lookup_vegetation_split(
    vegetation_class,
    table=None,
):
    table = table or DEFAULT_VEGETATION_SPLIT_TABLE
    key = normalize_vegetation_class(vegetation_class)
    if key in table:
        return table[key], key
    return table["default"], "default"


def lookup_heat_flux(
    vegetation_class,
    table=None,
):
    table = table or DEFAULT_HEAT_FLUX_TABLE
    key = normalize_vegetation_class(vegetation_class)
    if key in table:
        return table[key], key
    return table["default"], "default"


def validate_vegetation_split_table(
    table=None,
    tolerance=1.0e-9,
):
    table = table or DEFAULT_VEGETATION_SPLIT_TABLE
    for vegetation_class, split in table.items():
        total = split.smoldering_fraction + split.flaming_fraction
        if abs(total - 1.0) > tolerance:
            raise ValueError(
                f"vegetation split for {vegetation_class!r} sums to {total}, not 1"
            )
