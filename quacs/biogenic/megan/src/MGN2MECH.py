"""MEGAN3 19 emission classes -> 201 speciated species -> CBMZ.

Input/output units: nmol m-2 s-1 (no area/time conversion). Supports
single cell (19,), gridded (19, ny, nx), or arbitrary trailing dimensions.
The supplied tables are transcribed from SPC_NOCONVER.EXT and MAP_CV2CBMZ.EXT.
No automatic normalization or de-duplication of official coefficients.

Important: CBMZ coefficients describe mechanism-equivalent mole amounts,
not necessarily conserved molecule or carbon totals. This module only
maps gases; the special MEGAN soil-NO logic belongs upstream.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

INPUTS = Path(__file__).resolve().parent.parent / "inputs"


def _load_csv(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh))


class MGN2MECH:
    """Reusable, vectorized emission-speciation mapper.

    Parameters
    ----------
    input_dir : str or Path, optional
        Directory containing class_to_species.csv, species_to_CBMZ.csv, and species_catalog.csv.
    mapping_files : mapping, optional
        Explicit paths keyed by class_to_species, species_to_CBMZ, and species_catalog;
        overrides input_dir when supplied.
    Duplicate source-species/target-mechanism mappings always raise ValueError,
    regardless of conversion_factor.
    """

    def validate_dimensions(self, *, n_class, n_species):
        if self.n_classes != n_class:
            raise ValueError(f"n_class={n_class}, mapping has {self.n_classes} classes")
        if self.n_species != n_species:
            raise ValueError(f"n_species={n_species}, catalog has {self.n_species} species")

    def __init__(self, input_dir=None, mapping_files=None):
        d = Path(input_dir) if input_dir is not None else INPUTS
        if mapping_files is None:
            mapping_files = {
                name: d / f"{name}.csv"
                for name in ("class_to_species", "species_to_CBMZ", "species_catalog")
            }
        c = _load_csv(mapping_files["class_to_species"])
        m = _load_csv(mapping_files["species_to_CBMZ"])
        catalog = _load_csv(mapping_files["species_catalog"])
        ids = [int(row["species_id"]) for row in catalog]
        names = [row["species_name"] for row in catalog]
        if len(set(ids)) != len(ids) or len(set(names)) != len(names):
            raise ValueError("Duplicate species_id or species_name in species_catalog.csv")
        if sorted(ids) != list(range(1, len(ids) + 1)):
            raise ValueError("species_id must be contiguous 1..n_species")
        n_species = len(catalog)
        n_classes = max(int(row["class_id"]) for row in c)
        if not n_species or not n_classes:
            raise ValueError("Empty class or species catalog")
        self.n_species, self.n_classes = n_species, n_classes
        self.class_names = [None] * n_classes
        self.species_names = [None] * n_species
        self.species_mw = np.full(n_species, np.nan)
        self.class_to_species = np.zeros((n_species, n_classes), dtype=float)
        self.mechanism_names = {}
        self.species_to_mechanism = {}
        self.fraction_sums = defaultdict(float)
        source_ids = []
        for row in c:
            ci = int(row["class_id"]) - 1
            si = int(row["species_id"]) - 1
            if not 0 <= ci < n_classes or not 0 <= si < n_species:
                raise ValueError(f"Out-of-range class/species ID: {row}")
            cn = row["class_name"]
            sn = row["species_name"]
            if self.class_names[ci] not in (None, cn) or self.species_names[si] not in (None, sn):
                raise ValueError(f"Conflicting species/class name: {row}")
            self.class_names[ci] = cn
            self.species_names[si] = sn
            self.species_mw[si] = float(row["molecular_weight_g_mol"])
            frac = float(row["fraction"])
            self.class_to_species[si, ci] += frac
            self.fraction_sums[cn] += frac
            source_ids.append(si + 1)
        if any(x is None for x in self.class_names):
            raise ValueError("Incomplete class catalog")
        # Fill all 201 species names from the complete original catalog,
        # including the ID omitted from the class-to-species map.
        if len(catalog) != n_species:
            raise ValueError("Species catalog length changed unexpectedly")
        for row in catalog:
            si = int(row["species_id"]) - 1
            sn = row["species_name"]
            if self.species_names[si] not in (None, sn):
                raise ValueError(f"Catalog name conflict for ID {si + 1}")
            self.species_names[si] = sn
            self.species_mw[si] = float(row["molecular_weight_g_mol"])
        self.duplicate_species_ids = {k: v for k, v in Counter(source_ids).items() if v > 1}
        self.missing_species_ids = sorted(set(range(1, n_species + 1)) - set(source_ids))
        mech_names = {}
        pairs = []
        for row in m:
            si = int(row["species_id"]) - 1
            mi = int(row["mechanism_id"])
            name = row["mechanism_name"]
            fac = float(row["conversion_factor"])
            if not 0 <= si < n_species:
                raise ValueError(f"Invalid source species ID: {si + 1}")
            if mech_names.get(mi, name) != name:
                raise ValueError(f"Conflicting mechanism ID {mi}")
            mech_names[mi] = name
            pairs.append((si, mi, fac))
        # Reject duplicate source -> target mappings before building any matrix.
        # The conversion coefficient is deliberately NOT part of the key.
        pair_counts = Counter((si, mi) for si, mi, _ in pairs)
        duplicates = {
            (self.species_names[si], mech_names[mi]): {
                "count": count,
                "conversion_factors": [
                    fac for source, target, fac in pairs if source == si and target == mi
                ],
            }
            for (si, mi), count in pair_counts.items()
            if count > 1
        }
        if duplicates:
            raise ValueError(
                "Duplicate source-to-mechanism mappings are not allowed "
                "(even with different conversion factors): "
                f"{duplicates}"
            )
        self.duplicate_mechanism_rows = {}
        self.mechanism_ids = sorted(mech_names)
        self.mechanism_names = [mech_names[k] for k in self.mechanism_ids]
        self.species_to_mechanism = np.zeros((len(self.mechanism_ids), n_species))
        ix = {v: i for i, v in enumerate(self.mechanism_ids)}
        for si, mi, fac in pairs:
            self.species_to_mechanism[ix[mi], si] += fac

    @staticmethod
    def _map(flux, matrix, n_in, name):
        a = np.asarray(flux, dtype=float)
        if a.ndim < 1 or a.shape[0] != n_in:
            raise ValueError(f"{name} must have leading dimension {n_in}, got {a.shape}")
        return np.tensordot(matrix, a, axes=(1, 0))

    def to_species(self, class_flux):
        """Return 201 species in catalog-ID order, same flux units as input."""
        return self._map(class_flux, self.class_to_species, self.n_classes, "class_flux")

    def to_mechanism(self, species_flux):
        """Return CBMZ species in mechanism_ids order, same flux units."""
        return self._map(species_flux, self.species_to_mechanism, self.n_species, "species_flux")

    def convert(self, class_flux, mechanism="CBMZ", return_species=False):
        """Convert 19 classes to CBMZ; mechanism=None returns 201 species."""
        species = self.to_species(class_flux)
        if mechanism is None or str(mechanism).upper() in ("NONE", "201", "SPECIES"):
            return species
        if str(mechanism).upper() != "CBMZ":
            raise ValueError("Only CBMZ mapping was supplied; no MUSICA table is available")
        mapped = self.to_mechanism(species)
        return (mapped, species) if return_species else mapped

    def diagnostics(self):
        return {
            "fraction_sums": dict(self.fraction_sums),
            "duplicate_species_ids": dict(self.duplicate_species_ids),
            "missing_species_ids": list(self.missing_species_ids),
            "duplicate_mechanism_rows": dict(self.duplicate_mechanism_rows),
            "unmapped_species": [
                self.species_names[i]
                for i in range(self.n_species)
                if not np.any(self.species_to_mechanism[:, i])
            ],
        }


def convert_emissions(class_flux, mechanism="CBMZ", input_dir=None, return_species=False):
    """Convenience entry point for the QUACS MEGAN3 adapter."""
    mapper = MGN2MECH(input_dir=input_dir)
    return mapper.convert(class_flux, mechanism=mechanism, return_species=return_species)
