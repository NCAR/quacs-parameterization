"""Strict mapping validation using minimal standalone catalogs."""

import csv

import pytest

from quacs.biogenic.megan.src.MGN2MECH import MGN2MECH


def _write(path, filename, fields, rows):
    with (path / filename).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _mapper(tmp_path, pairs):
    _write(
        tmp_path,
        "species_catalog.csv",
        ["species_id", "species_name", "molecular_weight_g_mol"],
        [{"species_id": 1, "species_name": "a", "molecular_weight_g_mol": 100}],
    )
    _write(
        tmp_path,
        "class_to_species.csv",
        [
            "class_id",
            "class_name",
            "species_id",
            "species_name",
            "molecular_weight_g_mol",
            "fraction",
        ],
        [
            {
                "class_id": 1,
                "class_name": "X",
                "species_id": 1,
                "species_name": "a",
                "molecular_weight_g_mol": 100,
                "fraction": 1,
            }
        ],
    )
    _write(
        tmp_path,
        "species_to_CBMZ.csv",
        ["species_id", "mechanism_id", "mechanism_name", "conversion_factor"],
        [
            {"species_id": 1, "mechanism_id": mi, "mechanism_name": name, "conversion_factor": fac}
            for mi, name, fac in pairs
        ],
    )
    return MGN2MECH(input_dir=tmp_path)


@pytest.mark.parametrize("factors", [(1.0, 1.0), (1.0, 0.5)])
def test_duplicate_source_target_fails(tmp_path, factors):
    with pytest.raises(ValueError, match="Duplicate source-to-mechanism mappings"):
        _mapper(tmp_path, [(1, "SESQ", factors[0]), (1, "SESQ", factors[1])])


def test_distinct_targets_allowed(tmp_path):
    mapper = _mapper(tmp_path, [(1, "SESQ", 1.0), (2, "TERP", 0.5)])
    assert mapper.species_to_mechanism.shape == (2, 1)
