"""Tests for configuration-driven dimensions (no original large EF tables needed)."""

import csv
from pathlib import Path

import numpy as np

from quacs.biogenic.megan.src import MEGCAN, MEGVEA
from quacs.biogenic.megan.src.EFPreprocessor import process_grid
from quacs.biogenic.megan.src.MGN2MECH import MGN2MECH


def test_canopy_named_rows_extra_and_reordered(tmp_path):
    original = Path(MEGCAN.__file__).resolve().parents[1] / "inputs/canopy_characteristic.csv"
    rows = list(csv.DictReader(original.open()))
    columns = list(rows[0])
    rows.reverse()
    rows.append(dict(Characteristic="new_unused_characteristic", **{x: 2 for x in columns[1:]}))
    path = tmp_path / "canopy.csv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)
    loaded, names = MEGCAN.load_canopy_characteristics(path)
    assert len(loaded) == 18
    assert names == MEGCAN.CANOPY_TYPE_NAMES
    assert np.allclose(loaded["canopy_height_m"], MEGCAN.CANOPY_CHARACTERISTICS["canopy_height_m"])


def test_activity_extra_class(tmp_path):
    original = Path(MEGVEA.__file__).resolve().parents[1] / "inputs/activity_class_parameters.csv"
    rows = list(csv.DictReader(original.open()))
    rows.append(dict(rows[-1], class_name="NEWCLASS"))
    path = tmp_path / "activity.csv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    names, params = MEGVEA.load_activity_parameters(path)
    assert len(names) == 21 and len(params["beta"]) == 21


def test_dynamic_growthform_and_class_count(tmp_path):
    growth = tmp_path / "growth.csv"
    growth.write_text("ecotype_id,ecotype_fraction,moss_fraction\n1,1.0,1.0\n")
    ef = tmp_path / "ef.csv"
    ef.write_text("VegID,EF1,EF2,LDF3\n7,2.0,3.0,0.4\n")
    spec = tmp_path / "spec.csv"
    spec.write_text("EcoTypeID,VegID,MossSpecFrac\n1,7,1.0\n")
    raw = {
        "plant_ef": str(ef),
        "speciation": {"moss": {"file": str(spec), "specfrac_column": "MossSpecFrac"}},
    }
    result = process_grid(growth, n_classes=2, raw_ef_files=raw, ldf_classes=(3,))
    assert np.allclose(result.ef, [2, 3])
    assert np.allclose(result.ldf, [0.4])


def test_dynamic_mapper_sizes(tmp_path):
    source = Path(MEGCAN.__file__).resolve().parents[1] / "inputs"
    for name in ("class_to_species.csv", "species_catalog.csv", "species_to_CBMZ.csv"):
        (tmp_path / name).write_bytes((source / name).read_bytes())
    mapper = MGN2MECH(input_dir=tmp_path)
    assert mapper.n_classes == 19 and mapper.n_species == 201
    # Additional unmapped species can be represented in catalog without source-code edits.
    with (tmp_path / "species_catalog.csv").open("a") as f:
        f.write("202,synthetic_new_species,50\n")
    mapper = MGN2MECH(input_dir=tmp_path)
    assert mapper.n_species == 202
    assert mapper.to_species(np.ones(19)).shape == (202,)
