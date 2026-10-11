"""Regression guard for the MEGAN3.21 class labels and unchanged demo EF values."""
import csv
from pathlib import Path
from quacs.biogenic.megan.src.MGN2MECH import MGN2MECH

def test_ef_labels_and_values():
    base = Path(__file__).resolve().parents[1]
    with (base / 'inputs' / 'EF_LDF.csv').open(newline='') as fh:
        rows = list(csv.DictReader(fh))
    mapper = MGN2MECH(duplicate_policy='preserve')
    assert [r['Categories'] for r in rows] == mapper.class_names
    assert [float(r['Emission Factor (nmol/m-2/s-1)']) for r in rows] == [
        10, 0, 0.3, 0.09, 0.06, 0.12, 0.11, 0.075, 0.009,
        0.039, 1, 0.18, 0.18, 0.12, 0.6, 0.01, 0.01, 0.18, 0.003,
    ]
