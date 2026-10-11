"""Basic checks on the available grid-growth-form metadata."""
from pathlib import Path
import numpy as np
import pytest
from quacs.biogenic.megan.src.EFPreprocessor import read_grid_growthform

BASE = Path(__file__).resolve().parents[1]

def test_two_ecotypes():
    ecotypes, growth = read_grid_growthform(BASE/'inputs/grid_growthform.csv')
    assert len(ecotypes) == 2
    assert np.isclose(sum(frac for _, frac in ecotypes), 1)
    assert growth.shape == (4,)

def test_bad_ecotype_fraction(tmp_path):
    text=(BASE/'inputs/grid_growthform.csv').read_text().replace('1001,0.6,','1001,0.7,')
    path=tmp_path/'growth.csv';path.write_text(text)
    with pytest.raises(ValueError,match='sum to 1'):
        read_grid_growthform(path)
