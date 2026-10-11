"""Check original MEGAN class metadata against the new configurable tables."""
from quacs.biogenic.megan.src.MGN2MECH import MGN2MECH
from quacs.biogenic.megan.src.MEGVEA import ACTIVITY_CLASS_NAMES

def test_official_class_metadata():
    mapper = MGN2MECH()
    assert mapper.class_names == list(ACTIVITY_CLASS_NAMES[:19])
    assert len(mapper.species_names) == 201
