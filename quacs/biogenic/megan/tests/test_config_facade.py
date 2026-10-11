from pathlib import Path
import pytest
from quacs.biogenic.megan.src.config import MeganConfig
from quacs.biogenic.megan.biogenic_emission_megan_v3 import MeganModel


def test_config_reads_declared_dimensions():
    path = Path(__file__).resolve().parents[1] / 'configure.yaml'
    config = MeganConfig.from_yaml(path)
    assert config.dimensions['mechanism'] == 'CBMZ'
    assert config.dimensions['n_class'] == 19
    assert config.dimensions['NGFORM'] == 4
    assert len(config.raw_ef_files['speciation']) == config.dimensions['NGFORM']
    assert all(p.is_file() for p in config.mapping_files.values())


def test_config_rejects_bad_dimensions(tmp_path):
    path = tmp_path / 'configure.yaml'
    path.write_text('files: {}\nmodel: {n_class: 0, n_species: 201, NLayers: 5, NRTYP: 6, NGFORM: 4}\n')
    with pytest.raises(ValueError, match='n_class'):
        MeganConfig.from_yaml(path)


def test_facade_rejects_direct_injection():
    model = object.__new__(MeganModel)
    with pytest.raises(ValueError, match='managed by MeganModel'):
        model.run(emission_factor=[1])


def test_config_rejects_unsupported_mechanism(tmp_path):
    path = tmp_path / 'configure.yaml'
    path.write_text('files: {}\nmodel: {mechanism: MOZART, n_class: 19, n_species: 201, NLayers: 5, NRTYP: 6, NGFORM: 4}\n')
    with pytest.raises(ValueError, match='Unsupported model.mechanism'):
        MeganConfig.from_yaml(path)

def test_model_mechanism_uses_configuration():
    model = object.__new__(MeganModel)
    model.mechanism = 'SPECIES'
    class Mapper:
        def convert(self, flux, mechanism, return_species):
            return mechanism, return_species
    model.mapper = Mapper()
    assert model.to_mechanism([1], return_species=True) == ('SPECIES', True)
