from pathlib import Path
import pytest
from quacs.biogenic.megan.src.config import MeganConfig


def test_nested_files_and_mapping_paths():
    config = MeganConfig.from_yaml(Path(__file__).resolve().parents[1] / 'configure.yaml')
    assert set(config.files) == {'inputs', 'intermediate_outputs', 'outputs'}
    assert set(config.mapping_files) == {'species_catalog', 'class_to_species', 'species_to_CBMZ'}
    assert config.path(config.intermediate_outputs['ef_ldf_output']).parent.name == 'intermediate_outputs'
    assert config.dimensions['n_class'] == 19


def test_missing_mapping_path_rejected(tmp_path):
    path = Path(__file__).resolve().parents[1] / 'configure.yaml'
    content = path.read_text().replace('species_to_CBMZ: inputs/species_to_CBMZ.csv', 'species_to_CBMZ_bad: inputs/species_to_CBMZ.csv')
    new_path = tmp_path / 'configure.yaml'
    new_path.write_text(content)
    with pytest.raises(ValueError, match='mapping_files.species_to_CBMZ'):
        MeganConfig.from_yaml(new_path)
