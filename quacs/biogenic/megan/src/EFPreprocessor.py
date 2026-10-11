"""Single-grid MEGAN3.21 EF processor using EFv240120 and Speciation* tables.

Mirrors MEGEFP32 M3GEFP.build_interm_query/build_final_query: independently
weights species EF and LDF by species, growth-form and ecotype fractions.
Input/output EF units are preserved (no implicit conversion to nmol m-2 s-1).
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import numpy as np
import pandas as pd

def resolve_growth_forms(raw_ef_files):
    """Map arbitrary growth-form names to (speciation path, fraction column)."""
    spec = raw_ef_files['speciation']
    return {name: (Path(meta['file']), meta['specfrac_column']) for name, meta in spec.items()}


@dataclass(frozen=True)
class GridEF:
    ef: np.ndarray
    ldf: np.ndarray  # 4 entries: LDF3..LDF6, as in the EF processor output
    ecotypes: tuple
    growth_fractions: tuple


def _csv(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    # The distributed SpeciationTree file includes non-UTF8 characters.
    try:
        return pd.read_csv(path, encoding='utf-8-sig', low_memory=False)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding='latin1', low_memory=False)


def _fractions(values, label, total_one=False):
    a = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(a)) or np.any(a < 0) or np.any(a > 1):
        raise ValueError(f'{label}: fractions must be finite and between 0 and 1')
    if total_one and not np.isclose(a.sum(), 1, atol=1e-6):
        raise ValueError(f'{label}: fractions must sum to 1, got {a.sum()}')
    return a


def read_grid_growthform(path, growth_forms=None):
    """Read ecotype fractions and grid growth-form fractions for EF/LDF only."""
    frame = _csv(path)
    if growth_forms is None:
        growth_forms = ('tree', 'shrub', 'herb', 'crop')
    required = ['ecotype_id', 'ecotype_fraction'] + [f'{g}_fraction' for g in growth_forms]
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f'Grid growthform CSV missing columns: {sorted(missing)}')
    if frame.empty:
        raise ValueError('Grid growthform CSV is empty')
    ids = pd.to_numeric(frame.ecotype_id, errors='coerce').to_numpy(dtype=float)
    if not np.all(np.isfinite(ids)) or not np.all(ids == ids.astype(int)):
        raise ValueError('ecotype_id must contain integers')
    eco = ids.astype(int)
    if len(set(eco)) != len(eco):
        raise ValueError('ecotype_id must be unique')
    eco_frac = _fractions(frame.ecotype_fraction, 'ecotype', total_one=True)
    for key in required[2:]:
        if frame[key].nunique(dropna=False) != 1:
            raise ValueError(f'{key} must be identical on every ecotype row')
    growth = _fractions([frame.iloc[0][f'{g}_fraction'] for g in growth_forms], 'growth form')
    if growth.sum() > 1 + 1e-6 or growth.sum() <= 0:
        raise ValueError('Growth-form coverage must be >0 and <=1')
    return tuple(zip(eco.tolist(), eco_frac.tolist())), growth


def read_grid_canopy(path, *, n_canopy_types=6):
    """Read independent grid canopy fractions; IDs are zero-based (0..5)."""
    frame = _csv(path)
    required = ['canopy_type_id', 'canopy_fraction']
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f'Grid canopy CSV missing columns: {sorted(missing)}')
    if n_canopy_types != 6 or len(frame) != n_canopy_types:
        raise ValueError('Grid canopy CSV must have exactly 6 canopy types')
    ids = pd.to_numeric(frame.canopy_type_id, errors='coerce').to_numpy(dtype=float)
    if not np.all(np.isfinite(ids)) or not np.array_equal(np.sort(ids), np.arange(n_canopy_types)):
        raise ValueError('canopy_type_id must contain each integer 0..5 exactly once')
    canopy = np.zeros(n_canopy_types, dtype=float)
    canopy[ids.astype(int)] = _fractions(frame.canopy_fraction, 'canopy', total_one=True)
    return canopy


def process_grid(growthform_csv, input_dir=None, *, n_classes=19, raw_ef_files=None, ldf_classes=(3,4,5,6)):
    """Calculate grid EF1..EF19 and grid LDF3..LDF6, retaining source units."""
    if raw_ef_files is None:
        input_dir = Path(input_dir)
        raw_ef_files = {'plant_ef': str(input_dir / 'EFv240120.csv'), 'speciation': {g: {'file': str(input_dir / filename), 'specfrac_column': frac_col} for g,(filename,frac_col) in {'tree': ('SpeciationTree230921.csv','treeSpecfrac'), 'shrub': ('SpeciationShrub230822.csv','ShrubSpecFrac'), 'herb': ('SpeciationHerb230822.csv','HerbSpecFrac'), 'crop': ('SpeciationCrop230822.csv','CropSpecFrac')}.items()}}
    growth_forms = resolve_growth_forms(raw_ef_files)
    ecotypes, growth = read_grid_growthform(growthform_csv, growth_forms)
    ef_table = _csv(raw_ef_files['plant_ef'])
    columns = [f'EF{i}' for i in range(1, n_classes+1)] + [f'LDF{i}' for i in ldf_classes]
    if 'VegID' not in ef_table or any(c not in ef_table for c in columns):
        raise ValueError('EFv240120.csv missing VegID/EF1..EF19/LDF3..LDF6')
    if ef_table.VegID.duplicated().any():
        raise ValueError('Duplicate VegID in EFv240120.csv')
    ef_table = ef_table.set_index('VegID')
    results = np.zeros(len(columns), dtype=float)
    for g_idx, (growth_name, (filename, frac_col)) in enumerate(growth_forms.items()):
        if growth[g_idx] == 0:
            continue
        spec = _csv(filename)
        if not {'EcoTypeID', 'VegID', frac_col}.issubset(spec.columns):
            raise ValueError(f'{filename} missing required columns')
        for eco_id, eco_fraction in ecotypes:
            if eco_fraction == 0:
                continue
            subset = spec.loc[spec.EcoTypeID == eco_id]
            if subset.empty:
                raise ValueError(f'{filename}: EcoTypeID {eco_id} not found')
            fractions = _fractions(subset[frac_col], f'{filename} ecotype {eco_id}', total_one=True)
            missing_veg = set(subset.VegID) - set(ef_table.index)
            if missing_veg:
                raise ValueError(f'{filename}: VegID missing in EFv240120.csv: {sorted(missing_veg)[:10]}')
            vals = ef_table.loc[subset.VegID, columns].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
            if not np.all(np.isfinite(vals)) or np.any(vals[:, :n_classes] < 0) or np.any((vals[:, n_classes:] < 0) | (vals[:, n_classes:] > 1)):
                raise ValueError(f'Invalid EF/LDF values for ecotype {eco_id} in {filename}')
            results += eco_fraction * growth[g_idx] * (fractions @ vals)
    return GridEF(results[:n_classes], results[n_classes:], ecotypes, tuple(growth))


def write_grid_ef(path, result, *, class_names, ldf):
    """Write the original EF_LDF.csv layout (19 rows, one per MEGAN class).

    Output columns match adapter.load_species_parameters().
    """
    path = Path(path)
    ef_values = np.asarray(result.ef, dtype=float)
    ldf_values = np.asarray(ldf, dtype=float)
    if len(class_names) != len(ef_values) or ldf_values.shape != ef_values.shape:
        raise ValueError('EF_LDF output requires 19 class names, EF values, and LDF values')
    if not np.all(np.isfinite(ef_values)) or np.any(ef_values < 0):
        raise ValueError('Invalid EF values')
    if not np.all(np.isfinite(ldf_values)) or np.any((ldf_values < 0) | (ldf_values > 1)):
        raise ValueError('Invalid LDF values')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Categories', 'Emission Factor (nmol/m-2/s-1)', 'LDF'])
        for name, emission_factor, light_fraction in zip(class_names, ef_values, ldf_values):
            writer.writerow([name, float(emission_factor), float(light_fraction)])


def validate_growthform_count(raw_ef_files, expected):
    count = len(resolve_growth_forms(raw_ef_files))
    if count != expected:
        raise ValueError(f'NGFORM={expected}, configured speciation files={count}')

def assemble_ldf(class_names, activity_module, grid_result):
    """Match class-specific fixed/GRID LDF values without YAML class IDs."""
    grid_ids, ldf = activity_module.get_ldf_configuration(class_names)
    if len(grid_ids) != len(grid_result.ldf):
        raise ValueError('Grid LDF count differs from activity parameter GRID entries')
    for class_id, value in zip(grid_ids, grid_result.ldf):
        ldf[class_id - 1] = value
    if not np.isfinite(ldf).all() or np.any((ldf < 0) | (ldf > 1)):
        raise ValueError('Invalid combined LDF')
    return ldf
