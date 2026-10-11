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
import warnings
import numpy as np

INPUTS = Path(__file__).resolve().parent.parent / 'inputs'


def _load_csv(path):
    with Path(path).open(newline='') as fh:
        return list(csv.DictReader(fh))


class MGN2MECH:
    """Reusable, vectorized emission-speciation mapper.

    Parameters
    ----------
    input_dir : str or Path, optional
        Directory containing class_to_species.csv and species_to_CBMZ.csv.
    duplicate_policy : {'warn', 'error', 'preserve'}
        Whether to warn/error/preserve repeated identical source-target rows.
        'preserve' follows the supplied Fortran mapping literally.
    """
    def __init__(self, input_dir=None, duplicate_policy='warn'):
        if duplicate_policy not in ('warn','error','preserve'):
            raise ValueError('duplicate_policy must be warn, error, or preserve')
        d = Path(input_dir) if input_dir is not None else INPUTS
        c = _load_csv(d/'class_to_species.csv')
        m = _load_csv(d/'species_to_CBMZ.csv')
        if len(c)!=201 or len(m)!=148:
            raise ValueError('Unexpected mapping length; expected 201 and 148')
        self.class_names = [None]*19
        self.species_names = [None]*201
        self.species_mw = np.full(201,np.nan)
        self.class_to_species = np.zeros((201,19),dtype=float)
        self.mechanism_names = {}
        self.species_to_mechanism = {}
        self.fraction_sums = defaultdict(float)
        source_ids=[]
        for row in c:
            ci=int(row['class_id'])-1; si=int(row['species_id'])-1
            if not 0<=ci<19 or not 0<=si<201:
                raise ValueError(f'Out-of-range class/species ID: {row}')
            cn=row['class_name']; sn=row['species_name']
            if self.class_names[ci] not in (None,cn) or self.species_names[si] not in (None,sn):
                raise ValueError(f'Conflicting species/class name: {row}')
            self.class_names[ci]=cn; self.species_names[si]=sn
            self.species_mw[si]=float(row['molecular_weight_g_mol'])
            frac=float(row['fraction'])
            self.class_to_species[si,ci] += frac
            self.fraction_sums[cn] += frac
            source_ids.append(si+1)
        if any(x is None for x in self.class_names):
            raise ValueError('Incomplete 19-class catalog')
        # Fill all 201 species names from the complete original catalog,
        # including the ID omitted from the class-to-species map.
        catalog = _load_csv(d/'species_catalog.csv')
        if len(catalog)!=201:
            raise ValueError('Expected 201 species catalog entries')
        for row in catalog:
            si=int(row['species_id'])-1
            sn=row['species_name']
            if self.species_names[si] not in (None,sn):
                raise ValueError(f'Catalog name conflict for ID {si+1}')
            self.species_names[si]=sn
            self.species_mw[si]=float(row['molecular_weight_g_mol'])
        self.duplicate_species_ids={k:v for k,v in Counter(source_ids).items() if v>1}
        self.missing_species_ids=sorted(set(range(1,202))-set(source_ids))
        mech_names={}
        pairs=[]
        for row in m:
            si=int(row['species_id'])-1; mi=int(row['mechanism_id'])
            name=row['mechanism_name']; fac=float(row['conversion_factor'])
            if not 0<=si<201:
                raise ValueError(f'Invalid source species ID: {si+1}')
            if mech_names.get(mi,name)!=name:
                raise ValueError(f'Conflicting mechanism ID {mi}')
            mech_names[mi]=name
            pairs.append((si,mi,fac))
        self.mechanism_ids=sorted(mech_names)
        self.mechanism_names=[mech_names[k] for k in self.mechanism_ids]
        self.species_to_mechanism=np.zeros((len(self.mechanism_ids),201))
        ix={v:i for i,v in enumerate(self.mechanism_ids)}
        for si,mi,fac in pairs:
            self.species_to_mechanism[ix[mi],si]+=fac
        pair_counts=Counter((si,mi,fac) for si,mi,fac in pairs)
        self.duplicate_mechanism_rows={
            (self.species_names[si],mech_names[mi],fac):n
            for (si,mi,fac),n in pair_counts.items() if n>1
        }
        if self.duplicate_mechanism_rows:
            msg=f'Duplicate mechanism mapping rows (retained): {self.duplicate_mechanism_rows}'
            if duplicate_policy=='error': raise ValueError(msg)
            if duplicate_policy=='warn': warnings.warn(msg,UserWarning,stacklevel=2)

    @staticmethod
    def _map(flux, matrix, n_in, name):
        a=np.asarray(flux,dtype=float)
        if a.ndim<1 or a.shape[0]!=n_in:
            raise ValueError(f'{name} must have leading dimension {n_in}, got {a.shape}')
        return np.tensordot(matrix,a,axes=(1,0))

    def to_species(self, class_flux):
        """Return 201 species in catalog-ID order, same flux units as input."""
        return self._map(class_flux,self.class_to_species,19,'class_flux')

    def to_mechanism(self, species_flux):
        """Return CBMZ species in mechanism_ids order, same flux units."""
        return self._map(species_flux,self.species_to_mechanism,201,'species_flux')

    def convert(self, class_flux, mechanism='CBMZ', return_species=False):
        """Convert 19 classes to CBMZ; mechanism=None returns 201 species."""
        species=self.to_species(class_flux)
        if mechanism is None or str(mechanism).upper() in ('NONE','201','SPECIES'):
            return species
        if str(mechanism).upper()!='CBMZ':
            raise ValueError('Only CBMZ mapping was supplied; no MUSICA table is available')
        mapped=self.to_mechanism(species)
        return (mapped,species) if return_species else mapped

    def diagnostics(self):
        return {
            'fraction_sums':dict(self.fraction_sums),
            'duplicate_species_ids':dict(self.duplicate_species_ids),
            'missing_species_ids':list(self.missing_species_ids),
            'duplicate_mechanism_rows':dict(self.duplicate_mechanism_rows),
            'unmapped_species':[
                self.species_names[i] for i in range(201)
                if not np.any(self.species_to_mechanism[:,i])
            ],
        }


def convert_emissions(class_flux, mechanism='CBMZ', input_dir=None,
                      return_species=False, duplicate_policy='warn'):
    """Convenience entry point for the QUACS MEGAN3 adapter."""
    mapper=MGN2MECH(input_dir=input_dir,duplicate_policy=duplicate_policy)
    return mapper.convert(class_flux,mechanism=mechanism,return_species=return_species)
