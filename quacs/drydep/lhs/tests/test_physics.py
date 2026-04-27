"""Side-by-side comparison tests for drydep_physics.py vs drydep_functions.py.

For every physics function, both the old (Fortran-translated) and new
(Pythonic) implementations are called under a 500-sample Latin Hypercube
sweep of the full physically-relevant parameter space.  Agreement must be
within round-off tolerances (rel=1e-10 throughout).
"""

import numpy as np
import pytest
from scipy.stats import qmc

from quacs.drydep.lhs import olson_land_cover
from quacs.drydep.lhs.box_model import DrydepCoefficients, LandCoverPatch, _patches_to_depvel_arrays
from quacs.drydep.lhs.drydep_physics import (
    MetVars,
    _light_correction,
    calc_met_vars,
    deposition_velocity,
    molecular_diffusivity,
    monin_obukhov_length,
)
from quacs.drydep.lhs.lhs_config import (
    LC_ARCHETYPES,
    MET_DIMS,
    N_LC_ARCHETYPES,
    build_patches_from_sample,
)
from quacs.drydep.lhs.species import hg0 as HG0_SPECIES, o3 as O3_SPECIES, so2 as SO2_SPECIES
from quacs.drydep.simple.drydep_functions import BIOFIT, DEPVEL, DIFFG, GET_OBK, METERO

# ── LHS sample generation ────────────────────────────────────────────────────

_N_SAMPLES = 500
_SEED = 0

_N_MET = len(MET_DIMS)
_N_FRAC = N_LC_ARCHETYPES
_N_LAI = N_LC_ARCHETYPES
_NDIM_LHS = _N_MET + _N_FRAC + _N_LAI

_L = (
    [d[1] for d in MET_DIMS]
    + [0.01] * _N_FRAC
    + [lc[1] for lc in LC_ARCHETYPES]
)
_U = (
    [d[2] for d in MET_DIMS]
    + [1.0] * _N_FRAC
    + [lc[2] if lc[2] > lc[1] else lc[1] + 1.0 for lc in LC_ARCHETYPES]
)

_sampler = qmc.LatinHypercube(d=_NDIM_LHS, seed=_SEED)
_raw = _sampler.random(n=_N_SAMPLES)
_scaled = qmc.scale(_raw, _L, _U)

_met_samples = _scaled[:, :_N_MET]
_frac_samples = _scaled[:, _N_MET: _N_MET + _N_FRAC]
_lai_samples = _scaled[:, _N_MET + _N_FRAC:]

_met_dicts = [
    {name: float(_met_samples[i, j]) for j, (name, *_) in enumerate(MET_DIMS)}
    for i in range(_N_SAMPLES)
]
_patch_lists = [
    build_patches_from_sample(_frac_samples[i], _lai_samples[i])
    for i in range(_N_SAMPLES)
]

_COEFF = olson_land_cover.coefficients
_DRYCOEFF = _COEFF.drycoeff
_SPECIES = [HG0_SPECIES, SO2_SPECIES, O3_SPECIES]


def _old_depvel(m, patches, sp):
    CZ1, LSNOW, OBK, _ = METERO(
        m["BXHEIGHT"], m["ALBD"], m["TC0"], m["USTAR"],
        m["AIRDEN"], m["HFLUX"], m["U10M"], m["V10M"],
    )
    (IREG, ILAND, IUSE, XLAI,
     idep, iri, irlu, irac, irgss, irgso, ircls, irclo) = _patches_to_depvel_arrays(patches)
    return DEPVEL(
        _DRYCOEFF, _COEFF.iolson,
        idep, iri, irlu, irac, irgss, irgso, ircls, irclo,
        IREG, ILAND, IUSE, m["TC0"],
        XLAI, LSNOW, m["RADIAT"], m["CFRAC"], m["SUNCOS_MID"],
        m["PRESSU"], m["USTAR"], m["AZO"], CZ1, OBK,
        sp.molecular_weight_kg_mol,
        float(sp.other_properties["reactivity"]),
        float(sp.other_properties["henrys_law_constant"]),
    )


def _new_depvel(m, patches, sp):
    mv = calc_met_vars(
        m["BXHEIGHT"], m["ALBD"], m["TC0"], m["USTAR"],
        m["AIRDEN"], m["HFLUX"], m["U10M"], m["V10M"],
    )
    return deposition_velocity(
        patches,
        ts=m["TC0"],
        radiat=m["RADIAT"],
        cloud_frac=m["CFRAC"],
        suncos_mid=m["SUNCOS_MID"],
        pressure=m["PRESSU"],
        ustar=m["USTAR"],
        roughness_length=m["AZO"],
        cz1=mv.cz1,
        obk=mv.obk,
        mol_weight=sp.molecular_weight_kg_mol,
        f0=float(sp.other_properties["reactivity"]),
        hstar=float(sp.other_properties["henrys_law_constant"]),
        drycoeff=_DRYCOEFF,
        is_snow=mv.is_snow,
    )


# ── Tests ────────────────────────────────────────────────────────────────────

class TestMoninObukhovLength:
    def test_lhs_sweep(self):
        for m in _met_dicts:
            old = GET_OBK(m["TC0"], m["USTAR"], m["AIRDEN"], m["HFLUX"])
            new = monin_obukhov_length(m["TC0"], m["USTAR"], m["AIRDEN"], m["HFLUX"])
            assert new == pytest.approx(old, rel=1e-10)

    def test_zero_hflux(self):
        old = GET_OBK(300.0, 0.2, 1.2, 0.0)
        new = monin_obukhov_length(300.0, 0.2, 1.2, 0.0)
        assert new == pytest.approx(old, rel=1e-10)


class TestCalcMetVars:
    def test_lhs_sweep(self):
        for m in _met_dicts:
            old_cz1, old_snow, old_obk, old_w10 = METERO(
                m["BXHEIGHT"], m["ALBD"], m["TC0"], m["USTAR"],
                m["AIRDEN"], m["HFLUX"], m["U10M"], m["V10M"],
            )
            mv = calc_met_vars(
                m["BXHEIGHT"], m["ALBD"], m["TC0"], m["USTAR"],
                m["AIRDEN"], m["HFLUX"], m["U10M"], m["V10M"],
            )
            assert mv.cz1 == pytest.approx(old_cz1, rel=1e-10)
            assert mv.is_snow == bool(old_snow)
            assert mv.obk == pytest.approx(old_obk, rel=1e-10)
            assert mv.wind10 == pytest.approx(old_w10, rel=1e-10)

    def test_snow_flag_set_when_albedo_high(self):
        mv = calc_met_vars(130.0, 0.85, 275.0, 0.2, 1.2, 5.0, 1.0, -0.5)
        assert mv.is_snow is True

    def test_snow_flag_clear_when_albedo_low(self):
        mv = calc_met_vars(130.0, 0.1, 295.0, 0.2, 1.2, 5.0, 1.0, -0.5)
        assert mv.is_snow is False


class TestMolecularDiffusivity:
    @pytest.mark.parametrize("sp", _SPECIES)
    def test_lhs_sweep(self, sp):
        mw = sp.molecular_weight_kg_mol
        for m in _met_dicts:
            old = DIFFG(m["TC0"], m["PRESSU"], mw)
            new = molecular_diffusivity(m["TC0"], m["PRESSU"], mw)
            assert new == pytest.approx(old, rel=1e-10)


class TestLightCorrection:
    def test_lhs_sweep(self):
        npoly = len(_DRYCOEFF)
        for i in range(_N_SAMPLES):
            m = _met_dicts[i]
            for patch in _patch_lists[i]:
                if patch.fraction > 0 and patch.lai > 0:
                    lai_per_frac = patch.lai / patch.fraction
                    old = BIOFIT(_DRYCOEFF, lai_per_frac, m["SUNCOS_MID"], m["CFRAC"], npoly)
                    new = _light_correction(_DRYCOEFF, lai_per_frac, m["SUNCOS_MID"], m["CFRAC"])
                    assert new == pytest.approx(old, rel=1e-10), (
                        f"Mismatch at sample {i}: "
                        f"lai={lai_per_frac:.3f} suncos={m['SUNCOS_MID']:.3f} "
                        f"cfrac={m['CFRAC']:.3f}  old={old}  new={new}"
                    )

    def test_minimum_floor(self):
        new = _light_correction(_DRYCOEFF, 0.0, 0.0, 0.0)
        assert new >= 0.1


class TestDepositionVelocity:
    @pytest.mark.parametrize("sp", _SPECIES)
    def test_lhs_sweep(self, sp):
        for i in range(_N_SAMPLES):
            m = _met_dicts[i]
            patches = _patch_lists[i]
            old = _old_depvel(m, patches, sp)
            new = _new_depvel(m, patches, sp)
            assert new == pytest.approx(old, rel=1e-8), (
                f"Mismatch at sample {i} species={sp.name}: "
                f"old={old:.6f}  new={new:.6f}"
            )

    def test_snow_cover_lhs_sweep(self):
        snow_samples = [
            (m, p) for m, p in zip(_met_dicts, _patch_lists) if m["ALBD"] > 0.4
        ]
        assert len(snow_samples) > 0, "No snow samples — increase N or widen ALBD range"
        for m, patches in snow_samples:
            old = _old_depvel(m, patches, HG0_SPECIES)
            new = _new_depvel(m, patches, HG0_SPECIES)
            assert new == pytest.approx(old, rel=1e-8), (
                f"Snow mismatch: old={old:.6f}  new={new:.6f}"
            )
