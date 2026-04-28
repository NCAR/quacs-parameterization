"""Tests for quacs.drydep.lhs.box_model.

The canonical reference value (0.1287817554797227 cm/s) is produced by
running ex_drydep_simple.py with the same inputs against the unmodified
DEPVEL/METERO routines in drydep_functions.py.
"""

import musica.mechanism_configuration as mc
import numpy as np
import pytest

from quacs.drydep.lhs import olson_land_cover
from quacs.drydep.lhs.box_model import LandCoverPatch, compute_drydep_rate
from quacs.drydep.lhs.species import hg0 as HG0_SPECIES
from quacs.drydep.lhs.species import o3 as O3_SPECIES
from quacs.drydep.lhs.species import so2 as SO2_SPECIES

TEST_COEFFICIENTS = olson_land_cover.coefficients

TEST_MET = dict(
    TC0=302.0,
    CFRAC=0.8,
    RADIAT=606.0,
    AZO=1.5,
    USTAR=0.15,
    PRESSU=1e5,
    SUNCOS_MID=0.8,
    ALBD=0.1,
    BXHEIGHT=130.0,
    U10M=0.1,
    V10M=-0.05,
    AIRDEN=1.15,
    HFLUX=16.8,
)

TEST_LAND_COVER = [
    olson_land_cover.patch_from_olson(0, fraction=0.312, lai=0.10),
    olson_land_cover.patch_from_olson(29, fraction=0.200, lai=0.71),
    olson_land_cover.patch_from_olson(33, fraction=0.175, lai=0.61),
    olson_land_cover.patch_from_olson(43, fraction=0.313, lai=0.74),
]

REFERENCE_DVEL = 0.1287817554797227


class TestComputeDepdryRate:
    def test_default_matches_reference_value(self):
        """Must reproduce the result from ex_drydep_simple.py exactly."""
        dvel, _ = compute_drydep_rate(
            [HG0_SPECIES],
            box_height_m=1.0,
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert dvel == pytest.approx(REFERENCE_DVEL, rel=1e-10)

    def test_returns_positive_velocity(self):
        dvel, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert dvel > 0.0

    def test_k_formula(self):
        """k must equal v_d [cm/s] × 0.01 / H [m]."""
        box_height = 3.5
        dvel, k = compute_drydep_rate(
            [HG0_SPECIES],
            box_height_m=box_height,
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert k == pytest.approx(dvel * 0.01 / box_height, rel=1e-12)

    def test_box_height_scales_k(self):
        """Doubling box height must halve k (v_d is independent of box height)."""
        _, k1 = compute_drydep_rate(
            [HG0_SPECIES],
            box_height_m=1.0,
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        _, k2 = compute_drydep_rate(
            [HG0_SPECIES],
            box_height_m=2.0,
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert k1 == pytest.approx(2.0 * k2, rel=1e-10)

    def test_snow_reduces_velocity(self):
        """Snow-covered surface (albedo > 0.4) should lower deposition velocity."""
        dvel_normal, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        dvel_snow, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met={**TEST_MET, "ALBD": 0.8},
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert dvel_snow < dvel_normal

    def test_temperature_affects_velocity(self):
        dvel_ref, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        dvel_hot, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met={**TEST_MET, "TC0": 320.0},
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert dvel_ref != pytest.approx(dvel_hot, rel=1e-4)

    def test_met_does_not_mutate_caller_dict(self):
        caller_met = dict(TEST_MET)
        compute_drydep_rate(
            [HG0_SPECIES],
            met=caller_met,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert caller_met == TEST_MET


class TestSpeciesPresets:
    def test_presets_are_mc_species(self):
        for sp in [HG0_SPECIES, SO2_SPECIES, O3_SPECIES]:
            assert isinstance(sp, mc.Species)

    def test_hg0_preset_fields(self):
        assert HG0_SPECIES.name == "Hg0"
        assert HG0_SPECIES.molecular_weight_kg_mol == pytest.approx(201e-3)
        assert float(HG0_SPECIES.other_properties["henrys_law_constant"]) == pytest.approx(0.11)
        assert float(HG0_SPECIES.other_properties["reactivity"]) == pytest.approx(3e-5)

    def test_so2_preset_fields(self):
        assert SO2_SPECIES.name == "SO2"
        assert SO2_SPECIES.molecular_weight_kg_mol == pytest.approx(64e-3)
        assert float(SO2_SPECIES.other_properties["henrys_law_constant"]) == pytest.approx(1e5)
        assert float(SO2_SPECIES.other_properties["reactivity"]) == pytest.approx(0.0)

    def test_o3_preset_fields(self):
        assert O3_SPECIES.name == "O3"
        assert O3_SPECIES.molecular_weight_kg_mol == pytest.approx(48e-3)
        assert float(O3_SPECIES.other_properties["henrys_law_constant"]) == pytest.approx(1e-2)
        assert float(O3_SPECIES.other_properties["reactivity"]) == pytest.approx(1.0)

    def test_species_produce_different_velocities(self):
        dvel_hg, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        dvel_o3, _ = compute_drydep_rate(
            [O3_SPECIES],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        dvel_so2, _ = compute_drydep_rate(
            [SO2_SPECIES],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert dvel_hg != pytest.approx(dvel_o3, rel=1e-3)
        assert dvel_hg != pytest.approx(dvel_so2, rel=1e-3)

    def test_custom_species_other_properties(self):
        custom = mc.Species(
            name="TestGas",
            molecular_weight_kg_mol=30e-3,
            other_properties={"henrys_law_constant": "0.5", "reactivity": "0.5"},
        )
        dvel, k = compute_drydep_rate(
            [custom],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert dvel > 0.0


class TestOlsonLandCover:
    def test_patch_from_olson_returns_patch(self):
        patch = olson_land_cover.patch_from_olson(29, fraction=0.5, lai=0.71)
        assert isinstance(patch, LandCoverPatch)
        assert patch.fraction == pytest.approx(0.5)
        assert patch.lai == pytest.approx(0.71)

    def test_named_presets_return_patches(self):
        for constructor in [
            olson_land_cover.tropical_forest,
            olson_land_cover.boreal_forest,
            olson_land_cover.savanna,
            olson_land_cover.water,
            olson_land_cover.cropland,
        ]:
            patch = constructor()
            assert isinstance(patch, LandCoverPatch)
            assert 0.0 <= patch.fraction <= 1.0

    def test_custom_land_cover_changes_velocity(self):
        dvel_default, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met=TEST_MET,
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        water_patch = [olson_land_cover.water(fraction=1.0, lai=0.0)]
        dvel_water, _ = compute_drydep_rate(
            [HG0_SPECIES],
            met=TEST_MET,
            land_cover=water_patch,
            coefficients=TEST_COEFFICIENTS,
        )
        assert dvel_default != pytest.approx(dvel_water, rel=1e-3)

    def test_land_cover_fractions_sum_to_one(self):
        total = sum(p.fraction for p in TEST_LAND_COVER)
        assert total == pytest.approx(1.0, abs=1e-6)


class TestVectorizedMode:
    def test_vectorized_matches_scalar_per_cell(self):
        N = 4
        temps = np.linspace(280, 310, N)
        result = compute_drydep_rate(
            [HG0_SPECIES],
            met={**TEST_MET, "TC0": temps},
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        k_vec = result["Hg0"]
        for i, temp in enumerate(temps):
            _, k_scalar = compute_drydep_rate(
                [HG0_SPECIES],
                met={**TEST_MET, "TC0": float(temp)},
                land_cover=TEST_LAND_COVER,
                coefficients=TEST_COEFFICIENTS,
            )
            assert k_vec[i] == pytest.approx(k_scalar, rel=1e-10)

    def test_vectorized_returns_dict(self):
        result = compute_drydep_rate(
            [HG0_SPECIES, SO2_SPECIES],
            met={**TEST_MET, "TC0": np.array([295.0, 300.0])},
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        assert isinstance(result, dict)
        assert "Hg0" in result
        assert "SO2" in result
        assert len(result["Hg0"]) == 2

    def test_vectorized_multiple_species_shape(self):
        N = 10
        result = compute_drydep_rate(
            [HG0_SPECIES, SO2_SPECIES, O3_SPECIES],
            met={**TEST_MET, "TC0": np.linspace(275, 305, N)},
            land_cover=TEST_LAND_COVER,
            coefficients=TEST_COEFFICIENTS,
        )
        for sp in [HG0_SPECIES, SO2_SPECIES, O3_SPECIES]:
            assert result[sp.name].shape == (N,)
            assert np.all(result[sp.name] > 0)

    def test_per_cell_land_cover(self):
        patches_cell0 = [olson_land_cover.water(fraction=1.0)]
        patches_cell1 = [olson_land_cover.tropical_forest(fraction=1.0)]
        result = compute_drydep_rate(
            [HG0_SPECIES],
            met={**TEST_MET, "TC0": np.array([302.0, 302.0])},
            land_cover=[patches_cell0, patches_cell1],
            coefficients=TEST_COEFFICIENTS,
        )
        assert result["Hg0"][0] != pytest.approx(result["Hg0"][1], rel=1e-3)
