"""
Standard trace-gas species presets for dry deposition.

Each species carries the two properties required by ``compute_drydep_rate``:

- ``henrys_law_constant`` – effective Henry's law constant at pH 7 (M atm⁻¹)
- ``reactivity``          – biological reactivity factor F0 (dimensionless, 0–1)

Values are stored as strings because the MUSICA C++ binding requires
``Mapping[str, str]``; recover floats with ``float(sp.other_properties[key])``.
"""

import musica.mechanism_configuration as mc

hg0 = mc.Species(
    name="Hg0",
    molecular_weight_kg_mol=201e-3,
    other_properties={"henrys_law_constant": "0.11", "reactivity": "3e-5"},
)
so2 = mc.Species(
    name="SO2",
    molecular_weight_kg_mol=64e-3,
    other_properties={"henrys_law_constant": "1e5", "reactivity": "0.0"},
)
o3 = mc.Species(
    name="O3",
    molecular_weight_kg_mol=48e-3,
    other_properties={"henrys_law_constant": "1e-2", "reactivity": "1.0"},
)
