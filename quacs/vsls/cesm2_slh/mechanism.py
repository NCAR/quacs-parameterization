"""
MICM mechanism for the CESM2-SLH sea-salt bromine reactions.

Each reaction is a MICM ``UserDefined`` reaction whose rate constant is
supplied at run time from :func:`quacs.vsls.cesm2_slh.seasalt.seasalt_rate_constants`
via the state parameter ``USER.<reaction name>``. This follows the QUACS
pattern used in ``quacs.drydep``: Python computes the rate, MICM integrates.
"""

from typing import Dict, Tuple

import musica.mechanism_configuration as mc

from .seasalt import SEA_SALT_REACTIONS

MOLAR_MASS_KG_MOL = {
    "HOBr": 96.91e-3,
    "BrONO2": 141.91e-3,
    "BrNO2": 125.91e-3,
    "Br2": 159.81e-3,
    "BrCl": 115.36e-3,
}


def _species_first_tuples() -> bool:
    """Detect the (species, coefficient) tuple order accepted by this musica.

    Older musica releases expect ``(coefficient, Species)``; newer ones
    (e.g. 0.16.7) expect ``(Species, coefficient)``. Probe once instead of
    pinning a version.
    """
    probe = mc.Species(name="_probe")
    try:
        mc.UserDefined(name="_probe", products=[(probe, 1.0)])
        return True
    except (AttributeError, TypeError):
        return False


_SPECIES_FIRST = _species_first_tuples()


def _component(species: mc.Species, coefficient: float):
    """Product/reactant entry with a stoichiometric coefficient."""
    return (species, coefficient) if _SPECIES_FIRST else (coefficient, species)


def rate_parameter_key(reaction_name: str) -> str:
    """MICM state key for a UserDefined reaction's rate constant."""
    return f"USER.{reaction_name}"


def build_mechanism(
    name: str = "vsls_seasalt_bromine",
) -> Tuple[mc.Mechanism, Dict[str, mc.Species]]:
    """Build the sea-salt bromine mechanism.

    Returns
    -------
    mechanism : mc.Mechanism
    species : dict
        Species name -> ``mc.Species``.
    """
    species = {
        n: mc.Species(name=n, molecular_weight_kg_mol=mw) for n, mw in MOLAR_MASS_KG_MOL.items()
    }
    gas = mc.Phase(name="gas", species=list(species.values()))
    reactions = [
        mc.UserDefined(
            name=rxn.name,
            scaling_factor=1.0,
            reactants=[species[rxn.reactant]],
            products=[_component(species[p], y) for p, y in rxn.products],
            gas_phase=gas,
        )
        for rxn in SEA_SALT_REACTIONS
    ]
    mechanism = mc.Mechanism(
        name=name, species=list(species.values()), phases=[gas], reactions=reactions
    )
    return mechanism, species
