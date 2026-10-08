"""Tagged tracer prototype: CO tagged by source cell, with MUSICA and 1-D transport.

See README.md for the inputs and outputs.
"""

from .tagged_tracer_core import (
    as_1d_array,
    build_1d_grid,
    build_co_tracer_mechanism,
    build_musica_rate_parameters,
    build_single_cell_tracer_emissions,
    build_total_co_emissions,
    convert_surface_flux_to_concentration_tendency,
    derive_simple_met,
    history_to_dataframe,
    run_source_receptor,
    source_receptor_proof_table,
    transport_1d_advection_diffusion,
)

__all__ = [
    "as_1d_array",
    "build_1d_grid",
    "build_co_tracer_mechanism",
    "build_musica_rate_parameters",
    "build_single_cell_tracer_emissions",
    "build_total_co_emissions",
    "convert_surface_flux_to_concentration_tendency",
    "derive_simple_met",
    "history_to_dataframe",
    "run_source_receptor",
    "source_receptor_proof_table",
    "transport_1d_advection_diffusion",
]
