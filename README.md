# quacs-parameterization

This repository collects the target parameterizations for the QUACS project.
QUACS adds a scripting capability for chemistry and aerosol parameterizations to MUSICA and MPAS-A.
The parameterizations here guide the design of that capability.

Each parameterization has a README that tells its inputs, its outputs, and where each one comes from or goes to.
Most parameterizations also run in a MUSICA box model.

## Layout

Each parameterization is a scheme folder inside a process folder: `quacs/<process>/<scheme>/`.
See [CONTRIBUTING.md](CONTRIBUTING.md) for the full layout and the steps to add a scheme.

```
quacs/
└── <process>/
    ├── README.md          # table of the schemes in this process
    └── <scheme>/
        ├── README.md      # status, inputs, outputs, references
        ├── <scheme>.py    # the physics
        ├── examples/      # MUSICA box model scripts and notebooks
        └── tests/
```

## Target parameterizations

Status is one of: **implemented**, **in progress**, or **placeholder**.
A placeholder folder has only a README with the content of the target list.

| Target class | Process folder | Scheme | Status |
|---|---|---|---|
| Dust emission | [`dust`](quacs/dust/) | GOCART (Ginoux et al., 2001) | implemented ([PR #5](https://github.com/NCAR/quacs-parameterization/pull/5)) |
| Dust emission | [`dust`](quacs/dust/) | [LS99-FENGSHA](quacs/dust/fengsha/) | placeholder |
| Dust emission | [`dust`](quacs/dust/) | Kok et al. (2014) | implemented (planned PR) |
| Dust emission | [`dust`](quacs/dust/) | DEAD (Zender et al., 2003) | implemented (planned PR) |
| Dust emission | [`dust`](quacs/dust/) | Leung et al. (2023) | implemented (planned PR) |
| Dry deposition | [`drydep`](quacs/drydep/) | [Wesely (1989), from CAM](quacs/drydep/wesely/) | implemented |
| Dry deposition | [`drydep`](quacs/drydep/) | [GEOS-Chem, simple](quacs/drydep/simple/) | implemented |
| Dry deposition | [`drydep`](quacs/drydep/) | [GEOS-Chem, LHS](quacs/drydep/lhs/) | implemented |
| Aerosol wet deposition | [`wetdep`](quacs/wetdep/) | CAM wetdepa_v2 | implemented (planned PR) |
| Aerosol wet deposition | [`wetdep`](quacs/wetdep/) | [PNNL MOSAIC](quacs/wetdep/aerosol_mosaic/) | placeholder |
| Trace gas wet deposition | [`wetdep`](quacs/wetdep/) | [Neu et al. (2012)](quacs/wetdep/gas_neu/) | placeholder |
| MEGAN biogenic emissions | [`biogenic`](quacs/biogenic/) | MEGAN 3.0 | implemented ([PR #10](https://github.com/NCAR/quacs-parameterization/pull/10)) |
| pyroCb-aware plume rise | [`plumerise`](quacs/plumerise/) | Freitas et al. (2010) with pyroCb trigger | in progress ([PR #11](https://github.com/NCAR/quacs-parameterization/pull/11)) |
| VOC and NOx tagged tracers | [`tagged_tracers`](quacs/tagged_tracers/) | Lapaşcu & Butler (2019) | in progress (planned PR) |
| Lightning flash rate, parameterized convection | [`lightning`](quacs/lightning/) | [Price & Rind (1992)](quacs/lightning/flash_rate_price_rind/) | placeholder |
| Lightning flash rate, resolved convection | [`lightning`](quacs/lightning/) | [Cummings et al. (2024)](quacs/lightning/flash_rate_cummings/) | placeholder |
| Lightning NOx, parameterized convection | [`lightning`](quacs/lightning/) | [Ott et al. (2010)](quacs/lightning/lnox_ott/) | placeholder |
| Lightning NOx, resolved convection | [`lightning`](quacs/lightning/) | [DeCaria et al. (2000)](quacs/lightning/lnox_decaria/) | placeholder |
| Aerosol activation | [`aerosol_activation`](quacs/aerosol_activation/) | [GOCART-2G to Thompson](quacs/aerosol_activation/gocart2g_thompson/) | placeholder |
| New particle formation | [`nucleation`](quacs/nucleation/) | not selected | placeholder |
| Aerosol optics | [`aerosol_optics`](quacs/aerosol_optics/) | [Ghan & Zaveri (2007)](quacs/aerosol_optics/ghan_zaveri/) | placeholder |
| Aerosol and snow albedo | [`snow_albedo`](quacs/snow_albedo/) | [SNICAR-ADv3](quacs/snow_albedo/snicar_adv/) | placeholder |
| Heterogeneous chemistry | [`heterogeneous_chem`](quacs/heterogeneous_chem/) | [Uptake on dust](quacs/heterogeneous_chem/dust_uptake/) | placeholder |
| Halogenated VSLS | [`vsls`](quacs/vsls/) | CESM2-SLH sea-salt bromine | in progress (planned PR) |

## Installation

```bash
pip install -e ".[dev]"
```

## Running tests

```bash
pytest
```

This command runs the unit tests, each example script, each notebook, and the README checks.

## Style

[ruff](https://docs.astral.sh/ruff/) formats and lints the code:

```bash
ruff check .
ruff format .
```
