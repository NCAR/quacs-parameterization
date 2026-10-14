# WRF-Chem AFWA source provenance

Repository: https://github.com/wrf-model/WRF
Latest release: https://github.com/wrf-model/WRF/releases/tag/v4.8.0
Tag: v4.8.0
Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
Retrieved: 2026-10-11

## chem__module_gocart_dust_afwa.F

- Upstream path: chem/module_gocart_dust_afwa.F
- Tag: v4.8.0
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/module_gocart_dust_afwa.F
- Byte size: 27733
- sha256sum: 54f9fd0655848cc1f1a18b8abb164733edd56b8a3c539c63281d5f0e20c672fb
- Upstream Git blob SHA: 02ad19000d72818b639b04c631b5839162728002

## phys__module_data_gocart_dust.F

- Upstream path: phys/module_data_gocart_dust.F
- Tag: v4.8.0
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/phys/module_data_gocart_dust.F
- Byte size: 2169
- sha256sum: 8694100875645b54a58bcff38d462ad3c9fd115261fb39941358466c84b62329
- Upstream Git blob SHA: 24ebf6ed4eae9b82b0c9271aa46b9cde57ce378d

## chem__emissions_driver.F

- Upstream path: chem/emissions_driver.F
- Tag: v4.8.0
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/emissions_driver.F
- Byte size: 100045
- sha256sum: 8abc108df7b06bc2eba7b8d9a07f00091b5d8ab66d1d00372105139c7d409f9b
- Upstream Git blob SHA: 2394518d1ae6344f4619c378e1d48a13e1a7cd4a

## Registry__registry.chem

- Upstream path: Registry/registry.chem
- Tag: v4.8.0
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/Registry/registry.chem
- Byte size: 627854
- sha256sum: 6f614cccae4c0a151e707120d45d47f168f2427c6cae1853383470ee9216f312
- Upstream Git blob SHA: cac0a0db9e0308d8fe021c63517d992700e743f5

## Location checks

The Python docstring identifies source_dust and the dust bin/soil parameters. Both cited module paths are correct at this release, including phys/module_data_gocart_dust.F (not chem/module_data_gocart_dust.F).

chem/emissions_driver.F calls gocart_dust_afwa_driver at line 822. That routine calls source_dust in chem/module_gocart_dust_afwa.F at line 289. The candidate chem/dust_emissions_driver.F does not exist at this tag.

Registry/registry.chem lines 3902-3905 define dust_alpha, dust_gamma, dust_smtune, and dust_ustune, each with default 1. The folder holds the full Registry file.

All four source files match the upstream Git blob hashes and byte sizes. They keep the original line endings, including CRLF in phys/module_data_gocart_dust.F. The retrieval changed no file and ran no code. No requested source is missing.

## Registry__Registry.EM_COMMON

- Upstream path: Registry/Registry.EM_COMMON
- Tag: v4.8.0
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/Registry/Registry.EM_COMMON
- Retrieved: 2026-10-12
- Byte size: 447730
- sha256sum: 78aa665d0488da08b7441b8a7802c5c774cb21b24f64e99b8ef6cf805c84b0eb
- Why: it gives the `sf_surface_physics` numbers of the land-surface models that `dust_smois = 1` uses. `lsmscheme` is 2 (L3210), `ruclsmscheme` is 3 (L3211), `pxlsmscheme` is 7 (L3215).
