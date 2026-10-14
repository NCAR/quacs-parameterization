# ginoux: pinned Fortran sources

Upstream repository: https://github.com/wrf-model/WRF
Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
Tag used/verified: v4.8.0
Retrieved: 2026-10-11
Pin basis: commit named in the Python docstring
Status: MATCHED

## Pin selection

The Python docstring gives commit 06d4240ae989cc3e50af412bb472df3d9048783c. This is also the commit of WRF release v4.8.0. The line checks below give extra support but do not set the pin.

## Cited line checks

| Upstream file | Lines checked | Expected source content | Result |
| --- | --- | --- | --- |
| phys/module_data_gocart_dust.F | 2 | ndust=5, ndcls=3, ndsrc=1 | MATCH |
| phys/module_data_gocart_dust.F | 18-21 | den_dust, reff_dust, ipoint, frac_s | MATCH |
| chem/chemics_init.F | 1694 | ch_dust(:,:) = 0.8D-9 | MATCH |
| chem/chemics_init.F | 1723 | ch_dust(:,:) = 0.8D-9 | MATCH |
| chem/chemics_init.F | 1742 | ch_dust(:,:) = 0.8D-9 | MATCH |
| chem/chemics_init.F | 1842 | ch_dust(:,:) = 0.8D-9 | MATCH |
| chem/chemics_init.F | 2029 | ch_dust(:,:) = 0.8D-9 | MATCH |

## Files

### chem__module_gocart_dust.F

- Upstream repo: wrf-model/WRF
- Upstream path: chem/module_gocart_dust.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/module_gocart_dust.F
- Byte size: 8416
- sha256sum: df481177ab95ca91bf1155abf9690ade1b12c0a329d8a3914b48459c0c496fcb
- Upstream Git blob SHA: 1049227c677ae575164797cd7a6a19ac16aaee90
- Pin basis: commit named in the Python docstring

### phys__module_data_gocart_dust.F

- Upstream repo: wrf-model/WRF
- Upstream path: phys/module_data_gocart_dust.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/phys/module_data_gocart_dust.F
- Byte size: 2169
- sha256sum: 8694100875645b54a58bcff38d462ad3c9fd115261fb39941358466c84b62329
- Upstream Git blob SHA: 24ebf6ed4eae9b82b0c9271aa46b9cde57ce378d
- Pin basis: commit named in the Python docstring

### chem__chemics_init.F

- Upstream repo: wrf-model/WRF
- Upstream path: chem/chemics_init.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/chemics_init.F
- Byte size: 157790
- sha256sum: 6fc1ffbbfba806a86339d48d4523d91246cfb44cfae0583eb18a40c59df449bc
- Upstream Git blob SHA: a4dcda31f297e9c49c01a96d372d2c1a0a13cb22
- Pin basis: commit named in the Python docstring

### chem__emissions_driver.F

- Upstream repo: wrf-model/WRF
- Upstream path: chem/emissions_driver.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/emissions_driver.F
- Byte size: 100045
- sha256sum: 8abc108df7b06bc2eba7b8d9a07f00091b5d8ab66d1d00372105139c7d409f9b
- Upstream Git blob SHA: 2394518d1ae6344f4619c378e1d48a13e1a7cd4a
- Pin basis: commit named in the Python docstring

### Registry__registry.chem

- Upstream repo: wrf-model/WRF
- Upstream path: Registry/registry.chem
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/Registry/registry.chem
- Byte size: 627854
- sha256sum: 6f614cccae4c0a151e707120d45d47f168f2427c6cae1853383470ee9216f312
- Upstream Git blob SHA: cac0a0db9e0308d8fe021c63517d992700e743f5
- Pin basis: commit named in the Python docstring

## Paths and retrieval

All docstring file paths are correct at the selected commit. No requested files are missing.
The second WRF-Chem driver is chem/emissions_driver.F, which calls gocart_dust_driver at L814. The folder holds Registry/registry.chem in full at the same commit.

The files come from the SHA-pinned raw URLs. All bytes, including original line endings, match upstream Git blob hashes and sizes. The retrieval changed no file and ran no code.
