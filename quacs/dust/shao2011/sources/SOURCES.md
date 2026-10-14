# shao2011: pinned Fortran sources

Upstream repository: https://github.com/wrf-model/WRF
Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
Tag used/verified: v4.8.0
Retrieved: 2026-10-11
Pin basis: newest commit of the default branch
Status: MATCHED

## Pin selection

The docstring says master, which is a moving branch and not a fixed commit or release tag. Current master resolves to the commit above. Its cited lines match, so the pin needs no older commit. The newest change to module_qf03.F is commit 6273be20fb50ae819c5399bec038bf1a8d97d70e (2021-03-24).

Citation limitation: the docstring contains only two distinct line citations: L130 and L317-330. The table splits the latter into L317-323 and L324-330. This gives three checks of two original citations, not three independent citations.

## Cited line checks

| Upstream file | Lines checked | Expected source content | Result |
| --- | --- | --- | --- |
| chem/module_qf03.F | 130 | cys = 0.00001 | MATCH |
| chem/module_qf03.F | 317-323 | imod=3 branch with zeta, sigma_m and the start of the dust loop | MATCH |
| chem/module_qf03.F | 324-330 | cys*dpsdm, (1+sigma_m), qtotal*g/ustar**2, cover factor | MATCH |

## Files

### chem__module_uoc_dust.F

- Upstream repo: wrf-model/WRF
- Upstream path: chem/module_uoc_dust.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/module_uoc_dust.F
- Byte size: 15433
- sha256sum: 98a17e2feb8c81786d7ae8d694b124b50f8e4746cc7edc176305bc1ab4030ab8
- Upstream Git blob SHA: e8961e40eb1f8c85d3e37e3dc9f4664c54b2a375
- Pin basis: newest commit of the default branch

### chem__module_qf03.F

- Upstream repo: wrf-model/WRF
- Upstream path: chem/module_qf03.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/module_qf03.F
- Byte size: 18357
- sha256sum: 4c39512dac10ca3f37570e860e7fb475c883667772bb3f4eb743291e9881a63c
- Upstream Git blob SHA: 7fc24c8150652e0e990e8cf7529d53732d3076f2
- Pin basis: newest commit of the default branch

### chem__module_soilpsd.F

- Upstream repo: wrf-model/WRF
- Upstream path: chem/module_soilpsd.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/module_soilpsd.F
- Byte size: 8952
- sha256sum: 2ea2df2d7baa82af1399c9dca4e77fcd0a8dc788e799023fa5dc3b96c11e2609
- Upstream Git blob SHA: f36f6870dbf3662537d994c0953ee9abfd5eb2c7
- Pin basis: newest commit of the default branch

### chem__emissions_driver.F

- Upstream repo: wrf-model/WRF
- Upstream path: chem/emissions_driver.F
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/chem/emissions_driver.F
- Byte size: 100045
- sha256sum: 8abc108df7b06bc2eba7b8d9a07f00091b5d8ab66d1d00372105139c7d409f9b
- Upstream Git blob SHA: 2394518d1ae6344f4619c378e1d48a13e1a7cd4a
- Pin basis: newest commit of the default branch

### Registry__registry.chem

- Upstream repo: wrf-model/WRF
- Upstream path: Registry/registry.chem
- Commit SHA: 06d4240ae989cc3e50af412bb472df3d9048783c
- Tag: v4.8.0
- Raw URL: https://raw.githubusercontent.com/wrf-model/WRF/06d4240ae989cc3e50af412bb472df3d9048783c/Registry/registry.chem
- Byte size: 627854
- sha256sum: 6f614cccae4c0a151e707120d45d47f168f2427c6cae1853383470ee9216f312
- Upstream Git blob SHA: cac0a0db9e0308d8fe021c63517d992700e743f5
- Pin basis: newest commit of the default branch

## Paths and retrieval

All docstring file paths are correct at the selected commit. No requested files are missing.
The second WRF-Chem driver is chem/emissions_driver.F, which calls uoc_dust_driver at L843. The folder holds Registry/registry.chem in full at the same commit.

The files come from the SHA-pinned raw URLs. All bytes, including original line endings, match upstream Git blob hashes and sizes. The retrieval changed no file and ran no code.

## Land-model files for DRYSMC

The Shao driver reads `DRYSMC` from the land-model modules. Those modules fill it from `run/SOILPARM.TBL` at run time. All 4 files are at the same commit and come from raw.githubusercontent.com.

| File | Upstream path | Bytes | sha256 |
|---|---|---|---|
| `run__SOILPARM.TBL` | `run/SOILPARM.TBL` | 6557 | `1e2275a32d8cd3b48ca693d22c0816df0013f83b6594ac632716361db337d58f` |
| `phys__module_sf_noahlsm.F` | `phys/module_sf_noahlsm.F` | 206078 | `034bf1d5f48b0b734702f200e3bcefe67db441b170cdcec8592992eb7d358e8c` |
| `phys__module_sf_ruclsm.F` | `phys/module_sf_ruclsm.F` | 288509 | `3265f810d08dcbddfaf198371dc7f652e78e8d3a788f703a515c555a3bbb2a12` |
| `share__module_model_constants.F` | `share/module_model_constants.F` | 8280 | `5b80377fecdc18a5f0ad38d3b6c15cfc86ad5d76701adbbbb08a08698d0f7062` |
