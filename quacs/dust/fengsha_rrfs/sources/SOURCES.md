# fengsha_rrfs: pinned Fortran sources

Upstream repository: https://github.com/ufs-community/ccpp-physics
Commit SHA: 242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7
Tag used/verified: None. Pinned directly to commit
Retrieved: 2026-10-11
Pin basis: newest commit of the default branch
Status: MATCHED

## Pin selection

The docstring gives no fixed revision. The upstream default branch ufs/dev resolves to the commit above. All seven cited routine ranges have matching start and end definitions and the named source content. The newest change to dust_fengsha_mod.F90 is commit 02b3440378e5bb04ddc4ba50b44f2532eb7cab08 (2024-01-19).

The pin matches the cited routine ranges. It does not match every individual line annotation. The active horizontal-flux q assignment is at L485 (docstring L484). RHOSOIL is at L317 (docstring L318).

The values rhow and rhop are at L514 and L515 (docstring L513 and L514). Tracer remapping is at L369-L370 (docstring L370-L371). These are source-location observations only.

## Cited line checks

| Upstream file | Lines checked | Expected source content | Result |
| --- | --- | --- | --- |
| physics/smoke_dust/dust_fengsha_mod.F90 | 23-210 | gocart_dust_fengsha_driver definition through end | MATCH |
| physics/smoke_dust/dust_fengsha_mod.F90 | 213-372 | source_dust definition through end | MATCH |
| physics/smoke_dust/dust_fengsha_mod.F90 | 398-492 | DustEmissionFENGSHA definition through end | MATCH |
| physics/smoke_dust/dust_fengsha_mod.F90 | 494-527 | soilMoistureConvertVol2Grav definition through end | MATCH |
| physics/smoke_dust/dust_fengsha_mod.F90 | 529-563 | moistureCorrectionFecan definition through end | MATCH |
| physics/smoke_dust/dust_fengsha_mod.F90 | 565-594 | moistureCorrectionShao definition through end | MATCH |
| physics/smoke_dust/dust_fengsha_mod.F90 | 596-626 | DustFluxV2HRatioMB95 definition through end | MATCH |

## Files

### physics__smoke_dust__dust_fengsha_mod.F90

- Upstream repo: ufs-community/ccpp-physics
- Upstream path: physics/smoke_dust/dust_fengsha_mod.F90
- Commit SHA: 242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7
- Tag: None used
- Raw URL: https://raw.githubusercontent.com/ufs-community/ccpp-physics/242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7/physics/smoke_dust/dust_fengsha_mod.F90
- Byte size: 24255
- sha256sum: b387e6803846878e512aafbc0c0091e3759b703774a7f0ea5aa7c97e9b88199b
- Upstream Git blob SHA: 6ec8f8d4af6204d800e3db0f49ca8e90f0a6517a
- Pin basis: newest commit of the default branch

### physics__smoke_dust__dust_data_mod.F90

- Upstream repo: ufs-community/ccpp-physics
- Upstream path: physics/smoke_dust/dust_data_mod.F90
- Commit SHA: 242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7
- Tag: None used
- Raw URL: https://raw.githubusercontent.com/ufs-community/ccpp-physics/242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7/physics/smoke_dust/dust_data_mod.F90
- Byte size: 4245
- sha256sum: 7178cdfe1041bd367150bfe55b0985130c5e0636a7556dbc3178c24a304c3240
- Upstream Git blob SHA: eb809378d7cae204c9595257a1d86c06fdb050da
- Pin basis: newest commit of the default branch

## Paths and retrieval

All docstring file paths are correct at the selected commit. No requested files are missing.

The files come from the SHA-pinned raw URLs. All bytes, including original line endings, match upstream Git blob hashes and sizes. The retrieval changed no file. The comparison with the Python port uses these files compiled with gfortran from unchanged source.

## Dependency files

These files give four items. `rrfs_smoke_config.F90` gives the namelist defaults. `rrfs_smoke_wrapper.F90` (L419-421) and its `.meta` file (L553) give the call site and the gravity argument. `machine.F` gives `kind_phys`. `physcons.F90` (L55) gives the CCPP gravity `con_g = 9.80665`.

The CCPP host model sets `gravitational_acceleration`. `physcons.F90` is the standard CCPP source of that value, but the host mapping lives outside the pinned repository. Each file comes from the SHA-pinned raw URL via `curl`. Each Git blob SHA matches the GitHub tree API at commit `242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7`.

| Upstream path | Bytes | sha256 |
|---|---|---|
| `physics/smoke_dust/rrfs_smoke_config.F90` | 2616 | `541dd8c8903e088a8fa276dc4d2f9a8bb7838db80bc353369a6341437d32fb6c` |
| `physics/smoke_dust/rrfs_smoke_wrapper.F90` | 49465 | `9cffafd39af0897e0fd7f56d4277a6417e0fdcd4f97cb754694568c64e219a7d` |
| `physics/smoke_dust/rrfs_smoke_wrapper.meta` | 26765 | `5dac60ae712e2702ed5336c2931d2c7c355db9150e9529dc93306f1a1f8020f0` |
| `physics/hooks/machine.F` | 1455 | `bad679d0dd90e55265a0da8bfa4424187756c25acc06861de4af31dfc232e568` |
| `physics/hooks/physcons.F90` | 13335 | `9cca4be0183fd61dfe32128a27f3fb5876076f568ba3a0d990f124aa474128b6` |

## Older commit ab84c01a

Full commit SHA: ab84c01a110bbf13490fcf3243effe06cb26608d
Commit date: 2023-09-20T01:29:56Z (author and committer date)
Message line: "GF radar reflectivity, dust bug fix, and C3 updates"

The folder `older_ab84c01a/` holds both files. They come from the SHA-pinned raw URLs. Each file is the same as the upstream blob.

Raw URLs:
- https://raw.githubusercontent.com/ufs-community/ccpp-physics/ab84c01a110bbf13490fcf3243effe06cb26608d/physics/smoke_dust/dust_fengsha_mod.F90
- https://raw.githubusercontent.com/ufs-community/ccpp-physics/ab84c01a110bbf13490fcf3243effe06cb26608d/physics/smoke_dust/dust_data_mod.F90

| Path | Bytes | sha256 | Blob sha (API) | git hash-object | MATCH |
| --- | --- | --- | --- | --- | --- |
| `older_ab84c01a/physics__smoke_dust__dust_fengsha_mod.F90` | 24299 | 71c8fef4446ddad88c364c9d6f383300e8be9506e1cfed0cc1615aa312f7d8b3 | 1e24c89476b5b9c92f24d8211c52ab472eace7f9 | 1e24c89476b5b9c92f24d8211c52ab472eace7f9 | MATCH |
| `older_ab84c01a/physics__smoke_dust__dust_data_mod.F90` | 4245 | 7178cdfe1041bd367150bfe55b0985130c5e0636a7556dbc3178c24a304c3240 | eb809378d7cae204c9595257a1d86c06fdb050da | eb809378d7cae204c9595257a1d86c06fdb050da | MATCH |

Both paths exist at ab84c01a. No file was missing.

## UFS build flags

`ccpp-physics` at the pinned commit has no `CMakeLists.txt` and no CMake files. The GitHub tree API for `242a2fc0fce0fe71772e2b1a7dc76408fc96c6a7` lists only `.github`, `physics`, and `tools`. The other repository `ufs-community/ufs-weather-model` sets the compiler flags for UFS. The folder `ufs_weather_model_cd0c04c5/` holds two files from that repository.

- ufs-weather-model commit: `cd0c04c54dd879851bd2ec08e1d9113cf6b4f8ac` (committer date 2026-09-24T17:59:21Z). It is the newest commit of that repository at or before the committer date of the ccpp-physics pin (2026-10-01T13:41:52Z).
- This commit is not the one that a ufs-weather-model release links to the ccpp-physics pin. The link between the two commits is not verified.
- `cmake/GNU.cmake` was last changed at `de0da94e976b5c276391ad5dca047aacca62c18e` (2026-03-12). The file fetched at `cd0c04c5` is byte-identical to the file fetched at `de0da94e`.
- Each file comes from the SHA-pinned raw URL via `curl`. Each git blob SHA from `git hash-object` matches the GitHub contents API at that commit.

Raw URLs:
- https://raw.githubusercontent.com/ufs-community/ufs-weather-model/cd0c04c54dd879851bd2ec08e1d9113cf6b4f8ac/cmake/GNU.cmake
- https://raw.githubusercontent.com/ufs-community/ufs-weather-model/cd0c04c54dd879851bd2ec08e1d9113cf6b4f8ac/CMakeLists.txt

| Saved as | Upstream path | Bytes | sha256 | Git blob SHA (API) | git hash-object | MATCH |
| --- | --- | --- | --- | --- | --- | --- |
| `ufs_weather_model_cd0c04c5/cmake__GNU.cmake` | `cmake/GNU.cmake` | 1259 | eda95d9d9fe328f4b4fdb511d2e476b02b277f26c4e760ff7acc1220781fc73b | 29ccc20fab02d8d193a71435cd1ff1ecb280ce41 | 29ccc20fab02d8d193a71435cd1ff1ecb280ce41 | MATCH |
| `ufs_weather_model_cd0c04c5/CMakeLists.txt` | `CMakeLists.txt` | 16091 | 943a45b6fff0bc62f7e04ac533d2207fcaa17d7d937c732923305eaa6b797a3e | 7be7e0862528c4e34685d9ecf98aa7088b14d023 | 7be7e0862528c4e34685d9ecf98aa7088b14d023 | MATCH |

What the files show:
- `cmake/GNU.cmake` L8-10: `if(NOT 32BIT)` adds `-fdefault-real-8 -fdefault-double-8` to `CMAKE_Fortran_FLAGS`.
- `CMakeLists.txt` L70: `set(32BIT OFF CACHE BOOL ...)`. L71: `set(CCPP_32BIT OFF CACHE BOOL ...)`.
- Not shown: how `CCPP_32BIT` reaches the ccpp-physics build (for example a `SINGLE_PREC` define). That code is not in the pinned files.

## Host values of alpha and gamma

The wrapper copies `dust_alpha_in` and `dust_gamma_in` from the host at init (W:78-79). The host files below give these values. Each file comes from the SHA-pinned raw URL via `curl`. Each git blob SHA from `git hash-object` matches the GitHub contents API at that commit.

- ufsatm commit `6f461419f091c109d18b826ab88595da00ab336c` is the `UFSATM` submodule of ufs-weather-model `cd0c04c5` (contents API of the root tree at `cd0c04c5`).
- That ufsatm commit links to ccpp-physics `3e6660c6df54e95a0871e990c2294dd397ae3860`. The 7 ccpp-physics files pinned above are byte-identical at `3e6660c6` and `242a2fc` (`cmp`, 2026-10-12).

Raw URLs:
- https://raw.githubusercontent.com/NOAA-EMC/ufsatm/6f461419f091c109d18b826ab88595da00ab336c/ccpp/data/GFS_typedefs.F90
- https://raw.githubusercontent.com/ufs-community/ufs-weather-model/cd0c04c54dd879851bd2ec08e1d9113cf6b4f8ac/tests/parm/rrfs_conus13km_hrrr.nml.IN

| Saved as | Upstream path | Bytes | sha256 | Git blob SHA (API) | git hash-object | MATCH |
| --- | --- | --- | --- | --- | --- | --- |
| `ufsatm_6f461419/ccpp__data__GFS_typedefs.F90` | `ccpp/data/GFS_typedefs.F90` | 475433 | 484179cf21795d16cd33210dc68ada657447b787f8860c7054d023b8f027505c | 5318d5b3b4e314e119fd39ef9e7005edf0b7d2c5 | 5318d5b3b4e314e119fd39ef9e7005edf0b7d2c5 | MATCH |
| `ufs_weather_model_cd0c04c5/tests__parm__rrfs_conus13km_hrrr.nml.IN` | `tests/parm/rrfs_conus13km_hrrr.nml.IN` | 8544 | 93ecdc54a514c635020193924e6d0f2c0ab087f6d62e7d51a04c468c948ed065 | cd86be418377309cc71b68c1e0c8b1f6ed255434 | cd86be418377309cc71b68c1e0c8b1f6ed255434 | MATCH |

What the files show:
- `GFS_typedefs.F90:4163-4166,4169`: host defaults `dust_drylimit_factor = 1.0`, `dust_moist_correction = 1.0`, `dust_alpha = 0.`, `dust_gamma = 0.`, `dust_moist_opt = 1`. L4707-4711 copy them to `Model%`.
- `rrfs_conus13km_hrrr.nml.IN:159-160`: `dust_alpha = 0.01`, `dust_gamma = 1.3`. Of the test namelists at `cd0c04c5`, only this one sets them.
- Not shown: the values that operational RRFS uses. The RRFS workflow repository is not pinned.

## Real RRFS forecast namelist

The UFS test namelist above sets `dust_opt = 5`, so it does not run this dust code (W:418 needs `dust_opt == 1`). The RRFS workflow repository holds the namelist of real forecasts.

- Repository: NOAA-EMC/rrfs-workflow, default branch `release/v1.0`.
- Commit `d228600ef9f9174b1cbdff99a889078afdae3e14` (committer date 2026-10-07T20:33:25Z).
- Raw URL: https://raw.githubusercontent.com/NOAA-EMC/rrfs-workflow/d228600ef9f9174b1cbdff99a889078afdae3e14/parm/config/det/input.nml
- The same dust lines appear on the branches `main`, `rrfs-dev`, `rrfs-nco` and `dev-sci`, and in about 20 namelist variants.

| Saved as | Upstream path | Bytes | sha256 | Git blob SHA (API) | git hash-object | MATCH |
| --- | --- | --- | --- | --- | --- | --- |
| `rrfs_workflow_d228600/parm__config__det__input.nml` | `parm/config/det/input.nml` | 6778 | 12542e4d79694643c52c3df7670e249a30a352b1f74e42e4c157ab41b2ab6586 | 084199056d025ef7f5e0099b043a5e4d26d672d0 | 084199056d025ef7f5e0099b043a5e4d26d672d0 | MATCH |

What the file shows (`&gfs_physics_nml` starts at L155):
- L190 `dust_alpha = 10.0`, L191 `dust_drylimit_factor = 0.5`, L192 `dust_gamma = 1.3`, L193 `dust_moist_correction = 2.0`, L194 `dust_opt = 1`.
- L272 `rrfs_sd = .true.`
- `dust_moist_opt` is not set, so it keeps the host default 1 (`GFS_typedefs.F90:4169`).
- Not shown: which ccpp-physics commit the RRFS v1.0 release uses.

## Code version of real RRFS forecasts

The chain below links the real forecast workflow to its code version.

- rrfs-workflow `d228600` (`release/v1.0`), `sorc/Externals.cfg` L11-18: ufs-weather-model hash `d150b9a33a12ea8a783dc845a0d036e1f6f9d733` (branch `production/RRFS.v1`).
- ufs-weather-model `d150b9a`: submodule `FV3` is NOAA-EMC/fv3atm `957189f834c3b519c6d3f9445e486f8d4613da9e`.
- fv3atm `957189f`: submodule `ccpp/physics` is ufs-community/ccpp-physics `edd85a55502179611fdb7b902f0639c305cd1c3b` (2026-08-04).

Comparison of `edd85a5` with the pin `242a2fc` (raw download and `cmp`):
- `dust_fengsha_mod.F90`, `dust_data_mod.F90`, `rrfs_smoke_config.F90`: byte-identical.
- `rrfs_smoke_wrapper.F90`: differs only in arguments for fire smoke, optional flags, the MPI type and the no-fire land classes. No changed line touches a dust input, the dust call or the dust settings.
- `physcons.F90`: the pin adds constants for other schemes. `con_g = 9.80665` and `con_rd = 287.05` are the same in both.
- So real RRFS forecasts run the same dust code as the pin.
