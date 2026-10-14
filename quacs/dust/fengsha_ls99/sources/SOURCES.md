# fengsha_ls99: pinned Fortran sources

Upstream repository: https://github.com/USEPA/CMAQ
Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
Tag used/verified: None. Pinned directly to commit
Retrieved: 2026-10-11
Pin basis: commit named in the Python docstring
Status: MATCHED

## Pin selection

The docstring gives short commit 9bd3734, which the GitHub commit API resolves to the full SHA above. The docstring's CMAQv5.5 names a release family. The explicit commit sets this pin, and no release tag replaces it. The soil-remapping line check below gives extra support only.

## Cited line checks

| Upstream file | Lines checked | Expected source content | Result |
| --- | --- | --- | --- |
| CCTM/src/emis/emis/DUST_EMIS.F | 824-843 | soil remapping (informational) | MATCH |

## Files

### CCTM__src__emis__emis__DUST_EMIS.F

- Upstream repo: USEPA/CMAQ
- Upstream path: CCTM/src/emis/emis/DUST_EMIS.F
- Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
- Tag: None used
- Raw URL: https://raw.githubusercontent.com/USEPA/CMAQ/9bd3734176479c2e49139fea98e1d5e8a16170e3/CCTM/src/emis/emis/DUST_EMIS.F
- Byte size: 56253
- sha256sum: 45f830890c25d505e1880f49d4e06f26020c7a12d66e7d92c35b5f041d174110
- Upstream Git blob SHA: ab220f3c71527f3bdfe338780774424f3ddb33e7
- Pin basis: commit named in the Python docstring

### CCTM__src__emis__emis__lus_data_module.F

- Upstream repo: USEPA/CMAQ
- Upstream path: CCTM/src/emis/emis/lus_data_module.F
- Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
- Tag: None used
- Raw URL: https://raw.githubusercontent.com/USEPA/CMAQ/9bd3734176479c2e49139fea98e1d5e8a16170e3/CCTM/src/emis/emis/lus_data_module.F
- Byte size: 22706
- sha256sum: 0c23f6b9afe9fa21604dd9a0ff312e2ed2a8e9a5550b642b5aaaf9ba36f15be2
- Upstream Git blob SHA: 82b689a3b0ecebf34118e80b234c4db13ce5f242
- Pin basis: commit named in the Python docstring

## Paths and retrieval

All docstring file paths are correct at the selected commit. No requested files are missing.

The files come from the SHA-pinned raw URLs. All bytes, including original line endings, match upstream Git blob hashes and sizes. The retrieval changed no file and ran no code.

### CCTM__src__ICL__fixed__const__CONST.EXT

`dust_hflux` gets `grav` from this file through `include SUBST_CONST`.

- Upstream repo: USEPA/CMAQ
- Upstream path: CCTM/src/ICL/fixed/const/CONST.EXT
- Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
- Raw URL: https://raw.githubusercontent.com/USEPA/CMAQ/9bd3734176479c2e49139fea98e1d5e8a16170e3/CCTM/src/ICL/fixed/const/CONST.EXT
- Byte size: 5645
- sha256sum: 8a4de000b45b198a3eebd6791fc37e18ce2ddbfe532111f1dae9e846df3b1b98
- Git blob SHA (from `git hash-object`): 7b23f51e621c82fd91767f5eda65e973eaee8255
- Retrieved with `curl` from the SHA-pinned raw URL.

### CCTM__src__vdiff__acm2_m3dry__ASX_DATA_MOD.F

`get_dust_emis` uses `karman` at DUST_EMIS.F L972 through `use asx_data_mod`. This file defines `karman = 0.40` at L219. CMAQ has a second variant, `CCTM/src/vdiff/acm2_stage/ASX_DATA_MOD.F`, with the same value at L217. The default build script (`CCTM/scripts/bldit_cctm.csh` L113, `DepMod = m3dry`) selects the `acm2_m3dry` variant.

- Upstream repo: USEPA/CMAQ
- Upstream path: CCTM/src/vdiff/acm2_m3dry/ASX_DATA_MOD.F
- Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
- Raw URL: https://raw.githubusercontent.com/USEPA/CMAQ/9bd3734176479c2e49139fea98e1d5e8a16170e3/CCTM/src/vdiff/acm2_m3dry/ASX_DATA_MOD.F
- Byte size: 73051
- sha256sum: 358e6e19061eba65821f5d71a73866fa26c868732edc8653e54d157cfd019905
- Git blob SHA: 916b15f57f0f49a4fdd2d43145dcfc2934d84d1a (matches the GitHub tree API: 916b15f57f0f49a4fdd2d43145dcfc2934d84d1a)
- Retrieved with `curl` from the SHA-pinned raw URL.

### CCTM__src__emis__emis__LUS_DEFN.F

This file shows how CMAQ makes `ladut` from the land-use fractions. It scales `ladut` by 100 at L98 and merges land-use groups after that.

- Upstream repo: USEPA/CMAQ
- Upstream path: CCTM/src/emis/emis/LUS_DEFN.F
- Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
- Raw URL: https://raw.githubusercontent.com/USEPA/CMAQ/9bd3734176479c2e49139fea98e1d5e8a16170e3/CCTM/src/emis/emis/LUS_DEFN.F
- Byte size: 13904
- sha256sum: 3e994a9b3a80e69bbaa517d33b251c97142d40c51f3029bacc8548274e9fb4b8
- Git blob SHA: 996486531963c6020b67a3bdc12487fab8b44984 (matches the GitHub tree API)
- Retrieved with `curl` from the SHA-pinned raw URL.

### CCTM__src__cio__centralized_io_module.F

This file fills `ladut` from `LUFRAC` (L3363, L3423) or reads it by name (L3372-3378). It copies `dmsk` (L3256-3329), but no CMAQ code applies `dmsk`.

- Upstream repo: USEPA/CMAQ
- Upstream path: CCTM/src/cio/centralized_io_module.F
- Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
- Raw URL: https://raw.githubusercontent.com/USEPA/CMAQ/9bd3734176479c2e49139fea98e1d5e8a16170e3/CCTM/src/cio/centralized_io_module.F
- Byte size: 300675
- sha256sum: a91a47909a7ba20868d1e9345736431a174dfc57a6633fd6979530f4ee2eb4c9
- Git blob SHA: 53a57373dc5b725dbc5f8f351604cd6e3da00d95 (matches the GitHub tree API: 53a57373dc5b725dbc5f8f351604cd6e3da00d95)
- Retrieved with `curl` from the SHA-pinned raw URL.

### CCTM__src__aero__aero6__AERO_DATA.F

`get_dust_emis` takes `gpkg` (L737) and the dust densities from this module. `pi` is not from this file: `asx_data_mod` includes `SUBST_CONST` at module scope (ASX_DATA_MOD.F:95), so `pi = 3.14159265` (CONST.EXT:48), and `f6dpi = 6.0/pi` (DUST_EMIS.F:605). At this commit, `CCTM/src/aero/aero7/AERO_DATA.F` does not exist (HTTP 404), so aero6 is the only variant.

- Upstream repo: USEPA/CMAQ
- Upstream path: CCTM/src/aero/aero6/AERO_DATA.F
- Commit SHA: 9bd3734176479c2e49139fea98e1d5e8a16170e3
- Raw URL: https://raw.githubusercontent.com/USEPA/CMAQ/9bd3734176479c2e49139fea98e1d5e8a16170e3/CCTM/src/aero/aero6/AERO_DATA.F
- Byte size: 117354
- sha256sum: 2f631fe65651c06c78825748eee62c7f4a0d970becd1efebbfd8cc606167e424
- Git blob SHA: 16d7d41276ded77ae86a3c4371f8e6900a69cc4c (from `git hash-object`, not checked against the GitHub tree API)
- Retrieved with `curl` from the SHA-pinned raw URL.
