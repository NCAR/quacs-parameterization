# leung2023: pinned Fortran sources

- Upstream archive: Zenodo record 10621844, https://zenodo.org/records/10621844.
- Version: 1.0.
- DOI: 10.5281/zenodo.10621844.
- Archive: dmleung/dustmodelcesm2023-1.0.zip.
- Archive size: 82259 bytes.
- Archive published MD5: a388ad09ff1116ae2514c7f03d7793ba. Verified locally, matches.
- Archive sha256: 4e9a734b332d69abfa6cd30a04d5894cddafc6b6bcac24d2e3b328ccd110da17.
- Archive saved at: sources/zip/z.zip.
- GitHub source: https://github.com/dmleung/dustmodelcesm2023, commit 90771fefa1a36fb942eed6092a18e93896d01a4e.
- Archive root folder: dmleung-dustmodelcesm2023-90771fe/.
- Extracted to: sources/src/. The extraction drops the root folder and keeps member names.
- Retrieved: 2026-10-11.
- Pin basis: Zenodo record 10621844.
- Status: MATCHED. 54 of 54 citations MATCH, 0 OFFSET, 0 MISSING.

## Pin selection

The pin is the Zenodo record 10621844, version 1.0.

Verification run on 2026-10-11:

- The zip MD5 matches the Zenodo published value.
- The zip is 82259 bytes.
- src/DUSTMod.F90 and src/SoilStateInitTimeConstMod.F90 are byte-identical to their archive members. The cmp check passed.
- Their sha256 values match the member hashes.
- The GitHub commit 90771fe comes from the pin, with no separate fetch.

## Cited line checks

Source of citations: the Python scheme `leung2023.py`. The column "Cited at" gives the file and line in this folder (checked 2026-10-13).

Each cited line or range is one row. A Python line with two citations gives two rows.

| # | Cited at | Upstream file | Lines checked | Expected source content | Fortran at cited lines | Result |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | leung2023.py:26 | DUSTMod.F90 | 505-511 | lnd_frc_mbl bare/mobile fraction | 505-511 tlai test sets lnd_frc_mbl, 511 scales by frac_sno | MATCH |
| 2 | leung2023.py:27 | DUSTMod.F90 | 571-577 | Fecan moisture factor | 574 sqrt(1+1.21*(100*(gwc-gwc_thr))^0.68), 576 default 1 | MATCH |
| 3 | leung2023.py:28 | DUSTMod.F90 | 593-596 | Shao and Lu threshold, impact threshold | 593 tmp2 S&L00, 594 dry, 596 impact = B_it*tmp2 | MATCH |
| 4 | leung2023.py:29 | DUSTMod.F90 | 605-612 | standardized threshold, Cd, frag exponent | 605 std threshold, 606 Cd, 609-612 frag exponent with cap 3 | MATCH |
| 5 | leung2023.py:30 | DUSTMod.F90 | 623-641 | Okin/Pierre drag partition on wind | 623 condition, 630-635 K_length, ssr, frc_thr_rgh_fct, 639 else | MATCH |
| 6 | leung2023.py:31 | DUSTMod.F90 | 724 | K14 impact-threshold flux | 724 gate wnd_frc_slt > wnd_frc_thr_slt_it | MATCH |
| 7 | leung2023.py:31 | DUSTMod.F90 | 752 | K14 impact-threshold flux | 752 flx_mss_vrt_dst_ttl using wnd_frc_thr_slt_it | MATCH |
| 8 | leung2023.py:32 | DUSTMod.F90 | 755 | bare fraction, liqfrac, C_tune scaling | 755 multiplies by lnd_frc_mbl, C_tune, liqfrac | MATCH |
| 9 | leung2023.py:33 | DUSTMod.F90 | 778-811 | Comola intermittency factor | 778 u_mean_slt through 811 intrmtncy_fct | MATCH |
| 10 | leung2023.py:35 | SoilStateInitTimeConstMod.F90 | 675 | gwc_thr = 0.01*2*(0.17*clay + 0.0014*clay^2) | 675 same expression, a = 2 | MATCH |
| 11 | leung2023.py:36 | SoilStateInitTimeConstMod.F90 | 676 | mss_frc_cly_vld = min(clay/100, 0.20) | 676 min(clay*0.01, 0.20) | MATCH |
| 12 | leung2023.py:37 | SoilStateInitTimeConstMod.F90 | 329 | roughfct from F_eff | 329 ncd_io read of F_eff into roughfct2d | MATCH |
| 13 | leung2023.py:37 | SoilStateInitTimeConstMod.F90 | 343-345 | rockfrc A_r, vegefrc A_v, sparfrc A_s | 343 A_r, 344 A_v, 345 A_s | MATCH |
| 14 | leung2023.py:62 | DUSTMod.F90 | 596 | impact threshold uses dry threshold | 596 uses tmp2 (dry), not the wet factor | MATCH |
| 15 | leung2023.py:65 | DUSTMod.F90 | 593 | Dp = 130 um in code | 593 uses 130.0e-6 | MATCH |
| 16 | leung2023.py:67 | DUSTMod.F90 | 788 | stability hardcoded to 0 | 788 stblty(p) = 0 | MATCH |
| 17 | leung2023.py:786 | DUSTMod.F90 | 758 | coastal discard, *1e-3 where land frac < 0.5 | 758 condition frac(g) < 0.5. The *1e-3 is on 759 | MATCH (note 1) |
| 18 | leung2023.py:788 | DUSTMod.F90 | 829 | Zender bin partition, applied downstream | 829 loop over ndst, inside flx_mss_vrt_dst_prt block 826-834 | MATCH |
| 19 | leung2023.py:73 | DUSTMod.F90 | 382 | cst_slt = 2.61 declared, dead in this path | 382 declaration. The only other use, at 729, sits in a commented-out block | MATCH |
| 20 | leung2023.py:114 | DUSTMod.F90 | 593 | S&L00 threshold | 593 | MATCH |
| 21 | leung2023.py:118 | DUSTMod.F90 | 393 | DNS_SLT = 2650 | 393 dns_slt = 2650.0 | MATCH |
| 22 | leung2023.py:120 | SoilStateInitTimeConstMod.F90 | 675 | Fecan a = 2 | 675 0.01_r8*2_r8*(...) | MATCH |
| 23 | leung2023.py:123 | DUSTMod.F90 | 387-392 | Kok 2014 coefficients | 387 Cd0, 388 Ca, 389 Ce, 390 C_tune, 391 std min, 392 forc_rho_std | MATCH |
| 24 | leung2023.py:131 | DUSTMod.F90 | 386 | VAI_MBL_THR = 1.0. Comment group label says intermittency/impact | 386 vai_mbl_thr = 1.0 (VAI mobilization threshold) | MATCH (note 2) |
| 25 | leung2023.py:131 | DUSTMod.F90 | 395-399 | B_it, von Karman k, U_0, c_e | 395 B_it, 396 k, 398 U_0, 399 c_e | MATCH |
| 26 | leung2023.py:142-143 | DUSTMod.F90 | 790 | exponent 0.333 literal | 790 (12 - 0.5*stblty)**0.333_r8 | MATCH |
| 27 | leung2023.py:144 | DUSTMod.F90 | 806 | sqrt(2) literal 1.414 | 806 erf(... / (1.414_r8 * u_sd)) | MATCH |
| 28 | leung2023.py:144 | DUSTMod.F90 | 808 | sqrt(2) literal 1.414 | 808 same form for prb_crs_impct_thr | MATCH |
| 29 | leung2023.py:145 | DUSTMod.F90 | 633 | cube root written as **0.3333 | 633 (...)**(0.3333_r8) | MATCH |
| 30 | leung2023.py:271 | DUSTMod.F90 | 593-596 | S&L00 threshold and impact | 593-596 | MATCH |
| 31 | leung2023.py:293 | DUSTMod.F90 | 571-577 | Fecan moisture factor | 571-577 | MATCH |
| 32 | leung2023.py:293 | SoilStateInitTimeConstMod.F90 | 675 | threshold L675. File not named in comment, checked here | 675 gwc_thr_col | MATCH (note 3) |
| 33 | leung2023.py:313 | DUSTMod.F90 | 623-641 | drag partition | 623-641 | MATCH |
| 34 | leung2023.py:336 | DUSTMod.F90 | 778-811 | intermittency factor | 778-811 | MATCH |
| 35 | leung2023.py:342 | DUSTMod.F90 | 788 | stability = 0 | 788 | MATCH |
| 36 | leung2023.py:360 | DUSTMod.F90 | 791-792 | u_sd guard 0.001 | 791 else branch, 792 u_sd = 0.001 | MATCH |
| 37 | leung2023.py:166 | DUSTMod.F90 | 390 | C_tune default 0.05 | 390 C_tune = 0.05 | MATCH |
| 38 | leung2023.py:462 | SoilStateInitTimeConstMod.F90 | 676 | clay capped at 0.20 | 676 | MATCH |
| 39 | leung2023.py:465 | DUSTMod.F90 | 505-511 | bare/mobile fraction | 505-511 | MATCH |
| 40 | leung2023.py:473 | DUSTMod.F90 | 593-596 | thresholds | 593-596 | MATCH |
| 41 | leung2023.py:480 | DUSTMod.F90 | 605-612 | standardized threshold, Cd, frag exponent | 605-612 | MATCH |
| 42 | leung2023.py:484 | DUSTMod.F90 | 610-612 | frag exponent cap at 3 | 610-612 | MATCH |
| 43 | leung2023.py:486 | DUSTMod.F90 | 623-641 | drag partition on wind | 623-641 | MATCH |
| 44 | leung2023.py:489 | DUSTMod.F90 | 623 | Fortran condition for drag partition | 623 lnd_frc_mbl > 0 .AND. tlai_lu <= 1 | MATCH |
| 45 | leung2023.py:494 | DUSTMod.F90 | 724 | impact-threshold gate | 724 | MATCH |
| 46 | leung2023.py:494 | DUSTMod.F90 | 752 | impact-threshold flux | 752 | MATCH |
| 47 | leung2023.py:507 | DUSTMod.F90 | 755 | bare-fraction, tuning, liqfrac scaling | 755 | MATCH |
| 48 | leung2023.py:510 | DUSTMod.F90 | 778-811 | intermittency factor | 778-811 | MATCH |
| 49 | leung2023.py:510 | DUSTMod.F90 | 814-818 | NaN-safe multiply by intermittency | 814 NaN test, 817 multiply | MATCH |
| 50 | tests/test_leung2023.py:530 | DUSTMod.F90 | 593 | Dp = 130 um | 593 | MATCH |
| 51 | tests/test_leung2023.py:535 | SoilStateInitTimeConstMod.F90 | 675 | Fecan a must be 2 | 675 contains 0.01_r8*2_r8 | MATCH |
| 52 | tests/test_leung2023.py:551 | DUSTMod.F90 | 596 | u_it = B_it * dry threshold, moisture independent | 596 | MATCH |
| 53 | tests/test_leung2023.py:648 | DUSTMod.F90 | 790 | u_sd = u_s * 12^0.333 | 790 | MATCH |
| 54 | tests/test_leung2023.py:679 | DUSTMod.F90 | 755 | liqfrac enters linearly | 755 | MATCH |

Result counts: 54 checked. MATCH 54. OFFSET 0. MISSING 0.

Notes:

1. Row 17. The land-fraction test is on L758. The *1e-3 factor is on L759. Both sit in the same if block, so the citation holds.
2. Row 24. L386 holds vai_mbl_thr, the VAI threshold for dust mobilization. The value matches VAI_MBL_THR. The comment group label "Intermittency / impact threshold" does not describe L386. This is a label issue in the Python comment, not a line offset.
3. Row 32. The Python comment at line 172 names "threshold L675" without a file. It matches SoilStateInitTimeConstMod.F90.

Other checks on the docstring (not line numbers):

- Docstring line 25 names "subroutine DustEmission". It exists in DUSTMod.F90 at L320-L851.
- Docstring line 76 says cst_slt "appears only in the commented-out Z03 block". Confirmed. The only uses are L382 (declaration) and L729 (inside a commented block).
- Python line 613 names "DUSTMod.F90" with no line number. Not checked.

## Files

### DUSTMod.F90

- Archive member: dmleung-dustmodelcesm2023-90771fe/DUSTMod.F90.
- Saved path: sources/src/DUSTMod.F90.
- Byte size: 74689.
- Line count: 1309.
- sha256sum: dcfb4594d7edd2cadf37de13fb71667d9f4b5ea929bd7d0a1d67771905f2ae7f.
- Archive member sha256: dcfb4594d7edd2cadf37de13fb71667d9f4b5ea929bd7d0a1d67771905f2ae7f.
- Upstream Git blob SHA: not recorded.

### SoilStateInitTimeConstMod.F90

- Archive member: dmleung-dustmodelcesm2023-90771fe/SoilStateInitTimeConstMod.F90.
- Saved path: sources/src/SoilStateInitTimeConstMod.F90.
- Byte size: 33793.
- Line count: 692.
- sha256sum: 4199dfbfbcd22ef5fbc8ee05551873923ff74ff8bdc1fe8f88fb62328d862df4.
- Archive member sha256: 4199dfbfbcd22ef5fbc8ee05551873923ff74ff8bdc1fe8f88fb62328d862df4.
- Upstream Git blob SHA: not recorded.

### Other extracted members (not cited, not hashed in this file)

SoilStateType.F90, clm_varctl.F90, controlMod.F90, namelist_definition_clm4_5.xml, user_nl_cam, user_nl_clm, README.md, LICENSE.

## Paths and retrieval

- Zenodo archive saved as sources/zip/z.zip. Its local MD5 matches the Zenodo published MD5.
- The archive has one root folder, dmleung-dustmodelcesm2023-90771fe/, with flat files. Extracted copies drop that root folder. So sources/src/DUSTMod.F90 is archive member dmleung-dustmodelcesm2023-90771fe/DUSTMod.F90.
- The retrieval ran and built no code.

## CESM 2.1.5 base files

The Zenodo archive holds only the files that Leung changed. `DUSTMod.F90` also needs CTSM and CESM-share constants. The README says the code is for CESM 2.1. CESM `release-cesm2.1.5` (`Externals.cfg`) uses CTSM `release-clm5.0.37` and CIME `cime5.6.49`. The exact CESM 2.1.x that Leung used is not stated, so this choice is UNSURE.

The original `DUSTMod.F90` and `SoilStateInitTimeConstMod.F90` show what Leung changed. Folder: `base_cesm2.1.5/`. Each file comes from the tag-pinned raw URL via `curl`. Each Git blob SHA matches the GitHub contents API.

| Repo | Tag | Upstream path | Bytes | sha256 |
|---|---|---|---|---|
| ESCOMP/CTSM | release-clm5.0.37 | `src/main/clm_varcon.F90` | 18326 | `2cdd451915345c5918a04b8176fb831126321defa0944c1c58e660ccdf92a9b5` |
| ESCOMP/CTSM | release-clm5.0.37 | `src/main/clm_varpar.F90` | 9519 | `b8bd940e8b1148e6c7ae0d805e5bc8d1295c998abdc43d2ef6ed0a5f5b00c369` |
| ESCOMP/CTSM | release-clm5.0.37 | `src/biogeochem/DUSTMod.F90` | 44065 | `5ed79de55bc4ada72b9008b012051b07094394ae7e68e223654cf82d8d57401a` |
| ESCOMP/CTSM | release-clm5.0.37 | `src/biogeophys/SoilStateInitTimeConstMod.F90` | 29598 | `8fd8f889933ee5953e56bfaa840fe4906697ddff47818f463a122d8c1c486de7` |
| ESMCI/cime | cime5.6.49 | `src/share/util/shr_const_mod.F90` | 6055 | `1c229a98c2d75aa9a56a2dd5700469264b82a70c428100503df33bc129957b63` |

Gravity: `clm_varcon.F90` L52 sets `grav = SHR_CONST_G`. `shr_const_mod.F90` L22 sets `SHR_CONST_G = 9.80616`.

## Reference paper

Leung, D. M., et al. (2024), Atmos. Chem. Phys., 24, 2287-2318, https://doi.org/10.5194/acp-24-2287-2024, CC BY 4.0.

The PDF is at https://acp.copernicus.org/articles/24/2287/2024/acp-24-2287-2024.pdf. Section 2.3 names CESM 2.1 with no sub-version. The code-availability section cites Zenodo 10621844, the same record as the code pin.

## Extra CESM 2.1.5 base files

The compiled Fortran comparison needs three more upstream files. Each file comes from the tag-pinned raw URL (`https://raw.githubusercontent.com/<repo>/<tag>/<path>`) via `curl`. The local Git blob SHA (`sha1("blob <size>\0" + bytes)`) matches the `sha` field returned by the GitHub contents API (`https://api.github.com/repos/<repo>/contents/<path>?ref=<tag>`), and the byte size matches too.

| Repo | Tag | Upstream path | Saved as | Bytes | sha256 | Git blob SHA (API match) |
|---|---|---|---|---|---|---|
| ESMCI/cime | cime5.6.49 | `src/share/util/shr_kind_mod.F90` | `base_cesm2.1.5/src__share__util__shr_kind_mod.F90` | 1097 | `482493e1a14da351fe6ba17c7423cbc742e95a8edb7022dad0706059db2e7c0c` | `e9e7d170ca2506ac65c188a7b89230b6521d3751` |
| ESMCI/cime | cime5.6.49 | `src/share/util/shr_spfn_mod.F90` | `base_cesm2.1.5/src__share__util__shr_spfn_mod.F90` | 37293 | `44ad599f4045420b68d3fe9e876a872571ab0fa75eddecc13d0cb048222e58d3` | `7bc5a0ee858c18ed208427d93f4e22240db0aafe` |
| ESCOMP/CTSM | release-clm5.0.37 | `src/main/landunit_varcon.F90` | `base_cesm2.1.5/src__main__landunit_varcon.F90` | 4766 | `02e06d50f5eb1d19db73e52e734b12b5903cef0bf99ff71df4451b3c1b50ccd0` | `b6ddc7cf5c1305325c851833d435801336e5705e` |

Use:

- `shr_kind_mod.F90`: `SHR_KIND_R8 = selected_real_kind(12)` (L7). The comparison compiles this file as is for 8-byte reals. For 4-byte reals it changes only that line to `selected_real_kind( 6)`.
- `shr_spfn_mod.F90`: L14-17 define `HAVE_ERF_INTRINSICS` under `__GFORTRAN__`, and L183-185 (the r8 function `shr_spfn_erf_r8`, L167-198) make `shr_spfn_erf` call the Fortran `erf` intrinsic. `InitDustVars` (DUSTMod.F90 L1007) uses `shr_spfn_erf`, so the compiled Fortran calls the intrinsic `erf`. `DustEmission` has no `use` for `erf`, so it also resolves to the intrinsic. With Intel (`CPRINTEL`), CESM would use the file's own rational approximation instead. Which compiler Leung used is not stated (UNSURE).
- `landunit_varcon.F90`: L19-20 `istsoil = 1`, `istcrop = 2`.

## CAM files (dust_emis_fact)

- Repo: https://github.com/ESCOMP/CAM. Tag: cam_cesm2_1_rel_60 (annotated tag object 0b37e5fbda4d6a2b9756e6bcf9e5a85b648a6f58).
- Commit SHA: a03b84b7c4e34f965b115686f22a043b85739e56 (from GitHub API git/tags, and git rev-parse HEAD of the clone).
- Tag source: ESCOMP/CESM release-cesm2.1.5 Externals.cfg, [cam] tag = cam_cesm2_1_rel_60.
- Retrieved: 2026-10-12 by git clone --depth 1 --branch of the tag. The files come from git show HEAD:path (byte-exact blobs).
- Files found by grep -rI dust_emis_fact over the whole tag. Not saved: doc/ChangeLog, doc/ReleaseNotes (prose), use_cases/*.xml. In use_cases/*.xml, sd_waccm_ma_cam6.xml sets 0.7D0 for 1.9x2.5, and waccmxie_ma_2000_cam4.xml has a commented line.
- Saved at: sources/cam/.

| Saved file | URL | sha256 |
| --- | --- | --- |
| cam/src__chemistry__modal_aero__dust_model.F90 | https://raw.githubusercontent.com/ESCOMP/CAM/cam_cesm2_1_rel_60/src/chemistry/modal_aero/dust_model.F90 | 350dbc2896a25d71194e1cbb728b3717f959b89addbf9e5e330ed483fc82f399 |
| cam/src__chemistry__aerosol__soil_erod_mod.F90 | https://raw.githubusercontent.com/ESCOMP/CAM/cam_cesm2_1_rel_60/src/chemistry/aerosol/soil_erod_mod.F90 | 398d88cd3ae041bc04fabd0fa11f1c9e0d0a12839a2bfb21ce8635c175b38bf9 |
| cam/src__chemistry__bulk_aero__dust_model.F90 | https://raw.githubusercontent.com/ESCOMP/CAM/cam_cesm2_1_rel_60/src/chemistry/bulk_aero/dust_model.F90 | 1636a12fdc9cec65508c51e3d2134804b43b77e98cbd01b367db1ca691f29e92 |
| cam/bld__build-namelist | https://raw.githubusercontent.com/ESCOMP/CAM/cam_cesm2_1_rel_60/bld/build-namelist | b4f3ff8894c8ffed64c213d287fd5679faf0212373f8ffb4382b338f4c7b9a9c |
| cam/bld__namelist_files__namelist_defaults_cam.xml | https://raw.githubusercontent.com/ESCOMP/CAM/cam_cesm2_1_rel_60/bld/namelist_files/namelist_defaults_cam.xml | 6026cbe73b451dcd112a97c994309fcce6b0c5aa9293c2d0a934af138a51fc5a |
| cam/bld__namelist_files__namelist_definition.xml | https://raw.githubusercontent.com/ESCOMP/CAM/cam_cesm2_1_rel_60/bld/namelist_files/namelist_definition.xml | 8047fbd732e686cbec5577ff359f848014b68e5d49a37f4808e6a514980ad081 |

## Extra CTSM files for watsat(c,1)

- Repo: https://github.com/ESCOMP/CTSM. Tag: release-clm5.0.37, same as the base_cesm2.1.5 files above (CESM release-cesm2.1.5 Externals.cfg, as stated above).
- Commit SHA: 547ae9643e964b20dd66e3e95b55340c10f9720c (git rev-parse HEAD of the clone). The annotated tag object is 5cfa1ad60172476897eea19dc02b734db54ca3ec.
- Retrieved: 2026-10-12 by git clone --depth 1 --branch of the tag. The files come from git show HEAD:path (byte-exact blobs, cmp checked).
- Not saved but read: bld/namelist_files/namelist_defaults_clm4_5.xml L153 (soil_layerstruct clm5_0 = 20SL_8.5m).

| Repo | Tag | Upstream path | Saved as | Bytes | sha256 | URL |
|---|---|---|---|---|---|---|
| ESCOMP/CTSM | release-clm5.0.37 | `src/main/FuncPedotransferMod.F90` | `base_cesm2.1.5/src__main__FuncPedotransferMod.F90` | 5571 | `b9c9e8d2091ed23b0db2a42dc69da3a4a4bd0bdde483fc7baa4dfe11fab156e0` | https://raw.githubusercontent.com/ESCOMP/CTSM/release-clm5.0.37/src/main/FuncPedotransferMod.F90 |
| ESCOMP/CTSM | release-clm5.0.37 | `src/main/organicFileMod.F90` | `base_cesm2.1.5/src__main__organicFileMod.F90` | 3318 | `c9906dc8139fc305518b83a118a464814bc6901d56a774283612bd608571408d` | https://raw.githubusercontent.com/ESCOMP/CTSM/release-clm5.0.37/src/main/organicFileMod.F90 |
| ESCOMP/CTSM | release-clm5.0.37 | `src/main/initVerticalMod.F90` | `base_cesm2.1.5/src__main__initVerticalMod.F90` | 32183 | `ee8eb2227b56aa94a3f13b7ac2c328eda1fde92b224204615e0606c9afac5483` | https://raw.githubusercontent.com/ESCOMP/CTSM/release-clm5.0.37/src/main/initVerticalMod.F90 |

## CLM parameter file value (organic_max)

- File: `lnd/clm2/paramdata/clm5_params.c171117.nc`, the `phys="clm5_0"` default at `bld/namelist_files/namelist_defaults_clm4_5.xml:381` of release-clm5.0.37. Leung `user_nl_clm` sets no `paramfile`.
- Read on 2026-10-12 from `clm5_params.c171117.nc` in the CESM input data folder `lnd/clm2/paramdata/`, sha256 f255cceecfcdffbd44bc2a4e1fe6777193751baff30688389dd04c56029cbf0f.
- `ncdump -v organic_max`: `organic_max = 130` (kg m-3). The file is not copied into the repo (binary input data).
