# Dust emission

Schemes that compute the vertical flux of mineral dust from the surface.

| Scheme | Description | Status |
|---|---|---|
| [`ginoux`](ginoux/) | GOCART scheme of WRF-Chem (dust_opt=1), 5 bins | implemented |
| [`afwa`](afwa/) | AFWA scheme of WRF-Chem (dust_opt=3), 5 bins | implemented |
| [`shao2011`](shao2011/) | UoC Shao 2011 scheme of WRF-Chem (dust_opt=4), 5 bins | implemented |
| [`fengsha_ls99`](fengsha_ls99/) | FENGSHA of CMAQ (Lu and Shao 1999 form), 2 modes | implemented |
| [`fengsha_rrfs`](fengsha_rrfs/) | FENGSHA of RRFS-SD (Webb 2020 form), 5 bins | implemented |
| [`leung2023`](leung2023/) | Leung et al. 2023 scheme (CESM2.1 code), 4 bins | implemented |
| [`fengsha`](fengsha/) | LS99-FENGSHA, Ma et al. (2019), CMAQ 5.5 | placeholder |
| [`kok`](kok/) | Kok et al. (2014), CESM2/CLM5 | implemented |
| [`zender`](zender/) | DEAD, Zender et al. (2003), CLM4 | implemented |
| [`leung`](leung/) | Leung et al. (2023, 2024), CESM2/CLM5 | implemented |

Each scheme folder follows the layout in [CONTRIBUTING.md](../../CONTRIBUTING.md).
