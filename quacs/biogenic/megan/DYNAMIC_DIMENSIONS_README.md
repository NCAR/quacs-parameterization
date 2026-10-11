# MEGAN dynamic dimensions

Run from the repository root:

```bash
python -m pip install numpy pandas pyyaml
python -m quacs.biogenic.megan.examples.driver_megan_v3
```

`configure.yaml` controls model counts, input file paths, activity switches, LDF class indices, and the default per-class LDF array. The class names/order must match `inputs/activity_class_parameters.csv`, `inputs/class_to_species.csv`, and the `EF1...EFn` columns in the plant EF file. The activity CSV contains one row per activity class with all required response coefficients. For a new class, define its response coefficients, plant EF column, LDF and class-to-species mapping; a new physical response *formula* still requires new source code.

`inputs/canopy_characteristic.csv` is read by characteristic name; unused extra rows and row reordering are supported. Canopy type columns must be named identically and ordered consistently with `inputs/PFT_Fraction.csv` (whose fractions are percentages). New canopy types must have all required characteristic values. The current physics uses the same formulas for each type.

Growth forms are configured in `files.raw_ef_files.speciation`, each with `file` and `specfrac_column`. Add a `<name>_fraction` column to `grid_growthform.csv` and update `NGFORM`. The fraction column name is derived from the growth-form key. The ecotype-weighted growth-form fractions are assumed common across ecotypes, as in the original workflow.

`n_species` must match `species_catalog.csv`, which uses contiguous one-based IDs. The mapper supports species with no source class contribution (zero emissions), but a new species must have valid mapping data to contribute to a chemical mechanism.

The numeric parameters in YAML are validated against data and do not themselves invent physical parameters. For an arbitrary new activity class, supply the same parameter schema used by existing MEGAN formulas. `NLayers` is dynamically handled by the existing quadrature function, subject to valid positive values.

## Verification and limitations

- 33 tests passed for dynamic growth form, extra canopy characteristic, activity parameter row, species catalog size, and existing model behavior.
- A synthetic two-ecotype EF dataset successfully ran the full driver and wrote `outputs/EF_LDF.csv`; synthetic data were removed from this package.
- The archive supplied for this update does **not** contain the large original `EFv240120.csv` or `Speciation*.csv` datasets. Therefore no full original-data regression or numeric parity claim is made.
- Legacy tests referring to removed `inputs/EF_LDF.csv` and `grid_vegetation.csv` were updated to validate current input metadata.

### LDF and runtime configuration (updated)

`configure.yaml` retains file paths, model dimensions and activity switches.
Advanced physical constants live at the top of `biogenic_emission_megan_v3.py`.
`inputs/activity_class_parameters.csv` has a new `ldf` column: numeric values
are fixed LDF, while `GRID` means read the corresponding `LDF{class_id}`
column from the plant EF table and aggregate over ecotypes and growth forms.
Class IDs use the order of `class_to_species.csv`; the activity table must
match that order. The driver no longer constructs settings/switches objects.

Input checks are delegated to the adapter (canopy, PFT, meteorology),
`MGN2MECH` (species and class mapping), `MEGVEA` (activity classes and LDF),
and `EFPreprocessor` (growth forms, EF/LDF values).
