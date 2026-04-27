#!/usr/bin/env python3
"""
Latin Hypercube Sampling driver for multi-species dry deposition.

Generates N_CELLS sets of atmospheric conditions using LHS, runs a MUSICA
MICM box model with dry deposition for Hg0, SO2, and O3, then produces
sensitivity plots of deposition velocity and loss rate vs. each input
dimension.

Usage::

    python -m quacs.drydep.lhs.examples.lhs_driver
    python -m quacs.drydep.lhs.examples.lhs_driver --cells 200 --seed 42

Output files are written to ``output/lhs_results.csv`` and
``output/lhs_sensitivity_*.png``.
"""

import argparse
import os

import matplotlib.pyplot as plt
import musica
import musica.mechanism_configuration as mc
import numpy as np
import pandas as pd
from scipy.stats import qmc

from quacs.drydep.lhs import olson_land_cover
from quacs.drydep.lhs.box_model import compute_drydep_rate
from quacs.drydep.lhs.lhs_config import (
    LC_ARCHETYPES,
    MET_DIMS,
    N_LC_ARCHETYPES,
    build_patches_from_sample,
)
from quacs.drydep.lhs.species import hg0 as HG0_SPECIES, o3 as O3_SPECIES, so2 as SO2_SPECIES

# ── Configuration ─────────────────────────────────────────────────────────────

SPECIES_LIST = [HG0_SPECIES, SO2_SPECIES, O3_SPECIES]
BOX_HEIGHT_M = 1.0

TIME_STEP_S = 60.0
SIM_LENGTH_S = 3600.0

MET_BACKGROUND: dict = {}  # all met inputs sampled by LHS

N_MET = len(MET_DIMS)
N_LC_FRAC = N_LC_ARCHETYPES
N_LC_LAI = N_LC_ARCHETYPES
N_INIT_CONCS = len(SPECIES_LIST)

NDIM = N_MET + N_LC_FRAC + N_LC_LAI + N_INIT_CONCS

L_BOUNDS = (
    [d[1] for d in MET_DIMS]
    + [0.01] * N_LC_FRAC
    + [lc[1] for lc in LC_ARCHETYPES]
    + [0.1, 1.0, 1e6]  # init concs: Hg0 (ng/m³), SO2, O3
)
U_BOUNDS = (
    [d[2] for d in MET_DIMS]
    + [1.0] * N_LC_FRAC
    + [lc[2] if lc[2] > lc[1] else lc[1] + 1.0 for lc in LC_ARCHETYPES]
    + [5.0, 5e3, 1.2e8]
)


def run_lhs(n_cells: int = 100, seed: int = 0) -> pd.DataFrame:
    """Run the full LHS ensemble and return results as a DataFrame."""
    print(f"Generating {n_cells} LHS samples across {NDIM} dimensions …")
    sampler = qmc.LatinHypercube(d=NDIM, seed=seed)
    sample = sampler.random(n=n_cells)
    sample_scaled = qmc.scale(sample, L_BOUNDS, U_BOUNDS)

    col = 0
    met_sample = sample_scaled[:, col:col + N_MET]; col += N_MET
    frac_sample = sample_scaled[:, col:col + N_LC_FRAC]; col += N_LC_FRAC
    lai_sample = sample_scaled[:, col:col + N_LC_LAI]; col += N_LC_LAI
    init_sample = sample_scaled[:, col:col + N_INIT_CONCS]

    met_arrays: dict = {}
    for j, (name, *_) in enumerate(MET_DIMS):
        met_arrays[name] = met_sample[:, j]
    for key, val in MET_BACKGROUND.items():
        met_arrays[key] = np.full(n_cells, float(val))

    cell_land_cover = [
        build_patches_from_sample(frac_sample[i], lai_sample[i])
        for i in range(n_cells)
    ]

    print("Computing deposition velocities …")
    k_dict = compute_drydep_rate(
        SPECIES_LIST,
        box_height_m=BOX_HEIGHT_M,
        met=met_arrays,
        land_cover=cell_land_cover,
        coefficients=olson_land_cover.coefficients,
    )

    dvel_dict = {sp.name: k_dict[sp.name] * BOX_HEIGHT_M / 0.01 for sp in SPECIES_LIST}

    print("Setting up MUSICA MICM solver …")
    mc_species = list(SPECIES_LIST)
    gas = mc.Phase(name="gas", species=mc_species)
    reactions = [
        mc.FirstOrderLoss(name=f"{sp.name}_drydep", scaling_factor=1.0,
                          reactants=[mc_sp], gas_phase=gas)
        for sp, mc_sp in zip(SPECIES_LIST, mc_species)
    ]
    mechanism = mc.Mechanism(name="lhs_drydep", species=mc_species,
                             phases=[gas], reactions=reactions)
    solver = musica.MICM(mechanism=mechanism,
                         solver_type=musica.SolverType.rosenbrock_standard_order)

    state = solver.create_state(n_cells)
    state.set_conditions(temperatures=met_arrays["TC0"], pressures=met_arrays["PRESSU"])

    init_concs = {}
    for j, sp in enumerate(SPECIES_LIST):
        init_ng_m3 = init_sample[:, j]
        init_concs[sp.name] = init_ng_m3 * 1e-9 / (sp.molecular_weight_kg_mol * 1e3)

    state.set_concentrations(init_concs)
    state.set_user_defined_rate_parameters(
        {f"LOSS.{sp.name}_drydep": k_dict[sp.name] for sp in SPECIES_LIST}
    )

    print(f"Integrating {int(SIM_LENGTH_S / TIME_STEP_S)} steps of {TIME_STEP_S} s …")
    n_steps = int(SIM_LENGTH_S / TIME_STEP_S)
    curr_time = 0.0
    results_list = []

    concs_t0 = state.get_concentrations()
    for cell in range(n_cells):
        row = {"time_s": curr_time, "cell": cell}
        for sp in SPECIES_LIST:
            row[f"conc_{sp.name}_mol_m3"] = concs_t0[sp.name][cell]
        results_list.append(row)

    for _ in range(n_steps):
        solver.solve(state, TIME_STEP_S)
        curr_time += TIME_STEP_S
        cell_concs = state.get_concentrations()
        for cell in range(n_cells):
            row = {"time_s": curr_time, "cell": cell}
            for sp in SPECIES_LIST:
                row[f"conc_{sp.name}_mol_m3"] = cell_concs[sp.name][cell]
            results_list.append(row)

    df_concs = pd.DataFrame(results_list)

    rows = []
    for i in range(n_cells):
        row = {}
        for j, (name, *_) in enumerate(MET_DIMS):
            row[f"met_{name}"] = met_sample[i, j]
        fracs_norm = frac_sample[i] / frac_sample[i].sum()
        for j, (name, *_) in enumerate(LC_ARCHETYPES):
            row[f"lc_frac_{name}"] = float(fracs_norm[j])
            row[f"lc_lai_{name}"] = float(lai_sample[i, j])
        for j, sp in enumerate(SPECIES_LIST):
            row[f"init_{sp.name}_ng_m3"] = float(init_sample[i, j])
        for sp in SPECIES_LIST:
            row[f"vd_{sp.name}_cm_s"] = float(dvel_dict[sp.name][i])
            row[f"k_{sp.name}_s"] = float(k_dict[sp.name][i])
        df_final = df_concs[
            (df_concs["cell"] == i) & (df_concs["time_s"] == SIM_LENGTH_S)
        ]
        for sp in SPECIES_LIST:
            if not df_final.empty:
                row[f"conc_final_{sp.name}_mol_m3"] = float(
                    df_final[f"conc_{sp.name}_mol_m3"].values[0]
                )
        rows.append(row)

    return pd.DataFrame(rows), df_concs


def make_sensitivity_plots(df: pd.DataFrame, output_dir: str):
    """Generate sensitivity scatter plots of v_d vs. each input dimension."""
    os.makedirs(output_dir, exist_ok=True)

    met_names = [d[0] for d in MET_DIMS]
    met_labels = [d[3] for d in MET_DIMS]
    species_cols = [f"vd_{sp.name}_cm_s" for sp in SPECIES_LIST]
    species_names = [sp.name for sp in SPECIES_LIST]

    for sp_col, sp_name in zip(species_cols, species_names):
        ncols = 3
        nrows = int(np.ceil(len(met_names) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                                 constrained_layout=True)
        for ax, met_name, met_label in zip(axes.flatten(), met_names, met_labels):
            ax.scatter(df[f"met_{met_name}"], df[sp_col], s=10, alpha=0.5, color="steelblue")
            ax.set_xlabel(met_label, fontsize=9)
            ax.set_ylabel(f"$v_d$ {sp_name} (cm s⁻¹)", fontsize=9)
        for ax in axes.flatten()[len(met_names):]:
            ax.set_visible(False)
        fig.suptitle(f"Sensitivity of {sp_name} deposition velocity — met inputs\n"
                     f"(n = {len(df)} LHS samples)", fontsize=11)
        out_path = os.path.join(output_dir, f"lhs_sensitivity_{sp_name}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  Saved {out_path}")

    lc_names = [lc[0] for lc in LC_ARCHETYPES]
    for sp_col, sp_name in zip(species_cols, species_names):
        ncols = min(len(lc_names), 3)
        nrows = int(np.ceil(len(lc_names) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                                 constrained_layout=True)
        for ax, lc_name in zip(np.atleast_1d(axes).flatten(), lc_names):
            ax.scatter(df[f"lc_frac_{lc_name}"], df[sp_col], s=10, alpha=0.5,
                       color="darkorange")
            ax.set_xlabel(f"Fraction: {lc_name}", fontsize=9)
            ax.set_ylabel(f"$v_d$ {sp_name} (cm s⁻¹)", fontsize=9)
        for ax in np.atleast_1d(axes).flatten()[len(lc_names):]:
            ax.set_visible(False)
        fig.suptitle(f"Sensitivity of {sp_name} $v_d$ to land-cover fraction\n"
                     f"(n = {len(df)} LHS samples)", fontsize=11)
        out_path = os.path.join(output_dir, f"lhs_landcover_{sp_name}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  Saved {out_path}")


def make_timeseries_plot(df_concs: pd.DataFrame, output_dir: str):
    """Plot mean normalised concentration ± 1 std across cells over time."""
    os.makedirs(output_dir, exist_ok=True)
    times = sorted(df_concs["time_s"].unique())
    fig, axes = plt.subplots(1, len(SPECIES_LIST), figsize=(5 * len(SPECIES_LIST), 4),
                             constrained_layout=True)
    init_time = times[0]
    for ax, sp in zip(axes, SPECIES_LIST):
        col = f"conc_{sp.name}_mol_m3"
        init_vals = df_concs[df_concs["time_s"] == init_time].set_index("cell")[col]
        means, stds = [], []
        for t in times:
            normed = df_concs[df_concs["time_s"] == t].set_index("cell")[col] / init_vals
            means.append(normed.mean())
            stds.append(normed.std())
        times_min = np.array(times) / 60.0
        means, stds = np.array(means), np.array(stds)
        ax.plot(times_min, means, color="steelblue", label="mean")
        ax.fill_between(times_min, np.clip(means - stds, 0, None), means + stds,
                        alpha=0.3, color="steelblue", label="±1 std")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Normalised concentration")
        ax.set_title(sp.name)
        ax.legend(fontsize=8)
    fig.suptitle("Normalised concentration over time — all LHS cells", fontsize=11)
    out_path = os.path.join(output_dir, "lhs_timeseries.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run LHS ensemble for multi-species dry deposition."
    )
    parser.add_argument("--cells", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", default="output")
    args = parser.parse_args()

    df, df_concs = run_lhs(n_cells=args.cells, seed=args.seed)

    os.makedirs(args.output, exist_ok=True)
    csv_path = os.path.join(args.output, "lhs_results.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nResults written to {csv_path}")

    print("\nGenerating sensitivity plots …")
    make_sensitivity_plots(df, args.output)
    make_timeseries_plot(df_concs, args.output)

    print("\nDeposition velocity summary (across all LHS cells):")
    print(f"  {'Species':<6}  {'mean (cm/s)':>12}  {'min':>8}  {'max':>8}")
    print("  " + "-" * 42)
    for sp in SPECIES_LIST:
        col = df[f"vd_{sp.name}_cm_s"]
        print(f"  {sp.name:<6}  {col.mean():12.4f}  {col.min():8.4f}  {col.max():8.4f}")


if __name__ == "__main__":
    main()
