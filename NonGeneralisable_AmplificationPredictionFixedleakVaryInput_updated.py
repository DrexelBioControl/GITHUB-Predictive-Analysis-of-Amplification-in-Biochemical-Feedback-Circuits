# Non generalisatble code:
# use to get figure to match experimental amplification figure

# ============================================================
# PREDICT AMPLIFICATION FROM BEST-FIT MODEL
# ============================================================

import json
import sys
import numpy as np
import matplotlib.pyplot as plt
import xlrd

from pathlib import Path
from scipy.integrate import solve_ivp

# Make paths relative to this script, not the terminal folder
BASE_DIR = Path(__file__).resolve().parent
sys.path.append(str(BASE_DIR))

# Load model registry
from model_registry import get_model


# ============================================================
# 0) LOAD CONFIG
# ============================================================

CONFIG_FILE = BASE_DIR / "configs" / "Fig4H_after6C_editedfracparam" / "fit_exps12345_model3.json"

config = json.load(open(CONFIG_FILE))

RUN_NAME = config["run_name"]
MODEL_NAME = config["model_name"]

model = get_model(MODEL_NAME)

FIT_CONDITIONS = config["fit_conditions"]
ALL_CONDITIONS = config["all_conditions"]

# Identify held-out condition if one exists
held_out_list = [
    c for c in ALL_CONDITIONS
    if c not in FIT_CONDITIONS
]
HELD_OUT_CONDITION = (
    held_out_list[0]
    if len(held_out_list) > 0
    else None
)

FIT_PARAMS = config["fit_params"]

ALL_PARAM_INFO = {
    p: {"central": float(v)}
    for p, v in config["all_param_central"].items()
}

# ------------------------------------------------------------
# Load optional manually specified plateau-start rows
#
# Some configuration files store row_end_map as null because
# no condition-specific plateau cutoff rows were supplied.
# In that case, use an empty dictionary so the code can fall
# back to estimating the plateau from the final 20% of the trace.
# ------------------------------------------------------------

row_end_map_config = config.get("row_end_map") or {}

row_end_map = {
    int(column_index): int(row_end)
    for column_index, row_end in row_end_map_config.items()
}

if len(row_end_map) == 0:
    print(
        "No row_end_map was provided. "
        "The final 20% of each trace will be used to estimate DRL_0."
    )

EXCEL_FILE = BASE_DIR / config.get(
    "excel_file",
    "data/in_vitro_ctRSD_modeling_data.xls"
)

# ------------------------------------------------------------
# Excel worksheet used by this non-generalisable script
#
# This script was written specifically for the original Figure 4H
# worksheet layout, so override the worksheet name stored in the
# JSON configuration.
# ------------------------------------------------------------

EXCEL_SHEET = "Fig. 4H (SciAdv) fuel"

# ------------------------------------------------------------
# Locate saved fitting results
#
# Read both the result folder and the exact .npz filename from
# the JSON configuration. This avoids assuming that the saved
# result must follow the older "{RUN_NAME}_results.npz" format.
# ------------------------------------------------------------

RESULTS_DIR = BASE_DIR / config.get(
    "results_dir",
    "results"
)

FIT_RESULTS_FILENAME = config.get(
    "fit_results_file",
    f"{RUN_NAME}_results.npz"
)

FIT_RESULTS_FILE = RESULTS_DIR / FIT_RESULTS_FILENAME


# ------------------------------------------------------------
# Check that the saved fitting-result file exists
#
# Stop here with a clear error message rather than allowing
# np.load() to fail later with a less informative traceback.
# ------------------------------------------------------------

if not FIT_RESULTS_FILE.exists():

    raise FileNotFoundError(
        "\nCould not find the saved fitting-result file.\n\n"
        f"Expected file:\n"
        f"    {FIT_RESULTS_FILE}\n\n"
        "Check that:\n"
        "1. The fitting run has already been completed.\n"
        "2. The .npz file is inside the configured results folder.\n"
        "3. fit_results_file in the JSON matches the actual filename.\n"
    )


# ------------------------------------------------------------
# Create figure-output folder
#
# Use the configured figures_dir when one is supplied. Otherwise,
# fall back to the general figures folder.
# ------------------------------------------------------------

FIGURES_DIR = BASE_DIR / config.get(
    "figures_dir",
    "figures"
)

FIGURES_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 1) LOAD EXPERIMENTAL DATA
# ============================================================

wb = xlrd.open_workbook(str(EXCEL_FILE))
ws = wb.sheet_by_name(EXCEL_SHEET)

CONDITIONS = {
    "0in_0fuel":      (0, 0.0, 0.0),
    "1p25in_0fuel":  (1, 1.25, 0.0),
    "2p5in_0fuel":   (2, 2.5, 0.0),
    "0in_25fuel":    (3, 0.0, 25.0),
    "1p25in_25fuel": (4, 1.25, 25.0),
    "2p5in_25fuel":  (5, 2.5, 25.0),
}


# ============================================================
# 2) USER CHOICES
# ============================================================

IN_VALUES = [1.25, 2.5]

rate_floor = 1e-8   # nM/s


# ============================================================
# 3) LOAD BEST-FIT PARAMETERS
# ============================================================

def load_best_fit_params(results_file, all_param_info):

    # Load saved optimizer output
    data = np.load(results_file, allow_pickle=True)

    # Extract best-fit log10 parameter vector
    best_x = data["best_x"]

    # Retrieve fitted parameter names
    fit_params = list(data["fit_params"])

    # Start from central/default parameter values
    params = {
        p: all_param_info[p]["central"]
        for p in all_param_info
    }

    # Replace fitted parameters with optimized physical values
    for i, p in enumerate(fit_params):
        params[p] = 10**best_x[i]

    return params


best_params = load_best_fit_params(
    FIT_RESULTS_FILE,
    ALL_PARAM_INFO
)

print(f"Loaded best-fit parameters from: {FIT_RESULTS_FILE}")

print("\nFinal parameter values used for prediction:")

for p, v in best_params.items():
    print(f"{p:<15} = {v:.8e}")


# ============================================================
# 4) BUILD FULL DATASETS FOR PREDICTION
# ============================================================

def build_full_dataset(cond):

    # Look up column index and experimental concentrations
    col_idx, IN_conc, Fuel_conc = CONDITIONS[cond]

    # Lists for the full experimental timecourse
    t = []
    x = []


    # ------------------------------------------------------------
    # Load numerical trajectory from Excel
    #
    # The "All Data" worksheet can contain headings, labels, and
    # blank cells in addition to the numerical timecourse.
    #
    # Only retain rows in which both the time cell and the signal
    # cell contain numerical values. This prevents text values from
    # being converted into a string array by NumPy.
    # ------------------------------------------------------------

    for r in range(3, ws.nrows):

        # Time is stored in Excel column 0
        time_value = ws.cell_value(r, 0)

        # Signal column is shifted by one because time uses column 0
        signal_value = ws.cell_value(r, col_idx + 1)

        # Skip headings, labels, blanks, and other non-numerical rows
        if not isinstance(time_value, (int, float)):
            continue

        if not isinstance(signal_value, (int, float)):
            continue

        # Store numerical values explicitly as floats
        t.append(float(time_value))
        x.append(float(signal_value))


    # ------------------------------------------------------------
    # Check that numerical data were found
    # ------------------------------------------------------------

    if len(t) == 0:

        raise ValueError(
            f"No numerical data were found for condition '{cond}'. "
            f"Check EXCEL_SHEET and the column mapping in CONDITIONS."
        )


    # Convert time to a numerical NumPy array
    t = np.asarray(
        t,
        dtype=float
    )

    # Convert fraction reacted to nM reacted
    x = np.asarray(
        x,
        dtype=float
    ) * 500.0

    # Use configured plateau cutoff if available
    if col_idx in row_end_map:

        row_end = row_end_map[col_idx]

    # Otherwise use last 20% of trace as plateau estimate
    else:

        row_end = int(0.8 * ws.nrows)

    # Store numerical plateau-region signal values
    x_plateau = []


    # ------------------------------------------------------------
    # Read plateau region for reporter estimate
    #
    # As above, ignore blank cells and text labels so that only
    # numerical signal values are used to estimate DRL_0.
    # ------------------------------------------------------------

    for r in range(row_end, ws.nrows):

        signal_value = ws.cell_value(
            r,
            col_idx + 1
        )

        if not isinstance(signal_value, (int, float)):
            continue

        x_plateau.append(
            float(signal_value)
        )


    # Convert plateau signal to nM
    x_plateau = np.asarray(
        x_plateau,
        dtype=float
    ) * 500.0


    # Estimate initial reporter pool from the plateau.
    #
    # If no numerical plateau values remain after filtering,
    # use the final numerical point in the full trajectory.
    DRL_0 = (
        float(np.mean(x_plateau))
        if len(x_plateau) > 0
        else float(x[-1])
    )

    # Return dataset in same format expected by model.initial_conditions()
    return {
        "name": cond,
        "t_exp": t,
        "x_exp": x,
        "IN_conc": IN_conc,
        "Fuel_conc": Fuel_conc,
        "DRL_0": DRL_0,
        "RSD_temp": 25.0
    }


# Build all six experimental conditions needed for amplification
FULL_DATASETS = {
    cond: build_full_dataset(cond)
    for cond in CONDITIONS
}


# Print the number of numerical data points loaded for each condition
print("\nNumerical experimental data loaded:")

for cond, dataset in FULL_DATASETS.items():

    print(
        f"{cond:<18} : "
        f"{len(dataset['t_exp'])} points, "
        f"DRL_0 = {dataset['DRL_0']:.3f} nM"
    )


def get_full_dataset_by_condition(IN_conc, Fuel_conc):

    # Find dataset matching desired input/fuel condition
    for d in FULL_DATASETS.values():

        if (
            d["IN_conc"] == IN_conc and
            d["Fuel_conc"] == Fuel_conc
        ):
            return d

    raise ValueError(
        f"No dataset found for IN = {IN_conc}, Fuel = {Fuel_conc}"
    )


# ============================================================
# 5) SIMULATE MODEL AND EXTRACT dROL/dt FROM THE ODE
# ============================================================

def simulate_ROL_and_rate(d, params):

    # Read experimental time grid in minutes
    t_min = np.asarray(d["t_exp"], dtype=float)

    # Convert time grid to seconds for ODE solver
    t_sec = t_min * 60.0

    # Build model initial condition vector
    y0 = model.initial_conditions(d, params)

    # Solve model on experimental time grid
    sol = solve_ivp(

        lambda t, y: model.rhs(
            t,
            y,
            RSD_temp=d["RSD_temp"],
            IN_temp=d["IN_conc"],
            F_temp=d["Fuel_conc"],
            params=params
        ),

        (t_sec[0], t_sec[-1]),

        y0,

        t_eval=t_sec,

        method="LSODA",

        rtol=1e-6,
        atol=1e-6
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    # Extract simulated reporter output trajectory
    ROL = sol.y[model.output_index]

    # Allocate array for model-based reporter production rate
    dROL_dt = np.zeros_like(ROL)

    # Evaluate model RHS at every simulated time point
    for i in range(len(t_sec)):

        rhs_i = model.rhs(
            t_sec[i],
            sol.y[:, i],
            RSD_temp=d["RSD_temp"],
            IN_temp=d["IN_conc"],
            F_temp=d["Fuel_conc"],
            params=params
        )

        # Extract dROL/dt directly from ROL equation
        dROL_dt[i] = rhs_i[model.output_index]

    return ROL, dROL_dt


# ============================================================
# 6) COMPUTE PREDICTED AMPLIFICATION
# ============================================================

predicted_results = []

for IN_val in IN_VALUES:

    # Retrieve four conditions needed for normalized amplification
    d_IN_0  = get_full_dataset_by_condition(IN_val, 0.0)
    d_0_0   = get_full_dataset_by_condition(0.0,    0.0)
    d_IN_25 = get_full_dataset_by_condition(IN_val, 25.0)
    d_0_25  = get_full_dataset_by_condition(0.0,    25.0)

    # Use shared experimental time grid
    t_common_min = d_IN_0["t_exp"]

    if not (
        np.allclose(t_common_min, d_0_0["t_exp"]) and
        np.allclose(t_common_min, d_IN_25["t_exp"]) and
        np.allclose(t_common_min, d_0_25["t_exp"])
    ):

        raise ValueError("Time grids are not identical.")

    # Simulate four conditions and extract ODE-based production rates
    ROL_IN_0,  dIN0_dt  = simulate_ROL_and_rate(d_IN_0,  best_params)
    ROL_0_0,   d00_dt   = simulate_ROL_and_rate(d_0_0,   best_params)
    ROL_IN_25, dIN25_dt = simulate_ROL_and_rate(d_IN_25, best_params)
    ROL_0_25,  d025_dt  = simulate_ROL_and_rate(d_0_25,  best_params)

    # Allocate amplification array
    A_pred_new = np.full_like(
        t_common_min,
        np.nan,
        dtype=float
    )

    # Only compute amplification where all four rates are sufficiently nonzero
    valid = (
        (np.abs(dIN25_dt) >= rate_floor) &
        (np.abs(d025_dt)  >= rate_floor) &
        (np.abs(dIN0_dt)  >= rate_floor) &
        (np.abs(d00_dt)   >= rate_floor)
    )

    # Compute model-predicted normalized amplification
    A_pred_new[valid] = (

        dIN25_dt[valid] * d00_dt[valid]

    ) / (

        d025_dt[valid] * dIN0_dt[valid]

    )

    predicted_results.append({

        "IN_val": IN_val,

        "t_common_min": t_common_min,

        "ROL_IN_0": ROL_IN_0,
        "ROL_0_0": ROL_0_0,
        "ROL_IN_25": ROL_IN_25,
        "ROL_0_25": ROL_0_25,

        "dIN0_dt": dIN0_dt,
        "d00_dt": d00_dt,
        "dIN25_dt": dIN25_dt,
        "d025_dt": d025_dt,

        "A_pred_new": A_pred_new,

        "valid": valid
    })


print(
    f"\nComputed predicted amplification for "
    f"{len(predicted_results)} input conditions."
)


# ============================================================
# 7) PLOT PREDICTED ROL TRAJECTORIES AND PRODUCTION RATES
# ============================================================

for res in predicted_results:

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].plot(
        res["t_common_min"],
        res["ROL_IN_25"],
        label="IN, Fuel = 25"
    )

    axes[0].plot(
        res["t_common_min"],
        res["ROL_0_25"],
        label="OFF, Fuel = 25"
    )

    axes[0].set_xlabel("Time (min)")
    axes[0].set_ylabel("ROL reacted (nM)")
    axes[0].set_title(f"Predicted ROL | IN = {res['IN_val']} nM")
    axes[0].grid(True)
    axes[0].legend(fontsize=8)
    axes[0].set_xlim(0, 1000)

    axes[1].plot(
        res["t_common_min"],
        res["dIN25_dt"],
        label="IN rate"
    )

    axes[1].plot(
        res["t_common_min"],
        res["d025_dt"],
        label="OFF rate"
    )

    axes[1].set_xlabel("Time (min)")
    axes[1].set_ylabel("dROL/dt (nM/s)")
    axes[1].set_title("ODE-based production rate")
    axes[1].grid(True)
    axes[1].legend(fontsize=8)
    axes[1].set_xlim(0, 1000)

    A_plot = np.where(
        res["valid"],
        res["A_pred_new"],
        np.nan
    )

    axes[2].plot(
        res["t_common_min"],
        A_plot
    )

    axes[2].set_xlabel("Time (min)")
    axes[2].set_ylabel("Predicted amplification")
    axes[2].set_title("Normalized amplification")
    axes[2].grid(True)
    axes[2].set_xlim(0, 300)
    axes[2].set_ylim(-0.01, 2)

    fig.suptitle(
        f"Predicted amplification | {RUN_NAME} | IN = {res['IN_val']} nM"
    )

    plt.tight_layout()
    plt.show()


# ============================================================
# 8) OVERLAY PREDICTED AMPLIFICATION
# ============================================================

plt.figure(figsize=(10, 5))

for res in predicted_results:

    A_plot = np.where(
        res["valid"],
        res["A_pred_new"],
        np.nan
    )

    plt.plot(
        res["t_common_min"],
        A_plot,
        label=f'Input template = {res["IN_val"]} nM'
    )

plt.xlabel("Time (min)")
plt.ylabel("Predicted amplification")


# Build concise figure title
title_str = (
    f"Amplification prediction | "
    f"Model = {MODEL_NAME} | "
    f"Experiments fitted = {len(FIT_CONDITIONS)} | "
    f"Held-out = {HELD_OUT_CONDITION}"
)

plt.title(title_str)
plt.grid(True)
plt.legend()
plt.xlim(0, 300)
plt.ylim(-0.01, 2)
plt.tight_layout()

overlay_save_file = (
    FIGURES_DIR /
    f"{RUN_NAME}_predicted_amplification_overlay.png"
)

plt.savefig(
    overlay_save_file,
    dpi=300,
    bbox_inches="tight"
)

print(f"Saved figure: {overlay_save_file}")

plt.show()