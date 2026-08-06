# ============================================================
# MASTER AMPLIFICATION PREDICTION
# Heatmap + line plots for selected fitted model
# ============================================================

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from matplotlib.ticker import MaxNLocator


# ============================================================
# 0) PATHS AND USER CHOICES
# ============================================================

BASE_DIR = Path(__file__).resolve().parent
sys.path.append(str(BASE_DIR))

from model_registry import get_model


# ------------------------------------------------------------
# Choose fitted model/config here
# ------------------------------------------------------------
CONFIG_FILE = BASE_DIR / "configs" / "Fig4H_after6C_editedfracparam" / "fit_exps12345_model3.json"

# ------------------------------------------------------------
# Choose input concentration for ON condition
# ------------------------------------------------------------
IN_CONC = 2.5
IN_OFF = 0.0

# ------------------------------------------------------------
# Fuel sweep for heatmap
# ------------------------------------------------------------
FUELS = np.linspace(0.0, 50.0, 51)

# Fuel values to plot as line plots
FUEL_VALUES_TO_PLOT = [1,3,5,10,30,50]

# ------------------------------------------------------------
# Time grid for prediction
# ------------------------------------------------------------
T_END_MIN = 200
N_TIMEPOINTS = 1000

t_eval_sec = np.linspace(0.0, T_END_MIN * 60.0, N_TIMEPOINTS)
t_minutes = t_eval_sec / 60.0

# ------------------------------------------------------------
# Rate floor for stable amplification ratios
# ------------------------------------------------------------
RATE_FLOOR = 1e-8

# ------------------------------------------------------------
# Optional parameter sweep:
# For model3_misfoldgate_fuel, use "kmG_F"
# For model2_inactivegate_fuel, use "klk"
# For model1_activegate_fuel, use "klkg"
# For model0_transcriptional, use "basal_frac" or set to None
# ------------------------------------------------------------
SWEEP_PARAM = "kmG_F"
SWEEP_VALUES = np.logspace(np.log10(4.11276798e-09), np.log10(4.11276798e-01), 5)

# If you only want one row using the best-fit parameter value:
# SWEEP_PARAM = None
# SWEEP_VALUES = None


# ============================================================
# 1) LOAD CONFIG, MODEL, AND FITTED PARAMETERS
# ============================================================

config = json.load(open(CONFIG_FILE))

RUN_NAME = config["run_name"]
MODEL_NAME = config["model_name"]

model = get_model(MODEL_NAME)

FIT_PARAMS = config["fit_params"]

ALL_PARAM_INFO = {
    p: {"central": float(v)}
    for p, v in config["all_param_central"].items()
}

RESULTS_DIR = BASE_DIR / config.get("results_dir", "results")

# Try common result filenames
result_candidates = [
    RESULTS_DIR / f"{RUN_NAME}_results.npz",
    RESULTS_DIR / f"{RUN_NAME}.npz",
]

if "fit_results_file" in config:
    result_candidates.append(RESULTS_DIR / config["fit_results_file"])

FIT_RESULTS_FILE = None

for p in result_candidates:
    if p.exists():
        FIT_RESULTS_FILE = p
        break

if FIT_RESULTS_FILE is None:
    raise FileNotFoundError(
        "Could not find fitted results file. Tried:\n"
        + "\n".join(str(p) for p in result_candidates)
    )


def load_best_fit_params(results_file, all_param_info):
    """
    Load best-fit parameters from saved optimizer output.
    Assumes best_x is stored in log10 parameter space.
    """

    data = np.load(results_file, allow_pickle=True)

    best_x = data["best_x"]
    fit_params = list(data["fit_params"])

    params = {
        p: all_param_info[p]["central"]
        for p in all_param_info
    }

    for i, p in enumerate(fit_params):
        params[p] = float(10 ** best_x[i])

    return params


base_params = load_best_fit_params(
    FIT_RESULTS_FILE,
    ALL_PARAM_INFO
)

print(f"Loaded model: {MODEL_NAME}")
print(f"Loaded config: {CONFIG_FILE}")
print(f"Loaded fit: {FIT_RESULTS_FILE}")

print("\nBest-fit parameter values:")
for p, v in base_params.items():
    print(f"{p:<15} = {v:.8e}")


# ============================================================
# 2) SIMULATION HELPERS
# ============================================================

# Minimal dataset dictionary expected by model.initial_conditions()
def make_dataset(IN_value, Fuel_value):

    return {
        "name": f"IN_{IN_value}_Fuel_{Fuel_value}",
        "t_exp": t_minutes,
        "x_exp": np.zeros_like(t_minutes),
        "IN_conc": float(IN_value),
        "Fuel_conc": float(Fuel_value),
        "RSD_temp": 25.0,
    }


# Simulate model and extract ROL and dROL/dt directly from the ODE
def simulate_ROL_and_rate(IN_value, Fuel_value, params):

    d = make_dataset(IN_value, Fuel_value)

    y0 = model.initial_conditions(d, params)

    sol = solve_ivp(
        lambda t, y: model.rhs(
            t,
            y,
            RSD_temp=d["RSD_temp"],
            IN_temp=d["IN_conc"],
            F_temp=d["Fuel_conc"],
            params=params,
        ),
        (t_eval_sec[0], t_eval_sec[-1]),
        y0,
        t_eval=t_eval_sec,
        method="LSODA",
        rtol=1e-6,
        atol=1e-6,
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    ROL = sol.y[model.output_index]

    dROL_dt = np.zeros_like(ROL)

    for i in range(len(t_eval_sec)):

        rhs_i = model.rhs(
            t_eval_sec[i],
            sol.y[:, i],
            RSD_temp=d["RSD_temp"],
            IN_temp=d["IN_conc"],
            F_temp=d["Fuel_conc"],
            params=params,
        )

        dROL_dt[i] = rhs_i[model.output_index]

    return ROL, dROL_dt


def compute_amplification(params):
    """
    Compute normalized amplification:

        A0(t,f) = v_ON(t,f) / v_OFF(t,f)

        A(t,f) = A0(t,f) / A0(t,0)

                 [v_ON(t,f) * v_OFF(t,0)]
               = -------------------------
                 [v_OFF(t,f) * v_ON(t,0)]
    """

    # Zero-fuel baseline
    _, v_on_0  = simulate_ROL_and_rate(IN_CONC, 0.0, params)
    _, v_off_0 = simulate_ROL_and_rate(IN_OFF, 0.0, params)

    A = np.full(
        (len(FUELS), len(t_eval_sec)),
        np.nan,
        dtype=float,
    )

    # Compute amplification for each fuel value
    for i, fuel in enumerate(FUELS):

        _, v_on_f  = simulate_ROL_and_rate(IN_CONC, fuel, params)
        _, v_off_f = simulate_ROL_and_rate(IN_OFF, fuel, params)

        valid = (
            (np.abs(v_on_f) >= RATE_FLOOR) &
            (np.abs(v_off_f) >= RATE_FLOOR) &
            (np.abs(v_on_0) >= RATE_FLOOR) &
            (np.abs(v_off_0) >= RATE_FLOOR)
        )

        A[i, valid] = (
            v_on_f[valid] * v_off_0[valid]
        ) / (
            v_off_f[valid] * v_on_0[valid]
        )

    return A


# ============================================================
# 3) RUN PARAMETER SWEEP
# ============================================================

A_all = []
row_labels = []

if SWEEP_PARAM is None:

    A = compute_amplification(base_params.copy())
    A_all.append(A)
    row_labels.append("Best-fit parameters")

else:

    if SWEEP_PARAM not in base_params:
        raise ValueError(
            f"SWEEP_PARAM = {SWEEP_PARAM} is not in parameter set.\n"
            f"Available parameters are: {list(base_params.keys())}"
        )

    for val in SWEEP_VALUES:

        params = base_params.copy()
        params[SWEEP_PARAM] = float(val)

        print(f"\nComputing amplification for {SWEEP_PARAM} = {val:.3e}")

        A = compute_amplification(params)

        A_all.append(A)
        row_labels.append(rf"${SWEEP_PARAM} = {val:.0e}$")


A_stack = np.stack(A_all, axis=0)

print("\nAmplification calculation complete.")
print(f"A range: {np.nanmin(A_stack):.3g} to {np.nanmax(A_stack):.3g}")


# ============================================================
# 4) PLOT HORIZONTAL LINE PANELS WITH A SHARED Y AXIS
# ============================================================

ncols = len(A_all)

fig, axes = plt.subplots(
    1,
    ncols,
    figsize=(5.2 * ncols, 4.8),
    sharex=True,
    sharey=True,
    constrained_layout=True,
)

# Ensure axes is always a one-dimensional NumPy array,
# including when only one parameter value is plotted.
if ncols == 1:
    axes = np.array([axes])

# ------------------------------------------------------------
# Locate the closest simulated fuel concentration corresponding
# to each requested line-plot fuel concentration.
# ------------------------------------------------------------
fuel_indices = [
    int(np.argmin(np.abs(FUELS - fuel_value)))
    for fuel_value in FUEL_VALUES_TO_PLOT
]

# ------------------------------------------------------------
# Calculate one common y-axis range using only the curves that
# will actually appear in the line panels.
# ------------------------------------------------------------
selected_line_values = []

for A in A_all:
    for idx in fuel_indices:
        selected_line_values.append(A[idx, :])

selected_line_values = np.concatenate(selected_line_values)
finite_line_values = selected_line_values[
    np.isfinite(selected_line_values)
]

if finite_line_values.size == 0:
    raise RuntimeError(
        "No finite amplification values were available for plotting. "
        "Check RATE_FLOOR and the simulated production rates."
    )

y_min = min(0.0, float(np.nanmin(finite_line_values)))
y_max = max(1.0, float(np.nanmax(finite_line_values)))

y_range = y_max - y_min
if y_range == 0:
    y_range = 1.0

y_padding = 0.05 * y_range
shared_ylim = (
    y_min - y_padding,
    y_max + y_padding,
)

# ------------------------------------------------------------
# Plot one panel for each swept leak-parameter value.
# ------------------------------------------------------------
for panel_index, (A, ax) in enumerate(zip(A_all, axes)):

    # Zero-fuel normalized reference.
    ax.plot(
        t_minutes,
        np.ones_like(t_minutes),
        linestyle="--",
        linewidth=2.5,
        color="black",
        label="0 nM",
    )

    # Selected nonzero fuel conditions.
    for idx in fuel_indices:
        ax.plot(
            t_minutes,
            A[idx, :],
            linewidth=2.0,
            label=f"{FUELS[idx]:.0f} nM",
        )

    ax.set_title(
        row_labels[panel_index],
        fontsize=12,
    )

    ax.set_xlabel("Time (min)")
    ax.set_xlim(0, T_END_MIN)

    # Fixed shared y-axis from 0 to 6
    ax.set_ylim(-0.05, 4.5)
    ax.set_yticks(np.arange(0, 5, 1))

    ax.grid(
        True,
        alpha=0.3,
    )

    ax.xaxis.set_major_locator(
        MaxNLocator(integer=True)
    )

# A shared y axis means only the first panel needs the label.
axes[0].set_ylabel("Predicted amplification")

# Show a single legend for the full panel set.
axes[-1].legend(
    title="Fuel",
    loc="upper right",
    fontsize=9,
    title_fontsize=10,
    frameon=True,
)

fig.suptitle(
    f"Amplification prediction | {MODEL_NAME} | IN = {IN_CONC} nM",
    fontsize=16,
)

# Optional high-resolution export for PowerPoint.
OUTPUT_FILE = BASE_DIR / (
    f"amplification_timecourses_{MODEL_NAME}_"
    f"IN_{IN_CONC:g}nM.png"
)

fig.savefig(
    OUTPUT_FILE,
    dpi=400,
    bbox_inches="tight",
)

print(f"Saved figure: {OUTPUT_FILE}")

plt.show()
