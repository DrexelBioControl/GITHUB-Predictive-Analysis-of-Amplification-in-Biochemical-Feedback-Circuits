# ============================================================
# MASTER SBML AMPLIFICATION PREDICTION
# ============================================================
#
# PURPOSE OF THIS FILE
# --------------------
# This is the generalized SBML version of the previous Part 2 workflow.
# Part 1 asks, "Which candidate mechanism best explains the measured data?"
# Part 2 then asks, "What does that fitted mechanism predict outside the exact
# conditions used for fitting?"
#
# The information flow is:
#
#     selected_model.json from Part 1
#              +
#     selected SBML mechanism
#              +
#     best-fit parameter values
#              ↓
#     simulate ON and OFF states across fuel/design values
#              ↓
#     calculate normalized amplification
#              ↓
#     save prediction figures
#
# The current Figure 2D proof-of-concept candidates contain no fuel species,
# so this file is intentionally disabled in the JSON config. When genuine
# fuel-containing SBML models are added:
#   1. edit configs/amplification_prediction.json
#   2. set "enabled": true
#   3. provide the correct SBML species, reaction, and sweep IDs
# ============================================================

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from sbmltoodepy_solveivp_adapter import SBMLtoODEpySolveIVPModel


# ============================================================
# 0) PATHS AND USER CHOICES
# ============================================================
# As in the previous workflow, the scientific choices live in a config file.
# The Python below reads those choices and performs the same calculation for
# whichever model Part 1 selected.

BASE_DIR = Path(__file__).resolve().parent

CONFIG_FILE = (
    BASE_DIR
    / "configs"
    / "amplification_prediction.json"
)

with open(CONFIG_FILE, "r", encoding="utf-8") as f:
    config = json.load(f)

if not bool(config.get("enabled", False)):
    print(
        "Part 2 is currently disabled in:\n"
        f"{CONFIG_FILE}\n\n"
        "This is expected for the present Figure 2D proof-of-concept, "
        "because the two current SBML models do not contain a fuel species.\n\n"
        "After adding genuine fuel-containing SBML models, set "
        "\"enabled\": true and provide the correct SBML IDs."
    )
    raise SystemExit(0)


# ============================================================
# 1) LOAD SELECTED MODEL AND FITTED PARAMETERS
# ============================================================
# Part 1 writes a small selection record rather than forcing us to type the
# winning model name manually. This prevents Part 2 from accidentally loading
# a fit from one model and an SBML file from another.

def resolve_path(path_value):
    # Config files may contain either full paths or project-relative paths.
    # Relative paths are anchored to this script so the folder remains portable.
    path = Path(path_value)

    if path.is_absolute():
        return path

    return BASE_DIR / path


SELECTED_MODEL_FILE = resolve_path(
    config["selected_model_file"]
)

if not SELECTED_MODEL_FILE.exists():
    raise FileNotFoundError(
        "Part 1 selection record was not found:\n"
        f"{SELECTED_MODEL_FILE}\n\n"
        "Run Code1_MasterSBMLTimecourseFitting first."
    )

with open(
    SELECTED_MODEL_FILE,
    "r",
    encoding="utf-8",
) as f:
    selected = json.load(f)

MODEL_NAME = selected["model_name"]

MODEL_FILE = Path(
    selected["model_file"]
)

if not MODEL_FILE.is_absolute():
    MODEL_FILE = BASE_DIR / MODEL_FILE

if not MODEL_FILE.exists():
    fallback_model = (
        BASE_DIR
        / "results"
        / "selected_model.xml"
    )

    if fallback_model.exists():
        MODEL_FILE = fallback_model
    else:
        raise FileNotFoundError(
            f"Selected SBML file was not found:\n{MODEL_FILE}"
        )

FIT_RESULTS_FILE = Path(
    selected["fit_results_file"]
)

if not FIT_RESULTS_FILE.is_absolute():
    FIT_RESULTS_FILE = (
        BASE_DIR
        / FIT_RESULTS_FILE
    )

if not FIT_RESULTS_FILE.exists():
    raise FileNotFoundError(
        f"Selected fit-results file was not found:\n"
        f"{FIT_RESULTS_FILE}"
    )

saved_fit = np.load(
    FIT_RESULTS_FILE,
    allow_pickle=True,
)

best_x = np.asarray(
    saved_fit["best_x"],
    dtype=float,
)

fit_params = [
    str(name)
    for name in saved_fit["fit_params"]
]

# Part 1 fitted in log10 space, so convert the saved optimizer vector
# back into physical kinetic parameter values before changing the SBML model.
base_params = {
    parameter_name: float(
        10 ** best_x[index]
    )
    for index, parameter_name in enumerate(
        fit_params
    )
}

GENERATED_MODELS_DIR = BASE_DIR / config.get(
    "generated_models_dir",
    "generated_models",
)

rr = SBMLtoODEpySolveIVPModel(
    sbml_file=MODEL_FILE,
    generated_models_dir=GENERATED_MODELS_DIR,
)

print(f"Loaded selected model : {MODEL_NAME}")
print(f"Loaded SBML           : {MODEL_FILE}")
print(f"Loaded fit            : {FIT_RESULTS_FILE}")

print("\nBest-fit parameters:")

for p, value in base_params.items():
    print(f"{p:<20} = {value:.8e}")


# ============================================================
# 2) DESIGN-SPACE SETTINGS
# ============================================================
# These IDs tell the generic script which biological roles are played by which
# SBML objects. The code does not assume that every model uses names such as F,
# IN, or ROL; the mapping is made explicitly in the JSON configuration.

INPUT_SPECIES_ID = str(
    config["input_species_id"]
)

FUEL_SPECIES_ID = str(
    config["fuel_species_id"]
)

GATE_SPECIES_ID = str(
    config["gate_species_id"]
)

REPORTER_SPECIES_ID = str(
    config["reporter_species_id"]
)

OBSERVABLE_ID = str(
    config["observable_id"]
)

OUTPUT_RATE_REACTION_ID = str(
    config["output_rate_reaction_id"]
)

IN_CONC = float(
    config["input_on"]
)

IN_OFF = float(
    config["input_off"]
)

GATE_CONC = float(
    config["gate_concentration"]
)

REPORTER_CONC = float(
    config["reporter_concentration"]
)

fuel_config = config["fuels"]

FUELS = np.linspace(
    float(fuel_config["start"]),
    float(fuel_config["stop"]),
    int(fuel_config["points"]),
)

FUEL_VALUES_TO_PLOT = [
    float(value)
    for value in config[
        "fuel_values_to_plot"
    ]
]

T_END_MIN = float(
    config["time_end_min"]
)

N_TIMEPOINTS = int(
    config["time_points"]
)

MODEL_TIME_UNITS_PER_MINUTE = float(
    config.get(
        "model_time_units_per_minute",
        60.0,
    )
)

t_minutes = np.linspace(
    0.0,
    T_END_MIN,
    N_TIMEPOINTS,
)

# Keep minutes for plotting, but simulate in the SBML model's time units.
# For an SBML expressed in seconds, this multiplier will normally be 60.
t_model = (
    t_minutes
    * MODEL_TIME_UNITS_PER_MINUTE
)

RATE_FLOOR = float(
    config.get(
        "rate_floor",
        1e-8,
    )
)

SWEEP_PARAM = config.get(
    "sweep_parameter"
)

SWEEP_VALUES = config.get(
    "sweep_values"
)

if SWEEP_VALUES is not None:
    SWEEP_VALUES = np.asarray(
        SWEEP_VALUES,
        dtype=float,
    )


# ============================================================
# 3) VERIFY REQUIRED SBML IDS
# ============================================================
# Fail early with a readable message if the selected model does not contain the
# species, reporter reaction, or sweep parameter requested by the config. This
# is especially important because candidate SBML files may use different IDs.

species_ids = set(rr.species_ids)
parameter_ids = set(rr.parameter_ids)
reaction_ids = set(rr.reaction_ids)

required_species = {
    INPUT_SPECIES_ID,
    FUEL_SPECIES_ID,
    GATE_SPECIES_ID,
    REPORTER_SPECIES_ID,
    OBSERVABLE_ID,
}

missing_species = sorted(
    required_species.difference(
        species_ids
    )
)

if missing_species:
    raise ValueError(
        "The selected SBML model does not contain the required "
        "Part 2 species:\n"
        + "\n".join(
            f"  - {species}"
            for species in missing_species
        )
        + "\n\nUpdate amplification_prediction.json or use a "
        "fuel-containing candidate model."
    )

if (
    OUTPUT_RATE_REACTION_ID
    not in reaction_ids
):
    raise ValueError(
        f"Reaction {OUTPUT_RATE_REACTION_ID!r} was not found "
        f"in {MODEL_FILE.name}.\n"
        f"Available reactions: {sorted(reaction_ids)}"
    )

if (
    SWEEP_PARAM is not None
    and SWEEP_PARAM not in parameter_ids
):
    raise ValueError(
        f"Sweep parameter {SWEEP_PARAM!r} was not found.\n"
        f"Available global parameters: {sorted(parameter_ids)}"
    )


# ============================================================
# 4) SIMULATION HELPERS
# ============================================================
# Each ON/OFF/fuel combination must be an independent experiment. The model is
# therefore reset to its original SBML state before condition-specific initial
# concentrations and fitted parameters are applied.

def make_initial_values(input_concentration, fuel_concentration):
    """Build one independent ON/OFF/fuel initial-condition dictionary."""
    return {
        INPUT_SPECIES_ID: float(input_concentration),
        FUEL_SPECIES_ID: float(fuel_concentration),
        GATE_SPECIES_ID: GATE_CONC,
        REPORTER_SPECIES_ID: REPORTER_CONC,
    }


def simulate_output_and_rate(
    input_concentration,
    fuel_concentration,
    params,
):
    """Solve the generated SBML ODE and extract output plus d(output)/dt."""
    initial_values = make_initial_values(
        input_concentration=input_concentration,
        fuel_concentration=fuel_concentration,
    )

    output, output_rate = rr.simulate(
        initial_values=initial_values,
        parameters=params,
        fixed_parameters={},
        t_eval=t_model,
        observable_id=OBSERVABLE_ID,
        method="LSODA",
        rtol=float(config.get("rtol", 1e-6)),
        atol=float(config.get("atol", 1e-6)),
        return_output_rate=True,
    )

    return output, output_rate


def compute_amplification(params):
    """
    Calculate the normalized amplification used in the previous workflow.

    First compare ON and OFF reporter-production rates at one fuel value. Then
    divide that ratio by the corresponding zero-fuel ON/OFF ratio. Written as
    one expression:

        A(t,f) = [v_ON(t,f) * v_OFF(t,0)]
                 -------------------------
                 [v_OFF(t,f) * v_ON(t,0)]

    This normalization removes the baseline ON/OFF separation already present
    without fuel and asks specifically how fuel changes amplification.
    """

    _, v_on_0 = simulate_output_and_rate(
        IN_CONC,
        0.0,
        params,
    )

    _, v_off_0 = simulate_output_and_rate(
        IN_OFF,
        0.0,
        params,
    )

    amplification = np.full(
        (
            len(FUELS),
            len(t_model),
        ),
        np.nan,
        dtype=float,
    )

    for fuel_index, fuel in enumerate(
        FUELS
    ):
        _, v_on_f = simulate_output_and_rate(
            IN_CONC,
            fuel,
            params,
        )

        _, v_off_f = simulate_output_and_rate(
            IN_OFF,
            fuel,
            params,
        )

        # Rate ratios become numerically meaningless when any denominator
        # or reference rate is effectively zero. RATE_FLOOR masks those times
        # rather than allowing enormous artificial amplification values.
        valid = (
            (np.abs(v_on_f) >= RATE_FLOOR)
            & (np.abs(v_off_f) >= RATE_FLOOR)
            & (np.abs(v_on_0) >= RATE_FLOOR)
            & (np.abs(v_off_0) >= RATE_FLOOR)
        )

        amplification[
            fuel_index,
            valid,
        ] = (
            v_on_f[valid]
            * v_off_0[valid]
        ) / (
            v_off_f[valid]
            * v_on_0[valid]
        )

    return amplification


# ============================================================
# 5) RUN OPTIONAL PARAMETER SWEEP
# ============================================================
# With no sweep parameter, the prediction uses the fitted model exactly as it
# emerged from Part 1. With a sweep parameter, only that one parameter is varied
# while all other fitted values remain fixed. This reproduces the logic of the
# previous leak-rate design-space analysis.

A_all = []
row_labels = []

if SWEEP_PARAM is None:
    A_all.append(
        compute_amplification(
            base_params.copy()
        )
    )

    row_labels.append(
        "Best-fit parameters"
    )

else:
    if SWEEP_VALUES is None:
        raise ValueError(
            "sweep_parameter is set but sweep_values is null."
        )

    for value in SWEEP_VALUES:
        params = base_params.copy()
        params[SWEEP_PARAM] = float(
            value
        )

        print(
            f"Computing {SWEEP_PARAM} = "
            f"{value:.3e}"
        )

        A_all.append(
            compute_amplification(
                params
            )
        )

        row_labels.append(
            f"{SWEEP_PARAM} = {value:.2e}"
        )


# ============================================================
# 6) PLOT AMPLIFICATION TIME COURSES
# ============================================================
# Each panel corresponds to one value of the optional swept parameter. Within a
# panel, each curve is one fuel concentration. The zero-fuel reference is one by
# construction and is shown as a dashed horizontal line.

ncols = len(A_all)

fig, axes = plt.subplots(
    1,
    ncols,
    figsize=(
        5.2 * ncols,
        4.8,
    ),
    sharex=True,
    sharey=True,
    constrained_layout=True,
)

if ncols == 1:
    axes = np.asarray([axes])

# Requested plot values do not need to fall exactly on the full fuel grid.
# For each requested value, select the nearest simulated fuel concentration.
fuel_indices = [
    int(
        np.argmin(
            np.abs(
                FUELS
                - requested_fuel
            )
        )
    )
    for requested_fuel
    in FUEL_VALUES_TO_PLOT
]

for panel_index, (
    amplification,
    ax,
) in enumerate(
    zip(
        A_all,
        axes,
    )
):
    ax.plot(
        t_minutes,
        np.ones_like(t_minutes),
        linestyle="--",
        linewidth=2.5,
        color="black",
        label="0 nM",
    )

    for fuel_index in fuel_indices:
        ax.plot(
            t_minutes,
            amplification[
                fuel_index,
                :,
            ],
            linewidth=2.0,
            label=(
                f"{FUELS[fuel_index]:g} nM"
            ),
        )

    ax.set_title(
        row_labels[
            panel_index
        ]
    )

    ax.set_xlabel(
        "Time (min)"
    )

    ax.set_xlim(
        0,
        T_END_MIN,
    )

    ax.grid(
        True,
        alpha=0.3,
    )

axes[0].set_ylabel(
    "Predicted amplification"
)

axes[-1].legend(
    title="Fuel",
    loc="upper right",
)

fig.suptitle(
    f"Amplification prediction | "
    f"{MODEL_NAME} | "
    f"IN = {IN_CONC:g} nM"
)

output_file = resolve_path(
    config[
        "output_file"
    ]
)

output_file.parent.mkdir(
    parents=True,
    exist_ok=True,
)

fig.savefig(
    output_file,
    dpi=400,
    bbox_inches="tight",
)

print(
    f"Saved figure -> "
    f"{output_file}"
)

plt.show()
