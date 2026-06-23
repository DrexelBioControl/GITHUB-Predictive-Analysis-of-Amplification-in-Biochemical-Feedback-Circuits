# ============================================================
# Plot 0 input / 0 fuel condition with predictions from 3 models
# ============================================================

import json
import importlib.util
from pathlib import Path

import xlrd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ============================================================
# 0) PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parent

CONFIG_FILES = [
    BASE_DIR / "configs" / "fit_exps12345_model0.json",
    BASE_DIR / "configs" / "fit_exps12345_model1.json",
    BASE_DIR / "configs" / "fit_exps12345_model2.json",
]

MODEL_FILES = {
    "model0_transcriptional": BASE_DIR / "models" / "model0_transcriptional.py",
    "model1_activegate_fuel": BASE_DIR / "models" / "model1_activegate_fuel.py",
    "model2_inactivegate_fuel": BASE_DIR / "models" / "model2_inactivegate_fuel.py",
}

MODEL_LABELS = {
    "model0_transcriptional": "Model 0: transcriptional leak",
    "model1_activegate_fuel": "Model 1: active-gate fuel leak",
    "model2_inactivegate_fuel": "Model 2: inactive-gate fuel leak",
}


# ============================================================
# 1) HELPER FUNCTIONS
# ============================================================

def load_module_from_file(module_name, file_path):
    """Load a Python model file as a module."""
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolve_path(path_from_config):
    """Resolve paths written relative to the project root."""
    p = Path(path_from_config)
    if p.is_absolute():
        return p
    return BASE_DIR / p


def find_results_file(config):
    """Find saved fit result file."""
    results_dir = resolve_path(config.get("results_dir", "results"))
    run_name = config["run_name"]

    candidates = [
        results_dir / f"{run_name}_results.npz",
        results_dir / f"{run_name}.npz",
    ]

    if "fit_results_file" in config:
        candidates.append(results_dir / config["fit_results_file"])

    for p in candidates:
        if p.exists():
            return p

    raise FileNotFoundError(
        "Could not find results file. Tried:\n"
        + "\n".join(str(p) for p in candidates)
    )


def load_best_fit_params(config):
    """Load best-fit parameters from saved .npz file."""
    results_file = find_results_file(config)
    data = np.load(results_file, allow_pickle=True)

    best_x = data["best_x"]
    fit_params = list(config["fit_params"])

    params = {
        p: float(10 ** best_x[i])
        for i, p in enumerate(fit_params)
    }

    return params, results_file


def build_control_dataset(ws):
    """
    Build 0 input / 25 fuel dataset.

    Excel layout:
    column 0 = time
    column 4 = 0 input / 25 fuel normalized reporter signal
    """

    START_ROW = 3

    t_min = np.array(
        [ws.cell_value(r, 0) for r in range(START_ROW, ws.nrows)],
        dtype=float
    )

    y_norm = np.array(
        [ws.cell_value(r, 4) for r in range(START_ROW, ws.nrows)],
        dtype=float
    )

    # Convert normalized fraction reacted to nM using nominal reporter amount
    y_nM = 500.0 * y_norm

    # Estimate effective DRL_0 from last 20% of trace for visualization
    last_20_percent_start = int(0.8 * len(y_nM))
    DRL_0 = float(np.mean(y_nM[last_20_percent_start:]))

    return {
        "name": "0in_25fuel",
        "label": "Input 0 nM / Fuel 25 nM",
        "t_min": t_min,
        "x_nM": y_nM,
        "IN_conc": 0.0,
        "Fuel_conc": 25.0,
        "DRL_0": DRL_0,
        "RSD_temp": 25.0,
    }


def predict(model, dataset, params):
    """Simulate one model for the 0 input / 25 fuel condition."""
    t_eval_sec = dataset["t_min"] * 60.0

    sol = solve_ivp(
        lambda t, y: model.rhs(
            t,
            y,
            RSD_temp=dataset["RSD_temp"],
            IN_temp=dataset["IN_conc"],
            F_temp=dataset["Fuel_conc"],
            params=params,
        ),
        (t_eval_sec[0], t_eval_sec[-1]),
        model.initial_conditions(dataset, params),
        t_eval=t_eval_sec,
        method="LSODA",
        rtol=1e-6,
        atol=1e-6,
    )

    if not sol.success:
        print(f"WARNING: solver failed: {sol.message}")
        return np.full_like(dataset["t_min"], np.nan, dtype=float)

    return sol.y[model.output_index]


# ============================================================
# 2) LOAD MODELS AND BEST-FIT PARAMETERS
# ============================================================

model_runs = []

for config_file in CONFIG_FILES:

    with open(config_file, "r") as f:
        config = json.load(f)

    model_name = config["model_name"]
    model_file = MODEL_FILES[model_name]

    model = load_module_from_file(model_name, model_file)
    params, results_file = load_best_fit_params(config)

    model_runs.append({
        "model_name": model_name,
        "label": MODEL_LABELS[model_name],
        "model": model,
        "params": params,
        "results_file": results_file,
    })

    print(f"Loaded {model_name}")
    print(f"  Results file: {results_file}")
    for p, v in params.items():
        print(f"  {p:<12} = {v:.4e}")
    print()


# ============================================================
# 3) LOAD 0 INPUT / 0 FUEL EXPERIMENTAL DATA
# ============================================================

with open(CONFIG_FILES[0], "r") as f:
    base_config = json.load(f)

excel_file = resolve_path(base_config["excel_file"])
excel_sheet = base_config["excel_sheet"]

wb = xlrd.open_workbook(str(excel_file))
ws = wb.sheet_by_name(excel_sheet)

dataset = build_control_dataset(ws)

print(f"Control condition DRL_0 estimate = {dataset['DRL_0']:.2f} nM")


# ============================================================
# 4) PLOT EXPERIMENT AND MODEL PREDICTIONS
# ============================================================

plt.figure(figsize=(8, 5))

# Experimental data
plt.scatter(
    dataset["t_min"],
    dataset["x_nM"],
    s=12,
    alpha=0.55,
    color="gray",
    label="Experimental data"
)

# Model predictions
for run in model_runs:
    y_pred = predict(run["model"], dataset, run["params"])

    plt.plot(
        dataset["t_min"],
        y_pred,
        linewidth=2,
        label=run["label"]
    )

plt.xlabel("Time (min)")
plt.ylabel("Reacted reporter (nM)")
plt.title("0 input / 25 fuel control condition")
plt.grid(True)
plt.legend(fontsize=8)
plt.tight_layout()
plt.show()