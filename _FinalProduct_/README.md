# SBMLtoODEpy Fitting and Prediction Workflow

This project has two main parts:

1. `Code1_MasterSBMLTimecourseFitting.py` fits candidate SBML models to experimental time-course data and selects the best-performing model.
2. `Code2_MasterSBMLAmplificationPrediction.py` loads the selected model and fitted parameters for a design-space prediction.

The main settings are stored in JSON files so models, parameters, fitting ranges, and simulation settings can be changed without rewriting the Python scripts.

## Code 1: Time-course fitting

Code 1:

- Reads the experimental data from Excel.
- Finds the candidate SBML models in `models/`.
- Reads the experimental conditions and fitting settings from `configs/fit_all_models.json`.
- Automatically finds the transient region of each curve to fit.
- Reads the fitted parameters and their allowed ranges from the JSON file.
- Generates Sobol samples within the parameter ranges.
- Tests each sample as a starting point for `least_squares`.
- Keeps the fit with the lowest cost.
- Calculates RMSE values and checks whether fitted parameters are near their bounds.
- Saves plots, fitted parameters, error values, and model-specific result files.
- Compares the candidate models and saves the selected model for Code 2.

## Code 2: Amplification prediction

Code 2:

- Loads the model and fitted parameters selected by Code 1.
- Reads its settings from `configs/amplification_prediction.json`.
- Runs input-ON and input-OFF simulations.
- Tests a range of initial fuel concentrations.
- Calculates reporter-production rates and normalized amplification.
- Can repeat the calculation for several values of another parameter, such as a leak-rate parameter.
- Saves the resulting amplification plots.

The current Figure 2D models do not contain a genuine fuel-driven amplification mechanism. Placeholder species and parameters can be used to test the software workflow, but those results should not be interpreted biologically.

## Python compatibility

SBMLtoODEpy 1.0.4 is an older package that supports Python 3.5–3.7.

Create and activate the supplied environment:

```powershell
conda env create -f environment.yml
conda activate SBMLtoODEpyWorkflow
```

Confirm that Python 3.7 is active:

```powershell
python --version
```

## Numerical process

```text
SBML model
    |
    v
SBMLtoODEpy generates Python ODE code
    |
    v
The adapter changes _SolveReactions(y, t) to rhs(t, y)
    |
    v
solve_ivp simulates the model using LSODA
    |
    v
least_squares tests and fits parameter values
```

The adapter currently expects species controlled by reactions and does not yet support more complex features such as rate rules, algebraic rules, or events.

## Project structure

```text
_FinalProduct_/
|-- Code1_MasterSBMLTimecourseFitting.py
|-- Code2_MasterSBMLAmplificationPrediction.py
|-- sbmltoodepy_solveivp_adapter.py
|-- environment.yml
|
|-- configs/
|   |-- fit_all_models.json
|   `-- amplification_prediction.json
|
|-- data/
|   `-- testdata2D.xlsx
|
|-- models/
|   |-- Figure2D_base_model.xml
|   `-- Figure2D_base_model_decay.xml
|
|-- generated_models/
|-- results/
`-- figures/
```

## Running Code 1

If the SBML files have changed, delete the old generated models:

```powershell
Remove-Item -Recurse -Force .\generated_models
```

Then run:

```powershell
python .\Code1_MasterSBMLTimecourseFitting.py
```

Code 1 writes generated model files to `generated_models/`.

Do not edit these generated files manually. Make changes in the original SBML files and rerun Code 1.

## Running Code 2

After Code 1 successfully creates the selected-model files, run:

```powershell
python .\Code2_MasterSBMLAmplificationPrediction.py
```

If the amplification configuration contains:

```json
"enabled": false
```

Code 2 will print an explanation and exit without running the prediction.

## Important checks

Before using the results, confirm that:

- The fitted curves follow the experimental data.
- The RMSE values are finite.
- The program does not repeatedly return solver errors or penalty values.
- Fitted parameters are not all stuck at their search bounds.
- Generated models were recreated after SBML changes.
- Code 2 uses genuine fuel and leak components before its results are interpreted biologically.