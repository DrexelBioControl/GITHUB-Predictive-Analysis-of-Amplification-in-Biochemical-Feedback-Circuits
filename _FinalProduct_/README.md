# SBML Fitting and Amplification Prediction

This project fits several candidate SBML models to the same time-course data, compares their performance, and uses the selected fitted model to predict amplification.

## Table of contents

- [Big picture](#big-picture)
- [Setup and run](#setup-and-run)
- [Configure Code 1](#configure-code-1)
- [Configure Code 2](#configure-code-2)
- [Recommended workflow](#recommended-workflow)
- [Outputs](#outputs)
- [Important checks](#important-checks)

## Big picture

```text
Experimental data + candidate SBML models
                    |
                    v
       Code 1: fit and compare models
                    |
          selected model + fitted parameters
                    |
                    v
       Code 2: predict amplification
```

| File | Why it is needed |
|---|---|
| `models/*.xml` | Defines each candidate reaction mechanism and its parameters. |
| `data/testdata2D.xlsx` | Contains the experimental trajectories used for fitting and validation. |
| `configs/fit_all_models.json` | Tells Code 1 which data, conditions, species, parameters, and bounds to use. |
| `Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py` | Fits every candidate, predicts held-out conditions, compares models, and records the selected model. |
| `results/selected_model.json` | Connects the selected model and its fitted parameters to Code 2. It is created by Code 1. |
| `configs/amplification_prediction_general.json` | Describes each model's input, output, and amplification-control interface. |
| `Code2_MasterSBMLAmplificationPrediction_GeneralInterfaces.py` | Predicts how changing the amplification control affects the ON/OFF output-rate ratio. |
| `sbmltoodepy_solveivp_adapter.py` | Converts an SBML model into equations that the numerical solver can simulate. |

## Setup and run

The supplied environment uses Python 3.7 for compatibility with SBMLtoODEpy 1.0.4.

```powershell
conda env create -f environment.yml
conda activate SBMLtoODEpyWorkflow
python --version
```

Run the scripts from the `_FinalProduct_` folder:

```powershell
python .\Code1_MasterSBMLTimecourseFitting_WithHeldOutPredict.py
python .\Code2_MasterSBMLAmplificationPrediction_GeneralInterfaces.py
```

Run Code 2 only after Code 1 has produced `results/selected_model.json`.

## Configure Code 1

Edit `configs/fit_all_models.json`. JSON does not allow comments, so use the template below together with the field guide that follows it.

### Code 1 JSON template

This is a valid JSON example. Replace the example filenames, condition names, SBML IDs, and numerical values with values for your experiment.

```json
{
  "run_name": "example_model_comparison",
  "models_dir": "models",
  "excel_file": "data/example_data.xlsx",
  "excel_sheet": "Sheet1",
  "observable_id": "output",

  "fit_conditions": [
    "condition_1"
  ],
  "all_conditions": [
    "condition_1",
    "condition_2"
  ],
  "conditions": {
    "condition_1": {
      "data_column": "Condition 1",
      "initial_values": {
        "condition_species_1": 50.0,
        "fixed_species_1": 25.0,
        "fixed_species_2": 500.0
      }
    },
    "condition_2": {
      "data_column": "Condition 2",
      "initial_values": {
        "condition_species_1": 10.0,
        "fixed_species_1": 25.0,
        "fixed_species_2": 500.0
      }
    }
  },

  "fit_params": [
    "shared_rate"
  ],
  "all_param_central": {
    "shared_rate": 0.001
  },
  "custom_param_bounds": {
    "shared_rate": [0.000001, 0.1]
  },

  "model_fit_params": {
    "model_with_different_parameters": [
      "local_rate"
    ]
  },
  "model_param_central": {
    "model_with_different_parameters": {
      "local_rate": 0.002
    }
  },
  "model_param_bounds": {
    "model_with_different_parameters": {
      "local_rate": [0.000001, 0.1]
    }
  },
  "model_species_aliases": {
    "model_with_different_names": {
      "condition_species_1": "local_species_1",
      "fixed_species_1": "local_species_2",
      "fixed_species_2": "local_species_3"
    }
  },
  "model_observable_ids": {
    "model_with_different_names": "local_output"
  },

  "n_starts": 8,
  "random_seed": 20260806,
  "max_nfev": 400,
  "selection_metric": "held_out_normalized_rmse",
  "results_dir": "results",
  "figures_dir": "figures",
  "generated_models_dir": "generated_models"
}
```

### How to fill the Code 1 JSON

1. **Choose models and data.** Put candidate XML files in `models/`. Set `excel_file` and `excel_sheet`.
2. **Define every experimental condition.** Each key in `conditions` is your own condition label. `data_column` must match the experimental column. `initial_values` contains SBML species whose values Code 1 must set at time zero. These may be varied condition species, templates, pools, gates, reporters, fuels, or other fixed species; they are not automatically the analysis input or measured output.
3. **Choose fitted and held-out conditions.** Put fitting conditions in `fit_conditions`. Put fitting and prediction conditions in `all_conditions`. Anything only in `all_conditions` is predicted without refitting.
4. **Define shared parameters.** Put parameters used by most models in `fit_params`, their starting values in `all_param_central`, and positive lower/upper bounds in `custom_param_bounds`.
5. **Add model-specific parameter rules if needed.** Use `model_fit_params`, `model_param_central`, and `model_param_bounds` for different parameter IDs. `model_fixed_parameters` may set parameters that are not fitted.
6. **Map different species names.** The names in each condition are shared names. If one SBML uses different names, translate them with `model_species_aliases`. Set a different output with `model_observable_ids`. Do not rename the SBML file internally.
7. **Choose the ranking metric.** Use `held_out_normalized_rmse` for validation-based ranking or `global_normalized_rmse` for fitted-condition ranking.

Model-specific keys must exactly match the XML filename without `.xml`. For example, `model1_activegate_fuel.xml` uses `model1_activegate_fuel`.

The current configuration also contains advanced Excel-layout, cutoff, solver-tolerance, time-conversion, and output-directory settings. Keep their current values unless your data layout or model units require a change.

### Shared initial-species roles and model-specific names

Code 1 uses shared names for initial species so the same experimental condition can be applied to models with different ODE structures and SBML IDs.

| JSON section | Purpose |
|---|---|
| `conditions` | Assigns time-zero values to shared initial-species roles. |
| `model_species_aliases` | Translates each shared initial-species role into a model-specific SBML species ID. |
| `observable_id` | Gives the default output species. |
| `model_observable_ids` | Gives a model-specific output species when it differs from the default. |

```json
{
  "observable_id": "output",
  "conditions": {
    "condition_1": {
      "data_column": "Condition 1",
      "initial_values": {
        "condition_species_1": 50.0,
        "fixed_species_1": 25.0
      }
    },
    "condition_2": {
      "data_column": "Condition 2",
      "initial_values": {
        "condition_species_1": 10.0,
        "fixed_species_1": 25.0
      }
    }
  },
  "model_species_aliases": {
    "model_1": {
      "condition_species_1": "model1_species_name_1",
      "fixed_species_1": "model1_species_name_2"
    },
    "model_2": {
      "condition_species_1": "model2_species_name_1",
      "fixed_species_1": "model2_species_name_2"
    }
  },
  "model_observable_ids": {
    "model_1": "model1_output_name",
    "model_2": "model2_output_name"
  }
}
```

For `condition_2`, Code 1 sets `model1_species_name_1(0) = 10.0` in `model_1` and `model2_species_name_1(0) = 10.0` in `model_2`. It also sets each model's mapped fixed species to `25.0`.

The alias name on the left must exactly match a name in `initial_values`. The name on the right must exactly match an SBML species ID. The outer model name must exactly match the XML filename without `.xml`.

These entries set time-zero values for SBML species. For a dynamic species this is an ODE initial condition; for a constant species it represents a fixed model quantity. Biologically, the entries may be template or pool concentrations rather than the input signal or measured output. Only include species whose values Code 1 must override; other internal species should retain their SBML values.

'initial_values' contains shared names for SBML species whose time-zero values must be overridden. Use names that describe their experimental roles, then map them to exact model-specific SBML IDs with 'model_species_aliases'. Configure the measured output separately with 'observable_id' or 'model_observable_ids'.

Every shared initial-value role must be supported by each candidate model, either as the same SBML ID or through an alias. If only one candidate has an amplification variable, do not put that variable in the shared Code 1 conditions; define it later in that model's Code 2 interface.

### Fitting-effort settings

| Setting | Meaning |
|---|---|
| `n_starts` | Number of different starting parameter sets tried for each model. More starts reduce the chance of accepting a poor local optimum but increase runtime. |
| `random_seed` | Makes the generated starting sets reproducible. Keep it fixed when comparing models or repeating an analysis. The particular integer is arbitrary. |
| `max_nfev` | Maximum number of model evaluations allowed for each optimization start. Larger values give difficult fits more opportunity to converge but increase runtime. |

Runtime increases with both `n_starts` and `max_nfev`.

### Fit versus load

At the top of Code 1, use:

```python
MODE = "fit"
```

Use `fit` after changing data, models, parameters, or conditions. `load` is only for the exact same fitting setup. One run evaluates one held-out split; full leave-one-out validation requires repeated runs with a different condition omitted each time.

## Configure Code 2

Edit `configs/amplification_prediction_general.json` after the final Code 1 fit.

### Code 2 JSON template

```json
{
  "enabled": true,
  "selected_model_file": "results/selected_model.json",

  "model_interfaces": {
    "example_amplifier": {
      "input_control": {
        "type": "species_initial_concentration",
        "id": "input_species"
      },
      "observable_id": "output_species",
      "amplification_control": {
        "type": "species_initial_concentration",
        "id": "fuel_species"
      },
      "fixed_initial_species": {
        "reporter_species": 500.0
      },
      "fixed_parameters": {}
    }
  },

  "input_values": {
    "on_value": 50.0,
    "off_value": 0.0,
    "label": "Input",
    "units": "nM"
  },
  "amplification_values": {
    "reference_value": 0.0,
    "values": {
      "start": 0.0,
      "stop": 50.0,
      "points": 51
    },
    "values_to_plot": [1.0, 5.0, 10.0, 50.0],
    "label": "Fuel template",
    "units": "nM"
  },

  "time_end_min": 200.0,
  "time_points": 1000,
  "model_time_units_per_minute": 60.0,
  "sweep_parameter": null,
  "sweep_values": null,

  "output_file": "figures/amplification_prediction.png",
  "generated_models_dir": "generated_models"
}
```

### How to fill the Code 2 JSON

1. **Add one interface per possible selected model.** Use the XML filename without `.xml` as the key.
2. **Identify the three roles.** `input_control` defines ON/OFF, `observable_id` is the measured output, and `amplification_control` is the variable being tested.
3. **Choose each control type.** Use `species_initial_concentration` when changing a starting species concentration, or `parameter` when changing an SBML parameter.
4. **Set fixed values.** Put starting concentrations in `fixed_initial_species`; optional parameter overrides go in `fixed_parameters`.
5. **Define the experiment.** Set ON/OFF values, the amplification reference, the tested range, and the simulation duration. `amplification_values.values` may be a `start`/`stop`/`points` object or a direct list of numbers.
6. **Use the optional sensitivity sweep if needed.** JSON `null` disables it. Otherwise, set `sweep_parameter` to an exact SBML parameter ID and provide `sweep_values`. This tests different values during prediction; it does not refit the parameter.

A model without an amplification variable should use `"amplification_control": null`. Code 1 can still fit it, but Code 2 cannot vary amplification for it.

## Recommended workflow

1. Fit candidate models using a held-out condition and rank with `held_out_normalized_rmse`.
2. Inspect the fitted and held-out plots, error tables, and parameter-boundary checks.
3. Repeat with other held-out conditions if full leave-one-out validation is needed.
4. Select the model structure.
5. Run a final Code 1 fit using all conditions to obtain final parameter estimates.
6. Confirm that the selected model has a Code 2 interface and an amplification control.
7. Run Code 2.

## Outputs

### Code 1

| Output | Meaning |
|---|---|
| `results/model_comparison.csv` | Scores and ranks for all candidate models. |
| `results/selected_model.json` | Selected model and path to its fitted parameters. |
| `results/selected_model.xml` | Copy of the selected SBML structure. |
| `results/<model>/fitted_parameters.csv` | Readable fitted and fixed parameter values. |
| `results/<model>/<model>_results.npz` | Machine-readable fitted results used by Code 2. |
| `results/<model>/held_out_rmse.csv` | Prediction error for held-out conditions. |
| `figures/<model>/<model>_fit.png` | Fitted-condition plot. |
| `figures/<model>/<model>_held_out_prediction.png` | Held-out prediction plot. |

The copied `selected_model.xml` retains the SBML's original parameter values. The fitted values are stored in the `.npz` and CSV files.

### Code 2

Code 2 saves the amplification figure to the path in `output_file`, currently `figures/amplification_prediction_general_interface.png`.

## Important checks

- Species and parameter IDs in both JSON files must exactly match their SBML IDs after any configured alias mapping.
- All fitted parameter bounds must satisfy `0 < lower < upper` because fitting is performed in log space.
- A held-out plot is only produced when `all_conditions` contains a condition absent from `fit_conditions`.
- A parameter at its fitting bound may be poorly identified or may need a wider scientifically justified bound.
- When an XML file changes, delete only its matching file in `generated_models/` and rerun Code 1. Do not edit generated model files manually.
- The adapter supports the SBML features used here, but not arbitrary assignment rules, rate rules, algebraic rules, or events.
