# SBMLtoODEpy + solve_ivp workflow

This version preserves the previous two-part architecture:

1. `Code1_MasterSBMLTimecourseFitting.py` converts every `models/*.xml` file to generated Python ODE code, fits each candidate with Sobol multistart `least_squares`, compares normalized RMSE, and selects the best model.
2. `Code2_MasterSBMLAmplificationPrediction.py` loads the selected model and fit, then uses the generated ODE with `solve_ivp(method="LSODA")` for design-space prediction.

## Important compatibility note

SBMLtoODEpy 1.0.4 is an old package and officially lists Python 3.5-3.7 support. Use the supplied environment rather than the newer `LeakyAmplification` environment:

```powershell
conda env create -f environment.yml
conda activate SBMLtoODEpyWorkflow
python Code1_MasterSBMLTimecourseFitting.py
```

Generated model files are written automatically to `generated_models/`. Do not edit them manually; edit the SBML source and rerun instead.

## Numerical architecture

```text
SBML -> SBMLtoODEpy generated _SolveReactions(y,t)
     -> adapter changes argument order to rhs(t,y)
     -> scipy.integrate.solve_ivp(method="LSODA")
     -> scipy.optimize.least_squares
```

The adapter currently expects reaction-governed species without rate rules, events, or algebraic rules. This matches the two present Figure 2D candidate models. New SBML models should be validated against COPASI before being used for scientific conclusions.
