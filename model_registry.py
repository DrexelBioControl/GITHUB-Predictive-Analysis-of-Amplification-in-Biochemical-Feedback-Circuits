# ============================================================
# MODEL REGISTRY
#
# The config file gives a string, e.g.
#   "model_name": "model1_activegate_fuel"
#
# This file maps that string to the correct Python model module.
# ============================================================

from models import model_noleak
from models import model0_transcriptional
from models import model1_activegate_fuel
from models import model2_inactivegate_fuel
from models import model3_misfoldgate_fuel


MODEL_REGISTRY = {
    "model_noleak": model_noleak,
    "model0_transcriptional": model0_transcriptional,
    "model1_activegate_fuel": model1_activegate_fuel,
    "model2_inactivegate_fuel": model2_inactivegate_fuel,
    "model3_misfoldgate_fuel": model3_misfoldgate_fuel,
}


def get_model(model_name):
    if model_name not in MODEL_REGISTRY:
        raise ValueError(
            f"Unknown model_name: {model_name}. "
            f"Choose from: {list(MODEL_REGISTRY.keys())}"
        )

    return MODEL_REGISTRY[model_name]