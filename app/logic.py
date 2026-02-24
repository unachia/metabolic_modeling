import numpy as np
import pandas as pd

def run_metabolic_model(input_data: dict):
    """
    Core metabolic logic extracted from notebook
    """

    df = pd.DataFrame([input_data])

    score = df.mean(axis=1).iloc[0]

    return {
        "metabolic_score": float(score),
        "risk_category": "high" if score > 0.7 else "low"
    }
