from fastapi import FastAPI
from app.logic import run_metabolic_model

app = FastAPI(
    title="Metabolic Model API",
    description="API wrapper around Jupyter notebook logic",
    version="0.1.0"
)

@app.post("/run-model")
def run_model(input_data: dict):
    result = run_metabolic_model(input_data)
    return result