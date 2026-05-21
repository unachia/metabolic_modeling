from contextlib import asynccontextmanager
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from app.logic import init_community, run_metabolic_model


# ── Startup / shutdown ────────────────────────────────────────────────────────
@asynccontextmanager
async def lifespan(app: FastAPI):
    """Build the MICOM community once when the server starts."""
    print("Loading MICOM community... (this may take a moment)")
    init_community()
    yield
    print("Shutting down.")


app = FastAPI(
    title="Metabolic Model API",
    description="REST API for MICOM-based microbial community metabolic modelling.",
    version="1.0.0",
    lifespan=lifespan,
)


# ── Input schema ──────────────────────────────────────────────────────────────
class ModelRequest(BaseModel):
    low_mucin:   bool  = Field(False, description="Reduce mucin sugar availability")
    low_b12:     bool  = Field(False, description="Reduce B12/cobalt availability")
    mucin_level: float = Field(0.001, ge=0, description="Mucin exchange rate when low_mucin=True")
    b12_level:   float = Field(0.001, ge=0, description="B12 exchange rate when low_b12=True")
    fraction:    float = Field(0.3, ge=0, le=1, description="Cooperative tradeoff fraction")
    abundances:  dict  = Field(
        None,
        description=(
            "Custom species abundances (must sum to 1). "
            "Valid keys: A_muciniphila, F_prausnitzii, A_caccae, E_hallii. "
            "Example: {\"A_muciniphila\": 0.4, \"F_prausnitzii\": 0.3, \"A_caccae\": 0.2, \"E_hallii\": 0.1}"
        )
    )


# ── Endpoints ─────────────────────────────────────────────────────────────────
@app.get("/")
def root():
    return {"status": "running", "docs": "/docs"}


@app.post("/run-model")
def run_model(request: ModelRequest):
    """Run MICOM cooperative_tradeoff for a given condition."""
    try:
        return run_metabolic_model(
            low_mucin=request.low_mucin,
            low_b12=request.low_b12,
            mucin_level=request.mucin_level,
            b12_level=request.b12_level,
            fraction=request.fraction,
            abundances=request.abundances,
        )
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Model error: {str(e)}")


@app.post("/run-all-conditions")
def run_all_conditions():
    """Run all four standard conditions and return a comparison."""
    conditions = [
        {"label": "Baseline",            "low_mucin": False, "low_b12": False},
        {"label": "Low B12",             "low_mucin": False, "low_b12": True},
        {"label": "Low mucin",           "low_mucin": True,  "low_b12": False},
        {"label": "Low mucin + Low B12", "low_mucin": True,  "low_b12": True},
    ]
    results = []
    for c in conditions:
        try:
            r = run_metabolic_model(low_mucin=c["low_mucin"], low_b12=c["low_b12"])
            r["label"] = c["label"]
            results.append(r)
        except Exception as e:
            results.append({"label": c["label"], "error": str(e)})
    return {"conditions": results}
