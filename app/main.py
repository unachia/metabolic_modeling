from contextlib import asynccontextmanager
from typing import Optional
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field, model_validator
from app.logic import VALID_SPECIES, init_community, run_metabolic_model


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
    abundances:  Optional[dict[str, float]] = Field(
        None,
        description=(
            "Species relative abundances (fractions). "
            f"Valid keys: {sorted(VALID_SPECIES)}. "
            "Missing species default to 0. Values are renormalised to sum to 1. "
            "Omit to use default abundances."
        ),
    )
    low_mucin:   bool  = Field(False, description="Reduce mucin sugar availability")
    low_b12:     bool  = Field(False, description="Reduce B12/cobalt availability")
    mucin_level: float = Field(0.001, ge=0, description="Mucin exchange rate when low_mucin=True")
    b12_level:   float = Field(0.001, ge=0, description="B12 exchange rate when low_b12=True")
    fraction:    float = Field(0.3, ge=0, le=1, description="Cooperative tradeoff fraction")

    @model_validator(mode="after")
    def validate_abundances(self):
        if self.abundances is None:
            return self
        unknown = set(self.abundances) - VALID_SPECIES
        if unknown:
            raise ValueError(f"Unknown species: {sorted(unknown)}. Valid: {sorted(VALID_SPECIES)}")
        if any(v < 0 for v in self.abundances.values()):
            raise ValueError("Abundance values must be non-negative.")
        total = sum(self.abundances.values())
        if total == 0:
            raise ValueError("At least one abundance value must be > 0.")
        return self


# ── Endpoints ─────────────────────────────────────────────────────────────────
@app.get("/")
def root():
    return {"status": "running", "docs": "/docs"}


@app.post("/run-model")
def run_model(request: ModelRequest):
    """Run MICOM cooperative_tradeoff for a given condition."""
    try:
        return run_metabolic_model(
            abundances=request.abundances,
            low_mucin=request.low_mucin,
            low_b12=request.low_b12,
            mucin_level=request.mucin_level,
            b12_level=request.b12_level,
            fraction=request.fraction,
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
