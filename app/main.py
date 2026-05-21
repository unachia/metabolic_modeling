from contextlib import asynccontextmanager
from io import BytesIO

import matplotlib
matplotlib.use("Agg")  # non-interactive backend, safe for servers
import matplotlib.pyplot as plt
import numpy as np

from fastapi import FastAPI, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from app.logic import init_community, run_metabolic_model, METABOLITES


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


@app.post("/plot")
def plot_conditions():
    """Run all four conditions and return a bar chart as a PNG image."""
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
            raise HTTPException(status_code=500, detail=f"Error running {c['label']}: {str(e)}")

    # ── Build plot ────────────────────────────────────────────────────────────
    metabolite_names = list(METABOLITES.values())   # ["Acetate", "Butyrate", ...]
    condition_labels = [r["label"] for r in results]
    n_conditions     = len(results)
    n_metabolites    = len(metabolite_names)

    data = np.array([
        [r["net_community_flux"].get(name, 0) for name in metabolite_names]
        for r in results
    ])

    x      = np.arange(n_conditions)
    width  = 0.18
    colors = ["#378ADD", "#1D9E75", "#BA7517", "#D4537E"]

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (name, color) in enumerate(zip(metabolite_names, colors)):
        offset = (i - n_metabolites / 2 + 0.5) * width
        ax.bar(x + offset, data[:, i], width, label=name, color=color)

    ax.set_xticks(x)
    ax.set_xticklabels(condition_labels, rotation=15, ha="right")
    ax.set_ylabel("Net flux (mmol/gDW/h)")
    ax.set_title("Net community SCFA flux by condition")
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.legend()
    fig.tight_layout()

    # ── Return as PNG ─────────────────────────────────────────────────────────
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    buf.seek(0)

    return StreamingResponse(buf, media_type="image/png")
