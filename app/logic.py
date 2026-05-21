import pandas as pd
from micom import Community

# ── AGORA model paths ──────────────────────────────────────────────────────────
AGORA_MODELS = {
    "A_muciniphila": "agora_files/Akkermansia_muciniphila_ATCC_BAA_835.xml",
    "F_prausnitzii":  "agora_files/Faecalibacterium_prausnitzii_A2_165.xml",
    "A_caccae":       "agora_files/Anaerostipes_caccae_DSM_14662.xml",
    "E_hallii":       "agora_files/Eubacterium_hallii_DSM_3353.xml",
}

DEFAULT_ABUNDANCES = {
    "A_muciniphila": 0.30,
    "F_prausnitzii":  0.25,
    "A_caccae":       0.25,
    "E_hallii":       0.20,
}

METABOLITES = {
    "EX_ac_m":      "Acetate",
    "EX_but_m":     "Butyrate",
    "EX_ppa_m":     "Propionate",
    "EX_12ppd_S_m": "1,2-Propanediol",
}

MUCIN_PROXIES = ["EX_acgal_m", "EX_acgam_m", "EX_fuc_L_m", "EX_gal_m", "EX_man_m"]
B12_KEYWORDS  = ["b12", "cobal", "cbl", "vitb", "cobalamin", "cobalt"]
TRADEOFF      = 0.3

# ── Global state ───────────────────────────────────────────────────────────────
community       = None
baseline_medium = None


def init_community():
    """Build the MICOM community once at app startup."""
    global community, baseline_medium
    tax = pd.DataFrame({
        "id":        list(AGORA_MODELS.keys()),
        "file":      list(AGORA_MODELS.values()),
        "abundance": [DEFAULT_ABUNDANCES[s] for s in AGORA_MODELS.keys()],
    })
    community = Community(tax, id="crossfeeders")
    baseline_medium = community.medium.copy()
    print("MICOM community loaded successfully.")


# ── Helpers ───────────────────────────────────────────────────────────────────
def _build_medium(low_mucin: bool, low_b12: bool,
                  mucin_level: float, b12_level: float) -> dict:
    med = baseline_medium.copy()
    if low_mucin:
        for ex in MUCIN_PROXIES:
            if ex in med:
                med[ex] = mucin_level
    if low_b12:
        for k in list(med.keys()):
            if any(kw in k.lower() for kw in B12_KEYWORDS):
                med[k] = b12_level
    return med


def _extract_fluxes(fluxes: pd.DataFrame) -> dict:
    """Extract net SCFA fluxes from the medium row."""
    valid = [r for r in METABOLITES if r in fluxes.columns]
    if not valid:
        return {}
    if "medium" in fluxes.index:
        # medium row = net flux exchanged with the environment
        return fluxes.loc["medium", valid].rename(METABOLITES).to_dict()
    # fallback: sum species rows
    return fluxes[valid].sum(axis=0).rename(METABOLITES).to_dict()


def _apply_abundances(abundances: dict):
    """Validate and apply custom abundances to the community."""
    unknown = set(abundances) - set(AGORA_MODELS)
    if unknown:
        raise ValueError(f"Unknown species: {unknown}. Valid: {list(AGORA_MODELS.keys())}")
    total = sum(abundances.values())
    if not (0.99 <= total <= 1.01):
        raise ValueError(f"Abundances must sum to 1.0 (got {total:.3f}).")
    normalised = {k: v / total for k, v in abundances.items()}
    community.set_abundance(normalised)


# ── Main function called by the API ──────────────────────────────────────────
def run_metabolic_model(
    low_mucin: bool    = False,
    low_b12: bool      = False,
    mucin_level: float = 0.001,
    b12_level: float   = 0.001,
    fraction: float    = TRADEOFF,
    abundances: dict   = None,
) -> dict:
    if community is None:
        raise RuntimeError("Community not initialised.")

    if abundances:
        _apply_abundances(abundances)
    else:
        community.set_abundance(DEFAULT_ABUNDANCES)

    community.medium = _build_medium(low_mucin, low_b12, mucin_level, b12_level)

    sol = community.cooperative_tradeoff(fraction=fraction, pfba=True, fluxes=True)

    if sol is None:
        raise ValueError("MICOM returned no feasible solution for this condition.")

    net_flux = _extract_fluxes(sol.fluxes) if sol.fluxes is not None else {}

    return {
        "growth_rate":        sol.growth_rate,
        "net_community_flux": net_flux,
        "condition": {
            "low_mucin":   low_mucin,
            "low_b12":     low_b12,
            "mucin_level": mucin_level,
            "b12_level":   b12_level,
            "fraction":    fraction,
            "abundances":  abundances or DEFAULT_ABUNDANCES,
        }
    }
