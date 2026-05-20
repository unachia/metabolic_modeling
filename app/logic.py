from typing import Optional
import pandas as pd
from micom import Community

# ── AGORA model paths ──────────────────────────────────────────────────────────
AGORA_MODELS = {
    "A_muciniphila": "agora_files/Akkermansia_muciniphila_ATCC_BAA_835.xml",
    "F_prausnitzii": "agora_files/Faecalibacterium_prausnitzii_A2_165.xml",
    "A_caccae":      "agora_files/Anaerostipes_caccae_DSM_14662.xml",
    "E_hallii":      "agora_files/Eubacterium_hallii_DSM_3353.xml",
}

VALID_SPECIES = frozenset(AGORA_MODELS.keys())

DEFAULT_ABUNDANCES = {
    "A_muciniphila": 0.3,
    "F_prausnitzii": 0.25,
    "A_caccae":      0.25,
    "E_hallii":      0.2,
}

MUCIN_PROXIES = ["EX_acgal_m", "EX_acgam_m", "EX_fuc_L_m", "EX_gal_m", "EX_man_m"]
B12_KEYWORDS  = ["b12", "cobal", "cbl", "vitb", "cobalamin", "cobalt"]

EX_ACETATE     = "EX_ac_m"
EX_BUTYRATE    = "EX_but_m"
EX_PROPIONATE  = "EX_ppa_m"
EX_PROPANEDIOL = "EX_12ppd_S_m"

# ── Global state (populated at startup via lifespan) ──────────────────────────
community = None
baseline_medium = None


def init_community():
    """Build the MICOM community. Called once at app startup."""
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
    ids = [EX_ACETATE, EX_BUTYRATE, EX_PROPIONATE, EX_PROPANEDIOL]
    valid = [r for r in ids if r in fluxes.columns]
    if not valid:
        return {}
    df = fluxes[valid].copy()
    net = df.drop(index="medium", errors="ignore").sum(axis=0)
    return net.to_dict()


# ── Abundance helpers ─────────────────────────────────────────────────────────
def _resolve_abundances(user_input: Optional[dict[str, float]]) -> dict[str, float]:
    """Return a normalised abundance dict for all species.

    Missing species are filled with 0; values are renormalised to sum to 1.
    Falls back to DEFAULT_ABUNDANCES when user_input is None.
    """
    if user_input is None:
        return DEFAULT_ABUNDANCES.copy()

    raw = {sp: user_input.get(sp, 0.0) for sp in VALID_SPECIES}
    total = sum(raw.values())
    return {sp: v / total for sp, v in raw.items()}


# ── Main function called by the API ──────────────────────────────────────────
def run_metabolic_model(
    abundances: Optional[dict[str, float]] = None,
    low_mucin: bool    = False,
    low_b12: bool      = False,
    mucin_level: float = 0.001,
    b12_level: float   = 0.001,
    fraction: float    = 0.3,
) -> dict:
    if community is None:
        raise RuntimeError("Community not initialised. App may not have started correctly.")

    resolved = _resolve_abundances(abundances)
    community.abundances = pd.Series(resolved)

    med = _build_medium(low_mucin, low_b12, mucin_level, b12_level)
    community.medium = med

    sol = community.cooperative_tradeoff(fraction=fraction, pfba=True, fluxes=True)

    if sol is None:
        raise ValueError("MICOM returned no feasible solution for this condition.")

    net_flux = _extract_fluxes(sol.fluxes) if sol.fluxes is not None else {}

    return {
        "growth_rate":        sol.growth_rate,
        "net_community_flux": net_flux,
        "condition": {
            "abundances":  resolved,
            "low_mucin":   low_mucin,
            "low_b12":     low_b12,
            "mucin_level": mucin_level,
            "b12_level":   b12_level,
            "fraction":    fraction,
        }
    }
