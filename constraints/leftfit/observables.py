r"""Experimental inputs and the observable model.

Each observable maps a set of Wilson coefficients to the *actually measured*
quantity (e.g. the PDG reports P_mu*xi and xi*P_mu*delta/rho, not xi and delta
in isolation), so theory-side correlations between Michel parameters are handled
automatically.  We assume P_mu = 1.

Sources: thesis ``main.tex`` section "The Michel parameters" and
``particle2024review`` (PDG 2024).

  !! Values marked PLACEHOLDER must be replaced with the numbers/covariance you
     actually adopt in the thesis before quoting any bound.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from . import michel
from .coefficients import Coefficients

P_MU = 1.0  # assumed muon polarisation


@dataclass
class Observable:
    name: str
    model: Callable[[Coefficients], float]   # theory prediction
    exp: float                               # measured central value
    sigma: float                             # 1-sigma uncertainty
    sm: float                                # Standard-Model value
    note: str = ""


# ---------------------------------------------------------------------------
# G_F: theory-vs-experiment on (G_F/G_F,SM)^2.
#   G_F,exp = 1.1663787e-5   (muon lifetime, PDG; negligible error)
#   G_F,SM  = 1.1658(10)e-5  (thesis, {m_W,m_Z,alpha} scheme)
#   => (G_F,exp/G_F,SM)^2 = 1.00099, sigma ~ 2*(1e-3/1.1658) = 1.72e-3
# This reproduces the thesis band  -0.004 <= v^2 Re[L^VLL_2112] <= 0.002.
# ---------------------------------------------------------------------------
_GF_EXP = 1.1663787e-5
_GF_SM = 1.1658e-5
_GF_SM_ERR = 0.0010e-5
_GF_RATIO2 = (_GF_EXP / _GF_SM) ** 2
_GF_RATIO2_ERR = 2.0 * (_GF_SM_ERR / _GF_SM) * _GF_RATIO2


DEFAULT_OBSERVABLES = [
    Observable(
        "GF2_ratio", michel.gf_ratio_sq,
        exp=_GF_RATIO2, sigma=_GF_RATIO2_ERR, sm=1.0,
        note="(G_F,exp/G_F,SM)^2; sigma from G_F,SM error",
    ),
    Observable(
        "rho", michel.rho,
        exp=0.74979, sigma=0.00026, sm=0.75,
        note="PDG global fit",
    ),
    Observable(
        "eta", michel.eta,
        exp=0.057, sigma=0.034, sm=0.0,
        note="PLACEHOLDER sigma -- replace with adopted value",
    ),
    Observable(
        "Pmu_xi", lambda x: P_MU * michel.xi(x),
        exp=1.0009, sigma=0.0012, sm=1.0,
        note="P_mu * xi",
    ),
    Observable(
        "xi_Pmu_delta_over_rho",
        lambda x: P_MU * michel.xi(x) * michel.delta(x) / michel.rho(x),
        exp=1.00179, sigma=0.00114, sm=1.0,
        note="xi * P_mu * delta / rho",
    ),
]


def chi2(x: Coefficients, observables=None) -> float:
    """Diagonal (uncorrelated) chi^2.  Swap in a covariance matrix here once
    the PDG correlation coefficients are entered."""
    obs = observables if observables is not None else DEFAULT_OBSERVABLES
    total = 0.0
    for o in obs:
        total += ((o.model(x) - o.exp) / o.sigma) ** 2
    return total


def chi2_sm(observables=None) -> float:
    """chi^2 at the Standard-Model point (all coefficients zero)."""
    return chi2(Coefficients(), observables)
