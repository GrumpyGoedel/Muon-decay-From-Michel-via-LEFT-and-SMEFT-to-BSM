r"""Experimental inputs and the block-covariance chi^2.

Implements the likelihood of ``main.tex`` section "Experimental Michel-parameter
likelihood" (\label{sec:michel-likelihood}): the observables are grouped into
mutually independent blocks, with the experimental correlations *within* each
block retained.

    O_fit = (rho, delta, P_mu^pi xi)              [TWIST, 3x3 correlated]
          + (eta, eta'', alpha'/A, beta'/A)       [Danneberg, 4x4 correlated]
          + xi'                                    [Burkard, 1D]
          + xi''                                   [Prieels, 1D]

    chi2 = sum_b [O^th_b - O^exp_b]^T V_b^{-1} [O^th_b - O^exp_b]

Normalization: the Michel shape/polarization parameters do NOT constrain the
overall amplitude normalization.  Because we also study normalization-changing
coefficients (e.g. x^VLL_2112), we add a separate G_F block as the "equivalent
normalization constraint" the thesis calls for.  Use MICHEL_BLOCKS for the pure
shape/polarization likelihood without it.

P_mu^pi = 1 is assumed (a phenomenological choice, not part of the measurement).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List

import numpy as np
from scipy.linalg import block_diag

from . import michel
from .coefficients import Coefficients

P_MU = 1.0  # assumed muon polarisation (P_mu^pi = 1)


@dataclass
class Block:
    """A correlated block of observables with its own covariance matrix."""
    name: str
    labels: List[str]
    models: List[Callable[[Coefficients], float]]
    exp: np.ndarray
    cov: np.ndarray
    note: str = ""
    cov_inv: np.ndarray = field(init=False)

    def __post_init__(self):
        self.exp = np.asarray(self.exp, dtype=float)
        self.cov = np.atleast_2d(np.asarray(self.cov, dtype=float))
        self.cov_inv = np.linalg.inv(self.cov)

    def theory(self, x: Coefficients) -> np.ndarray:
        return np.array([m(x) for m in self.models], dtype=float)

    def residual(self, x: Coefficients) -> np.ndarray:
        return self.theory(x) - self.exp

    def chi2(self, x: Coefficients) -> float:
        r = self.residual(x)
        return float(r @ self.cov_inv @ r)


def _cov_from_corr(sigmas, corr):
    D = np.diag(np.asarray(sigmas, dtype=float))
    return D @ np.asarray(corr, dtype=float) @ D


# --------------------------------------------------------------------------
# TWIST block: (rho, delta, P_mu^pi xi)          main.tex eq:twist-*
# --------------------------------------------------------------------------
_sig_rho = np.hypot(0.00012, 0.00023)
_sig_delta = np.hypot(0.00021, 0.00027)
_sig_pxi = np.hypot(0.00029, (0.00165 + 0.00063) / 2)   # symmetrised
_R_TWIST = [[1.0, 0.19, 0.21],
            [0.19, 1.0, -0.72],
            [0.21, -0.72, 1.0]]

TWIST = Block(
    name="TWIST",
    labels=["rho", "delta", "Pmu_xi"],
    models=[michel.rho, michel.delta, lambda x: P_MU * michel.xi(x)],
    exp=[0.74977, 0.75049, 1.00084],
    cov=_cov_from_corr([_sig_rho, _sig_delta, _sig_pxi], _R_TWIST),
    note="native (rho, delta, P_mu^pi xi) basis; stat+syst in quadrature",
)

# --------------------------------------------------------------------------
# Danneberg block: (eta, eta'', alpha'/A, beta'/A)   main.tex eq:danneberg-*
# --------------------------------------------------------------------------
_D_DANN = [0.03734, 0.05235, 0.02186, 0.00800]
_R_DANN = [[1.0, 0.946, 0.0, 0.0],
           [0.946, 1.0, 0.0, 0.0],
           [0.0, 0.0, 1.0, -0.893],
           [0.0, 0.0, -0.893, 1.0]]

DANNEBERG = Block(
    name="Danneberg",
    labels=["eta", "eta''", "alpha'/A", "beta'/A"],
    models=[michel.eta, michel.eta_double_prime,
            michel.alpha_prime, michel.beta_prime],
    exp=[0.071, 0.105, -0.0034, -0.0005],
    cov=_cov_from_corr(_D_DANN, _R_DANN),
    note="eta-eta'' correlated (0.946); alpha'/A-beta'/A correlated (-0.893). "
         "beta'/A model uses the full-4-fermion sign (see michel.py).",
)

# --------------------------------------------------------------------------
# Single-observable blocks: xi' (Burkard), xi'' (Prieels)
# --------------------------------------------------------------------------
XI_PRIME = Block(
    name="Burkard", labels=["xi'"], models=[michel.xi_prime],
    exp=[0.998], cov=[[0.045 ** 2]], note="xi' = 0.998 +/- 0.045",
)
XI_DOUBLE_PRIME = Block(
    name="Prieels", labels=["xi''"], models=[michel.xi_double_prime],
    exp=[0.981], cov=[[0.0451 ** 2]], note="xi'' = 0.981 +/- 0.0451",
)

# --------------------------------------------------------------------------
# G_F normalization block (NOT part of the Michel likelihood; added because we
# include normalization-changing coefficients).   See main.tex L750-756.
#   (G_F,exp/G_F,SM)^2 = 1.00099, sigma from G_F,SM = 1.1658(10)e-5.
# --------------------------------------------------------------------------
_GF_EXP, _GF_SM, _GF_SM_ERR = 1.1663787e-5, 1.1658e-5, 0.0010e-5
_GF_RATIO2 = (_GF_EXP / _GF_SM) ** 2
_GF_RATIO2_ERR = 2.0 * (_GF_SM_ERR / _GF_SM) * _GF_RATIO2

GF_NORM = Block(
    name="GF_norm", labels=["GF2_ratio"], models=[michel.gf_ratio_sq],
    exp=[_GF_RATIO2], cov=[[_GF_RATIO2_ERR ** 2]],
    note="normalization constraint; (G_F,exp/G_F,SM)^2",
)

# Pure Michel shape/polarization likelihood (thesis eq:michel-chi2)
MICHEL_BLOCKS = [TWIST, DANNEBERG, XI_PRIME, XI_DOUBLE_PRIME]
# Default fit: Michel likelihood + normalization constraint
DEFAULT_BLOCKS = MICHEL_BLOCKS + [GF_NORM]


def chi2(x: Coefficients, blocks=None) -> float:
    blocks = DEFAULT_BLOCKS if blocks is None else blocks
    return sum(b.chi2(x) for b in blocks)


def chi2_sm(blocks=None) -> float:
    return chi2(Coefficients(), blocks)


# --- helpers for the Jacobian / flat-direction analysis --------------------

def flat_models(blocks=None):
    """[(name, model_callable), ...] over every scalar observable."""
    blocks = DEFAULT_BLOCKS if blocks is None else blocks
    out = []
    for b in blocks:
        for lab, m in zip(b.labels, b.models):
            out.append((f"{b.name}:{lab}", m))
    return out


def full_covariance(blocks=None) -> np.ndarray:
    """Block-diagonal covariance over all scalar observables."""
    blocks = DEFAULT_BLOCKS if blocks is None else blocks
    return block_diag(*[b.cov for b in blocks])
