"""Bookkeeping for the LEFT Wilson coefficients that enter muon decay.

We work throughout with the *dimensionless* coefficients

    x^X_{ab,st}  =  v^2 * c^{(6),X}_{nu e, ab st}          (v = 246 GeV)

so that x ~ (v/Lambda)^2 for an order-one UV coupling.  A coefficient is
labelled by

    X   in {'VLL','VLR','SLL','SLR','TLL'}   Dirac/chirality structure
    a,b in {1,2,3}                            neutrino generations (nu_a, nu_b)
    st  in {'12','21'}                        charged-lepton chirality slot
                                              ('12' = SM-like e/mu assignment)

Only the two vector structures 'VLL','VLR' conserve lepton number; the scalar
and tensor structures are lepton-number violating and therefore appear only
quadratically in the observables.
"""

from __future__ import annotations

STRUCTURES = ("VLL", "VLR", "SLL", "SLR", "TLL")
LNC = ("VLL", "VLR")
LNV = ("SLL", "SLR", "TLL")
FLAVORS = (1, 2, 3)
ST = ("12", "21")

VEW = 246.0  # electroweak vev in GeV, used only to convert |x| -> Lambda


class Coefficients(dict):
    """A sparse map (X, a, b, st) -> complex value; missing entries are 0."""

    def val(self, X: str, a: int, b: int, st: str) -> complex:
        return complex(self.get((X, a, b, st), 0.0))

    def set(self, X: str, a: int, b: int, st: str, value: complex) -> "Coefficients":
        self[(X, a, b, st)] = complex(value)
        return self

    @classmethod
    def single(cls, key, value: complex) -> "Coefficients":
        """A coefficient set with exactly one entry switched on."""
        c = cls()
        c[key] = complex(value)
        return c


def lambda_from_x(x_abs: float) -> float:
    """Convert a bound on |x| = v^2|c| to a new-physics scale Lambda (GeV),
    assuming an order-one coupling so |c| ~ 1/Lambda^2  =>  Lambda = v/sqrt|x|."""
    if x_abs <= 0:
        return float("inf")
    return VEW / (x_abs ** 0.5)


def all_keys(structures=STRUCTURES, flavors=FLAVORS, st=ST):
    """Enumerate every coefficient label in the chosen sub-space."""
    return [
        (X, a, b, s)
        for X in structures
        for a in flavors
        for b in flavors
        for s in st
    ]


def label(key) -> str:
    X, a, b, st = key
    return f"{X}_{a}{b}{st}"
