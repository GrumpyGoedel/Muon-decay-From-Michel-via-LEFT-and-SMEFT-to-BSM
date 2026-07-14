r"""Michel parameters as functions of the LEFT Wilson coefficients.

These are the *full four-fermion* results (all two-neutrino final states, any
generation, either lepton number), transcribed from the thesis
``main.tex`` -> chapter "The Low Energy Effective Field Theory (LEFT)",
section "Michel parameters in the LEFT", subsection "Full 4-fermion results".

    G_F^2/G_F,SM^2   main.tex ~L616-618
    rho              main.tex ~L624-628
    eta              main.tex ~L629-634
    xi               main.tex ~L636-641
    delta            main.tex ~L643-648
    xi'  xi''  eta''  alpha'/A  beta'/A     main.tex ~L655-673

--------------------------------------------------------------------------
!! VERIFY AGAINST YOUR OWN DERIVATION.  These are transcribed by hand from
   LaTeX; the sign/index conventions (in particular the "12<->21" swaps in
   xi and delta, and one apparently mixed-st term in xi flagged below) are
   exactly the places where a transcription slip would hide.  The validation
   harness in ``validate.py`` re-derives the analytic limits quoted in the
   thesis and will complain if a formula is off.
--------------------------------------------------------------------------

Notation used inside the functions:
    V(X,a,b,st)  -> x.val(X,a,b,st)           (a coefficient, complex)
    "12->21"     -> the same bracket summed over st in ('12','21')

We work to O(x^2) = O(v^4 c^2) = O(1/Lambda^4), matching the thesis.
"""

from __future__ import annotations

from .coefficients import FLAVORS, ST, Coefficients

_PAIRS = [(a, b) for a in FLAVORS for b in FLAVORS if a < b]   # a<b
_DIAG = [(a, a) for a in FLAVORS]                              # a=b
_FULL = [(a, b) for a in FLAVORS for b in FLAVORS]             # all a,b


def gf_ratio_sq(x: Coefficients) -> float:
    """G_F^2 / G_{F,SM}^2  (equals 1 in the SM)."""
    V = x.val
    r = 1.0 + 0j
    # linear, LNC vector, SM flavour only
    r += -0.5 * (V("VLL", 2, 1, "12") + V("VLL", 1, 2, "21"))
    # LNC vector quadratic, full flavour sum
    for a, b in _FULL:
        r += 0.25 * (
            V("VLL", a, b, "12") * V("VLL", b, a, "21")
            + V("VLR", a, b, "12") * V("VLR", b, a, "21")
        )
    # scalar, off-diagonal (a<b), both st slots
    for st in ST:
        for a, b in _PAIRS:
            r += (1 / 16) * (
                abs(V("SLL", a, b, st) + V("SLL", b, a, st)) ** 2
                + abs(V("SLR", a, b, st) + V("SLR", b, a, st)) ** 2
            )
    # scalar, diagonal (a=b), factor 2 for identical neutrinos, both st
    for st in ST:
        for a, _ in _DIAG:
            r += 2 * (1 / 16) * (
                abs(V("SLL", a, a, st)) ** 2 + abs(V("SLR", a, a, st)) ** 2
            )
    # tensor, off-diagonal, both st
    for st in ST:
        for a, b in _PAIRS:
            r += 3 * abs(V("TLL", a, b, st) - V("TLL", b, a, st)) ** 2
    return r.real


def rho(x: Coefficients) -> float:
    V = x.val
    r = 0.75
    for st in ST:
        for a, b in _PAIRS:
            r += -(3 / 64) * (
                abs(V("SLL", a, b, st) + V("SLL", b, a, st)) ** 2
                + abs(V("SLR", a, b, st) + V("SLR", b, a, st)) ** 2
            )
    for st in ST:
        for a, _ in _DIAG:
            r += -2 * (3 / 64) * (
                abs(V("SLL", a, a, st)) ** 2 + abs(V("SLR", a, a, st)) ** 2
            )
    for st in ST:
        for a, b in _PAIRS:
            r += (3 / 4) * abs(V("TLL", a, b, st) - V("TLL", b, a, st)) ** 2
    return r


def eta(x: Coefficients) -> float:
    V = x.val
    r = 0j
    # linear (SM flavour); c^VLR_1221 = (c^VLR_2112)* so this is (1/2)Re[..]
    r += 0.25 * (V("VLR", 2, 1, "12") + V("VLR", 1, 2, "21"))
    # vector-vector interference, SM flavour specific
    r += (1 / 8) * (
        V("VLL", 1, 2, "21") * V("VLR", 2, 1, "12")
        + V("VLL", 2, 1, "12") * V("VLR", 1, 2, "21")
    )
    r += (1 / 8) * (
        V("VLL", 2, 1, "12") * V("VLR", 2, 1, "12")
        + V("VLL", 1, 2, "21") * V("VLR", 1, 2, "21")
    )
    # vector-vector, full flavour sum
    for a, b in _FULL:
        r += -(1 / 8) * (
            V("VLL", b, a, "12") * V("VLR", a, b, "12")
            + V("VLL", a, b, "12") * V("VLR", b, a, "21")
        )
    # scalar-scalar interference (real part), off-diagonal, both st
    for st in ST:
        for a, b in _PAIRS:
            r += (1 / 8) * (
                (V("SLL", a, b, st) + V("SLL", b, a, st)).conjugate()
                * (V("SLR", a, b, st) + V("SLR", b, a, st))
            ).real
    # scalar-scalar interference, diagonal, both st
    for st in ST:
        for a, _ in _DIAG:
            r += (1 / 4) * (
                V("SLL", a, a, st).conjugate() * V("SLR", a, a, st)
            ).real
    return r.real


def xi(x: Coefficients) -> float:
    V = x.val
    r = 1.0 + 0j
    # LNC vector quadratic
    for a, b in _FULL:
        r += -0.5 * V("VLR", a, b, "12") * V("VLR", b, a, "21")
    # scalar off-diagonal.  NOTE: main.tex L637 writes the first modulus as
    # |c^SLR_ab12 + c^SLR_ba21| (mixed st) -- transcribed literally; VERIFY.
    for a, b in _PAIRS:
        r += (1 / 8) * (
            abs(V("SLR", a, b, "12") + V("SLR", b, a, "21")) ** 2
            - 2 * abs(V("SLL", a, b, "12") + V("SLL", b, a, "12")) ** 2
        )
    for a, b in _PAIRS:
        r += (1 / 8) * (
            -2 * abs(V("SLR", a, b, "21") + V("SLR", b, a, "21")) ** 2
            + abs(V("SLL", a, b, "21") + V("SLL", b, a, "21")) ** 2
        )
    for a, _ in _DIAG:
        r += (1 / 4) * (
            abs(V("SLL", a, a, "21")) ** 2 - 2 * abs(V("SLR", a, a, "21")) ** 2
            - 2 * abs(V("SLL", a, a, "12")) ** 2 + abs(V("SLR", a, a, "12")) ** 2
        )
    for a, b in _PAIRS:
        r += (1 / 8) * (
            32 * abs(V("TLL", a, b, "12") - V("TLL", b, a, "12")) ** 2
            - 80 * abs(V("TLL", a, b, "21") - V("TLL", b, a, "21")) ** 2
        )
    return r.real


def delta(x: Coefficients) -> float:
    V = x.val
    r = 0.75
    # scalar off-diagonal, antisymmetric under 12<->21
    for a, b in _PAIRS:
        base12 = (
            abs(V("SLL", a, b, "12") + V("SLL", b, a, "12")) ** 2
            - abs(V("SLR", a, b, "12") + V("SLR", b, a, "12")) ** 2
        )
        base21 = (
            abs(V("SLL", a, b, "21") + V("SLL", b, a, "21")) ** 2
            - abs(V("SLR", a, b, "21") + V("SLR", b, a, "21")) ** 2
        )
        r += (9 / 64) * (base12 - base21)
    # scalar diagonal, antisymmetric under 12<->21
    for a, _ in _DIAG:
        d12 = abs(V("SLL", a, a, "12")) ** 2 - abs(V("SLR", a, a, "12")) ** 2
        d21 = abs(V("SLL", a, a, "21")) ** 2 - abs(V("SLR", a, a, "21")) ** 2
        r += (9 / 32) * (d12 - d21)
    # tensor off-diagonal, antisymmetric under 12<->21
    for a, b in _PAIRS:
        t12 = abs(V("TLL", a, b, "12") - V("TLL", b, a, "12")) ** 2
        t21 = abs(V("TLL", a, b, "21") - V("TLL", b, a, "21")) ** 2
        r += -(9 / 4) * (t12 - t21)
    return r


# ---------------------------------------------------------------------------
# Electron-spin-sensitive parameters (main.tex L655-673).  Included for
# completeness; wire experimental values into observables.py to use them.
# ---------------------------------------------------------------------------

def xi_prime(x: Coefficients) -> float:
    V = x.val
    r = 1.0 + 0j
    for a, b in _FULL:
        r += -0.5 * V("VLR", a, b, "12") * V("VLR", b, a, "21")
    for a, b in _PAIRS:
        r += -(1 / 8) * (
            abs(V("SLL", a, b, "12") + V("SLL", b, a, "12")) ** 2
            + abs(V("SLR", a, b, "21") + V("SLR", b, a, "21")) ** 2
        )
    for a, _ in _DIAG:
        r += -(1 / 4) * (
            abs(V("SLL", a, a, "12")) ** 2 + abs(V("SLR", a, a, "21")) ** 2
        )
    for a, b in _PAIRS:
        r += -6 * abs(V("TLL", a, b, "12") - V("TLL", b, a, "12")) ** 2
    return r.real


def xi_double_prime(x: Coefficients) -> float:
    V = x.val
    r = 1.0 + 0j
    for st in ST:
        for a, b in _PAIRS:
            r += (1 / 8) * (
                abs(V("SLL", a, b, st) + V("SLL", b, a, st)) ** 2
                + abs(V("SLR", a, b, st) + V("SLR", b, a, st)) ** 2
            )
    for st in ST:
        for a, _ in _DIAG:
            r += (1 / 4) * (
                abs(V("SLL", a, a, st)) ** 2 + abs(V("SLR", a, a, st)) ** 2
            )
    for st in ST:
        for a, b in _PAIRS:
            r += -10 * abs(V("TLL", a, b, st) - V("TLL", b, a, st)) ** 2
    return r.real


def eta_double_prime(x: Coefficients) -> float:
    V = x.val
    r = 0j
    r += -0.25 * (V("VLR", 2, 1, "12") + V("VLR", 1, 2, "21"))
    r += -(1 / 8) * (
        V("VLL", 1, 2, "21") * V("VLR", 2, 1, "12")
        + V("VLL", 2, 1, "12") * V("VLR", 1, 2, "21")
    )
    r += -(1 / 8) * (
        V("VLL", 2, 1, "12") * V("VLR", 2, 1, "12")
        + V("VLL", 1, 2, "21") * V("VLR", 1, 2, "21")
    )
    for a, b in _FULL:
        r += (1 / 8) * (
            V("VLL", b, a, "21") * V("VLR", a, b, "12")
            + V("VLL", a, b, "12") * V("VLR", b, a, "21")
        )
    for st in ST:
        for a, b in _PAIRS:
            r += (3 / 8) * (
                (V("SLL", a, b, st) + V("SLL", b, a, st)).conjugate()
                * (V("SLR", a, b, st) + V("SLR", b, a, st))
            ).real
    for st in ST:
        for a, _ in _DIAG:
            r += (3 / 4) * (
                V("SLL", a, a, st).conjugate() * V("SLR", a, a, st)
            ).real
    return r.real


def alpha_prime(x: Coefficients) -> float:
    """alpha'/A -- T-odd, gets contributions only from LNV scalar interference."""
    V = x.val
    r = 0.0
    for st in ST:
        for a, b in _PAIRS:
            r += -(1 / 8) * (
                (V("SLL", a, b, st) + V("SLL", b, a, st)).conjugate()
                * (V("SLR", a, b, st) + V("SLR", b, a, st))
            ).imag
    for st in ST:
        for a, _ in _DIAG:
            r += -(1 / 4) * (
                V("SLL", a, a, st).conjugate() * V("SLR", a, a, st)
            ).imag
    return r


def beta_prime(x: Coefficients) -> float:
    """beta'/A -- fixes the imaginary part of the SM-flavour VLR coefficient."""
    V = x.val
    r = 0.0
    r += -0.25 * V("VLR", 2, 1, "12").imag
    r += -(1 / 8) * (V("VLL", 2, 1, "12").conjugate() * V("VLR", 2, 1, "12")).imag
    r += -(1 / 8) * (V("VLL", 2, 1, "12") * V("VLR", 2, 1, "12")).imag
    for a, b in _FULL:
        r += (1 / 8) * (
            V("VLL", a, b, "12").conjugate() * V("VLR", a, b, "12")
        ).imag
    return r


# Convenience registry
PARAMS = {
    "GF2_ratio": gf_ratio_sq,
    "rho": rho,
    "eta": eta,
    "xi": xi,
    "delta": delta,
    "xi_prime": xi_prime,
    "xi_double_prime": xi_double_prime,
    "eta_double_prime": eta_double_prime,
    "alpha_prime": alpha_prime,
    "beta_prime": beta_prime,
}
