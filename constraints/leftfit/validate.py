r"""Self-checks that re-derive the analytic results quoted in the thesis.

If a Michel-parameter formula was mis-transcribed, one of these will fail --
that is the whole point.  Run with:  python -m leftfit.validate
"""

from __future__ import annotations

from . import michel
from .coefficients import Coefficients, lambda_from_x
from .fit import build_single, single_operator_bound

_checks = []


def check(name, cond, detail=""):
    _checks.append((name, bool(cond), detail))


def run():
    # 1. Standard-Model point reproduces the SM Michel values.
    sm = Coefficients()
    check("SM rho = 3/4", abs(michel.rho(sm) - 0.75) < 1e-12)
    check("SM eta = 0", abs(michel.eta(sm)) < 1e-12)
    check("SM xi = 1", abs(michel.xi(sm) - 1.0) < 1e-12)
    check("SM delta = 3/4", abs(michel.delta(sm) - 0.75) < 1e-12)
    check("SM G_F ratio = 1", abs(michel.gf_ratio_sq(sm) - 1.0) < 1e-12)
    check("SM xi'' = 1", abs(michel.xi_double_prime(sm) - 1.0) < 1e-12)

    # 2. Linear G_F sensitivity to the SM-flavour VLL coefficient.
    #    With hermiticity, G_F^2/G_F,SM^2 = 1 - Re[x^VLL_2112] + O(x^2).
    key = ("VLL", 2, 1, "12")
    h = 1e-5
    slope = (michel.gf_ratio_sq(build_single(key, h)) - 1.0) / h
    check("G_F linear slope = -1", abs(slope + 1.0) < 1e-4,
          f"got {slope:.6f}")

    # 3. The G_F band on Re[x^VLL_2112] reproduces  -0.004 .. 0.002 (thesis).
    b = single_operator_bound(key, part="real")
    check("G_F band lower ~ -0.004", abs(b["t_lo"] - (-0.0044)) < 5e-4,
          f"t_lo = {b['t_lo']:.4f}")
    check("G_F band upper ~ +0.002", abs(b["t_hi"] - 0.0024) < 5e-4,
          f"t_hi = {b['t_hi']:.4f}")
    check("G_F -> Lambda ~ 5 TeV", 3.0 < b["Lambda_TeV"] < 8.0,
          f"Lambda = {b['Lambda_TeV']:.2f} TeV")

    # 4. eta fixes Re[x^VLR_2112] at leading order:  eta = (1/2) Re[x^VLR_2112].
    key2 = ("VLR", 2, 1, "12")
    c2 = build_single(key2, 0.02)     # real; partner 1221 = 0.02
    check("eta = (1/2) Re[x^VLR_2112]",
          abs(michel.eta(c2) - 0.5 * 0.02) < 1e-4,
          f"eta = {michel.eta(c2):.5f}")

    # 5. beta'/A fixes Im[x^VLR_2112] at leading order.
    #    NOTE: main.tex is internally inconsistent here -- the SM-like eq.
    #    (~L597) gives  beta'/A = +(1/4)Im[..], but the full-4-fermion eq.
    #    (~L671) gives  beta'/A = -(1/4)Im[..].  michel.py implements the
    #    full-4-fermion sign (-).  >>> RESOLVE THIS SIGN IN THE THESIS. <<<
    c3 = build_single(key2, 0.02j)    # imaginary; partner 1221 = -0.02j
    check("|beta'/A| = (1/4)|Im[x^VLR_2112]|  (sign per full-4f eq.)",
          abs(abs(michel.beta_prime(c3)) - 0.25 * 0.02) < 1e-4,
          f"beta'/A = {michel.beta_prime(c3):.5f} (expect -0.00500)")

    # 6. A scalar LNV operator enters rho quadratically and only downward.
    ks = ("SLR", 1, 1, "12")
    up = michel.rho(Coefficients.single(ks, 0.1))
    check("scalar SLR pushes rho down", up < 0.75,
          f"rho = {up:.6f}")
    #    ...while a tensor pushes rho up.
    kt = ("TLL", 1, 2, "12")
    upt = michel.rho(Coefficients.single(kt, 0.1))
    check("tensor TLL pushes rho up", upt > 0.75,
          f"rho = {upt:.6f}")

    # report
    ok = all(c for _, c, _ in _checks)
    width = max(len(n) for n, _, _ in _checks)
    print("=" * (width + 20))
    print("leftfit validation")
    print("=" * (width + 20))
    for name, passed, detail in _checks:
        tag = "PASS" if passed else "FAIL"
        line = f"[{tag}] {name.ljust(width)}"
        if detail and not passed:
            line += f"   {detail}"
        elif detail:
            line += f"   ({detail})"
        print(line)
    print("=" * (width + 20))
    print("ALL PASSED" if ok else ">>> SOME CHECKS FAILED -- inspect michel.py")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if run() else 1)
