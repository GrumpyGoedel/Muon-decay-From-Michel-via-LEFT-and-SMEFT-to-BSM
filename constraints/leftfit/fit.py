r"""Fitting machinery.

Three layers, in increasing sophistication:

1. single_operator_bound  -- switch on ONE coefficient (real or imaginary
   part), profile chi^2 over its value, return the confidence interval and the
   implied new-physics scale Lambda.  This is the "single-operator dominance"
   assumption.

2. global_fit / profile_bound -- turn on many coefficients at once, minimise
   the global chi^2, and profile a single (real) direction with all others
   floating.  This is the honest bound that accounts for cancellations.

3. flat_directions -- SVD of the observable Jacobian at the SM point plus the
   Fisher eigen-spectrum, to expose which coefficient combinations muon decay
   actually constrains and which are flat.

Key subtlety implemented here: for lepton-number-violating (scalar/tensor)
coefficients the observables depend only on |c|^2, so chi^2 is QUARTIC in the
real parameter and its curvature vanishes at c=0.  Gaussian/Fisher error
propagation is therefore invalid for those directions -- we scan the exact
(quartic) profile instead.  The Fisher analysis is reported only for the
linearly-sensitive (LNC) directions.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar

from . import michel
from .coefficients import Coefficients, label, lambda_from_x
from .observables import DEFAULT_OBSERVABLES, chi2

# Delta chi^2 thresholds
DCHI2_1PAR_95 = 3.8414588   # 1 dof, 95% CL
DCHI2_1PAR_68 = 1.0         # 1 dof, 68% CL


# ---------------------------------------------------------------------------
# 1. Single-operator bounds
# ---------------------------------------------------------------------------

def conjugate_key(key):
    """Hermitian-conjugate label: prst -> rpts, i.e. (X,a,b,st)->(X,b,a,rev st)."""
    X, a, b, st = key
    return (X, b, a, "21" if st == "12" else "12")


def is_hermitian_structure(X):
    """Only the LNC vector operators are Hermitian, so their coefficient and
    its conjugate partner are related; the LNV scalar/tensor coefficients are
    independent complex numbers."""
    return X in ("VLL", "VLR")


def build_single(key, value, hermitian=True):
    """A coefficient set with `key` switched on, and (for Hermitian LNC
    structures) its conjugate partner set to the conjugate value."""
    c = Coefficients.single(key, value)
    if hermitian and is_hermitian_structure(key[0]):
        pk = conjugate_key(key)
        if pk != key:
            c[pk] = complex(value).conjugate()
    return c


def single_operator_bound(key, part="real", cl=DCHI2_1PAR_95,
                          observables=None, span=5.0, hermitian=True):
    """Confidence interval on a single coefficient's real (or imaginary) part.

    Returns a dict with the chi^2 minimum, the interval [t_lo, t_hi] where
    Delta chi^2 <= cl, the tightest |t|, and the implied Lambda.
    """
    obs = observables if observables is not None else DEFAULT_OBSERVABLES

    def f(t):
        value = t if part == "real" else 1j * t
        return chi2(build_single(key, value, hermitian), obs)

    # locate minimum
    res = minimize_scalar(f, bounds=(-span, span), method="bounded")
    t_min, chi2_min = res.x, res.fun
    target = chi2_min + cl

    def g(t):
        return f(t) - target

    # bracket outward from the minimum on each side
    t_hi = _find_crossing(g, t_min, +1, span)
    t_lo = _find_crossing(g, t_min, -1, span)

    t_abs = max(abs(t_lo), abs(t_hi))
    return {
        "key": key,
        "label": label(key),
        "part": part,
        "t_min": t_min,
        "chi2_min": chi2_min,
        "t_lo": t_lo,
        "t_hi": t_hi,
        "t_abs_95": t_abs,
        "Lambda_TeV": lambda_from_x(t_abs) / 1000.0,
        "linear": linear_sensitivity(key, part, obs, hermitian) > 1e-6,
    }


def _find_crossing(g, t0, direction, span):
    """Find where g changes sign moving from t0 in the given direction."""
    step = span / 200.0
    t_prev = t0
    t = t0 + direction * step
    while abs(t) <= span:
        if g(t_prev) * g(t) <= 0:
            return brentq(g, min(t_prev, t), max(t_prev, t))
        t_prev, t = t, t + direction * step
    return direction * span  # no crossing within span -> unconstrained side


def linear_sensitivity(key, part, observables=None, hermitian=True, eps=1e-6):
    """Largest |dO_i/dt| at the SM point -- i.e. whether this coefficient
    (given the *observable set*) interferes with the SM at O(c).  Nonzero only
    for LNC vectors in the SM flavour slot that feed a measured parameter."""
    obs = observables if observables is not None else DEFAULT_OBSERVABLES
    slopes = []
    for o in obs:
        vp = eps if part == "real" else 1j * eps
        vm = -eps if part == "real" else -1j * eps
        op = o.model(build_single(key, vp, hermitian))
        om = o.model(build_single(key, vm, hermitian))
        slopes.append(abs(op - om) / (2 * eps))
    return max(slopes)


# ---------------------------------------------------------------------------
# 2. Global fit with profiling
# ---------------------------------------------------------------------------

def _pack(keys):
    """Parameter layout: [Re(k0), Im(k0), Re(k1), Im(k1), ...]."""
    return {k: (2 * i, 2 * i + 1) for i, k in enumerate(keys)}


def _coeffs_from_vec(keys, layout, theta):
    c = Coefficients()
    for k in keys:
        ir, ii = layout[k]
        c[k] = theta[ir] + 1j * theta[ii]
    return c


def global_fit(keys, observables=None, fixed=None, span=5.0, warm=None):
    """Minimise the global chi^2 over the real & imaginary parts of `keys`,
    inside the box [-span, span] (so flat directions cannot run to infinity).

    `fixed` maps {index_in_theta: value} to hold a direction during profiling.
    `warm`  is an optional starting vector (used to speed up profiling scans).
    """
    obs = observables if observables is not None else DEFAULT_OBSERVABLES
    layout = _pack(keys)
    n = 2 * len(keys)
    bounds = [(-span, span)] * n

    def objective(theta):
        if fixed:
            theta = np.array(theta, dtype=float)
            for idx, val in fixed.items():
                theta[idx] = val
        return chi2(_coeffs_from_vec(keys, layout, theta), obs)

    if warm is not None:
        seeds = [np.array(warm, dtype=float)]
    else:
        seeds = [np.zeros(n)] + [0.1 * np.sin(np.arange(n) + s) for s in range(1, 5)]
    best = None
    for x0 in seeds:
        x0 = np.array(x0, dtype=float)
        if fixed:
            for idx, val in fixed.items():
                x0[idx] = val
        r = minimize(objective, x0, method="L-BFGS-B", bounds=bounds,
                     options={"ftol": 1e-12, "gtol": 1e-10, "maxiter": 5000})
        if best is None or r.fun < best.fun:
            best = r
    return best, layout


def profile_bound(keys, target_key, part="real",
                  cl=DCHI2_1PAR_95, observables=None, span=5.0):
    """Profile chi^2 along one real direction with every other coefficient
    floating, via geometric root-finding on the profiled chi^2.  A direction
    with no crossing inside [-span, span] is reported as unconstrained (flat).
    """
    obs = observables if observables is not None else DEFAULT_OBSERVABLES
    layout = _pack(keys)
    idx = layout[target_key][0 if part == "real" else 1]

    base, _ = global_fit(keys, obs, span=span)
    chi2_min = base.fun
    t_best = float(base.x[idx])
    target = chi2_min + cl

    def g(t):
        r, _ = global_fit(keys, obs, fixed={idx: t}, span=span, warm=base.x)
        return r.fun - target

    hi = _crossing_expand(g, t_best, +1, span)
    lo = _crossing_expand(g, t_best, -1, span)
    unconstrained = (abs(hi) >= span - 1e-6) and (abs(lo) >= span - 1e-6)
    t_abs = max(abs(lo), abs(hi))
    return {
        "target": label(target_key),
        "part": part,
        "chi2_min": chi2_min,
        "t_lo": lo,
        "t_hi": hi,
        "t_abs_95": t_abs,
        "Lambda_TeV": float("inf") if unconstrained else lambda_from_x(t_abs) / 1000.0,
        "unconstrained": unconstrained,
    }


def _crossing_expand(g, t0, direction, span, step0=1e-4):
    """Find a sign change of g moving outward from t0 with geometrically
    growing steps, then bisect.  Returns direction*span if none is found."""
    step = step0
    t_prev = t0
    t = t0 + direction * step
    while abs(t) <= span:
        if g(t_prev) * g(t) <= 0:
            return brentq(g, min(t_prev, t), max(t_prev, t))
        step *= 1.5
        t_prev = t
        t = t + direction * step
    return direction * span


# ---------------------------------------------------------------------------
# 3. Flat-direction / sensitivity analysis
# ---------------------------------------------------------------------------

def observable_jacobian(keys, observables=None, eps=1e-6):
    """Numerical d O_i / d theta_j at the SM point (theta = 0)."""
    obs = observables if observables is not None else DEFAULT_OBSERVABLES
    layout = _pack(keys)
    n = 2 * len(keys)
    J = np.zeros((len(obs), n))
    for j in range(n):
        tp, tm = np.zeros(n), np.zeros(n)
        tp[j], tm[j] = eps, -eps
        op = np.array([o.model(_coeffs_from_vec(keys, layout, tp)) for o in obs])
        om = np.array([o.model(_coeffs_from_vec(keys, layout, tm)) for o in obs])
        J[:, j] = (op - om) / (2 * eps)
    return J, layout


def flat_directions(keys, observables=None, tol=1e-9):
    """SVD of the Jacobian at the SM point + Fisher eigen-spectrum.

    Directions with (near-)zero singular value receive NO linear sensitivity;
    for LNV coefficients that is *expected* -- they are bounded only at
    O(theta^2), via single_operator_bound / profile_bound.
    """
    obs = observables if observables is not None else DEFAULT_OBSERVABLES
    J, layout = observable_jacobian(keys, obs)
    sigma = np.array([o.sigma for o in obs])
    Jw = J / sigma[:, None]                       # weighted Jacobian
    U, s, Vt = np.linalg.svd(Jw, full_matrices=True)

    param_names = []
    for k in keys:
        param_names += [f"Re {label(k)}", f"Im {label(k)}"]

    constrained, flat = [], []
    for i in range(Vt.shape[0]):
        direction = Vt[i]
        sval = s[i] if i < len(s) else 0.0
        entry = {
            "singular_value": float(sval),
            "components": {param_names[j]: float(direction[j])
                           for j in range(len(param_names))
                           if abs(direction[j]) > 1e-6},
        }
        (constrained if sval > tol else flat).append(entry)

    fisher = Jw.T @ Jw
    evals = np.linalg.eigvalsh(fisher)
    return {
        "param_names": param_names,
        "singular_values": s.tolist(),
        "constrained": constrained,
        "flat": flat,
        "fisher_eigenvalues": evals.tolist(),
        "n_constrained": int((s > tol).sum()),
        "n_params": len(param_names),
    }
