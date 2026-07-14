# `leftfit` — constraining LEFT Wilson coefficients from muon-decay Michel data

A small, self-validating fit that turns the **full four-fermion Michel-parameter
results** of the thesis into experimental bounds on the LEFT Wilson coefficients,
and on the new-physics scale Λ. Built to feed the *"Constraints for the UV fields"*
chapter of `main.tex`.

## Quick start

```bash
cd constraints
pip install -r requirements.txt          # numpy, scipy
python run_constraints.py                # prints validation + all tables
python run_constraints.py --md           # also writes results.md
python -m leftfit.validate               # just the self-checks
```

## What it computes

1. **Validation** (`leftfit/validate.py`) — re-derives the analytic limits quoted
   in the thesis (SM Michel values; `G_F ⇒ Λ ≳ 5 TeV` band on `Re[x^{VLL}_{2112}]`;
   `η = ½ Re[x^{VLR}_{2112}]`; scalars push `ρ` down, tensors push it up). If a
   Michel formula was mis-transcribed, a check fails. **Run this first.**

2. **Single-operator bounds** (`fit.single_operator_bound`) — the single-operator
   dominance assumption: switch on one coefficient (real or imaginary part),
   profile χ², report the 95% CL interval and implied `Λ = v/√|x|`.

3. **Global profile** (`fit.profile_bound`) — turn on several coefficients at once,
   minimise the global χ², and profile one direction with the others floating.
   Exposes cancellations the single-operator bounds miss.

4. **Flat directions** (`fit.flat_directions`) — SVD of the observable Jacobian at
   the SM point + Fisher spectrum. Shows which coefficient combinations are
   linearly constrained and which are flat.

## The one physics subtlety the code is built around

Lepton-number-violating (scalar/tensor) coefficients, and every non-SM neutrino
flavour, enter the observables **only as `|c|²`**. So near the SM point:

- their χ² has a **quartic flat bottom** (zero curvature at `c = 0`);
- the Fisher matrix is **degenerate** for them, and Gaussian error propagation
  gives the wrong confidence interval;
- their reach scales as `σ^{1/4}`, not `σ^{1/2}` — hence the ~1 TeV bounds versus
  the tens-of-TeV linear (LNC-vector, SM-flavour) bounds.

The code therefore **scans the exact profile** for those directions and reports
Fisher/PCA only for the linearly-sensitive ones.

## Conventions

- Fit parameter is the **dimensionless** `x ≡ v² c ~ (v/Λ)²`, `v = 246 GeV`.
  `Λ` bounds assume an order-one UV coupling.
- Coefficient label `(X, a, b, st)`: structure `X ∈ {VLL,VLR,SLL,SLR,TLL}`,
  neutrino generations `a,b ∈ {1,2,3}`, charged-lepton slot `st ∈ {12,21}`.
- LNC vector operators are Hermitian, so `(a,b,st)` and its conjugate
  `(b,a, rev st)` are tied together (`build_single(..., hermitian=True)`).

## Caveats — read before quoting numbers

- **Formulas are hand-transcribed** from `main.tex` (`leftfit/michel.py`, with
  source line references). Verify against your own derivation; the validation
  harness only guards the limits it knows about.
- **Two things the harness already flags:**
  - `β'/A` sign is **inconsistent in the thesis** — the SM-like eq. (~L597) has
    `+¼ Im`, the full-4-fermion eq. (~L671) has `−¼ Im`. `michel.py` uses the
    full-4-fermion sign. *Resolve this in the thesis.*
  - `xi()` contains one apparently mixed-`st` term (`main.tex` L637,
    `|c^{SLR}_{ab12}+c^{SLR}_{ba21}|`), transcribed literally — check for a typo.
- **Experimental inputs** in `leftfit/observables.py` marked `PLACEHOLDER`
  (notably `σ_η`) must be replaced with your adopted values, and the χ² is
  currently **diagonal** — drop in the PDG correlation matrix in `chi2()` for a
  correct global fit.
- **EFT validity.** The numerical box allows `|x|` up to 5, i.e. `Λ < v`, which is
  unphysical. Meaningful bounds need `|x| ≲ 1`. Some "flat directions" in the
  global profile only open up via `|x| = O(4)` cancellations — artifacts of an
  invalid EFT regime; a perturbativity prior (`|x| ≲ 1`) closes them.

## Files

```
constraints/
  leftfit/
    coefficients.py   coefficient bookkeeping, Λ conversion
    michel.py         Michel parameters vs Wilson coefficients (full 4-fermion)
    observables.py    experimental inputs + χ²
    fit.py            single-op bounds, global profiling, flat-direction PCA
    validate.py       self-checks against thesis analytic limits
  run_constraints.py  driver → tables + results.md
  requirements.txt
```

## Extending

- **More observables:** add `Observable(...)` entries in `observables.py`
  (`ξ'`, `ξ''`, `η''`, `α'/A`, `β'/A` theory functions already exist in `michel.py`).
- **More operators:** just pass new `(X,a,b,st)` keys to the fit functions.
- **Correlations / Bayesian:** replace `chi2()` with a covariance-weighted form,
  or wrap the χ² in an MCMC for marginalised (rather than profiled) posteriors.
