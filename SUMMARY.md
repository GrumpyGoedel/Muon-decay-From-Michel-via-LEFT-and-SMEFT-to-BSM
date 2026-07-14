# Muon decay: From Michel via LEFT and SMEFT to BSM

Thesis by Carl Jakob Moritz (University of Vienna, supervisor: Tyler Corbett).
Source: `main.tex` (+ `literature.bib`).

## Goal

Trace how a heavy BSM field would imprint itself on muon decay, a process living
well below the electroweak scale. The strategy bridges three energy regimes:

```
BSM field  --(integrate out)-->  SMEFT (dim 6,7)  --(match/run)-->  LEFT (dim 5,6)  -->  Michel parameters
   ~10 TeV                            ~v = 246 GeV                       ~m_mu = 106 MeV        (experiment)
```

The Michel parameters of muon decay are perturbed by LEFT operators; from their
measured values one bounds the LEFT Wilson coefficients, runs/matches up to SMEFT,
and finally reads off which UV fields could be responsible.

## Chapter map

1. **Muon decay in the SM** — SM gauge structure, weak leptonic currents, integrating
   out the `W` to a 4-Fermi theory, total decay width, and the Michel parameters.
2. **The LEFT** — unbroken `SU(3)_c × U(1)_em` EFT below the EW scale. Operator basis,
   the Michel parameters expressed in LEFT Wilson coefficients (SM-like *and* full
   4-fermion / all-neutrino-output results), and LEFT RGEs up to `v`.
3. **The SMEFT** — dim-6 and dim-7 basis, tree-level matching onto the relevant LEFT
   operators, and RGE running up to the new-physics scale (`wilson`, `D7RGESolver`).
4. **BSM** — top-down vs. bottom-up (j-basis / on-shell) thinking; a UV dictionary
   (from Li et al.) mapping each SMEFT operator to candidate UV scalars/fermions/vectors.
5. **Outlook** — open questions (e.g. non-linear EWSB → HEFT).
- **Appendices** — kinematics & phase-space integrals, Michel "basis functions",
  full LEFT/SMEFT operator bases, RGEs, and the UV Lagrangian.

## Michel parameters and key experimental inputs

Standard set: `ρ` (isotropic shape), `δ` (anisotropic shape), `ξ` (parity violation),
`η` (scalar/tensor interference, needs mass insertion). Electron-spin-sensitive
extensions: `ξ'`, `ξ''`, `η''`, `α'/A`, `β'/A`. SM predicts `ρ = δ = 3/4`, `ξ = ξ' = ξ'' = 1`,
`η = η'' = α'/A = β'/A = 0`.

Current data used in the thesis:
- `ρ = 0.74979 ± 0.00026`
- `P_μ ξ = 1.0009 ± 0.0012`
- `ξ P_μ δ / ρ = 1.00179 ± 0.00114`
- `η ≈ 0.057` (low-energy spectrum end)
- `|g_RR^S| < 0.035` (PDG SPVAT bound, reused to constrain a LEFT coefficient)

## Relevant LEFT operators (dim 6, four-fermion `ν-e`)

`O^{V,LL}`, `O^{V,LR}` (lepton-number conserving), and `O^{S,LL}`, `O^{S,LR}`, `O^{T,LL}`
(lepton-number violating). Only the two `V` operators conserve lepton number; the
scalar/tensor ones enter only through their moduli-squared (no SM interference).

## Key results

- **SM-like neutrino output** (`μ⁻ → e⁻ ν̄_e ν_μ`): `ρ = δ = 3/4` unchanged at leading
  order; `η` and `β'/A` fix Re/Im of `L^{V,LR}_{νe,2112}`. Bounds obtained, e.g.
  `Λ ≳ 21 v ≈ 5.1 TeV` (95% CL) from `G_F`, and a joint `χ²` over `G_F, η, ξ`.
- **Full 4-fermion results**: closed-form Michel-parameter corrections summed over
  *all* two-neutrino final states (any generation `a,b`, either lepton number). LNV
  scalar/tensor coefficients appear quadratically (as `|c|²`) in `ρ, δ, ξ, ξ', ξ''`;
  the LNC vector coefficients enter linearly via interference. This full set is the
  basis for the discussion in the conversation with this repo.
- **LEFT RGEs**: the LNC muon-decay coefficients have vanishing anomalous dimension
  (`δ_st` kills them); scalar/tensor run only mildly (~2–3%) from `m_μ` to `v`.
- **SMEFT matching**: LNC `V` operators match at dim-6 onto `Q_ll`, `Q_le`, `Q_Hl^{(1,3)}`;
  the LNV `S`/`T` operators match only at dim-7 onto `O_LeHD`, `O_{ēLLLH}`.
- **UV dictionary**: e.g. `O_ll → S₁,S₆,V₁,V₄`; `O_Hl^{(3)} → F₁,F₂,F₅,F₆`;
  `O_le → S₄,V₁,V₃`; dim-7 operators need field *combinations* (e.g. `F₃+S₆`).
  Simplest example worked out: scalar `S₂ = (1,1)₁`, a single new Yukawa-like term.

## Status / notes

- Sections **"Constraints for the UV fields"** and **"Conclusion"** are currently
  stubs (empty), as is much of **Outlook**. The full-4-fermion → experiment
  constraint analysis is the active open thread.
- Some numbers are placeholders (e.g. `η = 0.057 ± 0`) pending final values.
- Tooling referenced: `wilson` (dim-6 running), `D7RGESolver` (dim-7 running).
