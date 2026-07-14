#!/usr/bin/env python3
r"""Driver: produce the constraint tables for the "Constraints for the UV
fields" chapter.

    python run_constraints.py            # print everything to stdout
    python run_constraints.py --md       # also write results.md

Outputs:
    1. validation of the transcribed Michel formulas
    2. single-operator 95% CL bounds table (one coefficient at a time)
    3. a small global profile over the LNC vector sector
    4. flat-direction / Fisher analysis
"""

from __future__ import annotations

import sys

from leftfit import fit, observables, validate
from leftfit.coefficients import label

# A representative set of single operators to bound.  Extend freely.
SINGLE_OPS = [
    # LNC vectors, SM flavour (linearly sensitive)
    (("VLL", 2, 1, "12"), "real"),
    (("VLR", 2, 1, "12"), "real"),
    (("VLR", 2, 1, "12"), "imag"),
    # LNC vectors, non-SM flavour (quadratic only)
    (("VLL", 3, 3, "12"), "real"),
    # LNV scalars / tensor (quadratic only)
    (("SLL", 1, 1, "12"), "real"),
    (("SLR", 1, 1, "12"), "real"),
    (("SLR", 1, 2, "12"), "real"),
    (("TLL", 1, 2, "12"), "real"),
]

# Global-fit sector: independent LNC coefficients (SM flavour, real & imag).
GLOBAL_KEYS = [("VLL", 2, 1, "12"), ("VLR", 2, 1, "12")]


def hline(c="-", n=72):
    return c * n


def section(title):
    print("\n" + hline("="))
    print(title)
    print(hline("="))


def single_operator_table(lines):
    section("2. SINGLE-OPERATOR 95% CL BOUNDS  (v^2 c, one coefficient at a time)")
    header = f"{'operator':<14}{'part':<6}{'interval [t_lo, t_hi]':<28}" \
             f"{'|t|_95':<10}{'Lambda [TeV]':<13}{'entry'}"
    print(header)
    print(hline())
    lines.append("## Single-operator 95% CL bounds (one coefficient at a time)")
    lines.append("")
    lines.append("| operator | part | interval | \\|t\\|₉₅ | Λ [TeV] | entry |")
    lines.append("|---|---|---|---|---|---|")
    for key, part in SINGLE_OPS:
        b = fit.single_operator_bound(key, part=part)
        kind = "linear" if b["linear"] else "quadratic"
        interval = f"[{b['t_lo']:+.4f}, {b['t_hi']:+.4f}]"
        print(f"{b['label']:<14}{part:<6}{interval:<28}"
              f"{b['t_abs_95']:<10.4f}{b['Lambda_TeV']:<13.2f}{kind}")
        lines.append(f"| `{b['label']}` | {part} | {interval} | "
                     f"{b['t_abs_95']:.4f} | {b['Lambda_TeV']:.2f} | {kind} |")


def global_section(lines):
    section("3. GLOBAL PROFILE  (LNC vector sector, all others floating)")
    lines.append("")
    lines.append("## Global profile (LNC vector sector, all others floating)")
    lines.append("")
    for key in GLOBAL_KEYS:
        for part in ("real", "imag"):
            r = fit.profile_bound(GLOBAL_KEYS, key, part=part)
            if r["unconstrained"]:
                verdict = "unconstrained (flat direction)"
                mdverdict = "unconstrained (flat)"
            else:
                verdict = (f"[{r['t_lo']:+.4f}, {r['t_hi']:+.4f}]  "
                           f"|t|_95 = {r['t_abs_95']:.4f}  "
                           f"Lambda = {r['Lambda_TeV']:.2f} TeV")
                mdverdict = (f"[{r['t_lo']:+.4f}, {r['t_hi']:+.4f}], "
                             f"Λ = {r['Lambda_TeV']:.2f} TeV")
            print(f"{label(key):<12} {part:<5}  profiled {verdict}")
            lines.append(f"- `{label(key)}` ({part}): profiled {mdverdict}")


def flat_section(lines):
    section("4. FLAT DIRECTIONS  (Jacobian SVD at the SM point)")
    # Use a broader sector so the flat structure is visible.
    keys = [("VLL", 2, 1, "12"), ("VLR", 2, 1, "12"),
            ("SLR", 1, 1, "12"), ("TLL", 1, 2, "12")]
    info = fit.flat_directions(keys)
    print(f"parameters: {info['n_params']}   "
          f"linearly constrained: {info['n_constrained']}   "
          f"flat: {info['n_params'] - info['n_constrained']}")
    print("\nsingular values:", [f"{s:.3g}" for s in info["singular_values"]])
    print("\nlinearly-sensitive directions:")
    for d in info["constrained"]:
        comps = ", ".join(f"{v:+.2f}*{n}" for n, v in d["components"].items())
        print(f"  sv={d['singular_value']:.3g}:  {comps}")
    print("\nflat directions (no linear sensitivity -- bounded only at O(c^2)):")
    for d in info["flat"]:
        comps = ", ".join(f"{v:+.2f}*{n}" for n, v in d["components"].items())
        print(f"  {comps}")
    lines.append("")
    lines.append(f"**Flat-direction analysis** — {info['n_constrained']} of "
                 f"{info['n_params']} directions are linearly constrained; the "
                 "rest (all LNV) are bounded only quadratically. The scalar/"
                 "tensor coefficients carry no linear sensitivity: their χ² has "
                 "a quartic flat bottom at the SM point, so their reach scales "
                 "as σ^(1/4) and must be read from the single-operator scan, "
                 "not from Gaussian error propagation.")


def main():
    write_md = "--md" in sys.argv
    lines = ["# Muon-decay constraints on LEFT Wilson coefficients",
             "",
             "_Auto-generated by `run_constraints.py`. Dimensionless "
             "coefficients `x = v² c ~ (v/Λ)²`; bounds at 95% CL._", ""]

    section("1. VALIDATION")
    ok = validate.run()
    if not ok:
        print("\n!! Validation failed -- bounds below are NOT trustworthy.\n")

    single_operator_table(lines)
    global_section(lines)
    flat_section(lines)

    if write_md:
        with open("results.md", "w") as f:
            f.write("\n".join(lines) + "\n")
        print("\nwrote results.md")


if __name__ == "__main__":
    main()
