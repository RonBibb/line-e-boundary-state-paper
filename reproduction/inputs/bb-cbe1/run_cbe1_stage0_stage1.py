"""BB-CBE-1 Stages 0-1 and homogeneous Stage 2, per the frozen program doc.

Program: Discussion/BB_COLLAPSE_BANG_SINGULAR_BOUNDARY_EQUIVALENCE_PROGRAM_2026-08-04.md
Units: 8*pi*G = 1 throughout. Declared action for Stage 1/2: GR + one massless
scalar (stiff fluid, w = 1), the program's Section-10 recommended first
realization.

Frozen before any comparison is computed (weights/tolerances clause, Sec. 5):
- Scored state-vector components on constant-invariant-clock surfaces:
  (Omega, Sigma^2, R2 = Ricci^2/H^4, w_eff, parity residual P = |lam*(H_col + H_exp)|).
  K3 = 0 (flat), W = 0 (FLRW conformally flat) are recorded identically-zero
  entries for the homogeneous control; anisotropic stress, currents, and
  constraint entries are identically satisfied by construction and recorded null.
- Weights: w_A = 1 for every scored component.
- Tolerance: D(lam) < 1.0e-6 at the fine integrator setting, and refinement
  ratio D_fine <= 0.5 * D_coarse unless both are below tolerance.
- Invariant clocks (Sec. C2): lam_rho = rho^{-1/2} and lam_R = Kretschmann^{-1/4};
  both must yield the same classification.
- Stage 2 frozen homogeneous perturbation set: Sigma^2_0 in {1e-4, 1e-2, 1e-1},
  integrated over N = ln a spanning 40 e-folds toward each boundary, stiff
  (w = 1) scored, radiation (w = 1/3) as contrast.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "CBE1_STAGE0_STAGE1.json"

TOL_D = 1.0e-6
REFINE_FACTOR = 0.5
SIGMA2_SET = (1.0e-4, 1.0e-2, 1.0e-1)


# ----------------------------------------------------------------- stage 0 --

def stage0_flrw_analytic():
    """FLRW stiff collapse vs expansion: exact time-even invariants, odd parity."""
    w = 1.0
    # exact scaling solution: H = 2/(3(1+w)tau) = 1/(3 tau)
    taus = np.geomspace(1.0e-8, 1.0, 25)
    rows = []
    for tau in taus:
        H_exp = 1.0 / (3.0 * tau)          # bang branch, tau > 0
        H_col = -1.0 / (3.0 * tau)         # crunch branch, mirrored clock
        rho = 3.0 * H_exp**2               # flat: Omega = 1
        omega = rho / (3.0 * H_exp**2)
        r2 = 9.0 * (1.0 + 3.0 * w**2)      # R_mn R^mn / H^4, exact
        rows.append({
            "tau": tau,
            "omega_diff": abs(omega - 1.0),
            "sigma2_diff": 0.0,
            "r2_diff": 0.0,
            "weyl_diff": 0.0,
            "parity_sum": abs(H_exp + H_col),   # exact odd parity
            "r2_value": r2,
        })
    D = max(
        np.sqrt(r["omega_diff"]**2 + r["sigma2_diff"]**2 + r["r2_diff"]**2
                + r["weyl_diff"]**2 + r["parity_sum"]**2)
        for r in rows
    )
    return {"pair": "FLRW stiff crunch vs FLRW stiff bang", "D_max": float(D),
            "passes": bool(D == 0.0)}


def stage0_schwarzschild_control():
    """Schwarzschild/KS interior vs FLRW: must early-terminate on Sec. 9 gates."""
    p = np.array([-1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0])   # interior Kasner limit
    # H_i = p_i / tau, H = sum(p)/3tau = 1/(3tau), sigma^2 = sum (H_i - H)^2
    coeff = np.sum((p - p.mean())**2)                   # = 2/3
    sigma2_limit = (coeff) / (6.0 * (1.0 / 3.0)**2)     # sigma^2 / 6H^2
    weyl_fraction = 1.0                                 # vacuum: Ricci = 0
    flrw = {"sigma2": 0.0, "weyl_fraction": 0.0}
    gates = {
        "shear_gate_C4": bool(abs(sigma2_limit - flrw["sigma2"]) > 0.1),
        "weyl_ricci_gate_C5": bool(abs(weyl_fraction - flrw["weyl_fraction"]) > 0.1),
        "vacuum_vs_matter_ricci": True,   # Ricci = 0 vs stiff rho != 0, no matter map
    }
    return {"pair": "Schwarzschild/Kantowski-Sachs interior vs FLRW",
            "sigma2_limit": float(sigma2_limit),
            "weyl_fraction": weyl_fraction,
            "early_termination_gates": gates,
            "early_terminated": bool(all(gates.values())),
            "classification_if_declared": "CBE-C — invariant mismatch",
            "passes": bool(all(gates.values()) and abs(sigma2_limit - 1.0) < 1.0e-12)}


# ----------------------------------------------------------------- stage 1 --

def flrw_scalar_rhs_free(t, y):
    a, H, phidot = y
    return [a * H, -0.5 * phidot**2, -3.0 * H * phidot]


def flrw_scalar_rhs_constrained(direction):
    """Hamiltonian-enforced flat-FLRW evolution: H = sign * sqrt(rho/3).

    Free evolution away from the stiff boundary amplifies constraint error
    like a^4 (rho ~ a^-6 dilutes faster than curvature-like error modes
    ~ a^-2); measured 2.1e-2 at tau=1 from an rtol=1e-11 start at tau=1e-7.
    That drift is a Sec.-C7 diagnostic, recorded in the output; scoring uses
    the constrained form, which is exact for flat FLRW up to the sign choice.
    """
    def rhs(t, y):
        a, phidot = y
        H = direction * abs(phidot) / np.sqrt(6.0)
        return [a * H, -3.0 * H * phidot]
    return rhs


def flrw_invariants(H, Hdot):
    ricci2 = 9.0 * (Hdot + H**2)**2 + 3.0 * (Hdot + 3.0 * H**2)**2
    kretsch = 12.0 * ((Hdot + H**2)**2 + H**4)
    return ricci2, kretsch


def run_flrw_branch(direction, rtol):
    """direction=+1: expanding from asymptotic prescription; -1: collapse toward boundary."""
    if direction > 0:
        tau0 = 1.0e-7
        H0 = 1.0 / (3.0 * tau0)
        a0 = tau0**(1.0 / 3.0)
        phidot0 = np.sqrt(6.0) * H0            # rho = 3H^2 exactly (flat)
        span, t_eval = (tau0, 1.0), np.geomspace(tau0, 1.0, 400)
    else:
        t_start = -1.0
        H0 = 1.0 / (3.0 * t_start)             # contracting: H < 0
        a0 = 1.0
        phidot0 = np.sqrt(6.0) * abs(H0)
        span = (t_start, -1.0e-7)
        t_eval = -np.geomspace(1.0, 1.0e-7, 400)
    direction = 1.0 if direction > 0 else -1.0
    free = solve_ivp(flrw_scalar_rhs_free, span, [a0, H0, phidot0],
                     t_eval=t_eval, rtol=rtol, atol=1.0e-14, method="DOP853")
    free_residual = float(np.max(
        np.abs(3.0 * free.y[1]**2 - 0.5 * free.y[2]**2) / (0.5 * free.y[2]**2)
    ))
    sol = solve_ivp(flrw_scalar_rhs_constrained(direction), span, [a0, phidot0],
                    t_eval=t_eval, rtol=rtol, atol=1.0e-14, method="DOP853")
    a, phidot = sol.y
    H = direction * np.abs(phidot) / np.sqrt(6.0)
    Hdot = -0.5 * phidot**2
    rho = 0.5 * phidot**2
    p_over_rho = np.ones_like(rho)             # stiff: p = rho
    ricci2, kretsch = flrw_invariants(H, Hdot)
    lam_rho = rho**-0.5
    lam_R = kretsch**-0.25
    return {
        "lam_rho": lam_rho, "lam_R": lam_R, "H": H,
        "omega": rho / (3.0 * H**2), "sigma2": np.zeros_like(H),
        "r2": ricci2 / H**4, "w_eff": p_over_rho,
        "constraint": 0.0, "free_evolution_residual": free_residual,
    }


def state_distance(col, exp, clock):
    lo = max(col[clock].min(), exp[clock].min())
    hi = min(col[clock].max(), exp[clock].max())
    grid = np.geomspace(lo * 1.001, hi * 0.999, 200)

    def onto(branch, key):
        order = np.argsort(branch[clock])
        return interp1d(branch[clock][order], branch[key][order])(grid)

    comps = [
        onto(col, "omega") - onto(exp, "omega"),
        onto(col, "sigma2") - onto(exp, "sigma2"),
        onto(col, "r2") / 36.0 - onto(exp, "r2") / 36.0,
        onto(col, "w_eff") - onto(exp, "w_eff"),
        grid * (onto(col, "H") + onto(exp, "H")),   # T-parity residual
    ]
    return float(np.max(np.sqrt(np.sum(np.square(comps), axis=0))))


def stage1():
    result = {}
    for label, rtol in (("coarse", 1.0e-8), ("fine", 1.0e-11)):
        col = run_flrw_branch(-1, rtol)
        exp = run_flrw_branch(+1, rtol)
        result[label] = {
            "D_lam_rho": state_distance(col, exp, "lam_rho"),
            "D_lam_R": state_distance(col, exp, "lam_R"),
            "constraint_residual_max": max(col["constraint"], exp["constraint"]),
            "free_evolution_residual_C7_diagnostic": max(
                col["free_evolution_residual"], exp["free_evolution_residual"]
            ),
        }
    fine, coarse = result["fine"], result["coarse"]
    clocks_agree = all(
        (fine[k] < TOL_D) == (fine["D_lam_rho"] < TOL_D)
        for k in ("D_lam_rho", "D_lam_R")
    )
    refined = all(
        fine[k] <= REFINE_FACTOR * coarse[k] or (fine[k] < TOL_D and coarse[k] < TOL_D)
        for k in ("D_lam_rho", "D_lam_R")
    )
    result["passes"] = bool(
        fine["D_lam_rho"] < TOL_D and fine["D_lam_R"] < TOL_D
        and clocks_agree and refined
    )
    result["clocks_agree"] = bool(clocks_agree)
    result["refinement_closes"] = bool(refined)
    return result


# ------------------------------------------------- stage 2 (homogeneous) ----

def bianchi_rhs(w):
    def rhs(N, y):
        (sigma,) = y
        omega = 1.0 - sigma**2
        q = 2.0 * sigma**2 + 0.5 * (1.0 + 3.0 * w) * omega
        return [(q - 2.0) * sigma]
    return rhs


def stage2_homogeneous():
    rows = []
    for sigma2_0 in SIGMA2_SET:
        s0 = np.sqrt(sigma2_0)
        for w, role in ((1.0, "scored-stiff"), (1.0 / 3.0, "contrast-radiation")):
            # toward the singular boundary: N decreasing (a -> 0)
            sol = solve_ivp(bianchi_rhs(w), (0.0, -40.0), [s0],
                            rtol=1.0e-12, atol=1.0e-16, dense_output=True)
            s_end = float(sol.y[0][-1])
            rows.append({
                "w": w, "role": role, "sigma2_initial": sigma2_0,
                "sigma2_at_boundary": s_end**2,
                "drift": abs(s_end**2 - sigma2_0),
                "isotropized": s_end**2 < 1.0e-8,
                "shear_dominated": s_end**2 > 0.9,
            })
    stiff = [r for r in rows if r["role"] == "scored-stiff"]
    stiff_conserved = all(r["drift"] < 1.0e-10 for r in stiff)
    radiation = [r for r in rows if r["role"] == "contrast-radiation"]
    radiation_dominated = all(r["shear_dominated"] for r in radiation)
    return {
        "rows": rows,
        "stiff_sigma2_conserved": bool(stiff_conserved),
        "radiation_shear_dominates_boundary": bool(radiation_dominated),
        "reading": (
            "Stiff/scalar collapse carries Sigma^2 to the boundary unchanged: the "
            "perturbed collapse does NOT return to the unperturbed expanding "
            "target's basin, but matches the expanding solution of the SAME "
            "Sigma^2 stratum. Open-basin toward a single target fails; pairwise "
            "type-equivalence survives stratum-by-stratum."
        ),
    }


def main():
    s0a = stage0_flrw_analytic()
    s0b = stage0_schwarzschild_control()
    s1 = stage1()
    s2 = stage2_homogeneous()
    stage0_pass = s0a["passes"] and s0b["passes"]
    classification = (
        "CBE-B — exact-control equivalence with Sigma^2-stratified basin "
        "(homogeneous scope; inhomogeneous Stage 2 and Stage 3 not executed)"
        if stage0_pass and s1["passes"] and s2["stiff_sigma2_conserved"]
        else "CBE-D — asymptotic numerical nonclosure"
    )
    result = {
        "program": "BB-CBE-1",
        "program_doc": "Discussion/BB_COLLAPSE_BANG_SINGULAR_BOUNDARY_EQUIVALENCE_PROGRAM_2026-08-04.md",
        "declared_action": "GR + one massless scalar (stiff, w=1), 8piG=1, flat homogeneous control",
        "frozen": {
            "weights": "unit weight on (Omega, Sigma^2, R2/36, w_eff, parity)",
            "tolerance_D": TOL_D, "refinement_factor": REFINE_FACTOR,
            "sigma2_perturbation_set": list(SIGMA2_SET),
            "clocks": ["lam_rho = rho^-1/2", "lam_R = Kretschmann^-1/4"],
        },
        "stage0": {"flrw_pair": s0a, "schwarzschild_control": s0b,
                   "passes": stage0_pass},
        "stage1": s1,
        "stage2_homogeneous": s2,
        "classification": classification,
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "stage0_passes": stage0_pass,
        "stage1_passes": s1["passes"],
        "stage1_D_fine": {k: s1["fine"][k] for k in ("D_lam_rho", "D_lam_R")},
        "stage2_stiff_sigma2_conserved": s2["stiff_sigma2_conserved"],
        "stage2_radiation_shear_dominates": s2["radiation_shear_dominates_boundary"],
        "classification": classification,
    }, indent=2))


if __name__ == "__main__":
    main()
