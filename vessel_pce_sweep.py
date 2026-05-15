#!/usr/bin/env python3
"""PCE-driven design of experiments for the parametric vessel generator.

Pipeline
--------
1. Draw N Sobol quasi-random samples in [0,1]^5.
2. Map each sample to a physical 5-parameter vessel (reference diameter,
   stenosis severity, stenosis length, shape sigma, centreline curvature).
3. Generate one watertight STL per sample via ``vessel_generator.generate_vessel``.
4. Evaluate an analytical proxy QoI (Young-Tsai-inspired pressure drop), which
   lets us validate the full PCE + Sobol pipeline before committing to real CFD.
5. Fit a sparse Polynomial Chaos Expansion (LARS + LOO-CV) on Legendre basis.
6. Report mean, variance, Sobol first-order and total indices per QoI.

The analytical QoI is a stand-in; replace ``evaluate_qoi`` with a CFD post-
processor that reads simulation output and returns the same dict keys.

Run
---
    python vessel_pce_sweep.py --n 64 --out sweep_out/
    python vessel_pce_sweep.py --n 64 --no-stl   # parameters + QoI only
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from vessel_generator import (
    Stenosis,
    arc_centreline,
    generate_vessel,
    straight_centreline,
)


# =============================================================================
# Parameter space  (5 active dimensions; Sudret hyperbolic truncation friendly)
# =============================================================================

# Each entry: (name, lo, hi, unit)
PARAM_SPEC = [
    ("d_ref",            4.0,   10.0,  "mm"),       # reference diameter
    ("severity",         0.0,   0.85,  ""),         # area reduction
    ("L_stenosis",       10.0,  40.0,  "mm"),       # stenosis axial length
    ("sigma_shape",      0.0,   1.0,   ""),         # shelf(0)..hourglass(1)
    ("curvature_inv",    0.0,   0.5,   "1/mm_norm"),  # 1/(R/L); 0=straight
]

L_TOTAL = 100.0   # vessel centreline length (mm); fixed for this demo


def unit_to_physical(u: np.ndarray) -> dict:
    """Map a single unit vector (5,) in [0,1]^5 to a parameter dict."""
    out = {}
    for v, (name, lo, hi, _unit) in zip(u, PARAM_SPEC):
        out[name] = lo + v * (hi - lo)
    return out


def build_case_geometry(params: dict) -> tuple[np.ndarray, float, Stenosis]:
    """Turn a parameter dict into (centreline, radius, stenosis)."""
    # Curvature: ``curvature_inv`` = L_total / R.  0 means straight.
    k = params["curvature_inv"]
    if k < 1e-4:
        centreline = straight_centreline(length=L_TOTAL, n=400)
    else:
        R = L_TOTAL / k
        arc_deg = math.degrees(L_TOTAL / R)
        centreline = arc_centreline(R=R, arc_deg=arc_deg, n=400)
    radius = 0.5 * params["d_ref"]
    stenosis = Stenosis(
        position=0.5,
        severity=params["severity"],
        length=params["L_stenosis"],
        shape="power_law",
        sigma=params["sigma_shape"],
    )
    return centreline, radius, stenosis


# =============================================================================
# Sobol experimental design
# =============================================================================

def sobol_sample(n: int, d: int, seed: int = 42) -> np.ndarray:
    """Draw n scrambled Sobol points in [0,1]^d."""
    from scipy.stats.qmc import Sobol
    engine = Sobol(d=d, scramble=True, seed=seed)
    return engine.random(n)


# =============================================================================
# Analytical proxy QoI  (replace with CFD post-processor for real work)
# =============================================================================

def evaluate_qoi(params: dict,
                 Q_inflow_ml_s: float = 5.0,
                 mu: float = 0.004,
                 rho: float = 1060.0) -> dict:
    """Young & Tsai (1973) flavoured pressure-drop estimate.

    ``dP = K_v * mu * Q / (pi * R0^3) + 0.5 * K_t * rho * (Q/A0)^2 * (A0/A_t - 1)^2``

    Augmented with small penalty terms for curvature and shape so that every
    parameter has non-zero sensitivity (useful for Sobol pipeline validation).

    Returns a dict of QoI scalars.
    """
    R0 = 0.5 * params["d_ref"] * 1e-3       # m
    alpha = params["severity"]
    L = params["L_stenosis"] * 1e-3         # m
    sigma = params["sigma_shape"]
    k = params["curvature_inv"]

    A0 = math.pi * R0 ** 2
    A_t = A0 * (1.0 - alpha)

    Q = Q_inflow_ml_s * 1e-6                # m^3/s
    K_v = 32.0 * L / (2.0 * R0)             # length-dependent viscous factor
    K_t = 1.52                              # Young-Tsai turbulent constant

    dP_visc = K_v * mu * Q / (math.pi * R0 ** 3)
    dP_iner = 0.5 * K_t * rho * (Q / A0) ** 2 * max(0.0, (A0 / A_t - 1.0)) ** 2

    # Shape sensitivity: shelf (sigma->0) has higher recovery losses
    shape_mult = 1.0 + 0.25 * (1.0 - sigma)
    dP_iner *= shape_mult

    # Curvature penalty (Dean number proxy: scales with curvature and Re)
    Re = 2.0 * R0 * (Q / A0) * rho / mu
    dP_curv = 0.1 * rho * (Q / A0) ** 2 * k * math.sqrt(Re) / 100.0

    dP_total = dP_visc + dP_iner + dP_curv

    # Convert to mmHg
    dP_mmHg = dP_total / 133.322

    # Peak WSS proxy at throat (Poiseuille)
    R_t = R0 * math.sqrt(1.0 - alpha)
    wss_peak = 4.0 * mu * Q / (math.pi * R_t ** 3)

    return {
        "dP_mmHg":   dP_mmHg,
        "wss_peak":  wss_peak,
        "Re":        Re,
    }


# =============================================================================
# Sparse PCE (Legendre basis, hyperbolic truncation, LARS + LOO)
# =============================================================================

def build_multi_index(M: int, p_max: int, q: float = 0.75) -> np.ndarray:
    def rec(current, result):
        if len(current) == M:
            a = np.array(current, dtype=float)
            if a.sum() == 0:
                result.append(list(current))
                return
            q_norm = np.sum(a ** q) ** (1.0 / q) if q < 1.0 else a.sum()
            if q_norm <= p_max + 1e-10:
                result.append(list(current))
            return
        for v in range(p_max + 1):
            rec(current + [v], result)
    out = []
    rec([], out)
    return np.array(sorted(out, key=lambda a: (sum(a), tuple(a))), dtype=int)


def legendre_basis(X: np.ndarray, alphas: np.ndarray) -> np.ndarray:
    """Evaluate normalised multivariate Legendre basis at rows of X in [-1,1]^M."""
    from numpy.polynomial.legendre import legval
    n, M = X.shape
    P = len(alphas)
    A = np.ones((n, P))
    for j, alpha in enumerate(alphas):
        for k in range(M):
            if alpha[k] > 0:
                c = np.zeros(alpha[k] + 1)
                c[-1] = 1.0
                vals = legval(X[:, k], c) * math.sqrt(2.0 * alpha[k] + 1.0)
                A[:, j] *= vals
    return A


def fit_sparse_pce(A: np.ndarray, y: np.ndarray):
    """Fit via LassoLarsCV on centred, scaled y."""
    from sklearn.linear_model import LassoLarsCV
    y_mean = float(np.mean(y))
    y_std = float(np.std(y))
    if y_std < 1e-15:
        coeffs = np.zeros(A.shape[1])
        coeffs[0] = y_mean
        return coeffs, np.array([0]), 0.0
    y_norm = (y - y_mean) / y_std
    try:
        model = LassoLarsCV(cv=min(10, len(y)), max_iter=500)
    except TypeError:
        model = LassoLarsCV(cv=min(10, len(y)))
    model.fit(A[:, 1:], y_norm)
    coeffs = np.zeros(A.shape[1])
    coeffs[0] = y_mean
    coeffs[1:] = model.coef_ * y_std
    active = np.where(np.abs(coeffs) > 1e-12)[0]
    y_pred = A @ coeffs
    loo = float(np.mean((y - y_pred) ** 2) / max(1e-30, np.var(y)))
    return coeffs, active, loo


def sobol_from_pce(coeffs: np.ndarray, alphas: np.ndarray, M: int):
    D = float(np.sum(coeffs[1:] ** 2))
    if D < 1e-20:
        return np.zeros(M), np.zeros(M), 0.0
    S1 = np.zeros(M)
    ST = np.zeros(M)
    for j in range(1, len(coeffs)):
        a = alphas[j]
        active = set(np.where(a > 0)[0])
        if len(active) == 1:
            i = next(iter(active))
            S1[i] += coeffs[j] ** 2
        for i in active:
            ST[i] += coeffs[j] ** 2
    return S1 / D, ST / D, D


# =============================================================================
# Driver
# =============================================================================

def run_sweep(n: int, out_dir: str, write_stl: bool, seed: int) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    M = len(PARAM_SPEC)

    print(f"Drawing {n} Sobol samples in {M}-D unit hypercube...")
    U = sobol_sample(n=n, d=M, seed=seed)

    print(f"Generating {'geometries + ' if write_stl else ''}analytical QoIs...")
    cases = []
    qoi_names = None
    for i, u in enumerate(U):
        params = unit_to_physical(u)
        case_id = f"case_{i:03d}"
        stl_path = None
        if write_stl:
            centreline, radius, stenosis = build_case_geometry(params)
            stl_path = os.path.join(out_dir, f"{case_id}.stl")
            generate_vessel(centreline, radius=radius, stenosis=stenosis,
                            n_sectors=48, out=stl_path)
        qois = evaluate_qoi(params)
        if qoi_names is None:
            qoi_names = list(qois.keys())
        cases.append({
            "case_id": case_id,
            "unit_coords": u.tolist(),
            "params": params,
            "qois": qois,
            "stl": stl_path,
        })

    manifest = {
        "n_samples": n,
        "M": M,
        "params": [{"name": p[0], "lo": p[1], "hi": p[2], "unit": p[3]} for p in PARAM_SPEC],
        "qoi_names": qoi_names,
        "cases": cases,
    }
    manifest_path = os.path.join(out_dir, "manifest.json")
    Path(manifest_path).write_text(json.dumps(manifest, indent=2))
    print(f"  wrote {manifest_path}")
    return manifest


def build_and_report_pce(manifest: dict, out_dir: str, p_max: int = 4) -> dict:
    U = np.array([c["unit_coords"] for c in manifest["cases"]])
    X = 2.0 * U - 1.0     # map [0,1] -> [-1,1] for Legendre
    M = manifest["M"]
    n = len(manifest["cases"])

    results = {}
    for qoi in manifest["qoi_names"]:
        y = np.array([c["qois"][qoi] for c in manifest["cases"]])
        print()
        print(f"=== QoI: {qoi} ===")
        print(f"  y range [{y.min():.3g}, {y.max():.3g}], mean {y.mean():.3g}, std {y.std():.3g}")

        best = None
        for p in range(1, p_max + 1):
            alphas = build_multi_index(M, p, q=0.75)
            A = legendre_basis(X, alphas)
            coeffs, active, loo = fit_sparse_pce(A, y)
            print(f"  p={p}: P={len(alphas)} terms, {len(active)} active, "
                  f"LOO={loo:.3e}, n/P={n/len(alphas):.2f}")
            if best is None or loo < best[3]:
                best = (p, alphas, coeffs, loo, active)

        p_best, alphas, coeffs, loo, active = best
        S1, ST, D = sobol_from_pce(coeffs, alphas, M)

        print(f"  best p={p_best}, LOO={loo:.3e}, total variance D={D:.3e}")
        print(f"  {'parameter':<15s} {'S1':>8s} {'ST':>8s}")
        for i, (name, *_rest) in enumerate(PARAM_SPEC):
            bar = "#" * int(S1[i] * 40)
            print(f"  {name:<15s} {S1[i]:>8.4f} {ST[i]:>8.4f}  {bar}")

        results[qoi] = {
            "best_p": p_best,
            "loo": loo,
            "D": D,
            "mean": float(coeffs[0]),
            "S1": {p[0]: float(v) for p, v in zip(PARAM_SPEC, S1)},
            "ST": {p[0]: float(v) for p, v in zip(PARAM_SPEC, ST)},
            "n_active": int(len(active)),
        }

    results_path = os.path.join(out_dir, "pce_results.json")
    Path(results_path).write_text(json.dumps(results, indent=2))
    print()
    print(f"wrote {results_path}")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=64,
                        help="number of Sobol samples (power of 2 recommended)")
    parser.add_argument("--out", default="vessel_sweep_out",
                        help="output directory")
    parser.add_argument("--no-stl", action="store_true",
                        help="skip STL generation (parameters + QoI only)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--p-max", type=int, default=4,
                        help="max polynomial degree for PCE")
    args = parser.parse_args()

    manifest = run_sweep(n=args.n, out_dir=args.out,
                         write_stl=not args.no_stl, seed=args.seed)
    build_and_report_pce(manifest, args.out, p_max=args.p_max)


if __name__ == "__main__":
    main()
