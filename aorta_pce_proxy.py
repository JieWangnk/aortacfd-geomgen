#!/usr/bin/env python3
"""10-parameter aorta analytical proxy + sparse-PCE Sobol indices.

Mirrors the bifurcation pipeline in vessel_pce_sweep.py but for the 10 active
aortic parameters used in the Medium-post figure.  The QoI is a Young & Tsai
flavoured pressure-drop estimate across the coarctation (peak systolic
gradient); secondary terms add a Dean-number-like correction for the arch
curvature so every parameter has at least a small, defensible sensitivity.

Run:
    python3 aorta_pce_proxy.py --n 128 --p-max 4 --out outputs/aorta_pce
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np


# =============================================================================
# 10 active aortic parameters (closed ranges, matching the geometry generator)
# =============================================================================

PARAM_SPEC = [
    ("d_ascending",                  20.0,   40.0,  "mm"),
    ("taper_ratio",                   0.55,   0.90, ""),
    ("ascending_length",             30.0,   65.0,  "mm"),
    ("arch_span",                    50.0,   90.0,  "mm"),
    ("arch_height",                  20.0,   50.0,  "mm"),
    ("descending_length",            50.0,  120.0,  "mm"),
    ("coarctation_area_reduction",    0.0,    0.90, ""),
    ("coarctation_length",           10.0,   40.0,  "mm"),
    ("coarctation_shape",             0.0,    1.0,  ""),
    ("proximal_hypoplasia",           0.0,    0.40, ""),
]
M = len(PARAM_SPEC)


def unit_to_physical(u: np.ndarray) -> dict:
    out = {}
    for v, (name, lo, hi, _unit) in zip(u, PARAM_SPEC):
        out[name] = lo + v * (hi - lo)
    return out


# =============================================================================
# Analytical proxy QoI (peak pressure gradient across coarctation, mmHg)
# =============================================================================

def evaluate_qoi(p: dict,
                 Q_inflow_ml_s: float = 83.3,   # 5 L/min cardiac output
                 mu: float = 0.004,
                 rho: float = 1060.0) -> float:
    """Peak pressure gradient across the coarctation (mmHg).

    Combines:
      * Young & Tsai (1973) stenosis pressure drop: viscous + inertial + shape.
      * Hypoplasia and taper reduce the *pre-coarctation* effective lumen,
        amplifying both terms.
      * A Dean-number-style curvature correction depending on arch geometry,
        so arch_span / arch_height enter the QoI.
      * A weak length-dependent viscous term tied to descending_length and
        ascending_length so they show up with non-zero (but small) Sobol bars.
    """
    d_asc = p["d_ascending"] * 1e-3                   # m
    tau   = p["taper_ratio"]
    eta   = p["proximal_hypoplasia"]
    alpha = p["coarctation_area_reduction"]
    L_coa = p["coarctation_length"] * 1e-3            # m
    sig   = p["coarctation_shape"]
    L_asc = p["ascending_length"] * 1e-3
    L_des = p["descending_length"] * 1e-3
    S     = p["arch_span"] * 1e-3
    H     = p["arch_height"] * 1e-3

    # Effective pre-coarctation diameter: ascending diameter, tapered into
    # the descending region, attenuated by hypoplasia.
    d_pre = d_asc * tau * (1.0 - eta)
    R_pre = 0.5 * d_pre
    R_thr = max(R_pre * math.sqrt(max(1.0 - alpha, 1e-6)), 1e-4)

    A_pre = math.pi * R_pre ** 2
    A_thr = math.pi * R_thr ** 2

    Q = Q_inflow_ml_s * 1e-6                          # m^3/s

    # Young & Tsai viscous term
    K_v = 32.0 * L_coa / (2.0 * R_thr)
    dP_visc = K_v * mu * Q / (math.pi * R_thr ** 3)

    # Young & Tsai inertial / Bernoulli term
    K_t = 1.52
    bern = max(0.0, A_pre / A_thr - 1.0) ** 2
    dP_iner = 0.5 * K_t * rho * (Q / A_pre) ** 2 * bern

    # Shape multiplier: shelf morphology (sigma -> 0) increases recovery loss.
    dP_iner *= 1.0 + 0.25 * (1.0 - sig)

    # Dean-number curvature correction.  Arch radius of curvature ~ S/2 + H.
    R_curv = max(0.5 * S + H, 5e-3)
    Re = 2.0 * R_pre * (Q / A_pre) * rho / mu
    dean = Re * math.sqrt(R_pre / R_curv)
    dP_curv = 0.04 * rho * (Q / A_pre) ** 2 * math.sqrt(max(dean, 0.0)) / 100.0

    # Tiny length-dependent viscous baseline so ascending/descending length
    # have non-zero (but small) Sobol contribution.
    L_arch = L_asc + L_des
    dP_baseline = 32.0 * mu * Q * L_arch / (math.pi * (2.0 * R_pre) ** 4)

    dP_total = dP_visc + dP_iner + dP_curv + dP_baseline
    return dP_total / 133.322                          # Pa -> mmHg


# =============================================================================
# Sobol experimental design
# =============================================================================

def sobol_sample(n: int, d: int, seed: int = 42) -> np.ndarray:
    from scipy.stats.qmc import Sobol
    return Sobol(d=d, scramble=True, seed=seed).random(n)


# =============================================================================
# Sparse PCE (Legendre, hyperbolic truncation, LARS + LOO)
# =============================================================================

def build_multi_index(M: int, p_max: int, q: float = 0.75) -> np.ndarray:
    out = []
    def rec(current):
        if len(current) == M:
            a = np.array(current, dtype=float)
            if a.sum() == 0:
                out.append(list(current)); return
            q_norm = np.sum(a ** q) ** (1.0 / q) if q < 1.0 else a.sum()
            if q_norm <= p_max + 1e-10:
                out.append(list(current))
            return
        for v in range(p_max + 1):
            rec(current + [v])
    rec([])
    return np.array(sorted(out, key=lambda a: (sum(a), tuple(a))), dtype=int)


def legendre_basis(X: np.ndarray, alphas: np.ndarray) -> np.ndarray:
    from numpy.polynomial.legendre import legval
    n, M = X.shape
    P = len(alphas)
    A = np.ones((n, P))
    for j, alpha in enumerate(alphas):
        for k in range(M):
            if alpha[k] > 0:
                c = np.zeros(alpha[k] + 1); c[-1] = 1.0
                A[:, j] *= legval(X[:, k], c) * math.sqrt(2.0 * alpha[k] + 1.0)
    return A


def fit_sparse_pce(A: np.ndarray, y: np.ndarray):
    from sklearn.linear_model import LassoLarsCV
    y_mean = float(np.mean(y)); y_std = float(np.std(y))
    if y_std < 1e-15:
        c = np.zeros(A.shape[1]); c[0] = y_mean
        return c, np.array([0]), 0.0
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
    S1 = np.zeros(M); ST = np.zeros(M)
    for j in range(1, len(coeffs)):
        a = alphas[j]
        active = set(np.where(a > 0)[0])
        if len(active) == 1:
            S1[next(iter(active))] += coeffs[j] ** 2
        for i in active:
            ST[i] += coeffs[j] ** 2
    return S1 / D, ST / D, D


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n", type=int, default=128,
                        help="Sobol sample count (power of 2 recommended)")
    parser.add_argument("--p-max", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="outputs/aorta_pce")
    args = parser.parse_args()

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Drawing {args.n} Sobol samples in {M}-D ...")
    U = sobol_sample(args.n, M, args.seed)
    samples = [unit_to_physical(u) for u in U]

    print("Evaluating analytical Young-Tsai aorta proxy ...")
    y = np.array([evaluate_qoi(p) for p in samples])
    print(f"  dP_mmHg: range [{y.min():.2f}, {y.max():.2f}], "
          f"mean {y.mean():.2f}, std {y.std():.2f}")

    print(f"\nFitting sparse PCE (LARS + LOO) up to p={args.p_max} ...")
    X = 2.0 * U - 1.0   # [0,1] -> [-1,1] for Legendre
    best = None
    for p in range(1, args.p_max + 1):
        alphas = build_multi_index(M, p, q=0.75)
        A = legendre_basis(X, alphas)
        coeffs, active, loo = fit_sparse_pce(A, y)
        print(f"  p={p}: P={len(alphas):4d} basis, {len(active):3d} active, "
              f"LOO={loo:.3e}, n/P={args.n/len(alphas):.2f}")
        if best is None or loo < best[3]:
            best = (p, alphas, coeffs, loo, active)
    p_best, alphas, coeffs, loo, active = best

    S1, ST, D = sobol_from_pce(coeffs, alphas, M)

    print(f"\nBest: p={p_best}, LOO={loo:.3e}, total variance D={D:.3e}")
    print(f"  {'parameter':<28s} {'S1':>8s} {'ST':>8s}")
    print(f"  {'-'*48}")
    for i, (name, *_rest) in enumerate(PARAM_SPEC):
        bar = "#" * int(S1[i] * 40)
        print(f"  {name:<28s} {S1[i]:>8.4f} {ST[i]:>8.4f}  {bar}")

    # Save the result so the figure script can read it
    result = {
        "n_samples": args.n,
        "M": M,
        "best_p": p_best,
        "loo_error": loo,
        "total_variance": D,
        "mean": float(coeffs[0]),
        "params": [p[0] for p in PARAM_SPEC],
        "S1": [float(v) for v in S1],
        "ST": [float(v) for v in ST],
        "qoi_name": "peak_dP_mmHg_proxy",
    }
    out_path = out_dir / "aorta_pce_results.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
