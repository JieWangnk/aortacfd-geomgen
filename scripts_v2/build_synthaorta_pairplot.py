#!/usr/bin/env python3
"""Joint-distribution pairplot for the cli_v2 Sobol cube.

Companion to ``build_synthaorta_compat_figure.py`` — the marginal violin
shows our 1-D distributions match the SynthAorta Table I distributions,
but it tells us nothing about *pairwise* coverage. A scrambled Sobol
sequence in 7-D could in principle introduce phantom correlations
between dimensions that the violin would miss.

This script builds a 7×7 scatter matrix:

  - **diagonal**: marginal histogram of v2 samples + reference PDF curve
  - **lower triangle**: v2 scatter (black dots) over a 2-D reference
    contour (faint yellow band — KDE of independent Monte-Carlo from
    the same marginals)
  - **upper triangle**: Pearson r between the pair, colour-coded by
    |r|. Threshold for "consistent with independence" at N=100 samples
    is roughly |r| < 2/√N ≈ 0.20 (two-sigma).

Validation criteria — for the SynthAorta parametrisation the seven
parameters are INDEPENDENT by construction, so:

  * All off-diagonal Pearson |r| should be < ~0.20
  * Lower-triangle scatters should look like the reference contour:
    no diagonals, no banding, no holes

A summary JSON records the worst-case correlation and whether the
"all independent" hypothesis holds for the supplied spec + N.

Usage::

    python scripts_v2/build_synthaorta_pairplot.py \\
        --spec specs_v2/sample_sobol_synthaorta_v2.json \\
        --out figures/synthaorta_pairplot.png
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = Path(__file__).resolve().parent
for p in (HERE, SCRIPTS_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from cli_v2 import PARAMETERS, sample_cases  # noqa: E402
from build_synthaorta_compat_figure import reference_samples  # noqa: E402


# Independence threshold for Pearson r at sample size N — two-sigma rule of thumb
def independence_threshold(n: int) -> float:
    return 2.0 / np.sqrt(max(n, 4))


def correlation_colour(r: float, threshold: float):
    """Background colour for the upper-triangle Pearson-r cells.

    Green: |r| < 0.5*threshold (clearly independent)
    Yellow: |r| in [0.5*threshold, threshold] (borderline)
    Red: |r| > threshold (suspicious — investigate)
    """
    a = abs(r)
    if a < 0.5 * threshold:
        return "#d4edda"
    if a < threshold:
        return "#fff3cd"
    return "#f8d7da"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spec", type=Path,
                   default=HERE / "specs_v2" / "sample_sobol_synthaorta_v2.json")
    p.add_argument("--out", type=Path,
                   default=HERE / "figures" / "synthaorta_pairplot.png")
    p.add_argument("--n-ref", type=int, default=10_000)
    p.add_argument("--ref-seed", type=int, default=2027)
    args = p.parse_args()

    print(f"Loading v2 samples from {args.spec.name} …")
    spec = json.loads(args.spec.read_text())
    cases = sample_cases(spec)
    n_v2 = len(cases)

    # Parameter order from the spec — preserves the order columns appear in cli_v2
    names = list(cases[0].keys())
    n_params = len(names)
    print(f"  → N={n_v2}, dim={n_params}: {names}")

    # Collect v2 samples + per-parameter reference samples
    v2 = np.column_stack([np.asarray([c[name] for c in cases], dtype=float)
                          for name in names])

    overrides = spec.get("distribution_overrides", {})

    def _effective_dist(name: str):
        return overrides.get(name) or PARAMETERS[name]["default_dist"]

    ref_cols = []
    for j, name in enumerate(names):
        ref_cols.append(reference_samples(_effective_dist(name), n=args.n_ref,
                                           seed=args.ref_seed + j))
    ref = np.column_stack(ref_cols)

    # Pearson r matrix on the v2 samples (off-diagonal is what we care about)
    r_mat = np.corrcoef(v2.T)
    threshold = independence_threshold(n_v2)
    print(f"\nIndependence threshold for N={n_v2}: |r| < {threshold:.3f}")
    worst_pair = None
    worst_r = 0.0
    for i in range(n_params):
        for j in range(i + 1, n_params):
            if abs(r_mat[i, j]) > abs(worst_r):
                worst_r = r_mat[i, j]
                worst_pair = (names[i], names[j])

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(n_params, n_params, figsize=(2.0 * n_params, 2.0 * n_params))
    rng = np.random.default_rng(0)

    for i in range(n_params):
        for j in range(n_params):
            ax = axes[i, j]

            if i == j:
                # Diagonal: marginal v2 hist + reference PDF curve
                vals_v2 = v2[:, i]
                vals_ref = ref[:, i]
                ax.hist(vals_v2, bins=12, density=True, color="#1f77b4",
                        edgecolor="black", alpha=0.65, label="v2")
                # Reference KDE as a smooth line
                kde = stats.gaussian_kde(vals_ref)
                xs = np.linspace(vals_ref.min(), vals_ref.max(), 200)
                ax.plot(xs, kde(xs), color="#bcbd22", lw=2, label="ref")
                ax.set_yticks([])
                if i == 0:
                    ax.legend(fontsize=7, loc="upper right")

            elif i > j:
                # Lower triangle: v2 scatter over reference 2-D contour
                x_ref = ref[:, j]
                y_ref = ref[:, i]
                try:
                    kde2 = stats.gaussian_kde(np.vstack([x_ref, y_ref]))
                    # Evaluate on a coarse grid for contour
                    pad_x = 0.05 * (x_ref.max() - x_ref.min())
                    pad_y = 0.05 * (y_ref.max() - y_ref.min())
                    xs = np.linspace(x_ref.min() - pad_x, x_ref.max() + pad_x, 40)
                    ys = np.linspace(y_ref.min() - pad_y, y_ref.max() + pad_y, 40)
                    X, Y = np.meshgrid(xs, ys)
                    Z = kde2(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
                    ax.contourf(X, Y, Z, levels=6, cmap="YlOrBr", alpha=0.45)
                except np.linalg.LinAlgError:
                    pass  # singular covariance — skip contour
                ax.scatter(v2[:, j], v2[:, i], s=10, c="black", alpha=0.7, zorder=3)

            else:  # i < j: upper triangle, correlation
                r = r_mat[i, j]
                bg = correlation_colour(r, threshold)
                ax.set_facecolor(bg)
                ax.text(0.5, 0.55, f"r = {r:+.3f}", ha="center", va="center",
                        fontsize=12, fontweight="bold", transform=ax.transAxes)
                badge = "indep ✓" if abs(r) < 0.5 * threshold else (
                        "borderline" if abs(r) < threshold else "SUSPECT")
                ax.text(0.5, 0.30, badge, ha="center", va="center",
                        fontsize=9, transform=ax.transAxes,
                        color="#155724" if "✓" in badge else
                              ("#856404" if "borderline" in badge else "#721c24"))
                ax.set_xticks([])
                ax.set_yticks([])

            # Axis labels only on the bottom row and leftmost column
            if i == n_params - 1:
                ax.set_xlabel(names[j], fontsize=8, rotation=20, ha="right")
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel(names[i], fontsize=8)
            else:
                ax.set_yticklabels([])
            ax.tick_params(axis="both", labelsize=7)

    # Summary title
    n_pass = sum(1 for i in range(n_params) for j in range(i + 1, n_params)
                 if abs(r_mat[i, j]) < threshold)
    n_pairs = n_params * (n_params - 1) // 2

    verdict = "ALL INDEPENDENT" if n_pass == n_pairs else (
              f"{n_pairs - n_pass}/{n_pairs} pair(s) above |r|={threshold:.2f}")
    fig.suptitle(
        "SynthAorta pairplot — cli_v2 Sobol joint-distribution check\n"
        f"v2 Sobol N={n_v2} • independence threshold |r|<{threshold:.3f} "
        f"(2/√N) • worst pair: {worst_pair[0]}↔{worst_pair[1]} r={worst_r:+.3f} "
        f"• verdict: {verdict}",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=130)
    print(f"\nWrote {args.out}")

    # JSON summary
    summary = {
        "spec": str(args.spec.name),
        "v2_sample_size": n_v2,
        "n_parameters": n_params,
        "parameters": names,
        "independence_threshold": float(threshold),
        "pearson_r_matrix": r_mat.tolist(),
        "worst_offdiag_pair": list(worst_pair),
        "worst_offdiag_r": float(worst_r),
        "n_pairs_within_threshold": int(n_pass),
        "n_pairs_total": int(n_pairs),
        "verdict": verdict,
    }
    summary_path = args.out.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_path}")

    # Console summary
    print("\nPearson r matrix (off-diagonal):")
    for i in range(n_params):
        for j in range(i + 1, n_params):
            badge = ("OK" if abs(r_mat[i, j]) < 0.5 * threshold else
                     "borderline" if abs(r_mat[i, j]) < threshold else "SUSPECT")
            print(f"  {names[i]:20s} ↔ {names[j]:20s}  r={r_mat[i, j]:+.3f}  [{badge}]")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
