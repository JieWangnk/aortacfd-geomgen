#!/usr/bin/env python3
"""Build the SynthAorta-comparability violin figure for cli_v2.

For each of the four SynthAorta-validated parameters (r_ascending, r_arch,
r_descending, arch_R_c), plots a side-by-side violin comparing:

  - "v2 Sobol (N=100)"        — what cli_v2 actually produces from
                                 specs_v2/sample_sobol_synthaorta_v2.json
  - "Published dist (N=10000)" — fresh Monte-Carlo reference drawn from
                                 the same paper Table I distribution with
                                 numpy default_rng (independent of our
                                 Sobol mapping)

Plus a fifth subplot for arch_angle_deg labelled as an *extension*
parameter (no SynthAorta equivalent) so the validation gap is visible.

A two-sample Kolmogorov-Smirnov test statistic and its p-value annotate
each subplot. The interpretation: a large p-value (> 0.05) means we
cannot reject the hypothesis that the two samples come from the same
distribution — i.e., our Sobol mapping reproduces the published
distribution.

Usage::

    python scripts_v2/build_synthaorta_compat_figure.py \\
        --spec specs_v2/sample_sobol_synthaorta_v2.json \\
        --out figures/synthaorta_compat_violin.png

Pure Python — no Blender, no STL generation. Uses only the sampler +
distribution mapping inside cli_v2 (which is what we want to validate).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent.parent  # repo root
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cli_v2 import PARAMETERS, sample_cases  # noqa: E402


# Parameters with a SynthAorta-published distribution (the ones we can validate)
SYNTHAORTA_PARAMS = ["r_ascending", "r_arch", "r_descending", "arch_R_c"]
# Extension parameter (no SynthAorta equivalent — gets a marked subplot)
EXTENSION_PARAMS = ["arch_angle_deg"]


def reference_samples(dist: dict, n: int, seed: int) -> np.ndarray:
    """Independent reference draw from the published distribution.

    Uses scipy's frozen-distribution .rvs() with a different RNG to
    SciPy's Sobol — i.e. a fully independent sampling path. Truncation
    is applied by rejection sampling so the truncated mass is correctly
    re-normalised (unlike clip-based truncation).
    """
    rng = np.random.default_rng(seed)
    low = dist.get("low", -np.inf)
    high = dist.get("high", np.inf)

    if dist["type"] == "normal":
        rv = stats.norm(loc=dist["mean"], scale=dist["std"])
    elif dist["type"] == "gumbel":
        rv = stats.gumbel_r(loc=dist["loc"], scale=dist["scale"])
    elif dist["type"] == "uniform":
        rv = stats.uniform(loc=dist["low"], scale=dist["high"] - dist["low"])
    else:
        raise ValueError(f"Unknown dist type {dist['type']!r}")

    out = []
    # Generous oversample then rejection-truncate
    chunk = max(2 * n, 4096)
    while len(out) < n:
        raw = rv.rvs(size=chunk, random_state=rng)
        keep = raw[(raw >= low) & (raw <= high)]
        out.extend(keep.tolist())
    return np.asarray(out[:n], dtype=float)


def collect_v2_samples(spec_path: Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Run cli_v2.sample_cases on the spec; return (per-param samples, spec)."""
    spec = json.loads(spec_path.read_text())
    cases = sample_cases(spec)
    out: dict[str, np.ndarray] = {}
    for name in cases[0]:
        out[name] = np.asarray([c[name] for c in cases], dtype=float)
    return out, spec


def _split_violin(ax, data_left, data_right, label_left, label_right,
                  colour_left, colour_right) -> None:
    """Draw a split (two-half) violin at x=0 with the two datasets."""
    parts_l = ax.violinplot([data_left], positions=[0], widths=0.85,
                            showmeans=False, showmedians=False, showextrema=False)
    parts_r = ax.violinplot([data_right], positions=[0], widths=0.85,
                            showmeans=False, showmedians=False, showextrema=False)
    for b in parts_l["bodies"]:
        # Keep only the LEFT half
        verts = b.get_paths()[0].vertices
        verts[:, 0] = np.minimum(verts[:, 0], 0)
        b.set_facecolor(colour_left)
        b.set_edgecolor("black")
        b.set_alpha(0.85)
    for b in parts_r["bodies"]:
        verts = b.get_paths()[0].vertices
        verts[:, 0] = np.maximum(verts[:, 0], 0)
        b.set_facecolor(colour_right)
        b.set_edgecolor("black")
        b.set_alpha(0.85)

    # Overlay individual v2 samples as horizontal jittered ticks
    ax.scatter(np.random.default_rng(0).uniform(-0.05, -0.005, size=len(data_left)),
               data_left, s=6, c="black", alpha=0.5, zorder=3)

    # Mean / quartile markers
    for d, side, c in [(data_left, -1, colour_left), (data_right, +1, colour_right)]:
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax.plot([side * 0.07], [med], "o", color="white", markeredgecolor="black",
                markersize=6, zorder=4)
        ax.vlines(side * 0.07, q1, q3, color="black", lw=2.5)

    # Legend handles (custom)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=colour_left, edgecolor="black", label=label_left),
                       Patch(facecolor=colour_right, edgecolor="black", label=label_right)],
              fontsize=8, loc="upper right")
    ax.set_xticks([])
    ax.set_xlim(-0.6, 0.6)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--spec", type=Path,
                   default=HERE / "specs_v2" / "sample_sobol_synthaorta_v2.json",
                   help="v2 sample-mode spec to take v2 samples from.")
    p.add_argument("--out", type=Path,
                   default=HERE / "figures" / "synthaorta_compat_violin.png",
                   help="Output PNG path.")
    p.add_argument("--n-ref", type=int, default=10_000,
                   help="Reference sample size from the published distribution.")
    p.add_argument("--ref-seed", type=int, default=2026,
                   help="RNG seed for the reference sample (independent of Sobol).")
    args = p.parse_args()

    print(f"Loading v2 samples from {args.spec.name} …")
    v2_samples, spec = collect_v2_samples(args.spec)
    print(f"  → {len(next(iter(v2_samples.values())))} samples × {len(v2_samples)} params")

    # Resolve effective distribution per validated parameter (default_dist unless overridden)
    overrides = spec.get("distribution_overrides", {})

    def _effective_dist(name: str) -> dict:
        if name in overrides:
            return overrides[name]
        return PARAMETERS[name]["default_dist"]

    # Build subplots: 2x3 grid (4 SynthAorta + 1 extension + 1 summary panel)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.5))
    axes_flat = axes.ravel()

    COL_V2 = "#1f77b4"
    COL_REF = "#bcbd22"
    COL_EXT = "#7f7f7f"

    cell_idx = 0
    for name in SYNTHAORTA_PARAMS:
        ax = axes_flat[cell_idx]
        cell_idx += 1
        v2 = v2_samples[name]
        dist = _effective_dist(name)
        ref = reference_samples(dist, n=args.n_ref, seed=args.ref_seed + cell_idx)

        _split_violin(ax, v2, ref,
                      label_left=f"v2 Sobol (N={len(v2)})",
                      label_right=f"Published dist (N={args.n_ref})",
                      colour_left=COL_V2, colour_right=COL_REF)

        # Two-sample KS for distribution-shape agreement
        ks, pval = stats.ks_2samp(v2, ref)
        verdict = "OK" if pval > 0.05 else ("borderline" if pval > 0.01 else "MISMATCH")

        ax.set_title(
            f"{name}\n"
            f"v2: μ={v2.mean():.2f} σ={v2.std():.2f}    "
            f"ref: μ={ref.mean():.2f} σ={ref.std():.2f}\n"
            f"KS={ks:.3f}, p={pval:.3f} [{verdict}]",
            fontsize=9,
        )
        info = PARAMETERS[name]
        unit = "mm" if "radius" in info["description"].lower() else ""
        ax.set_ylabel(f"{name} [{unit}]" if unit else name)
        ax.grid(alpha=0.3, axis="y")

    # Extension parameter subplot
    for name in EXTENSION_PARAMS:
        if name not in v2_samples:
            continue
        ax = axes_flat[cell_idx]
        cell_idx += 1
        v2 = v2_samples[name]
        parts = ax.violinplot([v2], positions=[0], widths=0.85,
                              showmeans=False, showmedians=True, showextrema=True)
        for b in parts["bodies"]:
            b.set_facecolor(COL_EXT)
            b.set_alpha(0.85)
            b.set_edgecolor("black")
        ax.scatter(np.random.default_rng(0).uniform(-0.08, 0.08, size=len(v2)),
                   v2, s=8, c="black", alpha=0.5, zorder=3)
        ax.set_title(
            f"{name}\nEXTENSION (no SynthAorta equivalent)\n"
            f"v2: μ={v2.mean():.2f} σ={v2.std():.2f}",
            fontsize=9,
        )
        ax.set_ylabel(f"{name} [deg]")
        ax.set_xticks([])
        ax.set_xlim(-0.6, 0.6)
        ax.grid(alpha=0.3, axis="y")
        # Mark engineering range as a band
        info = PARAMETERS[name]
        ax.axhspan(info["min"], info["max"], alpha=0.08, color=COL_EXT,
                   label=f"workshop range [{info['min']}-{info['max']}]")
        ax.legend(loc="upper right", fontsize=8)

    # Summary text in the last cell
    ax_summary = axes_flat[-1]
    ax_summary.axis("off")
    ks_results = []
    for name in SYNTHAORTA_PARAMS:
        dist = _effective_dist(name)
        ref = reference_samples(dist, n=args.n_ref, seed=args.ref_seed + 100)
        ks, pval = stats.ks_2samp(v2_samples[name], ref)
        ks_results.append((name, ks, pval))

    n_v2 = len(next(iter(v2_samples.values())))
    n_pass = sum(1 for _, _, p in ks_results if p > 0.05)
    summary_lines = [
        "SynthAorta-comparability summary",
        "",
        f"v2 sampler  : Sobol scrambled, N={n_v2}",
        f"reference   : Monte-Carlo from Table I dist, N={args.n_ref}",
        "",
        f"KS p > 0.05 (same distribution):  {n_pass}/{len(ks_results)}",
        "",
        "Per parameter (validated):",
    ]
    for name, ks, pval in ks_results:
        verdict = "OK" if pval > 0.05 else ("borderline" if pval > 0.01 else "MISMATCH")
        summary_lines.append(f"  {name:14s} KS={ks:.3f}  p={pval:.3f}  [{verdict}]")
    summary_lines += [
        "",
        "Per parameter (extension, no reference):",
    ]
    for name in EXTENSION_PARAMS:
        if name in v2_samples:
            v2 = v2_samples[name]
            summary_lines.append(f"  {name:14s} μ={v2.mean():.2f} σ={v2.std():.2f}")

    ax_summary.text(0.05, 0.95, "\n".join(summary_lines),
                    family="monospace", fontsize=9, verticalalignment="top",
                    transform=ax_summary.transAxes)

    fig.suptitle(
        "SynthAorta-comparability check — cli_v2 sampler vs published Table I distributions",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"\nWrote {args.out}")

    # Also dump a small machine-readable summary alongside
    summary_path = args.out.with_suffix(".json")
    summary_payload = {
        "spec": str(args.spec.name),
        "v2_sample_size": n_v2,
        "reference_sample_size": args.n_ref,
        "ks_per_parameter": [
            {"name": name, "ks_stat": float(ks), "p_value": float(pval),
             "verdict": "OK" if pval > 0.05
                       else "borderline" if pval > 0.01 else "MISMATCH"}
            for name, ks, pval in ks_results
        ],
        "extension_parameters": [
            {"name": name, "mean": float(v2_samples[name].mean()),
             "std": float(v2_samples[name].std()),
             "min": float(v2_samples[name].min()),
             "max": float(v2_samples[name].max())}
            for name in EXTENSION_PARAMS if name in v2_samples
        ],
    }
    summary_path.write_text(json.dumps(summary_payload, indent=2) + "\n")
    print(f"Wrote {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
