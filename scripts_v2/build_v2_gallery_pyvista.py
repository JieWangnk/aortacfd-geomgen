#!/usr/bin/env python3
"""PyVista 3D figures for the v2 healthy-aorta generator.

Replaces the matplotlib statistical violin/pairplot figures with
actual 3D renders of generated STL geometries — what the v2 cohort
*looks like*, not just what its parameter marginals look like.

Two figures are produced:

  1. ``figures/v2_cohort_diversity_gallery.png`` — N×N grid of cases
     sampled at random from a cohort directory. Same viewpoint for
     every panel so you can compare radii, lengths, arch shape, and
     out-of-plane wobble across cases.

  2. ``figures/v2_planar_vs_nonplanar.png`` — top-down + oblique view
     of the same parameter set rendered first with δ_3=δ_4=0 (planar)
     and then with δ_3=δ_4=1 (SynthAorta nominal). Shows that the
     Fourier displacement actually does something.

Run with the project virtualenv that has pyvista::

    /home/mchi4jw4/GitHub/.venv/bin/python \\
        scripts_v2/build_v2_gallery_pyvista.py \\
        --planar-cohort outputs/v2_sobol_100 \\
        --nonplanar-cohort outputs/v2_sobol_nonplanar_demo

If a cohort path doesn't exist, that figure is skipped with a notice.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

# Headless rendering — works on machines without an active display
pv.OFF_SCREEN = True
pv.global_theme.background = "white"
pv.global_theme.font.color = "black"


HERE = Path(__file__).resolve().parent.parent  # repo root


def load_case_wall(case_dir: Path) -> pv.PolyData | None:
    """Read ``wall_aorta.stl`` from a v2 case directory."""
    stl_path = case_dir / "wall_aorta.stl"
    if not stl_path.exists():
        return None
    return pv.read(str(stl_path))


def list_cases(cohort_dir: Path) -> list[Path]:
    if not cohort_dir.exists():
        return []
    return sorted(p for p in cohort_dir.iterdir()
                  if p.is_dir() and (p / "wall_aorta.stl").exists())


def setup_subplot(p: pv.Plotter, mesh: pv.PolyData, *, view: str = "iso",
                  title: str | None = None, zoom: float = 0.85) -> None:
    """Add the mesh to the current subplot with a consistent style + view.

    Uses PyVista's built-in auto-fit views (which size the camera based on
    the mesh's actual bounding box) followed by ``reset_camera`` and a
    multiplicative ``zoom`` so the geometry sits comfortably inside the
    panel with margin. zoom < 1.0 → camera moves further back.
    """
    p.add_mesh(mesh, color="#c44e52", smooth_shading=True,
               specular=0.2, opacity=1.0, show_edges=False)
    if title:
        p.add_text(title, font_size=8, position="upper_left",
                   color="black", shadow=False)

    # Auto-fit camera: PyVista picks an appropriate distance from the
    # mesh's bounding box. Then nudge the zoom for margin.
    if view == "iso":
        p.view_isometric()
    elif view == "front":
        # Look along -y axis (mesh's xz extent fills the panel)
        p.view_xz()
    elif view == "side":
        # Look along -x axis (mesh's yz extent fills the panel)
        p.view_yz()
    elif view == "top":
        # Look down -z axis (sees the y-wobble most clearly)
        p.view_xy()
    else:
        p.view_isometric()

    p.reset_camera()
    p.camera.zoom(zoom)


# ── Figure 1: cohort diversity gallery ────────────────────────────────────


def build_diversity_gallery(cohort_dir: Path, out_path: Path,
                            grid: tuple[int, int] = (5, 5),
                            seed: int = 0) -> bool:
    """Render an n_rows × n_cols gallery from a v2 cohort directory."""
    cases = list_cases(cohort_dir)
    if not cases:
        print(f"  [skip] no cases found under {cohort_dir}")
        return False

    n_rows, n_cols = grid
    n_panels = n_rows * n_cols
    rng = random.Random(seed)
    pick = rng.sample(cases, k=min(n_panels, len(cases)))

    print(f"  cohort       : {cohort_dir}  ({len(cases)} cases)")
    print(f"  picking      : {len(pick)} cases for {n_rows}x{n_cols} grid")

    p = pv.Plotter(shape=(n_rows, n_cols), off_screen=True,
                   window_size=(320 * n_cols, 360 * n_rows),
                   border=True, border_color="lightgray")

    for k in range(n_panels):
        r, c = divmod(k, n_cols)
        p.subplot(r, c)
        if k >= len(pick):
            p.add_text("(empty)", font_size=8, color="gray")
            continue
        case_dir = pick[k]
        mesh = load_case_wall(case_dir)
        if mesh is None:
            p.add_text("(no STL)", font_size=8, color="gray")
            continue
        setup_subplot(p, mesh, view="iso", title=case_dir.name, zoom=0.75)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    p.screenshot(str(out_path), return_img=False)
    p.close()
    print(f"  → wrote {out_path}")
    return True


# ── Figure 2: planar vs non-planar comparison ─────────────────────────────


def build_planar_vs_nonplanar(planar_cohort: Path, nonplanar_cohort: Path,
                              out_path: Path, n_pairs: int = 4) -> bool:
    """Side-by-side: same case-index from planar + nonplanar cohorts.

    Top row = planar (δ=0). Bottom row = nonplanar (δ ≈ 1). Each pair
    uses the same Sobol index, so the underlying parameter realisations
    are matched apart from the δ_3, δ_4 dimensions.

    Two views per geometry: oblique (top-half) + top-down (bottom-half)
    so the y-wobble is visible.
    """
    planar_cases = list_cases(planar_cohort)
    np_cases = list_cases(nonplanar_cohort)
    if not planar_cases:
        print(f"  [skip] planar cohort {planar_cohort} is empty")
        return False
    if not np_cases:
        print(f"  [skip] nonplanar cohort {nonplanar_cohort} is empty")
        return False

    n_pairs = min(n_pairs, len(planar_cases), len(np_cases))
    print(f"  planar       : {planar_cohort}  ({len(planar_cases)} cases)")
    print(f"  nonplanar    : {nonplanar_cohort}  ({len(np_cases)} cases)")
    print(f"  comparing    : first {n_pairs} cases of each, two views each")

    n_rows = 4  # 2 views (oblique, top-down) × 2 cohorts
    n_cols = n_pairs

    p = pv.Plotter(shape=(n_rows, n_cols), off_screen=True,
                   window_size=(340 * n_cols, 360 * n_rows),
                   border=True, border_color="lightgray")

    row_labels = [
        ("Planar — oblique view", "iso", planar_cases, 0.75),
        ("Planar — top-down view (y-axis horizontal)", "top", planar_cases, 0.85),
        ("Non-planar (δ_3, δ_4 ~ 1) — oblique view", "iso", np_cases, 0.75),
        ("Non-planar (δ_3, δ_4 ~ 1) — top-down view", "top", np_cases, 0.85),
    ]
    for r, (label, view, cohort, zoom) in enumerate(row_labels):
        for c in range(n_cols):
            p.subplot(r, c)
            mesh = load_case_wall(cohort[c])
            if mesh is None:
                continue
            title = f"{cohort[c].name}\n{label}" if c == 0 else cohort[c].name
            setup_subplot(p, mesh, view=view, title=title, zoom=zoom)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    p.screenshot(str(out_path), return_img=False)
    p.close()
    print(f"  → wrote {out_path}")
    return True


# ── Figure 3: single hero render (multi-view) ─────────────────────────────


def build_single_hero(case_dir: Path, out_path: Path) -> bool:
    mesh = load_case_wall(case_dir)
    if mesh is None:
        print(f"  [skip] no wall_aorta.stl in {case_dir}")
        return False

    print(f"  single case  : {case_dir}")

    p = pv.Plotter(shape=(1, 3), off_screen=True,
                   window_size=(1500, 600),
                   border=True, border_color="lightgray")
    for c, (view, title, zoom) in enumerate(
        [("iso", "Oblique", 0.75),
         ("front", "Front (y-axis into page)", 0.80),
         ("top", "Top-down (y-axis horizontal)", 0.85)]
    ):
        p.subplot(0, c)
        setup_subplot(p, mesh, view=view, title=title, zoom=zoom)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    p.screenshot(str(out_path), return_img=False)
    p.close()
    print(f"  → wrote {out_path}")
    return True


# ── Main ──────────────────────────────────────────────────────────────────


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--planar-cohort", type=Path,
                        default=HERE / "outputs" / "v2_sobol_100",
                        help="Directory of generated planar v2 cases.")
    parser.add_argument("--nonplanar-cohort", type=Path,
                        default=HERE / "outputs" / "v2_sobol_nonplanar_demo",
                        help="Directory of generated non-planar v2 cases.")
    parser.add_argument("--hero-case", type=Path, default=None,
                        help="Single case dir for the hero render. If omitted, "
                             "uses the first case of --planar-cohort.")
    parser.add_argument("--figures-dir", type=Path,
                        default=HERE / "figures",
                        help="Where to write the PNGs.")
    parser.add_argument("--gallery-grid", type=str, default="auto",
                        help="Diversity gallery grid as ROWSxCOLS (e.g. 4x2). "
                             "'auto' picks the smallest grid that holds the available cases "
                             "(prefers near-square or 2:1 layouts).")
    parser.add_argument("--n-compare-pairs", type=int, default=4,
                        help="Number of planar-vs-nonplanar pairs to show.")
    parser.add_argument("--seed", type=int, default=0,
                        help="RNG seed for case picking in the gallery.")
    args = parser.parse_args()

    if args.gallery_grid.lower() == "auto":
        # Pick smallest near-square grid that holds the available planar cases
        n_cases = len(list_cases(args.planar_cohort))
        if n_cases == 0:
            rows, cols = 1, 1
        else:
            cols = max(1, int(math.ceil(math.sqrt(n_cases))))
            rows = max(1, int(math.ceil(n_cases / cols)))
        print(f"[gallery] auto-grid: {rows}x{cols} for {n_cases} cases")
    else:
        rows, cols = (int(x) for x in args.gallery_grid.lower().split("x"))

    print("=" * 72)
    print("PyVista figures for v2 healthy-aorta cohort")
    print("=" * 72)

    print("\n[1/3] diversity gallery")
    build_diversity_gallery(
        cohort_dir=args.planar_cohort,
        out_path=args.figures_dir / "v2_cohort_diversity_gallery.png",
        grid=(rows, cols),
        seed=args.seed,
    )

    print("\n[2/3] planar vs nonplanar comparison")
    build_planar_vs_nonplanar(
        planar_cohort=args.planar_cohort,
        nonplanar_cohort=args.nonplanar_cohort,
        out_path=args.figures_dir / "v2_planar_vs_nonplanar.png",
        n_pairs=args.n_compare_pairs,
    )

    print("\n[3/3] single hero render")
    hero_case = args.hero_case
    if hero_case is None:
        cases = list_cases(args.planar_cohort)
        hero_case = cases[0] if cases else None
    if hero_case is None or not hero_case.exists():
        print(f"  [skip] no hero case available (planar cohort empty?)")
    else:
        build_single_hero(
            case_dir=hero_case,
            out_path=args.figures_dir / "v2_single_hero.png",
        )

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
