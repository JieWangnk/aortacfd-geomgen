#!/usr/bin/env python3
"""PyVista 4x4 gallery for the 16-case v3 Sobol cohort with gradual twist.

Two output figures:
  - figures/v3_sobol_gallery_with_twist.png       (4x4 oblique view)
  - figures/v3_sobol_gallery_topdown.png          (4x4 top-down view —
                                                   shows the twist hooks)
"""

from __future__ import annotations

import sys
from pathlib import Path

import pyvista as pv

pv.OFF_SCREEN = True
pv.global_theme.background = "white"
pv.global_theme.font.color = "black"

HERE = Path(__file__).resolve().parent.parent
COHORT = HERE / "outputs" / "v3_sobol_gallery"


def render_grid(view: str, out_path: Path, zoom: float) -> None:
    cases = sorted(p for p in COHORT.iterdir()
                   if p.is_dir() and (p / "wall_aorta.stl").exists())
    if not cases:
        print(f"No cases in {COHORT}")
        return
    n = len(cases)
    print(f"  {n} cases → {view} view")

    rows, cols = 4, 4
    p = pv.Plotter(shape=(rows, cols), off_screen=True,
                   window_size=(280 * cols, 320 * rows),
                   border=True, border_color="lightgray")
    for k in range(rows * cols):
        r, c = divmod(k, cols)
        p.subplot(r, c)
        if k >= n:
            p.add_text("(empty)", font_size=8, color="gray")
            continue
        case_dir = cases[k]
        mesh = pv.read(str(case_dir / "wall_aorta.stl"))
        p.add_mesh(mesh, color="#c44e52", smooth_shading=True, specular=0.2)
        if view == "iso":
            p.view_isometric()
        else:
            p.view_xy()
        p.reset_camera()
        p.camera.zoom(zoom)
        p.add_text(case_dir.name, font_size=8, position="upper_left", color="black")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    p.screenshot(str(out_path), return_img=False)
    p.close()
    print(f"  → {out_path}")


def main() -> int:
    print("[v3 Sobol gallery] oblique view")
    render_grid("iso", HERE / "figures" / "v3_sobol_gallery_with_twist.png", zoom=0.75)

    print("\n[v3 Sobol gallery] top-down view")
    render_grid("top", HERE / "figures" / "v3_sobol_gallery_topdown.png", zoom=0.85)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
