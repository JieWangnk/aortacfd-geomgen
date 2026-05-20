#!/usr/bin/env python3
"""Side-by-side comparison of rigid torsion vs gradual twist at the same angle.

Two cases, both at 30°:
  - LEFT  : torsion_deg=30 (rigid rotation of arch+descending around z-axis;
            arch stays planar in a rotated plane)
  - RIGHT : twist_deg=30 (gradual twist along the arch; arch becomes a
            non-planar 3D curve)

Two views per case (top row = oblique, bottom row = top-down) so the
difference between rigid rotation and gradual twist is clearly visible.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pyvista as pv

pv.OFF_SCREEN = True
pv.global_theme.background = "white"
pv.global_theme.font.color = "black"

HERE = Path(__file__).resolve().parent.parent
RIGID = HERE / "outputs" / "v3_compare_rigid" / "baseline_v3"
GRADUAL = HERE / "outputs" / "v3_compare_gradual" / "baseline_v3"


def main() -> int:
    cases = [(RIGID, "torsion_deg = 30°  (RIGID rotation)"),
             (GRADUAL, "twist_deg = 30°  (GRADUAL twist)")]
    for d, _ in cases:
        if not (d / "wall_aorta.stl").exists():
            print(f"Missing: {d}")
            return 1

    p = pv.Plotter(shape=(2, 2), off_screen=True,
                   window_size=(720, 760),
                   border=True, border_color="lightgray")
    for r, (view, label) in enumerate([("iso", "Oblique"),
                                        ("top", "Top-down (y-axis horizontal)")]):
        for c, (case_dir, title) in enumerate(cases):
            p.subplot(r, c)
            mesh = pv.read(str(case_dir / "wall_aorta.stl"))
            p.add_mesh(mesh, color="#c44e52", smooth_shading=True, specular=0.2)
            if view == "iso":
                p.view_isometric()
            else:
                p.view_xy()
            p.reset_camera()
            p.camera.zoom(0.7 if view == "iso" else 0.8)
            heading = f"{title}\n{label}" if r == 0 else label
            p.add_text(heading, font_size=10, position="upper_left", color="black")

    out_path = HERE / "figures" / "v3_torsion_vs_twist.png"
    p.screenshot(str(out_path), return_img=False)
    p.close()
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
