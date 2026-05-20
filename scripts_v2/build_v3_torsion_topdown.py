#!/usr/bin/env python3
"""Top-down PyVista figure for the v3 torsion sweep.

Shows that 'torsion_deg' actually rotates the arch+descending around
the inlet z-axis — invisible in a side view (planar tube), but clearly
visible from above as the descending tube swings around the inlet.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pyvista as pv

pv.OFF_SCREEN = True
pv.global_theme.background = "white"
pv.global_theme.font.color = "black"

HERE = Path(__file__).resolve().parent.parent

# Pick 6 cases evenly across tors_001 .. tors_010
COHORT = HERE / "outputs" / "v3_torsion_sweep"
PICK_INDICES = [1, 3, 5, 6, 8, 10]  # roughly evenly across -20° to +20°


def main() -> int:
    cases = []
    for i in PICK_INDICES:
        d = COHORT / f"tors_{i:03d}"
        if (d / "wall_aorta.stl").exists():
            cases.append(d)
    if not cases:
        print(f"No cases found under {COHORT}")
        return 1

    n_cols = len(cases)
    p = pv.Plotter(shape=(2, n_cols), off_screen=True,
                   window_size=(280 * n_cols, 500),
                   border=True, border_color="lightgray")

    for r, (view, label) in enumerate([("iso", "Oblique"),
                                        ("top", "Top-down (y-axis horizontal)")]):
        for c, case_dir in enumerate(cases):
            p.subplot(r, c)
            mesh = pv.read(str(case_dir / "wall_aorta.stl"))
            p.add_mesh(mesh, color="#c44e52", smooth_shading=True,
                       specular=0.2, opacity=1.0)
            if view == "iso":
                p.view_isometric()
            else:
                p.view_xy()
            p.reset_camera()
            p.camera.zoom(0.75 if view == "iso" else 0.85)
            title = f"{case_dir.name}\n{label}" if c == 0 else case_dir.name
            p.add_text(title, font_size=8, position="upper_left", color="black")

    out_path = HERE / "figures" / "v3_torsion_oblique_vs_topdown.png"
    p.screenshot(str(out_path), return_img=False)
    p.close()
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
