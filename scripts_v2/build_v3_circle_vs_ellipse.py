#!/usr/bin/env python3
"""Side-by-side comparison of arch_shape='circle' vs 'ellipse'.

Three pairs at increasingly extreme W/H ratios:
  - Pair 1: W=H=45 (square arch, OK in both modes)
  - Pair 2: W=90, H=45 (canonical U-arch, ratio 2.0)
  - Pair 3: W=30, H=80 (tall narrow — REJECTED in circle mode, OK in ellipse)

The circle column shows what the user would get if they had stuck with
the default circle mode (or, for pair 3, the error they would see).
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pyvista as pv

pv.OFF_SCREEN = True
pv.global_theme.background = "white"
pv.global_theme.font.color = "black"

HERE = Path(__file__).resolve().parent.parent
OUT = HERE / "outputs" / "v3_circle_vs_ellipse"


PAIRS = [
    ("square_W45_H45", 45.0, 45.0),
    ("U_arch_W90_H45", 90.0, 45.0),
    ("tall_W30_H80",   30.0, 80.0),
]


def run_one(shape: str, w: float, h: float, label: str) -> Path | None:
    case_dir = OUT / f"{shape}_{label}"
    if case_dir.exists():
        shutil.rmtree(case_dir.parent / case_dir.name, ignore_errors=True)
    cmd = [
        "python3", str(HERE / "cli_v3.py"),
        "--spec", str(HERE / "specs_v3" / "single_baseline_v3.json"),
        "--output", str(OUT / f"{shape}_{label}"),
        "--yes",
        "--param", f"arch_shape={shape}",
        "--param", f"arch_width_mm={w}",
        "--param", f"arch_height_mm={h}",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  [{shape}/{label}] FAIL: {res.stderr.strip().splitlines()[-1] if res.stderr.strip() else 'unknown'}")
        return None
    inner = OUT / f"{shape}_{label}" / "baseline_v3"
    return inner if (inner / "wall_aorta.stl").exists() else None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    results = {}  # (label, shape) → Path or None
    for label, w, h in PAIRS:
        for shape in ("circle", "ellipse"):
            results[(label, shape)] = run_one(shape, w, h, label)

    # 2 rows × 3 cols (rows: circle / ellipse; cols: pairs)
    rows, cols = 2, 3
    p = pv.Plotter(shape=(rows, cols), off_screen=True,
                   window_size=(360 * cols, 380 * rows),
                   border=True, border_color="lightgray")

    shapes = ["circle", "ellipse"]
    for r, shape in enumerate(shapes):
        for c, (label, w, h) in enumerate(PAIRS):
            p.subplot(r, c)
            case = results[(label, shape)]
            if case is None:
                p.add_text(f"{shape}\n{label}\n(REJECTED)",
                           font_size=11, position="upper_left", color="#aa0000")
                continue
            mesh = pv.read(str(case / "wall_aorta.stl"))
            p.add_mesh(mesh, color="#c44e52", smooth_shading=True, specular=0.2)
            p.view_xz()  # front view — best for seeing arch profile
            p.reset_camera()
            p.camera.zoom(0.75)
            title = f"{shape}: W={w:.0f}, H={h:.0f}"
            p.add_text(title, font_size=10, position="upper_left", color="black")

    out_path = HERE / "figures" / "v3_circle_vs_ellipse.png"
    p.screenshot(str(out_path), return_img=False)
    p.close()
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
