# Aorta geometry generator — v3

Parametric pipe U-bend (aortic-arch-like) STL generator. Dial a handful
of knobs, get a watertight tube split into `inlet.stl`, `outlet1.stl`,
`wall_aorta.stl`.

![v3 baseline](figures/v3_baseline_hero.png)

## Install

```bash
pip install -r requirements.txt        # numpy, scipy, numpy-stl + pyvista/matplotlib
```

Plus **Blender 3.x+** on your PATH (`blender --version` should work). See
the [main README](./README.md) for per-OS Blender install commands if needed.

## Quick start

```bash
# 1. Baseline pipe
python3 cli_v3.py --spec specs_v3/single_baseline_v3.json --output /tmp/v3 --yes

# 2. Discover all knobs
python3 cli_v3.py --list-params

# 3. Override on the fly
python3 cli_v3.py --spec specs_v3/single_baseline_v3.json --output /tmp/v3_custom --yes \
    --param r_inlet=16 --param arch_R_c_mm=40 --param twist_deg=15
```

## Output per case

```
<output>/<case_id>/
  inlet.stl              # ascending cap (radius = r_inlet)
  outlet1.stl            # descending cap (radius = r_outlet)
  wall_aorta.stl         # vessel wall
  geometry.meta.json     # provenance (v3 inputs + v2 translation)
```

## The knobs

| Knob | Default | What it is |
|---|---|---|
| `r_inlet` | 14 mm | Tube radius at the inlet (ascending) |
| `r_outlet` | 10 mm | Tube radius at the outlet (descending) |
| `arch_radius_mm` | 0 (auto) | Tube radius at the arch segment. 0 = midpoint of inlet/outlet. |
| `taper_mode` | `smoothstep` | Lumen transition: `smoothstep` / `linear` / `piecewise` |
| `arch_width_mm` | 90 mm | Horizontal arch extent |
| `arch_height_mm` | 45 mm | Arch peak height above ascending top |
| `arch_R_c_mm` | 0 (off) | Centerline curvature shortcut. > 0 auto-sets W=2R, H=R. |
| `arch_shape` | `circle` | `circle` (constraint H ≤ W ≤ 2H) or `ellipse` (any positive W, H) |
| `torsion_deg` | 0° | Rigid tilt of arch+descending around inlet z-axis |
| `twist_deg` | 0° | Gradual twist along arch (non-planar 3D curve) |
| `ascending_length` | 50 mm | Straight ascending length before arch |
| `descending_length` | 200 mm | Straight descending length after arch |

> **Two "radius" knobs — easy to confuse:**
> - `arch_radius_mm` = **tube** cross-section radius at the arch (like `r_inlet`, `r_outlet`)
> - `arch_R_c_mm` = **centerline curvature** radius (how sharply the path bends)

## Example specs

| Spec | Mode | Demonstrates |
|---|---|---|
| `single_baseline_v3.json` | single | Canonical U-arch (the default in the image above) |
| `single_R_c_v3.json` | single | Same geometry specified by curvature radius `arch_R_c_mm=45` |
| `single_ellipse_v3.json` | single | Tall narrow arch (W=30, H=80) — only possible in ellipse mode |
| `single_linear_taper_v3.json` | single | Linear taper between segment radii instead of smoothstep |
| `sweep_torsion_v3.json` | sweep | 10-step torsion sweep, -20° → +20° (rigid plane tilt) |
| `sweep_twist_v3.json` | sweep | 10-step gradual-twist sweep, -30° → +30° (non-planar curl) |

Run any of them:
```bash
python3 cli_v3.py --spec specs_v3/<spec>.json --output outputs/<name> --yes
```

## Sobol sampling (all 8 v3 knobs at once)

v3 itself only does single + sweep modes. For Sobol over the full v3
parameter space, use the pre-built v2 spec:

```bash
python3 cli_v2.py --spec specs_v2/sample_sobol_v3_8d_gradual_twist.json \
    --output outputs/v3_sobol_gallery --yes
# 16 cases, ~50 s, all with twist_deg ∈ [10°, 30°] for visible non-planarity
```

![16-case Sobol gallery](figures/v3_sobol_gallery_with_twist.png)

For 256+ cases (PCE sensitivity), edit `n_cases` in the spec — see the
[main README](./README.md) for sample-size guidance.

## Render a gallery from any cohort

```bash
/home/mchi4jw4/GitHub/.venv/bin/python scripts_v2/build_v2_gallery_pyvista.py \
    --planar-cohort outputs/<your_cohort> \
    --hero-case outputs/<your_cohort>/<one_case>
# → figures/v2_cohort_diversity_gallery.png + v2_single_hero.png
```

## When v3 isn't enough

Drop down to [`cli_v2.py`](./README_v2.md) for:
- Direct `arch_R_c` + `arch_angle_deg` control (any angle, no W/H constraints)
- SynthAorta non-planar Fourier multipliers (`δ_3`, `δ_4`)
- Junction blend width tuning, mesh resolution knobs
- Sobol / LHS / random sample modes

v3 keeps the interface minimal; v2 exposes everything.

## Files

| File | Purpose |
|---|---|
| `cli_v3.py` | 12-knob orchestrator (wraps `blender_aorta_v2.py` internally) |
| `specs_v3/` | Example specs (one per scenario above) |
| `tests/test_v3.py` | 41 tests, no Blender required |
