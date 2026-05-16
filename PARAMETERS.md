# Aorta geometry generator — parameter reference

Generated from `cli.py --list-params --markdown`. Every parameter the orchestrator passes to the Blender script, with workshop-sensible ranges. Ranges are *suggestions*, not hard limits — the generator will still run outside them but the result may not be physiologically realistic.

## Main tube

| Parameter | Type | Default | Workshop range | Description |
|---|---|---|---|---|
| `diameter` | float | `24.0` | 18.0–40.0 | Main lumen diameter [mm] |
| `ascending_length` | float | `45.0` | 30.0–60.0 | Ascending aorta length [mm] |
| `arch_span` | float | `70.0` | 55.0–90.0 | Arch span ascending→descending [mm] |
| `arch_height` | float | `35.0` | 20.0–50.0 | Arch rise above ascending top [mm] |
| `descending_length` | float | `80.0` | 60.0–250.0 | Descending aorta length [mm] |

## Branches

| Parameter | Type | Default | Workshop range | Description |
|---|---|---|---|---|
| `branch_count` | int | `3` | 1–3 | Number of supra-aortic branches (1, 2, or 3) |
| `branch_diameter_ratio` | float | `0.42` | 0.3–0.55 | Branch diameter / main diameter |
| `branch_length` | float | `35.0` | 20.0–60.0 | Branch outlet extension length [mm] |
| `branch_spacing` | float | `14.0` | 8.0–20.0 | Spacing of branch origins along arch [mm] |
| `branch_tilt_deg` | float | `60.0` | 30.0–80.0 | Branch take-off tilt from +x toward +z [deg] |
| `branch_splay_deg` | float | `22.0` | 0.0–40.0 | Out-of-plane branch splay magnitude [deg] |

## Coarctation

| Parameter | Type | Default | Workshop range | Description |
|---|---|---|---|---|
| `coarctation_area_reduction` | float | `0.65` | 0.0–0.9 | Area reduction at the throat (0=none, 0.9=critical) |
| `coarctation_length` | float | `30.0` | 10.0–40.0 | Axial length of the smooth coarctation [mm] |
| `coarctation_centre_fraction` | float | `0.72` | 0.5–0.9 | Centre location along centreline arc [0-1] |
| `proximal_hypoplasia` | float | `0.0` | 0.0–0.25 | Smooth proximal arch diameter reduction [0-0.25] |

## Roughness

| Parameter | Type | Default | Workshop range | Description |
|---|---|---|---|---|
| `roughness` | bool | `False` | true / false | Enable mild radial surface roughness |
| `noise_amplitude` | float | `0.35` | 0.1–1.0 | Roughness amplitude [mm] |
| `noise_scale` | float | `0.08` | 0.05–0.5 | Roughness spatial scale |

## Geometry mesh

| Parameter | Type | Default | Workshop range | Description |
|---|---|---|---|---|
| `segments_radial` | int | `64` | 32–128 | Circumferential ring vertices |
| `curve_samples` | int | `220` | 100–400 | Main centreline sample count |
| `cylinder_vertices` | int | `64` | 32–128 | Cylinder vertices for branches |

