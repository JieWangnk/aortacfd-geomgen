# Aorta geometry generator v2 — parameter reference

Generated from `cli_v2.py --list-params --markdown`. Defaults and default sample-mode distributions are from the SynthAorta paper (Bošnjak et al. 2025) Table I unless noted.

## Radii

| Parameter | Type | Default | Workshop range | Default distribution | Citation | Description |
|---|---|---|---|---|---|---|
| `r_ascending` | float | `13.7` | 8.0–22.0 | normal(μ=13.7, σ=2.3) | Schäfer 2018; Wolak 2008 | Ascending aorta radius [mm] |
| `r_arch` | float | `13.0` | 8.0–20.0 | normal(μ=13.0, σ=2.0) | Marrocco-Trischitta; Saitta 2022 | Arch radius [mm] |
| `r_descending` | float | `12.2` | 8.0–20.0 | normal(μ=12.2, σ=2.3) | Bouti 2017; Schäfer 2018 | Descending aorta radius [mm] |
| `taper_mode` | str | `smoothstep` | piecewise / linear / smoothstep | (fixed) | — | Radius blending across segment boundaries |

## Lengths

| Parameter | Type | Default | Workshop range | Default distribution | Citation | Description |
|---|---|---|---|---|---|---|
| `ascending_length` | float | `50.0` | 40.0–90.0 | uniform(40.0-90.0) | Mills 1970; Bouti 2017 | Ascending aorta length [mm] |
| `descending_length` | float | `200.0` | 150.0–300.0 | uniform(150.0-300.0) | anatomy textbooks | Descending aorta length [mm] |

## Arch curvature

| Parameter | Type | Default | Workshop range | Default distribution | Citation | Description |
|---|---|---|---|---|---|---|
| `arch_R_c` | float | `40.4` | 25.0–60.0 | gumbel(loc=40.4, scale=2.4) | Choi 2017; Saitta 2022 (SynthAorta Table I) | Arch radius of curvature [mm] |
| `arch_angle_deg` | float | `180.0` | 120.0–200.0 | normal(μ=180.0, σ=15.0) | engineering default (Madhwal arch-type classification context) | Subtended angle of the arch arc [deg] |
| `arch_tilt_deg` | float | `0.0` | -30.0–30.0 | normal(μ=0.0, σ=8.0) | anatomy textbooks — typical leftward tilt 5-15° | RIGID rotation of arch+descending around inlet z-axis [deg] (arch stays planar) |
| `arch_twist_deg` | float | `0.0` | -45.0–45.0 | normal(μ=0.0, σ=10.0) | engineering default — physiologically plausible helical descending | GRADUAL twist around z-axis along the arch [deg] (arch becomes non-planar) |
| `junction_blend_mm` | float | `12.0` | 0.0–40.0 | (fixed) | — | Cubic-Bezier blend width at each arch junction [mm] (0 = sharp circular-arc corners) |

## Arch curvature (alt direct)

| Parameter | Type | Default | Workshop range | Default distribution | Citation | Description |
|---|---|---|---|---|---|---|
| `arch_span_mm` | float | `80.8` | 30.0–120.0 | (fixed) | — | Arch horizontal extent ascending→descending [mm]. Use WITH arch_height_mm to override arch_R_c+angle. |
| `arch_height_mm` | float | `40.4` | 20.0–60.0 | (fixed) | — | Arch PEAK height above ascending top [mm]. Use WITH arch_span_mm to override arch_R_c+angle. |

## Non-planar Fourier

| Parameter | Type | Default | Workshop range | Default distribution | Citation | Description |
|---|---|---|---|---|---|---|
| `delta_3` | float | `0.0` | 0.0–1.5 | normal(μ=1.0, σ=0.09) | Bošnjak et al. 2025 Table I (δ_3) | SynthAorta δ_3: cos(2w·||x||) multiplier [0=planar, 1=SynthAorta nominal] |
| `delta_4` | float | `0.0` | 0.0–1.5 | normal(μ=1.0, σ=0.09) | Bošnjak et al. 2025 Table I (δ_4) | SynthAorta δ_4: sin(2w·||x||) multiplier [0=planar, 1=SynthAorta nominal] |

## Geometry mesh

| Parameter | Type | Default | Workshop range | Default distribution | Citation | Description |
|---|---|---|---|---|---|---|
| `segments_radial` | int | `96` | 32–192 | (fixed) | — | Circumferential ring vertices |
| `curve_samples` | int | `300` | 100–600 | (fixed) | — | Total centreline sample count |

