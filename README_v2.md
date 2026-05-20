# Aorta geometry generator — v2 (healthy aorta, no branches, no coarctation)

Parallel "healthy aorta" parameter-sweep tool that lives alongside the
production v1 (`cli.py`). Designed for SynthAorta-paper-comparable
sensitivity studies and ML training-data generation on a simpler topology.

**What's the same as v1:**

- JSON spec schema (single / sweep / sample / grid modes)
- `--param key=value` overrides, `--dry-run`, `--limit`, `--list-params`,
  validator "did you mean" hints, sweep-manifest CSV, cost warning, etc.
- Sobol / LHS / random samplers (via `sampler.py`)
- STL patch splitter (via `split_patches.py`)

**What's different from v1:**

| | v1 (`cli.py`) | v2 (`cli_v2.py`) |
|---|---|---|
| Topology | Aortic arch + 1-3 supra-aortic branches | Aortic arch only, no branches |
| Pathology | Coarctation + proximal hypoplasia | Healthy only |
| Main radius | Single `diameter` (constant along the tube) | Three radii (ascending/arch/descending) with smoothstep taper |
| Arch curvature | Euclidean `arch_height` + `arch_span` | Clinical `arch_R_c` + `arch_angle_deg` |
| Sample mode | Uniform between `low`/`high` per parameter | Normal / Gumbel / Uniform via `distribution_overrides`; defaults from SynthAorta Table I |
| Output patches | `inlet.stl`, `outlet1..N.stl`, `wall_aorta.stl` | `inlet.stl`, `outlet1.stl`, `wall_aorta.stl` |
| Parameter count | 21 | 12 (10 planar + 2 non-planar Fourier) |
| Generator | `blender_aorta_like_generator.py` | `blender_aorta_v2.py` |

---

## Design overview

### Geometry construction (`blender_aorta_v2.py`)

```
   inlet                   peak (z = asc_length + R_c)
     |                          ◜◝
     |                       ◜      ◝
ascending  ─────────────►  arch          ──────────► descending
  (z-axis,                 (circular arc in xz                  outlet
   y=0)                     plane, radius R_c,
                            subtended angle θ)
```

**The pipeline, top to bottom:**

1. **Build the centreline** as three concatenated segments
   (`build_centreline`):
   - **Ascending**: straight line from (0, 0, 0) up the z-axis to
     (0, 0, `ascending_length`).
   - **Arch**: a circular arc of radius `arch_R_c` subtending angle
     `arch_angle_deg`, lying in the xz-plane. Centre at
     (`arch_R_c`, 0, `ascending_length`). The arc starts heading +z
     (tangent to the ascending segment) and rotates over to head into
     the descending direction.
   - **Descending**: straight line from the arch endpoint, in the
     tangent direction at the arc's end (flipped if needed so it
     heads toward -z). Length `descending_length`.
   - Total ~370 mm of centreline at the SynthAorta-mean defaults.
   - Samples distributed proportionally across the three segments
     based on `curve_samples`.

2. **Apply non-planar Fourier displacement** (optional, gated on
   `delta_3 != 0 or delta_4 != 0`). Adds an out-of-plane y-component
   to each centreline point per SynthAorta Eq 13 (Bošnjak et al. 2025):

   ```
   d_y(x) = A₁·cos(w·||x||) + B₁·sin(w·||x||)
          + A₂·δ₃·cos(2w·||x||) + B₂·δ₄·sin(2w·||x||)
   ```

   with constants `A₁=-0.798, B₁=-0.453, A₂=1.517, B₂=2.699, w=0.027`
   from SynthAorta Table II (fitted to their published patient
   centreline). `γ₁=γ₂=1` is hardcoded (the paper's sensitivity
   analysis identifies them as non-influential).

3. **Compute per-point radii** (`radius_along_arc`) using `taper_mode`:
   - `piecewise`: r = `r_ascending` / `r_arch` / `r_descending` per segment
   - `linear`: piecewise-linear blend across a 15 mm window at each
     segment boundary
   - `smoothstep` (default): 3t²-2t³ cosine-Hermite blend across the
     same window — C¹-continuous, monotonic between segment values

4. **Build the tube** (`build_tube_mesh`):
   - Rotation-minimising frames (Wang 2008 double-reflection) along
     the centreline so the cross-section doesn't twist artificially
   - Ring of `segments_radial` vertices at each centreline point
   - Quad faces between adjacent rings
   - Triangulated cap fans at the inlet (s=0) and outlet (s=1)

5. **Export STL** + sidecar JSON. The JSON carries
   `inlet_xyz_mm` / `outlet_xyz_mm` / `outlet_normal_xyz` etc., which
   `split_patches.py` uses to flood-fill the cap regions into
   `inlet.stl`, `outlet1.stl`, `wall_aorta.stl`.

### Orchestration (`cli_v2.py`)

- Reads a JSON spec, validates it against the `PARAMETERS` schema
  (with "did you mean" hints on typos), expands it into per-case
  parameter dicts (`expand_cases`), and runs them through Blender
  one at a time.
- **Sample mode** does *not* use uniform [low, high] mapping like v1.
  Instead, each parameter has a `default_dist` (Normal / Gumbel /
  Uniform with optional truncation bounds), and the Sobol unit cube
  is mapped through that distribution's inverse CDF
  (`scipy.stats.truncnorm.ppf` / `gumbel_r.ppf`). This is what makes
  the v2 sampling **statistically equivalent to SynthAorta** instead
  of just uniformly bounded.
- Reuses `sampler.py` and `split_patches.py` from v1; reuses
  validator helpers (`_suggest`, `warn_if_expensive`,
  `write_sweep_manifest`) by `from cli import …`.
- Writes a `sweep_manifest.csv` (one row per case with status, mode,
  sampler, sample_index, seed, full parameter vector) and a
  per-case `geometry.meta.json` with blake2b STL checksums for
  provenance.

### Why the parameter choices

The full design rationale, including the comparison with the
SynthAorta paper that motivated the v2 redesign, is documented
inline in `cli_v2.py`'s `PARAMETERS` dict. Short version:

- **Per-segment radii** (`r_ascending`, `r_arch`, `r_descending`)
  replace v1's single `diameter` so real-aorta tapering is captured;
  defaults from Schäfer 2018, Wolak 2008 (JACC), Bouti 2017.
- **Clinical `arch_R_c`** replaces v1's `arch_height` + `arch_span`
  because cardiologists measure radius of curvature, not Euclidean
  height. Default Gumbel(40.4, 2.4) from Choi 2017 / Saitta 2022.
- **`arch_angle_deg`** added so we can produce shallow-V and
  over-arched morphologies — not in SynthAorta (they use one base
  patient) but motivated by Madhwal et al's three arch-type clinical
  classification.
- **`delta_3` / `delta_4`** (non-planar Fourier) added to break the
  strict in-plane behaviour that the first v2 release had. Defaults
  to 0 (planar, backwards-compat); sample around 1.0 ± 0.09 for
  SynthAorta-equivalent y-wobble.

---

## Install

### Python (always needed)

```bash
pip install numpy numpy-stl scipy
# Optional, only for the validation figure scripts:
pip install matplotlib
```

Tested with Python 3.10 / 3.12.

### Blender (needed to actually generate STLs)

`cli_v2.py` launches Blender headless to build each geometry. Without
Blender installed, the CLI still works for `--list-params`,
`--dry-run`, validator checks, and the test suite — but no STL is
produced.

| OS | Install |
|---|---|
| **Linux** (Debian/Ubuntu) | `sudo apt install blender` (3.x in 22.04+, 4.x via [snap](https://snapcraft.io/blender) or [download](https://www.blender.org/download/)) |
| **Linux** (Fedora/Arch) | `sudo dnf install blender` / `sudo pacman -S blender` |
| **macOS** | `brew install --cask blender` or [download .dmg](https://www.blender.org/download/) |
| **Windows** | [Installer from blender.org](https://www.blender.org/download/) — add the install dir to PATH |

Verify with `blender --version` — should show 3.0 or newer.

If Blender lives somewhere other than the system PATH:

```bash
python3 cli_v2.py --blender /opt/blender-5.1/blender --spec ... --output ...
# or
export BLENDER=/opt/blender-5.1/blender
python3 cli_v2.py --spec ... --output ...
```

### What's actually used inside Blender

We use Blender purely as a mesh kernel — `bpy.data.meshes.new`,
`mesh.from_pydata`, the triangulate modifier, and the STL export
operator. No materials, lighting, animation, or Cycles rendering.
That means **Blender installs without GPU drivers** are fine.

## Quick start

```bash
# Discover parameters
python cli_v2.py --list-params

# Smoke test: one geometry at SynthAorta means
python cli_v2.py --spec specs_v2/single_baseline_v2.json --output /tmp/v2_single

# Sweep ascending-aorta radius from 10 → 18 mm
python cli_v2.py --spec specs_v2/sweep_r_ascending_v2.json --output /tmp/v2_sweep_ra

# Sweep arch radius of curvature from 25 → 60 mm
python cli_v2.py --spec specs_v2/sweep_R_c_v2.json --output /tmp/v2_sweep_Rc

# Sweep arch subtended angle from 120° → 200°
python cli_v2.py --spec specs_v2/sweep_arch_angle_v2.json --output /tmp/v2_sweep_ang

# 256-case Sobol over a 7-D SynthAorta-matched cube (PLANAR centrelines) — ~13 min
python cli_v2.py --spec specs_v2/sample_sobol_synthaorta_v2.json \
    --output outputs/v2_sobol_planar_256 --yes

# 256-case Sobol over a 9-D cube WITH SynthAorta δ_3, δ_4 non-planar Fourier — ~13 min
# centrelines acquire ~few-mm y-axis wobble, breaking strict planarity.
python cli_v2.py --spec specs_v2/sample_sobol_synthaorta_nonplanar_v2.json \
    --output outputs/v2_sobol_nonplanar_256 --yes
```

For a quick smoke run before committing to 13 minutes, add `--limit 8`
to either Sobol command — produces the first 8 cases (~25 s).

### Non-planar centrelines (δ_3, δ_4)

By default the v2 centreline is strictly planar (everything in y=0). To
introduce the SynthAorta-style out-of-plane wobble, set `delta_3` and
`delta_4` non-zero — the multipliers gate the Eq 13 Fourier displacement.

```bash
# Single non-planar geometry at SynthAorta nominal (δ_3 = δ_4 = 1.0)
python cli_v2.py --spec specs_v2/single_baseline_v2.json --output /tmp/v2_nonplanar \
    --param delta_3=1.0 --param delta_4=1.0
```

When sampled (default distribution Normal(1, 0.09)) you get the
statistically-validated SynthAorta variability. When fixed at 0 (default)
you get the legacy planar geometry — useful for planar-vs-non-planar
ablation studies.

---

## Sample size — how many cases do I need?

The shipped Sobol specs default to **N = 256** (Sobol-native power of 2
just above the 30·dim rule of thumb at dim=7). For the 9-D non-planar
cube N=256 is borderline for stable Sobol indices — bump to 512 for
that use case.

### What the 9-D nonplanar cube actually samples

The `sample_sobol_synthaorta_nonplanar_v2.json` spec varies **all 9
parameters simultaneously** — it's not a separate "non-planar only"
study. Every Sobol point varies:

1. `r_ascending` (Normal 13.7 ± 2.3 mm)
2. `r_arch` (Normal 13.0 ± 2.0 mm)
3. `r_descending` (Normal 12.2 ± 2.3 mm)
4. `arch_R_c` (Gumbel 40.4, 2.4 mm)
5. `arch_angle_deg` (Normal 180 ± 15°)
6. `ascending_length` (Uniform 40-90 mm)
7. `descending_length` (Uniform 150-300 mm)
8. `delta_3` (Normal 1.0 ± 0.09)
9. `delta_4` (Normal 1.0 ± 0.09)

`taper_mode`, `segments_radial`, `curve_samples` are fixed in `fixed:`
so the sampled cube stays purely continuous.

### N=256 verdict by downstream use

| Downstream use | N=256 verdict |
|---|---|
| Visual diversity / "what does the design space look like" | ✅ Plenty |
| Marginal-distribution validation (KS test vs published) | ✅ Plenty |
| Joint-distribution independence check (Pearson \|r\| < 2/√N ≈ 0.125) | ✅ Plenty |
| ML surrogate training (small MLP, ≲10⁴ weights) | ✅ Comfortable |
| ML surrogate training (GNN, neural operator) | ⚠️ Tight — 500-2000 typical |
| Sparse-PCE Sobol-index sensitivity at 7-D (planar) | ✅ Above rule-of-thumb (30·dim = 210) |
| Sparse-PCE Sobol-index sensitivity at 9-D (non-planar) | ⚠️ Borderline (30·dim = 270); bump to 512 |
| Full-quadratic PCE at 9-D | ⚠️ Bump to 512 or 1024 |

Compute scales linearly at **~3 s/case** on a modest workstation:

| N | Wall-clock |
|---|---|
| 100 | ~5 min |
| 256 (default) | ~13 min |
| 512 | ~26 min |
| 1024 | ~52 min |

To bump, edit the spec and re-run:

```bash
# In specs_v2/sample_sobol_synthaorta_nonplanar_v2.json:
#   "n_cases": 512
python cli_v2.py --spec specs_v2/sample_sobol_synthaorta_nonplanar_v2.json \
    --output outputs/v2_sobol_512_nonplanar --yes
```

## Generating figures

PyVista-based 3D renders of generated cohorts are produced by
`scripts_v2/build_v2_gallery_pyvista.py`. PyVista is **not** in the
system Python; use the project venv:

```bash
# Three figures land in figures/v2_*.png:
#   1. v2_cohort_diversity_gallery.png   — N×N grid of random picks
#   2. v2_planar_vs_nonplanar.png        — same-row planar vs δ_3,δ_4≈1
#   3. v2_single_hero.png                — three views of one hero case

/home/mchi4jw4/GitHub/.venv/bin/python scripts_v2/build_v2_gallery_pyvista.py \
    --planar-cohort   outputs/v2_sobol_planar_256 \
    --nonplanar-cohort outputs/v2_sobol_nonplanar_256
```

The script auto-skips any cohort that doesn't exist on disk and picks
a near-square gallery grid based on how many cases are available
(`--gallery-grid auto`, or set explicitly like `5x5`). It uses
PyVista's auto-fit `view_isometric` / `view_xy` / `view_xz` with a
margin zoom so the full geometry is visible regardless of arch height
or descending length.

Run end-to-end (generate cohorts → render figures):

```bash
# 1. Generate the 256-case planar cohort  (~13 min)
python cli_v2.py --spec specs_v2/sample_sobol_synthaorta_v2.json \
    --output outputs/v2_sobol_planar_256 --yes

# 2. Generate the 256-case nonplanar cohort  (~13 min)
python cli_v2.py --spec specs_v2/sample_sobol_synthaorta_nonplanar_v2.json \
    --output outputs/v2_sobol_nonplanar_256 --yes

# 3. Build all three PyVista figures
/home/mchi4jw4/GitHub/.venv/bin/python scripts_v2/build_v2_gallery_pyvista.py \
    --planar-cohort   outputs/v2_sobol_planar_256 \
    --nonplanar-cohort outputs/v2_sobol_nonplanar_256
```

## Output per case

```
<output>/<case_id>/
  inlet.stl              # ascending aorta inlet patch
  outlet1.stl            # descending aorta outlet patch
  wall_aorta.stl         # everything else (the vessel wall)
  geometry.meta.json     # parameters, derived geometry, patch checksums
```

The case-folder layout matches v1 so `AortaCFD-app/scripts/package_cases.py`
can stamp configs onto v2 output identically.

## When to use v2 vs v1

| If you want… | Use |
|---|---|
| ML training data for healthy-aorta surrogates | v2 |
| SynthAorta-comparable statistics | v2 |
| Coarctation severity sweep, branch position sensitivity, hypoplasia studies | v1 |
| Workshop demo of pressure-drop vs stenosis | v1 |
| Cleaner / simpler parameter space for sensitivity analysis | v2 |
| Cohort with realistic 3-vessel arch topology | v1 |

v1 specs (`geometry: "arch_branched_coarctation"`) are rejected by `cli_v2.py`
and vice versa — the two pipelines do not share spec files.

## Spec schema (v2)

```json
{
  "schema_version": "2.0",
  "name": "your_experiment_name",
  "mode": "single|sweep|sample|grid",
  "geometry": "healthy_arch_v2",

  "params": { ... },     // single mode: explicit values
                         // sample mode: {paramname: {}} (use default dist)
                         //              {paramname: {"low": ..., "high": ...}} (uniform)

  "distribution_overrides": {   // sample mode only — override per-parameter
    "r_ascending": {"type": "normal", "mean": 13.7, "std": 2.3,
                    "low": 8, "high": 22},
    "arch_R_c":    {"type": "gumbel", "loc": 40.4, "scale": 2.4,
                    "low": 25, "high": 60}
  },

  "sweep": { "param": "...", "low": 0, "high": 1, "n": 10 },
  "grid":  { "params": {"name1": [v1, v2], "name2": [v3, v4]} },

  "sampler": "sobol|lhs|random",
  "n_cases": 256,
  "seed": 42,

  "fixed": { ... }       // parameters held constant across all cases
}
```

When a sample-mode param entry is `{}` (empty dict), the sampler uses the
default distribution declared in `cli_v2.py:PARAMETERS[name]["default_dist"]`
— the SynthAorta-aligned distribution. This is the recommended default.

## Files

| File | Purpose |
|---|---|
| `cli_v2.py` | The Python orchestrator. |
| `blender_aorta_v2.py` | The Blender script: 3-segment tube, RMF, smoothstep taper, optional SynthAorta Eq 13 Fourier displacement. |
| `specs_v2/*.json` | Example experiment specs (1 single, 3 sweeps, 2 Sobol — planar + non-planar). |
| `PARAMETERS_v2.md` | Auto-generated parameter reference with literature citations. |
| `scripts_v2/build_v2_gallery_pyvista.py` | PyVista 3D figure builder (gallery + planar-vs-nonplanar + hero). |
| `tests/test_v2_*.py` | 4 test modules, ~67 tests, no Blender required. |
| `figures/v2_*.png` | Shipped PyVista renders (cohort diversity, planar-vs-nonplanar, hero). |

See `PARAMETERS_v2.md` for the full list of parameters and their citations.
