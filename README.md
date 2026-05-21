# Aorta geometry generator (Blender)

Parametric synthetic aorta STLs for CFD studies. Three flavours of CLI in
increasing complexity:

| Variant | Topology | Knobs | Use when |
|---|---|---|---|
| [**v3**](./README_v3.md) — `cli_v3.py` | Pipe U-bend, no branches | 12 | You just want a clean pipe U-bend STL |
| [**v2**](./README_v2.md) — `cli_v2.py` | Healthy arch, no branches | 18 + sample modes | ML training data, SynthAorta-comparable studies, sensitivity analysis |
| **v1** — `cli.py` | Aortic arch + 1-3 branches + coarctation | 21 | Pathology studies, branch sweeps |

All three produce the same case-folder layout (`inlet.stl`,
`outlet1..N.stl`, `wall_aorta.stl`, `geometry.meta.json`) so they hand off
to `AortaCFD-app/scripts/package_cases.py` identically.

## Install

```bash
pip install -r requirements.txt        # numpy, scipy, numpy-stl + pyvista/matplotlib
```

Plus **Blender 3.x+** on your PATH:

| OS | Install |
|---|---|
| Windows | [Installer from blender.org](https://www.blender.org/download/) — check "Add to PATH" during install |
| macOS | `brew install --cask blender` |
| Linux | `sudo apt install blender` (or `dnf`, `pacman`, `snap`) |

Verify with `blender --version`.

If you hit `error: externally-managed-environment` on Ubuntu 24.04+, use
a venv: `python3 -m venv .venv && source .venv/bin/activate`.

## Quick start (v3 — recommended for new users)

```bash
# Baseline U-bend
python3 cli_v3.py --spec specs_v3/single_baseline_v3.json --output /tmp/v3 --yes

# Discover all 12 knobs
python3 cli_v3.py --list-params

# Tweak any knob
python3 cli_v3.py --spec specs_v3/single_baseline_v3.json --output /tmp/v3 --yes \
    --param r_inlet=16 --param twist_deg=20
```

See [`README_v3.md`](./README_v3.md) for the full walkthrough.

## v3 in one image

![v3 baseline](figures/v3_baseline_hero.png)

*v3 baseline render — `cli_v3.py --spec specs_v3/single_baseline_v3.json`.
Three views of the canonical U-bend.*

![16-case Sobol gallery](figures/v3_sobol_gallery_with_twist.png)

*16 Sobol samples over the 8 v3 knobs with gradual twist forced into
`[10°, 30°]` — every case shows visible non-planarity.*

## v2 in one paragraph

v2 has the full SynthAorta-compatible 18-knob parameter space (per-segment
radii, arch curvature, lengths, non-planar Fourier multipliers, mesh
resolution) plus Sobol / LHS / random sampling. Use it when you need ML
training data, sensitivity analysis, or distributions grounded in clinical
literature. See [`README_v2.md`](./README_v2.md).

## v1 in one paragraph

v1 is the original generator with supra-aortic branches and coarctation
for pathology studies. Documentation below.

---

## v1 (legacy)

![v1 baseline — branched arch with coarctation](figures/synthetic_aorta_geometry.png)

```bash
python3 cli.py --spec specs/single_baseline.json --output /tmp/v1 --yes
python3 cli.py --list-params         # 21 knobs
```

| Mode | Spec example |
|---|---|
| `single` | `specs/single_baseline.json` |
| `sweep` | `specs/sweep_severity.json` — 10-step coarctation sweep |
| `sample` | `specs/sample_sobol_50.json` — 50-case Sobol cohort |

See `PARAMETERS.md` for the full v1 knob list with workshop-sensible ranges
and `figures/severity_sweep_demo.png` for the canonical pathology-sweep
output.

## Block A — what this is in the bigger picture

This repo is **Block A** of the four-block aortic CFD workflow:

1. **A. Geometry** — this repo (v1/v2/v3) → STL cases
2. **B. Packaging** — `AortaCFD-app/scripts/package_cases.py` stamps config on each case
3. **C. CFD** — `AortaCFD-app/run_batch.py` runs OpenFOAM
4. **D. Aggregation** — `AortaCFD-app/scripts/compare_cohort.py` joins results

A and B are decoupled — hand-off is by filesystem (a folder of cases), not
Python imports.

## Branches

- `main` / `v3` — latest (everything: v1 + v2 + v3)
- `v2` — frozen snapshot when v2 was complete (v1 + v2 only)
- `v1` — frozen snapshot at v1-complete (v1 only)
