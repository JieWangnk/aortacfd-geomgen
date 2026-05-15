# Aorta geometry generator (Blender)

Parametric synthetic aortic-arch geometries for CFD studies. Standalone
Blender script + a thin Python CLI that supports three modes:

- **single** — one parameter set, one STL
- **sweep** — vary one named parameter from low to high in N linear steps
- **sample** — Sobol / LHS / uniform-random sampling for ML training data
  or global sensitivity studies

Output is per-case folders with already-split inlet/outlet/wall STL
patches, ready to be packaged into runnable AortaCFD-app cases via
[`AortaCFD-app/scripts/package_cases.py`](../AortaCFD-app/scripts/package_cases.py).

This is **Block A** of the four-block parametric-study workflow
(generation → packaging → running → aggregating). It is independent of
AortaCFD-app — hand-off is by filesystem (a folder of cases), not by
Python imports.

## Install

You need:

- **Blender 3.x+** (tested with 5.1) on your PATH (`blender --version`)
- **Python 3.10+** with `numpy`, `numpy-stl`, `scipy`

```bash
pip install numpy numpy-stl scipy
```

## Three modes

### 1. Single case (smoke test)

```bash
python cli.py --spec specs/single_baseline.json --output /tmp/gen_single
```

→ `/tmp/gen_single/baseline/{inlet,outlet1,outlet2,outlet3,wall_aorta}.stl + geometry.meta.json`

### 2. Linear sweep (one parameter)

```bash
python cli.py --spec specs/sweep_severity.json --output /tmp/gen_sweep
```

→ 10 cases varying `coarctation_area_reduction` from 0 → 0.9.

### 3. Statistical sample (multi-parameter)

```bash
python cli.py --spec specs/sample_sobol_50.json --output /tmp/gen_sobol
```

→ 50 cases sampled from a 6-D parameter space using a scrambled Sobol
sequence.

## Spec schema

```json
{
  "schema_version": "1.0",
  "name": "your_experiment_name",
  "mode": "single|sweep|sample",
  "geometry": "arch_branched_coarctation",
  "params": { ... },     // single mode: explicit values
  "sweep": { ... },      // sweep mode: {param, low, high, n}
  "params": { ... },     // sample mode: {paramname: {low, high}}
  "sampler": "sobol|lhs|random",
  "n_cases": 50,
  "seed": 42,
  "fixed": { ... }       // params held constant for all cases
}
```

See `specs/*.json` for working examples. New parameter names go through
`DIRECT_FLAGS` in [`cli.py`](cli.py) — add a mapping if you want to
expose a new Blender CLI flag.

## Output per case

```
<output>/<case_id>/
  inlet.stl              # ascending aorta inlet patch
  outlet1.stl            # descending aorta outlet patch
  outlet2.stl ...        # one per supra-aortic branch
  wall_aorta.stl         # everything else (the vessel wall)
  geometry.meta.json     # parameters, seed, sampler, patch checksums
```

Plus at the output root:

- `sweep_manifest.csv` — one row per case (case_id, status, params, error if failed)
- A copy of the spec JSON for reproducibility

## What's in `geometry.meta.json`

```json
{
  "schema_version": "1.0",
  "case_id": "sobol_007",
  "generated_utc": "2026-05-15T13:42:01+00:00",
  "generator": "blender_aorta_like_generator",
  "mode": "sample",
  "sampler": "sobol",
  "sample_index": 7,
  "seed": 42,
  "spec_name": "sample_sobol_50",
  "geometry": "arch_branched_coarctation",
  "params": { ... },
  "patch_checksums": {
    "inlet.stl": "1234abcd...",
    "wall_aorta.stl": "5678efef..."
  }
}
```

The case packager (`AortaCFD-app/scripts/package_cases.py`) reads this
file to produce a unified `case.meta.json` that downstream tooling
(`compare_cohort.py`) joins on `case_id` when aggregating CFD results.

## Composing with the AortaCFD-app workflow

```bash
# 1. Block A: generate the geometries
python cli.py --spec specs/sample_sobol_50.json --output /tmp/gen_sobol

# 2. Block B: stamp configs (in AortaCFD-app)
cd ~/GitHub/AortaCFD-app
python -m scripts.package_cases /tmp/gen_sobol \
    --config-template examples/templates/config_workshop_quick.json \
    --output cases_input/sobol_demo

# 3. Block C: run them
python run_batch.py --cases-dir cases_input/sobol_demo --workers 2 --quick

# 4. Block D: aggregate
python -m scripts.compare_cohort output/
```

## Files in this directory

| File | Purpose |
|---|---|
| `blender_aorta_like_generator.py` | The Blender script that builds one aorta from CLI flags. Standalone — runs inside `blender -b -P ...`. |
| `cli.py` | The Python orchestrator. Reads a spec, spawns Blender per case, calls the splitter. |
| `sampler.py` | Sobol / LHS / uniform random samplers. Pure numpy + scipy, no Blender. |
| `split_patches.py` | Flood-fill STL splitter. Reads monolithic STL + sidecar JSON, writes split patches. |
| `specs/*.json` | Example experiment specs (single / sweep / sample). |
| `tests/` | Unit tests for the samplers and spec expansion. |

## Tests

```bash
pytest tests/
```

(no Blender required; tests cover sampler + spec expansion only.)
