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

## How to customise your sweep

Three workflows, pick whichever fits the change you want to make.

### Recipe A — change values in an existing spec

```bash
cp specs/sweep_severity.json specs/my_sweep.json
$EDITOR specs/my_sweep.json     # change "sweep.high": 0.5, or any param value
python cli.py --spec specs/my_sweep.json --output /tmp/my_run
```

Best when you want to vary the *same* parameter as one of the example specs
but over a different range (e.g. severity 0→0.5 instead of 0→0.9).

### Recipe B — override at the command line, no file edits

```bash
python cli.py --spec specs/sweep_severity.json \
    --param diameter=28 \
    --param arch_height=40 \
    --output /tmp/my_run
```

`--param` is repeatable. In `single` mode it overrides `params`; in `sweep`,
`sample`, or `grid` mode it overrides `fixed` (so the swept axis is preserved).
Best for demo-time live tweaks ("watch what happens with a bigger diameter").

### Recipe C — author a new sweep / sample / grid from scratch

```bash
python cli.py --list-params           # find the parameter name you want
```

Then copy whichever existing spec is structurally closest, edit its params
and ranges, and run. The four modes are documented in
[`PARAMETERS.md`](PARAMETERS.md), one entry per Blender parameter, with
workshop-sensible ranges.

### Spec schema (reference)

```json
{
  "schema_version": "1.0",
  "name": "your_experiment_name",
  "mode": "single|sweep|sample|grid",
  "geometry": "arch_branched_coarctation",

  "params": { ... },     // single mode: explicit values
                         // sample mode: {paramname: {low, high}}

  "sweep": {             // sweep mode only
    "param": "...", "low": 0, "high": 1, "n": 10
  },

  "grid": {              // grid mode only
    "params": {"name1": [v1, v2], "name2": [v3, v4]}
  },

  "sampler": "sobol|lhs|random",   // sample mode only
  "n_cases": 50,
  "seed": 42,

  "fixed": { ... }       // params held constant across all cases
}
```

See `specs/*.json` for one working example per mode.

### How not to break things

The validator catches the easy mistakes:

- **Parameter typo** — `"diametr"` → `Error: Unknown parameter 'diametr'.
  Did you mean 'diameter'?`
- **Bad range** — `"low": 30, "high": 20` → `Error: low must be < high`
- **Too few sample cases** — `n_cases: 2` → `Error: n_cases must be at least 4`
- **Cost runaway** — specs that would produce > 30 cases print a
  warning with an estimated wall-clock before launching Blender

`PARAMETERS.md` and `python cli.py --list-params` are the single source of
truth for what parameter names exist and what ranges are sensible.

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
