"""Block A (v2) — CLI orchestrator for the healthy-aorta geometry generator.

Companion to ``blender_aorta_v2.py``. Parameter sweeps over a simpler
topology (no branches, no coarctation) with literature-grounded
distributions matching SynthAorta (Bošnjak et al. 2025).

Differences vs v1 (``cli.py``):

  - Different ``PARAMETERS`` table — v2 set of 10 knobs.
  - Different generator: ``blender_aorta_v2.py``.
  - Spec ``geometry`` field accepts only ``"healthy_arch_v2"``.
  - Sample mode supports ``"distribution_overrides"`` per parameter:
      ``{"type": "normal|gumbel|uniform", ...}`` with optional truncation
      ``{"low", "high"}``. Default distributions (when only ``params: {param: {}}``
      is given) come from the SynthAorta paper Table I.

The orchestration scaffolding (sweep manifest, ``--param`` overrides, dry-run,
``--list-params``, validator "did you mean" hints) is identical to v1; we
import those helpers from ``cli.py`` directly so they stay in one place.

Usage::

    python cli_v2.py --spec specs_v2/single_baseline_v2.json --output /tmp/v2_single
    python cli_v2.py --spec specs_v2/sample_sobol_synthaorta_v2.json --output /tmp/v2_sobol
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

# Sibling-module importability when invoked directly.
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

# Reuse v1 helpers that aren't parameter-schema-specific.
from cli import (  # noqa: E402
    _suggest,
    estimate_case_count,
    warn_if_expensive,
    write_sweep_manifest,
)

logger = logging.getLogger("blender_cli_v2")


# ---------------------------------------------------------------------------
# v2 parameter schema
# ---------------------------------------------------------------------------
#
# Each entry: name -> {type, default, min, max, group, description,
#                      default_dist (optional)}
# default_dist is the SynthAorta-grounded sampling distribution used when a
# sample-mode spec lists the parameter under ``params`` with no explicit
# ``distribution_overrides`` entry. Workshop min/max are not hard validators
# — sample-mode truncates to these bounds when sampling.

PARAMETERS: dict[str, dict[str, Any]] = {
    # ── Radii ────────────────────────────────────────────────────────────────
    "r_ascending": {
        "type": "float", "default": 13.7, "min": 8.0, "max": 22.0,
        "group": "Radii", "description": "Ascending aorta radius [mm]",
        "default_dist": {"type": "normal", "mean": 13.7, "std": 2.3,
                         "low": 8.0, "high": 22.0},
        "citation": "Schäfer 2018; Wolak 2008",
    },
    "r_arch": {
        "type": "float", "default": 13.0, "min": 8.0, "max": 20.0,
        "group": "Radii", "description": "Arch radius [mm]",
        "default_dist": {"type": "normal", "mean": 13.0, "std": 2.0,
                         "low": 8.0, "high": 20.0},
        "citation": "Marrocco-Trischitta; Saitta 2022",
    },
    "r_descending": {
        "type": "float", "default": 12.2, "min": 8.0, "max": 20.0,
        "group": "Radii", "description": "Descending aorta radius [mm]",
        "default_dist": {"type": "normal", "mean": 12.2, "std": 2.3,
                         "low": 8.0, "high": 20.0},
        "citation": "Bouti 2017; Schäfer 2018",
    },
    "taper_mode": {
        "type": "str", "default": "smoothstep", "min": None, "max": None,
        "group": "Radii", "description": "Radius blending across segment boundaries",
        "choices": ["piecewise", "linear", "smoothstep"],
    },
    # ── Lengths ──────────────────────────────────────────────────────────────
    "ascending_length": {
        "type": "float", "default": 50.0, "min": 40.0, "max": 90.0,
        "group": "Lengths", "description": "Ascending aorta length [mm]",
        "default_dist": {"type": "uniform", "low": 40.0, "high": 90.0},
        "citation": "Mills 1970; Bouti 2017",
    },
    "descending_length": {
        "type": "float", "default": 200.0, "min": 150.0, "max": 300.0,
        "group": "Lengths", "description": "Descending aorta length [mm]",
        "default_dist": {"type": "uniform", "low": 150.0, "high": 300.0},
        "citation": "anatomy textbooks",
    },
    # ── Arch curvature ───────────────────────────────────────────────────────
    "arch_R_c": {
        "type": "float", "default": 40.4, "min": 25.0, "max": 60.0,
        "group": "Arch curvature", "description": "Arch radius of curvature [mm]",
        "default_dist": {"type": "gumbel", "loc": 40.4, "scale": 2.4,
                         "low": 25.0, "high": 60.0},
        "citation": "Choi 2017; Saitta 2022 (SynthAorta Table I)",
    },
    "arch_angle_deg": {
        "type": "float", "default": 180.0, "min": 120.0, "max": 200.0,
        "group": "Arch curvature", "description": "Subtended angle of the arch arc [deg]",
        "default_dist": {"type": "normal", "mean": 180.0, "std": 15.0,
                         "low": 120.0, "high": 200.0},
        "citation": "engineering default (Madhwal arch-type classification context)",
    },
    "arch_tilt_deg": {
        "type": "float", "default": 0.0, "min": -30.0, "max": 30.0,
        "group": "Arch curvature",
        "description": "Rotation of the arch+descending around the inlet z-axis [deg]",
        "default_dist": {"type": "normal", "mean": 0.0, "std": 8.0,
                         "low": -25.0, "high": 25.0},
        "citation": "anatomy textbooks — typical leftward tilt 5-15°",
    },
    "junction_blend_mm": {
        "type": "float", "default": 12.0, "min": 0.0, "max": 40.0,
        "group": "Arch curvature",
        "description": "Cubic-Bezier blend width at each arch junction [mm] "
                       "(0 = sharp circular-arc corners)",
    },
    # ── Non-planar Fourier (SynthAorta Eq 13) ────────────────────────────────
    # When δ_3 = δ_4 = 0 (the scalar defaults), the centreline stays in the
    # xz-plane (backwards-compat with v2 pre-2026-05-20). When sampled around
    # 1.0 ± 0.09 (the default distribution), the centreline acquires the
    # ~few-mm y-axis wobble that matches SynthAorta's published non-planar
    # statistics. γ_1, γ_2 are hardcoded to 1.0 inside the generator (paper
    # sensitivity analysis: non-influential).
    "delta_3": {
        "type": "float", "default": 0.0, "min": 0.0, "max": 1.5,
        "group": "Non-planar Fourier",
        "description": "SynthAorta δ_3: cos(2w·||x||) multiplier [0=planar, 1=SynthAorta nominal]",
        "default_dist": {"type": "normal", "mean": 1.0, "std": 0.09,
                         "low": 0.5, "high": 1.5},
        "citation": "Bošnjak et al. 2025 Table I (δ_3)",
    },
    "delta_4": {
        "type": "float", "default": 0.0, "min": 0.0, "max": 1.5,
        "group": "Non-planar Fourier",
        "description": "SynthAorta δ_4: sin(2w·||x||) multiplier [0=planar, 1=SynthAorta nominal]",
        "default_dist": {"type": "normal", "mean": 1.0, "std": 0.09,
                         "low": 0.5, "high": 1.5},
        "citation": "Bošnjak et al. 2025 Table I (δ_4)",
    },
    # ── Mesh resolution (NOT the CFD mesh; the Blender output mesh) ──────────
    "segments_radial": {
        "type": "int", "default": 96, "min": 32, "max": 192,
        "group": "Geometry mesh", "description": "Circumferential ring vertices",
    },
    "curve_samples": {
        "type": "int", "default": 300, "min": 100, "max": 600,
        "group": "Geometry mesh", "description": "Total centreline sample count",
    },
}


# Param name -> Blender flag
DIRECT_FLAGS: dict[str, str] = {
    "r_ascending": "--r_ascending",
    "r_arch": "--r_arch",
    "r_descending": "--r_descending",
    "taper_mode": "--taper_mode",
    "ascending_length": "--ascending_length",
    "descending_length": "--descending_length",
    "arch_R_c": "--arch_R_c",
    "arch_angle_deg": "--arch_angle_deg",
    "arch_tilt_deg": "--arch_tilt_deg",
    "junction_blend_mm": "--junction_blend_mm",
    "delta_3": "--delta_3",
    "delta_4": "--delta_4",
    "segments_radial": "--segments_radial",
    "curve_samples": "--curve_samples",
}

INT_FLAGS = {"segments_radial", "curve_samples"}


# ---------------------------------------------------------------------------
# Discovery + validation
# ---------------------------------------------------------------------------


def _format_params_table(markdown: bool = False) -> str:
    groups: dict[str, list[tuple[str, dict[str, Any]]]] = {}
    for name, info in PARAMETERS.items():
        groups.setdefault(info["group"], []).append((name, info))

    if markdown:
        lines = [
            "# Aorta geometry generator v2 — parameter reference",
            "",
            "Generated from `cli_v2.py --list-params --markdown`. Defaults and "
            "default sample-mode distributions are from the SynthAorta paper "
            "(Bošnjak et al. 2025) Table I unless noted.",
            "",
        ]
        for group, entries in groups.items():
            lines += [f"## {group}", ""]
            lines += [
                "| Parameter | Type | Default | Workshop range | Default distribution | Citation | Description |",
                "|---|---|---|---|---|---|---|",
            ]
            for name, info in entries:
                if info["type"] == "str":
                    rng = " / ".join(info.get("choices", []))
                else:
                    rng = f"{info['min']}–{info['max']}"
                dist = info.get("default_dist")
                if dist is None:
                    dist_str = "(fixed)"
                else:
                    dist_str = (
                        f"{dist['type']}"
                        + (f"(μ={dist['mean']}, σ={dist['std']})"
                           if dist["type"] == "normal" else "")
                        + (f"(loc={dist['loc']}, scale={dist['scale']})"
                           if dist["type"] == "gumbel" else "")
                        + (f"({dist['low']}-{dist['high']})"
                           if dist["type"] == "uniform" else "")
                    )
                cite = info.get("citation", "—")
                lines.append(
                    f"| `{name}` | {info['type']} | `{info['default']}` | {rng} "
                    f"| {dist_str} | {cite} | {info['description']} |"
                )
            lines.append("")
        return "\n".join(lines)

    lines = [f"Available v2 parameters ({len(PARAMETERS)}):", ""]
    for group, entries in groups.items():
        lines.append(f"  {group}")
        for name, info in entries:
            if info["type"] == "str":
                rng = f"choices={'/'.join(info.get('choices', []))}"
            else:
                rng = f"range={info['min']}–{info['max']}"
            default = f"default={info['default']}"
            cite = info.get("citation", "")
            lines.append(
                f"    {name:20s} {info['type']:6s} {default:18s} {rng:28s} "
                f"{info['description']}"
                + (f"  [{cite}]" if cite else "")
            )
        lines.append("")
    return "\n".join(lines)


def _check_param_name(name: str, source: str) -> None:
    if name in PARAMETERS:
        return
    suggestions = _suggest(name, list(PARAMETERS.keys()))
    hint = f"  Did you mean {' or '.join(repr(s) for s in suggestions)}?" if suggestions else ""
    raise ValueError(
        f"Unknown parameter {name!r} in {source}.{hint}\n"
        f"  Run `python cli_v2.py --list-params` or see PARAMETERS_v2.md for the full list."
    )


def validate_spec(payload: dict[str, Any], *, source: str = "spec") -> None:
    if not isinstance(payload, dict):
        raise ValueError(f"Spec must be a JSON object, got {type(payload).__name__}")
    mode = payload.get("mode")
    if mode not in {"single", "sweep", "sample", "grid"}:
        raise ValueError(f"{source}.mode must be one of single/sweep/sample/grid, got {mode!r}")

    geom = payload.get("geometry")
    if geom is not None and geom != "healthy_arch_v2":
        raise ValueError(
            f"{source}.geometry must be 'healthy_arch_v2' for cli_v2.py, got {geom!r}. "
            f"v1 specs (geometry='arch_branched_coarctation') are not compatible — use cli.py."
        )

    for section in ("params", "fixed"):
        block = payload.get(section)
        if block is None:
            continue
        if not isinstance(block, dict):
            raise ValueError(f"{source}.{section} must be an object")
        for k in block:
            _check_param_name(k, f"{source}.{section}")

    if mode == "sweep":
        sw = payload.get("sweep")
        if not isinstance(sw, dict):
            raise ValueError(f"{source}.sweep must be {{'param','low','high','n'}}")
        for key in ("param", "low", "high", "n"):
            if key not in sw:
                raise ValueError(f"{source}.sweep is missing required key {key!r}")
        _check_param_name(sw["param"], f"{source}.sweep.param")
        if float(sw["low"]) >= float(sw["high"]):
            raise ValueError(f"{source}.sweep.low ({sw['low']}) must be < high ({sw['high']})")
        if int(sw["n"]) < 2:
            raise ValueError(f"{source}.sweep.n must be at least 2, got {sw['n']}")

    if mode == "sample":
        ps = payload.get("params")
        if not isinstance(ps, dict) or not ps:
            raise ValueError(f"{source}.params must have at least one parameter in sample mode")
        for pname, prange in ps.items():
            _check_param_name(pname, f"{source}.params")
            if not isinstance(prange, dict):
                raise ValueError(
                    f"{source}.params[{pname!r}] must be a dict (empty {{}} to use the "
                    f"default distribution, or {{'low','high'}} to override range)"
                )
            if "low" in prange and "high" in prange:
                if float(prange["low"]) >= float(prange["high"]):
                    raise ValueError(f"{source}.params[{pname!r}].low must be < high")
        if int(payload.get("n_cases", 10)) < 4:
            raise ValueError(f"{source}.n_cases must be at least 4 in sample mode")

        # Validate distribution_overrides if present
        overrides = payload.get("distribution_overrides", {})
        if not isinstance(overrides, dict):
            raise ValueError(f"{source}.distribution_overrides must be a dict")
        for pname, dist in overrides.items():
            _check_param_name(pname, f"{source}.distribution_overrides")
            _validate_distribution(dist, source=f"{source}.distribution_overrides[{pname!r}]")

    if mode == "grid":
        grid = payload.get("grid")
        if not isinstance(grid, dict) or "params" not in grid:
            raise ValueError(f"{source}.grid must be {{'params': {{name: [v1,v2,...], ...}}}}")
        gps = grid["params"]
        if not isinstance(gps, dict) or not gps:
            raise ValueError(f"{source}.grid.params must have at least one parameter")
        for pname, values in gps.items():
            _check_param_name(pname, f"{source}.grid.params")
            if not isinstance(values, list) or len(values) < 2:
                raise ValueError(f"{source}.grid.params[{pname!r}] must be a list of >= 2 values")


def _validate_distribution(dist: Any, *, source: str) -> None:
    if not isinstance(dist, dict) or "type" not in dist:
        raise ValueError(f"{source} must be {{'type': 'normal|gumbel|uniform', ...}}")
    dtype = dist["type"]
    if dtype == "normal":
        for k in ("mean", "std"):
            if k not in dist:
                raise ValueError(f"{source}.{k} required for type='normal'")
        if float(dist["std"]) <= 0:
            raise ValueError(f"{source}.std must be > 0")
    elif dtype == "gumbel":
        for k in ("loc", "scale"):
            if k not in dist:
                raise ValueError(f"{source}.{k} required for type='gumbel'")
        if float(dist["scale"]) <= 0:
            raise ValueError(f"{source}.scale must be > 0")
    elif dtype == "uniform":
        for k in ("low", "high"):
            if k not in dist:
                raise ValueError(f"{source}.{k} required for type='uniform'")
        if float(dist["low"]) >= float(dist["high"]):
            raise ValueError(f"{source}.low must be < high for type='uniform'")
    else:
        raise ValueError(f"{source}.type must be one of normal|gumbel|uniform, got {dtype!r}")


def _parse_param_override(s: str) -> tuple[str, Any]:
    if "=" not in s:
        raise ValueError(f"--param expects key=value, got {s!r}")
    key, _, value = s.partition("=")
    key, value = key.strip(), value.strip()
    _check_param_name(key, "--param")
    t = PARAMETERS[key]["type"]
    if t == "bool":
        return key, value.lower() in {"true", "1", "yes", "on"}
    if t == "int":
        return key, int(value)
    if t == "str":
        return key, value
    return key, float(value)


def apply_param_overrides(payload: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    if not overrides:
        return payload
    out = json.loads(json.dumps(payload))
    section = "params" if out.get("mode") == "single" else "fixed"
    out.setdefault(section, {}).update(overrides)
    return out


# ---------------------------------------------------------------------------
# Sample-mode distribution mapping (ppf-based)
# ---------------------------------------------------------------------------


def _map_unit_to_distribution(u: np.ndarray, dist: dict[str, Any]) -> np.ndarray:
    """Map [0,1] samples through the inverse CDF of ``dist`` (truncated to [low, high])."""
    from scipy import stats

    dtype = dist["type"]
    low = float(dist.get("low", -np.inf))
    high = float(dist.get("high", np.inf))

    if dtype == "normal":
        mu, sigma = float(dist["mean"]), float(dist["std"])
        if np.isfinite(low) or np.isfinite(high):
            a = (low - mu) / sigma if np.isfinite(low) else -np.inf
            b = (high - mu) / sigma if np.isfinite(high) else np.inf
            return stats.truncnorm.ppf(u, a, b, loc=mu, scale=sigma)
        return stats.norm.ppf(u, loc=mu, scale=sigma)

    if dtype == "gumbel":
        loc, scale = float(dist["loc"]), float(dist["scale"])
        v = stats.gumbel_r.ppf(u, loc=loc, scale=scale)
        if np.isfinite(low) or np.isfinite(high):
            v = np.clip(v, low if np.isfinite(low) else v.min(),
                          high if np.isfinite(high) else v.max())
        return v

    if dtype == "uniform":
        return low + (high - low) * u

    raise ValueError(f"Unknown distribution type {dtype!r}")


def sample_cases(spec: dict[str, Any]) -> list[dict[str, Any]]:
    """Sample N cases from the per-parameter distributions in a v2 sample-mode spec."""
    from sampler import get_sampler

    params_spec = spec["params"]
    sampler_name = spec.get("sampler", "sobol")
    n_cases = int(spec.get("n_cases", 10))
    seed = int(spec.get("seed", 0))
    overrides = spec.get("distribution_overrides", {})

    # Resolve effective distribution per parameter
    dists: dict[str, dict[str, Any]] = {}
    for pname, p_entry in params_spec.items():
        if pname in overrides:
            dists[pname] = overrides[pname]
        elif "low" in p_entry and "high" in p_entry:
            dists[pname] = {"type": "uniform",
                            "low": float(p_entry["low"]),
                            "high": float(p_entry["high"])}
        elif "default_dist" in PARAMETERS[pname]:
            dists[pname] = PARAMETERS[pname]["default_dist"]
        else:
            raise ValueError(
                f"Parameter {pname!r} has no default distribution and no explicit "
                f"override or low/high in spec.params"
            )

    names = list(dists.keys())
    sampler_obj = get_sampler(sampler_name)
    unit = sampler_obj.sample(n_cases, len(names), seed=seed)

    cases: list[dict[str, Any]] = []
    for row in unit:
        case: dict[str, Any] = {}
        for j, name in enumerate(names):
            val = float(_map_unit_to_distribution(np.asarray([row[j]]), dists[name])[0])
            case[name] = val
        cases.append(case)
    return cases


# ---------------------------------------------------------------------------
# Spec parsing + expansion
# ---------------------------------------------------------------------------


def load_spec(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    validate_spec(payload, source=f"{path.name}")
    return payload


def expand_cases(spec: dict[str, Any]) -> list[dict[str, Any]]:
    mode = spec.get("mode")
    fixed = dict(spec.get("fixed") or {})
    name_prefix = spec.get("case_prefix", "case")
    cases: list[dict[str, Any]] = []

    if mode == "single":
        params = dict(fixed)
        params.update(spec.get("params", {}))
        case_id = spec.get("case_id", f"{name_prefix}_001")
        cases.append({"case_id": case_id, "params": params, "mode": "single",
                      "sampler": None, "sample_index": 0, "seed": None})

    elif mode == "sweep":
        from sampler import linear_sweep

        sw = spec["sweep"]
        values = linear_sweep(float(sw["low"]), float(sw["high"]), int(sw["n"]))
        for i, v in enumerate(values, start=1):
            params = dict(fixed)
            params[sw["param"]] = float(v)
            cases.append({
                "case_id": f"{name_prefix}_{i:03d}",
                "params": params, "mode": "sweep",
                "sampler": None, "sample_index": i, "seed": None,
            })

    elif mode == "sample":
        sampled = sample_cases(spec)
        sampler_name = spec.get("sampler", "sobol")
        seed = int(spec.get("seed", 0))
        for i, s in enumerate(sampled, start=1):
            params = dict(fixed)
            params.update(s)
            cases.append({
                "case_id": f"{name_prefix}_{i:03d}",
                "params": params, "mode": "sample",
                "sampler": sampler_name, "sample_index": i, "seed": seed,
            })

    elif mode == "grid":
        from sampler import grid_product

        gps = spec["grid"]["params"]
        for i, sampled in enumerate(grid_product(gps), start=1):
            params = dict(fixed)
            params.update(sampled)
            cases.append({
                "case_id": f"{name_prefix}_{i:03d}",
                "params": params, "mode": "grid",
                "sampler": None, "sample_index": i, "seed": None,
            })

    else:
        raise ValueError(f"spec.mode must be one of single/sweep/sample/grid, got {mode!r}")

    return cases


# ---------------------------------------------------------------------------
# Blender invocation
# ---------------------------------------------------------------------------


def build_blender_cmd(blender_path: str | Path, generator: Path,
                      params: dict[str, Any], out_stl: Path,
                      save_blend: bool = False) -> list[str]:
    cmd: list[str] = [
        str(blender_path), "-b", "-P", str(generator), "--",
        "--output", str(out_stl), "--metadata", "--triangulate",
    ]
    if save_blend:
        cmd.append("--save_blend")
    for name, flag in DIRECT_FLAGS.items():
        if name not in params:
            continue
        v = params[name]
        if name in INT_FLAGS:
            v = int(round(float(v)))
        cmd += [flag, str(v)]
    return cmd


def run_blender_and_split(blender_path: str | Path, generator: Path,
                          params: dict[str, Any], case_dir: Path,
                          save_blend: bool = False) -> Path:
    case_dir.mkdir(parents=True, exist_ok=True)
    stl_path = case_dir / f"{case_dir.name}.stl"
    cmd = build_blender_cmd(blender_path, generator, params, stl_path, save_blend=save_blend)
    logger.debug("Blender cmd: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)  # noqa: S603
    if result.returncode != 0:
        raise RuntimeError(
            f"Blender failed for {case_dir.name} (rc={result.returncode}).\n"
            f"STDOUT (tail):\n{result.stdout[-2000:]}\n"
            f"STDERR (tail):\n{result.stderr[-2000:]}"
        )
    from split_patches import split_stl  # noqa: E402

    split_stl(stl_path, case_dir)
    return case_dir


# ---------------------------------------------------------------------------
# Per-case provenance
# ---------------------------------------------------------------------------


def write_geometry_meta(case_dir: Path, spec: dict[str, Any], params: dict[str, Any],
                        mode: str, sampler: str | None, sample_index: int,
                        seed: int | None, generator_path: Path) -> None:
    patches: dict[str, str] = {}
    for p in sorted(case_dir.glob("*.stl")):
        h = hashlib.blake2b(p.read_bytes(), digest_size=8).hexdigest()
        patches[p.name] = h

    meta = {
        "schema_version": "2.0",
        "case_id": case_dir.name,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "blender_aorta_v2",
        "generator_path": str(generator_path),
        "mode": mode,
        "sampler": sampler,
        "sample_index": sample_index,
        "seed": seed,
        "spec_name": spec.get("name"),
        "geometry": spec.get("geometry", "healthy_arch_v2"),
        "params": params,
        "patch_checksums": patches,
    }
    (case_dir / "geometry.meta.json").write_text(json.dumps(meta, indent=2) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="blender_cli_v2", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spec", type=Path, help="Path to a v2 JSON spec.")
    parser.add_argument("--output", "-o", type=Path, help="Output dir for case folders.")
    parser.add_argument("--blender", default=os.environ.get("BLENDER", "blender"),
                        help="Path to Blender executable (default: env $BLENDER or `blender`).")
    parser.add_argument("--generator", type=Path, default=HERE / "blender_aorta_v2.py",
                        help="Path to blender_aorta_v2.py (default: sibling of this script).")
    parser.add_argument("--param", action="append", default=[], metavar="KEY=VALUE",
                        help="Override a parameter from the spec (repeatable).")
    parser.add_argument("--save-blend", action="store_true",
                        help="Save .blend alongside STL (debugging).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the cases that would be generated without running Blender.")
    parser.add_argument("--limit", type=int, default=None, help="Stop after N cases.")
    parser.add_argument("--list-params", action="store_true",
                        help="Print every v2 parameter (name, type, default, range, default distribution).")
    parser.add_argument("--markdown", action="store_true",
                        help="With --list-params, emit Markdown (to regenerate PARAMETERS_v2.md).")
    parser.add_argument("--yes", "-y", action="store_true",
                        help="Skip the cost-estimate warning for large sweeps.")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(levelname)s %(message)s")

    if args.list_params:
        print(_format_params_table(markdown=args.markdown))
        return 0

    if args.spec is None or args.output is None:
        parser.error("--spec and --output are required (unless --list-params)")

    spec = load_spec(args.spec)

    overrides: dict[str, Any] = {}
    for s in args.param:
        k, v = _parse_param_override(s)
        overrides[k] = v
    if overrides:
        spec = apply_param_overrides(spec, overrides)

    cases = expand_cases(spec)
    if args.limit is not None:
        cases = cases[: args.limit]

    args.output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.spec, args.output / args.spec.name)

    print(f"Spec       : {args.spec}")
    print(f"Mode       : {spec['mode']}")
    print(f"Cases      : {len(cases)}")
    if overrides:
        print(f"Overrides  : {overrides}")
    print(f"Output     : {args.output}")
    print(f"Generator  : {args.generator}")
    print(f"Blender    : {args.blender}")

    if not args.dry_run and not args.yes:
        warn_if_expensive(len(cases))

    if args.dry_run:
        print("\n[dry-run] Cases that WOULD be generated:")
        for c in cases:
            short = {k: (round(v, 3) if isinstance(v, float) else v)
                     for k, v in c["params"].items()}
            print(f"  {c['case_id']}: {short}")
        return 0

    rows: list[dict[str, Any]] = []
    for i, c in enumerate(cases, start=1):
        case_dir = args.output / c["case_id"]
        print(f"\n[{i}/{len(cases)}] {c['case_id']}")
        try:
            run_blender_and_split(
                blender_path=args.blender,
                generator=args.generator,
                params=c["params"],
                case_dir=case_dir,
                save_blend=args.save_blend,
            )
            write_geometry_meta(
                case_dir=case_dir, spec=spec, params=c["params"],
                mode=c["mode"], sampler=c["sampler"],
                sample_index=c["sample_index"], seed=c["seed"],
                generator_path=args.generator,
            )
            rows.append({
                "case_id": c["case_id"], "status": "ok",
                "mode": c["mode"], "sampler": c["sampler"],
                "sample_index": c["sample_index"], "seed": c["seed"],
                **c["params"],
            })
            print(f"  OK -> {case_dir}")
        except Exception as e:
            rows.append({
                "case_id": c["case_id"], "status": "failed", "error": str(e),
                "mode": c["mode"], "sampler": c["sampler"],
                "sample_index": c["sample_index"], "seed": c["seed"],
                **c["params"],
            })
            print(f"  FAIL: {e}")

    manifest_path = write_sweep_manifest(args.output, rows)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_fail = len(rows) - n_ok
    print(f"\nSweep manifest -> {manifest_path}")
    print(f"Done. {n_ok} succeeded, {n_fail} failed.")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as e:
        sys.stderr.write(f"\nError: {e}\n")
        raise SystemExit(2)
