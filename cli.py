"""Block A — CLI orchestrator for the aorta geometry generator.

Three modes, all driven by a single JSON spec file:

  - **single**: one parameter set → one case folder
  - **sweep**: vary one named parameter from low to high in N linear steps
  - **sample**: Sobol / LHS / uniform random sampling over multiple parameters

For each generated case, the orchestrator (a) runs the Blender
generator (``blender_aorta_like_generator.py``) as a subprocess to
produce a monolithic STL + sidecar JSON, then (b) calls
``split_patches.split_stl`` to break the STL into AortaCFD-app-style
patches (``inlet.stl``, ``outlet1..N.stl``, ``wall_aorta.stl``), then
(c) writes ``geometry.meta.json`` for provenance.

Output layout::

    <output>/
      <case_id>/
        inlet.stl
        outlet1.stl ... outletN.stl
        wall_aorta.stl
        geometry.meta.json
      sweep_manifest.csv
      <spec_filename>          # copy of the spec for reproducibility

Hand this output directory off to AortaCFD-app's ``scripts/package_cases.py``
to stamp a ``config.json`` on every case.

Usage::

    blender -b -P cli.py -- --spec specs/sample_sobol_50.json --output /tmp/gen_sobol
    # or as a regular Python script (no Blender on PATH; spawns blender):
    python cli.py --spec specs/single_baseline.json --output /tmp/gen_single
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

# Make sibling modules importable when this script is invoked directly
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


logger = logging.getLogger("blender_cli")


# ---------------------------------------------------------------------------
# Parameter <-> Blender CLI flag mapping
# ---------------------------------------------------------------------------

# Spec param name -> Blender generator CLI flag
DIRECT_FLAGS: dict[str, str] = {
    "diameter": "--diameter",
    "ascending_length": "--ascending_length",
    "arch_span": "--arch_span",
    "arch_height": "--arch_height",
    "descending_length": "--descending_length",
    "branch_count": "--branch_count",
    "branch_diameter_ratio": "--branch_diameter_ratio",
    "branch_length": "--branch_length",
    "branch_spacing": "--branch_spacing",
    "branch_tilt_deg": "--branch_tilt_deg",
    "branch_splay_deg": "--branch_splay_deg",
    "coarctation_area_reduction": "--coarctation_area_reduction",
    "coarctation_length": "--coarctation_length",
    "coarctation_centre_fraction": "--coarctation_centre_fraction",
    "proximal_hypoplasia": "--proximal_hypoplasia",
    "noise_amplitude": "--noise_amplitude",
    "noise_scale": "--noise_scale",
    "segments_radial": "--segments_radial",
    "curve_samples": "--curve_samples",
    "cylinder_vertices": "--cylinder_vertices",
}

# Params whose value should be cast to int before becoming a CLI string
INT_FLAGS = {"branch_count", "segments_radial", "curve_samples", "cylinder_vertices"}

# Truthy boolean flags (--coarctation, --roughness)
BOOL_FLAGS = {"coarctation", "roughness"}


# ---------------------------------------------------------------------------
# Parameter reference — single source of truth for discoverability + validation
# ---------------------------------------------------------------------------
#
# Each entry: name -> {type, default, min, max, group, description}
# `min` and `max` are *workshop-sensible* bounds, NOT hard validators —
# the generator will accept values outside these but the result may not
# be physiologically realistic. Used by `--list-params`, the spec validator
# (typo detection + "did you mean"), and the auto-generated PARAMETERS.md.

PARAMETERS: dict[str, dict[str, Any]] = {
    # ── Main tube ────────────────────────────────────────────────────────────
    "diameter": {"type": "float", "default": 24.0, "min": 18.0, "max": 40.0,
                 "group": "Main tube", "description": "Main lumen diameter [mm]"},
    "ascending_length": {"type": "float", "default": 45.0, "min": 30.0, "max": 60.0,
                         "group": "Main tube", "description": "Ascending aorta length [mm]"},
    "arch_span": {"type": "float", "default": 70.0, "min": 55.0, "max": 90.0,
                  "group": "Main tube", "description": "Arch span ascending→descending [mm]"},
    "arch_height": {"type": "float", "default": 35.0, "min": 20.0, "max": 50.0,
                    "group": "Main tube", "description": "Arch rise above ascending top [mm]"},
    "descending_length": {"type": "float", "default": 80.0, "min": 60.0, "max": 250.0,
                          "group": "Main tube", "description": "Descending aorta length [mm]"},
    # ── Branches ─────────────────────────────────────────────────────────────
    "branch_count": {"type": "int", "default": 3, "min": 1, "max": 3,
                     "group": "Branches", "description": "Number of supra-aortic branches (1, 2, or 3)"},
    "branch_diameter_ratio": {"type": "float", "default": 0.42, "min": 0.30, "max": 0.55,
                              "group": "Branches", "description": "Branch diameter / main diameter"},
    "branch_length": {"type": "float", "default": 35.0, "min": 20.0, "max": 60.0,
                      "group": "Branches", "description": "Branch outlet extension length [mm]"},
    "branch_spacing": {"type": "float", "default": 14.0, "min": 8.0, "max": 20.0,
                       "group": "Branches", "description": "Spacing of branch origins along arch [mm]"},
    "branch_tilt_deg": {"type": "float", "default": 60.0, "min": 30.0, "max": 80.0,
                        "group": "Branches", "description": "Branch take-off tilt from +x toward +z [deg]"},
    "branch_splay_deg": {"type": "float", "default": 22.0, "min": 0.0, "max": 40.0,
                         "group": "Branches", "description": "Out-of-plane branch splay magnitude [deg]"},
    # ── Coarctation ──────────────────────────────────────────────────────────
    "coarctation_area_reduction": {"type": "float", "default": 0.65, "min": 0.0, "max": 0.9,
                                   "group": "Coarctation", "description": "Area reduction at the throat (0=none, 0.9=critical)"},
    "coarctation_length": {"type": "float", "default": 30.0, "min": 10.0, "max": 40.0,
                           "group": "Coarctation", "description": "Axial length of the smooth coarctation [mm]"},
    "coarctation_centre_fraction": {"type": "float", "default": 0.72, "min": 0.5, "max": 0.9,
                                    "group": "Coarctation", "description": "Centre location along centreline arc [0-1]"},
    "proximal_hypoplasia": {"type": "float", "default": 0.0, "min": 0.0, "max": 0.25,
                            "group": "Coarctation", "description": "Smooth proximal arch diameter reduction [0-0.25]"},
    # ── Roughness ────────────────────────────────────────────────────────────
    "roughness": {"type": "bool", "default": False, "min": False, "max": True,
                  "group": "Roughness", "description": "Enable mild radial surface roughness"},
    "noise_amplitude": {"type": "float", "default": 0.35, "min": 0.1, "max": 1.0,
                        "group": "Roughness", "description": "Roughness amplitude [mm]"},
    "noise_scale": {"type": "float", "default": 0.08, "min": 0.05, "max": 0.5,
                    "group": "Roughness", "description": "Roughness spatial scale"},
    # ── Geometry mesh resolution (NOT the CFD mesh) ──────────────────────────
    "segments_radial": {"type": "int", "default": 64, "min": 32, "max": 128,
                        "group": "Geometry mesh", "description": "Circumferential ring vertices"},
    "curve_samples": {"type": "int", "default": 220, "min": 100, "max": 400,
                      "group": "Geometry mesh", "description": "Main centreline sample count"},
    "cylinder_vertices": {"type": "int", "default": 64, "min": 32, "max": 128,
                          "group": "Geometry mesh", "description": "Cylinder vertices for branches"},
}


def _format_params_table(markdown: bool = False) -> str:
    """Format the PARAMETERS dict as a human-readable table (terminal or markdown)."""
    groups: dict[str, list[tuple[str, dict[str, Any]]]] = {}
    for name, info in PARAMETERS.items():
        groups.setdefault(info["group"], []).append((name, info))

    if markdown:
        lines = [
            "# Aorta geometry generator — parameter reference",
            "",
            "Generated from `cli.py --list-params --markdown`. Every parameter the "
            "orchestrator passes to the Blender script, with workshop-sensible "
            "ranges. Ranges are *suggestions*, not hard limits — the generator will "
            "still run outside them but the result may not be physiologically realistic.",
            "",
        ]
        for group, entries in groups.items():
            lines += [f"## {group}", ""]
            lines += [
                "| Parameter | Type | Default | Workshop range | Description |",
                "|---|---|---|---|---|",
            ]
            for name, info in entries:
                rng = f"{info['min']}–{info['max']}" if info["type"] != "bool" else "true / false"
                lines.append(f"| `{name}` | {info['type']} | `{info['default']}` | {rng} | {info['description']} |")
            lines.append("")
        return "\n".join(lines)

    lines = [f"Available parameters ({len(PARAMETERS)}):", ""]
    for group, entries in groups.items():
        lines.append(f"  {group}")
        for name, info in entries:
            rng = f"range={info['min']}–{info['max']}" if info["type"] != "bool" else ""
            default = f"default={info['default']}"
            lines.append(f"    {name:30s} {info['type']:6s} {default:18s} {rng:22s} {info['description']}")
        lines.append("")
    return "\n".join(lines)


def _suggest(name: str, candidates: list[str], n: int = 3) -> list[str]:
    import difflib
    return difflib.get_close_matches(name, candidates, n=n, cutoff=0.6)


def _check_param_name(name: str, source: str) -> None:
    """Raise ValueError if name not in PARAMETERS, with a 'did you mean' hint."""
    if name in PARAMETERS:
        return
    suggestions = _suggest(name, list(PARAMETERS.keys()))
    hint = f"  Did you mean {' or '.join(repr(s) for s in suggestions)}?" if suggestions else ""
    raise ValueError(
        f"Unknown parameter {name!r} in {source}.{hint}\n"
        f"  Run `python cli.py --list-params` or see PARAMETERS.md for the full list."
    )


def validate_spec(payload: dict[str, Any], *, source: str = "spec") -> None:
    """Validate a spec payload. Raises ValueError with a helpful message on failure."""
    if not isinstance(payload, dict):
        raise ValueError(f"Spec must be a JSON object, got {type(payload).__name__}")
    mode = payload.get("mode")
    if mode not in {"single", "sweep", "sample", "grid"}:
        raise ValueError(
            f"{source}.mode must be one of single/sweep/sample/grid, got {mode!r}"
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
            raise ValueError(
                f"{source}.sweep.low ({sw['low']}) must be < high ({sw['high']})"
            )
        if int(sw["n"]) < 2:
            raise ValueError(f"{source}.sweep.n must be at least 2, got {sw['n']}")

    if mode == "sample":
        ps = payload.get("params")
        if not isinstance(ps, dict) or not ps:
            raise ValueError(f"{source}.params must have at least one parameter in sample mode")
        for pname, prange in ps.items():
            if not isinstance(prange, dict) or "low" not in prange or "high" not in prange:
                raise ValueError(
                    f"{source}.params[{pname!r}] must be {{'low': ..., 'high': ...}}"
                )
            if float(prange["low"]) >= float(prange["high"]):
                raise ValueError(f"{source}.params[{pname!r}].low must be < high")
        if int(payload.get("n_cases", 10)) < 4:
            raise ValueError(f"{source}.n_cases must be at least 4 in sample mode")

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
                raise ValueError(
                    f"{source}.grid.params[{pname!r}] must be a list of >= 2 values"
                )


def estimate_case_count(payload: dict[str, Any]) -> int:
    """Cheap estimate of how many cases a spec will produce."""
    mode = payload.get("mode")
    if mode == "single":
        return 1
    if mode == "sweep":
        return int(payload.get("sweep", {}).get("n", 0))
    if mode == "sample":
        return int(payload.get("n_cases", 0))
    if mode == "grid":
        from functools import reduce
        from operator import mul
        gps = payload.get("grid", {}).get("params", {})
        return reduce(mul, (len(v) for v in gps.values()), 1) if gps else 0
    return 0


_BLENDER_SECONDS_PER_CASE = 30  # workshop_quick profile, Blender + split (no CFD)


def warn_if_expensive(n_cases: int, threshold: int = 30, stream=sys.stderr) -> None:
    if n_cases > threshold:
        minutes = (n_cases * _BLENDER_SECONDS_PER_CASE) / 60.0
        stream.write(
            f"\n[notice] This spec will produce {n_cases} cases.\n"
            f"         Block A (geometry only) will take ~{minutes:.0f} minutes "
            f"at ~{_BLENDER_SECONDS_PER_CASE}s/case.\n"
            f"         CFD on top adds ~1–2 min/case (workshop_quick template).\n\n"
        )


def _parse_param_override(s: str) -> tuple[str, Any]:
    """Parse a 'key=value' string and cast value based on the PARAMETERS type."""
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
    return key, float(value)


def apply_param_overrides(payload: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    """Merge ``--param key=value`` overrides into the spec.

    In ``single`` mode the overrides go to ``params``.
    In every other mode they go to ``fixed`` so swept ranges are not clobbered.
    """
    if not overrides:
        return payload
    out = json.loads(json.dumps(payload))  # deep copy
    section = "params" if out.get("mode") == "single" else "fixed"
    out.setdefault(section, {}).update(overrides)
    return out


# ---------------------------------------------------------------------------
# Spec parsing
# ---------------------------------------------------------------------------


def load_spec(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    validate_spec(payload, source=f"{path.name}")
    return payload


def expand_cases(spec: dict[str, Any]) -> list[dict[str, Any]]:
    """Return a list of per-case dicts: ``{case_id, params, mode, sampler, sample_index, seed}``."""
    mode = spec.get("mode")
    if mode not in {"single", "sweep", "sample", "grid"}:
        raise ValueError(f"spec.mode must be one of single/sweep/sample/grid, got {mode!r}")
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
                "params": params,
                "mode": "sweep",
                "sampler": None,
                "sample_index": i,
                "seed": None,
            })

    elif mode == "sample":
        from sampler import get_sampler, map_unit_to_physical

        params_spec = spec["params"]
        sampler_name = spec.get("sampler", "sobol")
        n_cases = int(spec.get("n_cases", 10))
        seed = int(spec.get("seed", 0))
        sampler_obj = get_sampler(sampler_name)
        unit = sampler_obj.sample(n_cases, len(params_spec), seed=seed)
        physical = map_unit_to_physical(unit, params_spec)
        for i, sampled in enumerate(physical, start=1):
            params = dict(fixed)
            params.update(sampled)
            cases.append({
                "case_id": f"{name_prefix}_{i:03d}",
                "params": params,
                "mode": "sample",
                "sampler": sampler_name,
                "sample_index": i,
                "seed": seed,
            })

    elif mode == "grid":
        from sampler import grid_product

        gps = spec["grid"]["params"]
        for i, sampled in enumerate(grid_product(gps), start=1):
            params = dict(fixed)
            params.update(sampled)
            cases.append({
                "case_id": f"{name_prefix}_{i:03d}",
                "params": params,
                "mode": "grid",
                "sampler": None,
                "sample_index": i,
                "seed": None,
            })

    return cases


# ---------------------------------------------------------------------------
# Blender invocation
# ---------------------------------------------------------------------------


def build_blender_cmd(
    blender_path: str | Path,
    generator: Path,
    geometry: str,
    params: dict[str, Any],
    out_stl: Path,
    save_blend: bool = False,
) -> list[str]:
    cmd: list[str] = [
        str(blender_path),
        "-b",
        "-P",
        str(generator),
        "--",
        "--geometry",
        geometry,
        "--output",
        str(out_stl),
        "--metadata",
        "--triangulate",
    ]
    # Boolean flags
    if params.get("coarctation_area_reduction", 0.0) > 0.01:
        cmd.append("--coarctation")
    if params.get("roughness", False):
        cmd.append("--roughness")
    if save_blend:
        cmd.append("--save_blend")
    # Direct value flags
    for name, flag in DIRECT_FLAGS.items():
        if name not in params:
            continue
        v = params[name]
        if name in INT_FLAGS:
            v = int(round(float(v)))
        cmd += [flag, str(v)]
    return cmd


def run_blender_and_split(
    blender_path: str | Path,
    generator: Path,
    geometry: str,
    params: dict[str, Any],
    case_dir: Path,
    save_blend: bool = False,
) -> Path:
    case_dir.mkdir(parents=True, exist_ok=True)
    stl_path = case_dir / f"{case_dir.name}.stl"
    cmd = build_blender_cmd(blender_path, generator, geometry, params, stl_path, save_blend=save_blend)
    logger.debug("Blender cmd: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)  # noqa: S603 — args are validated
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


def write_geometry_meta(
    case_dir: Path,
    spec: dict[str, Any],
    params: dict[str, Any],
    mode: str,
    sampler: str | None,
    sample_index: int,
    seed: int | None,
    generator_path: Path,
) -> None:
    patches: dict[str, str] = {}
    for p in sorted(case_dir.glob("*.stl")):
        h = hashlib.blake2b(p.read_bytes(), digest_size=8).hexdigest()
        patches[p.name] = h

    meta = {
        "schema_version": "1.0",
        "case_id": case_dir.name,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "blender_aorta_like_generator",
        "generator_path": str(generator_path),
        "mode": mode,
        "sampler": sampler,
        "sample_index": sample_index,
        "seed": seed,
        "spec_name": spec.get("name"),
        "geometry": spec.get("geometry"),
        "params": params,
        "patch_checksums": patches,
    }
    (case_dir / "geometry.meta.json").write_text(json.dumps(meta, indent=2) + "\n")


# ---------------------------------------------------------------------------
# Manifest CSV
# ---------------------------------------------------------------------------


def write_sweep_manifest(output_dir: Path, rows: list[dict[str, Any]]) -> Path:
    path = output_dir / "sweep_manifest.csv"
    if not rows:
        path.write_text("")
        return path
    fields: list[str] = []
    seen: set[str] = set()
    preferred = ["case_id", "status", "mode", "sampler", "sample_index", "seed", "error"]
    for k in preferred:
        if any(k in r for r in rows) and k not in seen:
            fields.append(k)
            seen.add(k)
    for r in rows:
        for k in r:
            if k not in seen:
                fields.append(k)
                seen.add(k)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="blender_cli", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spec", type=Path, help="Path to a JSON spec defining mode + params.")
    parser.add_argument("--output", "-o", type=Path, help="Output dir for case folders.")
    parser.add_argument("--blender", default=os.environ.get("BLENDER", "blender"),
                        help="Path to Blender executable (default: env $BLENDER or `blender`).")
    parser.add_argument("--generator", type=Path, default=HERE / "blender_aorta_like_generator.py",
                        help="Path to blender_aorta_like_generator.py (default: sibling of this script).")
    parser.add_argument("--param", action="append", default=[], metavar="KEY=VALUE",
                        help="Override a parameter from the spec (repeatable). "
                             "In single mode applies to params; otherwise to fixed. "
                             "e.g. --param diameter=28 --param arch_height=40.")
    parser.add_argument("--save-blend", action="store_true", help="Save .blend alongside STL (debugging).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the cases that would be generated without running Blender.")
    parser.add_argument("--limit", type=int, default=None, help="Stop after N cases (smoke testing).")
    parser.add_argument("--list-params", action="store_true",
                        help="Print every available parameter (name, type, default, workshop range) and exit.")
    parser.add_argument("--markdown", action="store_true",
                        help="With --list-params, emit Markdown (use to regenerate PARAMETERS.md).")
    parser.add_argument("--yes", "-y", action="store_true",
                        help="Skip the cost-estimate warning for large sweeps.")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s %(message)s")

    # Pseudo-subcommand: --list-params is a no-spec discovery flag.
    if args.list_params:
        print(_format_params_table(markdown=args.markdown))
        return 0

    if args.spec is None or args.output is None:
        parser.error("--spec and --output are required (unless --list-params)")

    spec = load_spec(args.spec)

    # Apply --param key=value overrides (after load, before expand)
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
    geometry = spec.get("geometry", "arch_branched")

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
            print(f"  {c['case_id']}: {c['params']}")
        return 0

    rows: list[dict[str, Any]] = []
    for i, c in enumerate(cases, start=1):
        case_dir = args.output / c["case_id"]
        print(f"\n[{i}/{len(cases)}] {c['case_id']}")
        try:
            run_blender_and_split(
                blender_path=args.blender,
                generator=args.generator,
                geometry=geometry,
                params=c["params"],
                case_dir=case_dir,
                save_blend=args.save_blend,
            )
            write_geometry_meta(
                case_dir=case_dir,
                spec=spec,
                params=c["params"],
                mode=c["mode"],
                sampler=c["sampler"],
                sample_index=c["sample_index"],
                seed=c["seed"],
                generator_path=args.generator,
            )
            rows.append({
                "case_id": c["case_id"],
                "status": "ok",
                "mode": c["mode"],
                "sampler": c["sampler"],
                "sample_index": c["sample_index"],
                "seed": c["seed"],
                **c["params"],
            })
            print(f"  OK -> {case_dir}")
        except Exception as e:
            rows.append({
                "case_id": c["case_id"],
                "status": "failed",
                "error": str(e),
                "mode": c["mode"],
                "sampler": c["sampler"],
                "sample_index": c["sample_index"],
                "seed": c["seed"],
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
        # Spec / validator errors — emit cleanly without a Python traceback.
        sys.stderr.write(f"\nError: {e}\n")
        raise SystemExit(2)
