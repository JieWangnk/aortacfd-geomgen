"""Block A (v3) — minimal 5-knob healthy-aorta CLI.

Thin wrapper over v2. Exposes just the parameters a clinician or
engineer typically dials by hand:

    r_inlet           — inlet (ascending) radius [mm]
    r_outlet          — outlet (descending) radius [mm]
    arch_width_mm     — arch horizontal extent (ascending→descending) [mm]
    arch_height_mm    — arch peak height above ascending top [mm]
    torsion_deg       — rotation of arch+descending around inlet z-axis [deg]

Optional secondary knobs (sensible defaults from the workshop baseline):
    ascending_length  — straight ascending length before the arch [mm]
    descending_length — straight descending length after the arch [mm]

Everything else (taper mode, junction blend, non-planar Fourier, mesh
resolution, …) is fixed inside this script — v3 is the "I just want
five knobs that match the picture" interface. Use cli_v2.py if you
need finer control.

Internally the 5 v3 knobs are translated to v2 names and passed to
``blender_aorta_v2.py`` via the same orchestration pipeline cli_v2.py
uses (so split_patches.py, manifest CSV, geometry.meta.json behave
identically).

Usage::

    python cli_v3.py --spec specs_v3/single_baseline_v3.json --output /tmp/v3_one
    python cli_v3.py --list-params
    python cli_v3.py --spec specs_v3/single_baseline_v3.json --output /tmp/v3 \\
        --param torsion_deg=15 --param arch_height_mm=50
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

# Reuse the v2 + v1 plumbing — only the parameter schema differs.
from cli import _suggest, warn_if_expensive, write_sweep_manifest  # noqa: E402
from cli_v2 import (  # noqa: E402
    _resolve_arch_params,
    run_blender_and_split,
    write_geometry_meta,
)
from sampler import linear_sweep  # noqa: E402

logger = logging.getLogger("blender_cli_v3")


# ---------------------------------------------------------------------------
# v3 parameter schema — 5 primary + 2 optional knobs only
# ---------------------------------------------------------------------------

PARAMETERS: dict[str, dict[str, Any]] = {
    "r_inlet": {
        "type": "float", "default": 14.0, "min": 8.0, "max": 22.0,
        "group": "Radii", "description": "Inlet (ascending) radius [mm]",
    },
    "r_outlet": {
        "type": "float", "default": 10.0, "min": 8.0, "max": 20.0,
        "group": "Radii", "description": "Outlet (descending) radius [mm]",
    },
    "arch_width_mm": {
        "type": "float", "default": 90.0, "min": 30.0, "max": 120.0,
        "group": "Arch", "description": "Arch width (ascending→descending horizontal extent) [mm]",
    },
    "arch_height_mm": {
        "type": "float", "default": 45.0, "min": 20.0, "max": 60.0,
        "group": "Arch", "description": "Arch peak height above ascending top [mm]",
    },
    "torsion_deg": {
        "type": "float", "default": 0.0, "min": -30.0, "max": 30.0,
        "group": "Arch",
        "description": "RIGID arch tilt around inlet z-axis [deg] (arch stays planar)",
    },
    "twist_deg": {
        "type": "float", "default": 0.0, "min": -45.0, "max": 45.0,
        "group": "Arch",
        "description": "GRADUAL twist along arch [deg] (arch becomes non-planar 3D curve)",
    },
    "arch_shape": {
        "type": "str", "default": "circle", "min": None, "max": None,
        "group": "Arch",
        "description": "'circle' = constrained to H ≤ W ≤ 2H; 'ellipse' = independent W + H",
        "choices": ["circle", "ellipse"],
    },
    "arch_radius_mm": {
        "type": "float", "default": 0.0, "min": 0.0, "max": 22.0,
        "group": "Radii",
        "description": "TUBE cross-section radius at the ARCH segment [mm] (between "
                       "r_inlet and r_outlet — middle lumen radius). Default 0 = "
                       "auto-derive as midpoint (r_inlet + r_outlet) / 2. Set > 0 "
                       "to dial it independently. NOT to be confused with arch_R_c_mm "
                       "(which is centerline curvature).",
    },
    "taper_mode": {
        "type": "str", "default": "smoothstep", "min": None, "max": None,
        "group": "Radii",
        "description": "How the lumen radius transitions between r_inlet, arch_radius_mm, "
                       "and r_outlet across segment boundaries: 'smoothstep' "
                       "(cosine-Hermite blend — default, smoothest), 'linear' "
                       "(piecewise-linear blend across 30 mm window), 'piecewise' "
                       "(constant radius per segment, hard step at boundaries).",
        "choices": ["piecewise", "linear", "smoothstep"],
    },
    "arch_R_c_mm": {
        "type": "float", "default": 0.0, "min": 0.0, "max": 100.0,
        "group": "Arch",
        "description": "CENTERLINE curvature radius shortcut [mm]. When > 0, derives "
                       "arch_width_mm and arch_height_mm automatically. Circle mode: "
                       "sets W=2·R, H=R (canonical U-arch). Ellipse mode: treated as "
                       "R_peak (curvature at top); combine with W or H to derive the "
                       "third (alone gives degenerate ellipse = circle). NOT the tube "
                       "radius — that's arch_radius_mm.",
    },
    "ascending_length": {
        "type": "float", "default": 50.0, "min": 40.0, "max": 90.0,
        "group": "Lengths (optional)",
        "description": "Straight ascending segment length before arch [mm]",
    },
    "descending_length": {
        "type": "float", "default": 200.0, "min": 60.0, "max": 300.0,
        "group": "Lengths (optional)",
        "description": "Straight descending segment length after arch [mm]",
    },
}

# v3 name → v2 name
V3_TO_V2: dict[str, str] = {
    "r_inlet": "r_ascending",
    "r_outlet": "r_descending",
    "arch_radius_mm": "r_arch",       # tube radius at arch segment (was midpoint-auto)
    "arch_width_mm": "arch_span_mm",
    "arch_height_mm": "arch_height_mm",
    "taper_mode": "taper_mode",       # lumen taper between segment radii
    "torsion_deg": "arch_tilt_deg",
    "twist_deg": "arch_twist_deg",
    "arch_shape": "arch_shape",
    "ascending_length": "ascending_length",
    "descending_length": "descending_length",
}

# v2 params hard-wired (NOT exposed in v3) — workshop-quality defaults
# taper_mode is now a v3 knob (was here) — when v3 user doesn't set it,
# blender_aorta_v2.py's --taper_mode default "smoothstep" applies.
V2_FIXED: dict[str, Any] = {
    "delta_3": 0.0,
    "delta_4": 0.0,
    "junction_blend_mm": 12.0,
    "segments_radial": 96,
    "curve_samples": 300,
}


# ---------------------------------------------------------------------------
# v3 → v2 translation
# ---------------------------------------------------------------------------


def translate_v3_to_v2(v3_params: dict[str, Any]) -> dict[str, Any]:
    """Map v3 parameter names to v2 names + inject fixed defaults.

    Two preprocessing steps:

      1. arch_R_c_mm (CURVATURE shortcut): when > 0, expands to
         arch_width_mm + arch_height_mm based on arch_shape:
           circle  → W = 2·R, H = R  (canonical U-arch)
           ellipse → R is R_peak; combine with W or H to derive the third
                     (alone → degenerate U-arch a=b=R).
         Mutually exclusive with W+H in circle mode; in ellipse mode the
         R+W+H combination is over-determined.

      2. arch_radius_mm (TUBE radius at arch segment): when > 0, sets
         v2.r_arch directly. When 0 (default), auto-derived as the midpoint
         of r_inlet and r_outlet so the main lumen tapers smoothly.

    Calls ``cli_v2._resolve_arch_params`` so v2 sees ``arch_R_c`` and
    ``arch_angle_deg`` (not ``arch_span_mm`` / ``arch_height_mm``) when in
    circle mode.
    """
    v3_params = dict(v3_params)  # copy so we can mutate

    # ── arch_R_c_mm curvature shortcut (NOT in V3_TO_V2 — handled here) ────
    R_c_shortcut = float(v3_params.pop("arch_R_c_mm", 0.0))
    if R_c_shortcut > 0:
        shape = v3_params.get("arch_shape", "circle")
        has_w = "arch_width_mm" in v3_params
        has_h = "arch_height_mm" in v3_params

        if shape == "circle":
            if has_w or has_h:
                raise ValueError(
                    "arch_R_c_mm cannot be combined with arch_width_mm or "
                    "arch_height_mm in circle mode — they would conflict "
                    "(R_c uniquely determines both W and H of the U-arch)."
                )
            v3_params["arch_width_mm"] = 2.0 * R_c_shortcut
            v3_params["arch_height_mm"] = R_c_shortcut
        else:  # ellipse — arch_R_c_mm = R_peak (curvature at top)
            if has_w and has_h:
                raise ValueError(
                    "arch_R_c_mm + arch_width_mm + arch_height_mm is "
                    "over-determined in ellipse mode (R_peak is fully "
                    "determined by W and H alone — drop one of them)."
                )
            if has_w:
                W = float(v3_params["arch_width_mm"])
                v3_params["arch_height_mm"] = (W * W) / (4.0 * R_c_shortcut)
            elif has_h:
                H = float(v3_params["arch_height_mm"])
                v3_params["arch_width_mm"] = 2.0 * math.sqrt(R_c_shortcut * H)
            else:
                v3_params["arch_width_mm"] = 2.0 * R_c_shortcut
                v3_params["arch_height_mm"] = R_c_shortcut

    # ── Standard v3 → v2 name translation ──────────────────────────────────
    v2: dict[str, Any] = dict(V2_FIXED)
    for v3_key, value in v3_params.items():
        if v3_key not in V3_TO_V2:
            raise ValueError(f"Unknown v3 parameter {v3_key!r}")
        v2[V3_TO_V2[v3_key]] = value

    # Auto-derive r_arch (tube radius at arch segment) when not set explicitly.
    r_in = float(v2.get("r_ascending", PARAMETERS["r_inlet"]["default"]))
    r_out = float(v2.get("r_descending", PARAMETERS["r_outlet"]["default"]))
    if float(v2.get("r_arch", 0.0)) <= 0.0:
        v2["r_arch"] = 0.5 * (r_in + r_out)
    return _resolve_arch_params(v2)


# ---------------------------------------------------------------------------
# Discovery + validation
# ---------------------------------------------------------------------------


def _format_params_table(markdown: bool = False) -> str:
    groups: dict[str, list[tuple[str, dict[str, Any]]]] = {}
    for name, info in PARAMETERS.items():
        groups.setdefault(info["group"], []).append((name, info))
    if markdown:
        lines = [
            "# Aorta geometry generator v3 — parameter reference",
            "",
            "Five-knob minimal interface. Translates to v2 internally; everything "
            "not listed here is fixed at the v2 workshop defaults "
            f"(taper={V2_FIXED['taper_mode']}, junction_blend_mm={V2_FIXED['junction_blend_mm']}, "
            f"segments_radial={V2_FIXED['segments_radial']}, curve_samples={V2_FIXED['curve_samples']}, "
            f"non-planar Fourier OFF).",
            "",
        ]
        for group, entries in groups.items():
            lines += [f"## {group}", ""]
            lines += [
                "| Parameter | Type | Default | Workshop range | Description |",
                "|---|---|---|---|---|",
            ]
            for name, info in entries:
                rng = f"{info['min']}–{info['max']}"
                lines.append(
                    f"| `{name}` | {info['type']} | `{info['default']}` | {rng} | {info['description']} |"
                )
            lines.append("")
        return "\n".join(lines)

    lines = [f"v3 minimal parameters ({len(PARAMETERS)}):", ""]
    for group, entries in groups.items():
        lines.append(f"  {group}")
        for name, info in entries:
            lines.append(
                f"    {name:18s} {info['type']:6s} default={info['default']:<8} "
                f"range={info['min']}–{info['max']:<6}  {info['description']}"
            )
        lines.append("")
    lines += [
        f"v2 fixed defaults (not exposed in v3): {V2_FIXED}",
        "Use cli_v2.py for finer control.",
    ]
    return "\n".join(lines)


def _check_param_name(name: str, source: str) -> None:
    if name in PARAMETERS:
        return
    suggestions = _suggest(name, list(PARAMETERS.keys()))
    hint = f"  Did you mean {' or '.join(repr(s) for s in suggestions)}?" if suggestions else ""
    raise ValueError(
        f"Unknown v3 parameter {name!r} in {source}.{hint}\n"
        f"  Run `python cli_v3.py --list-params`."
    )


def validate_spec(payload: dict[str, Any], *, source: str = "spec") -> None:
    if not isinstance(payload, dict):
        raise ValueError(f"Spec must be a JSON object, got {type(payload).__name__}")
    mode = payload.get("mode")
    if mode not in {"single", "sweep"}:
        raise ValueError(
            f"{source}.mode must be 'single' or 'sweep' (v3 keeps things minimal). "
            f"Got {mode!r}. Use cli_v2.py for sample/grid modes."
        )
    geom = payload.get("geometry")
    if geom is not None and geom != "healthy_arch_v3":
        raise ValueError(
            f"{source}.geometry must be 'healthy_arch_v3', got {geom!r}."
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
                raise ValueError(f"{source}.sweep missing key {key!r}")
        _check_param_name(sw["param"], f"{source}.sweep.param")
        if float(sw["low"]) >= float(sw["high"]):
            raise ValueError(f"{source}.sweep.low must be < high")
        if int(sw["n"]) < 2:
            raise ValueError(f"{source}.sweep.n must be ≥ 2")


def _parse_param_override(s: str) -> tuple[str, Any]:
    if "=" not in s:
        raise ValueError(f"--param expects key=value, got {s!r}")
    key, _, value = s.partition("=")
    key, value = key.strip(), value.strip()
    _check_param_name(key, "--param")
    t = PARAMETERS[key]["type"]
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
# Spec expansion (single / sweep only)
# ---------------------------------------------------------------------------


def load_spec(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    validate_spec(payload, source=f"{path.name}")
    return payload


def expand_cases(spec: dict[str, Any]) -> list[dict[str, Any]]:
    mode = spec["mode"]
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
    return cases


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="blender_cli_v3", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spec", type=Path, help="Path to a v3 JSON spec.")
    parser.add_argument("--output", "-o", type=Path, help="Output dir for case folders.")
    parser.add_argument("--blender", default=os.environ.get("BLENDER", "blender"),
                        help="Path to Blender executable.")
    parser.add_argument("--generator", type=Path, default=HERE / "blender_aorta_v2.py",
                        help="Path to blender_aorta_v2.py (the v2 generator is reused for v3).")
    parser.add_argument("--param", action="append", default=[], metavar="KEY=VALUE",
                        help="Override a v3 parameter (repeatable).")
    parser.add_argument("--save-blend", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--list-params", action="store_true",
                        help="Print the 5+2 v3 parameters and exit.")
    parser.add_argument("--markdown", action="store_true")
    parser.add_argument("--yes", "-y", action="store_true")
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
    print(f"Output     : {args.output}")
    print(f"Generator  : {args.generator}")
    print(f"Blender    : {args.blender}")

    if not args.dry_run and not args.yes:
        warn_if_expensive(len(cases))

    if args.dry_run:
        print("\n[dry-run] v3 cases:")
        for c in cases:
            v3_params = c["params"]
            v2_params = translate_v3_to_v2(v3_params)
            print(f"  {c['case_id']}: v3 = {v3_params}")
            print(f"           → v2 = {v2_params}")
        return 0

    rows: list[dict[str, Any]] = []
    for i, c in enumerate(cases, start=1):
        case_dir = args.output / c["case_id"]
        print(f"\n[{i}/{len(cases)}] {c['case_id']}")
        try:
            v3_params = c["params"]
            v2_params = translate_v3_to_v2(v3_params)
            run_blender_and_split(
                blender_path=args.blender,
                generator=args.generator,
                params=v2_params,
                case_dir=case_dir,
                save_blend=args.save_blend,
            )
            # Record provenance in v3-flavoured terms but keep v2 details for traceability.
            write_geometry_meta(
                case_dir=case_dir, spec=spec, params={**v3_params,
                                                        "_translated_to_v2": v2_params},
                mode=c["mode"], sampler=c["sampler"],
                sample_index=c["sample_index"], seed=c["seed"],
                generator_path=args.generator,
            )
            rows.append({
                "case_id": c["case_id"], "status": "ok",
                "mode": c["mode"], **v3_params,
            })
            print(f"  OK -> {case_dir}")
        except Exception as e:
            rows.append({"case_id": c["case_id"], "status": "failed",
                         "error": str(e), **c["params"]})
            print(f"  FAIL: {e}")

    manifest_path = write_sweep_manifest(args.output, rows)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"\nSweep manifest -> {manifest_path}")
    print(f"Done. {n_ok}/{len(rows)} succeeded.")
    return 0 if n_ok == len(rows) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as e:
        sys.stderr.write(f"\nError: {e}\n")
        raise SystemExit(2)
