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
# Spec parsing
# ---------------------------------------------------------------------------


def load_spec(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    mode = payload.get("mode")
    if mode not in {"single", "sweep", "sample"}:
        raise ValueError(f"spec.mode must be one of single/sweep/sample, got {mode!r}")
    return payload


def expand_cases(spec: dict[str, Any]) -> list[dict[str, Any]]:
    """Return a list of per-case dicts: ``{case_id, params, mode, sampler, sample_index, seed}``."""
    mode = spec.get("mode")
    if mode not in {"single", "sweep", "sample"}:
        raise ValueError(f"spec.mode must be one of single/sweep/sample, got {mode!r}")
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
    parser.add_argument("--spec", type=Path, required=True, help="Path to a JSON spec defining mode + params.")
    parser.add_argument("--output", "-o", type=Path, required=True, help="Output dir for case folders.")
    parser.add_argument("--blender", default=os.environ.get("BLENDER", "blender"),
                        help="Path to Blender executable (default: env $BLENDER or `blender`).")
    parser.add_argument("--generator", type=Path, default=HERE / "blender_aorta_like_generator.py",
                        help="Path to blender_aorta_like_generator.py (default: sibling of this script).")
    parser.add_argument("--save-blend", action="store_true", help="Save .blend alongside STL (debugging).")
    parser.add_argument("--dry-run", action="store_true", help="Print the cases that would be generated without running Blender.")
    parser.add_argument("--limit", type=int, default=None, help="Stop after N cases (smoke testing).")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s %(message)s")

    spec = load_spec(args.spec)
    cases = expand_cases(spec)
    if args.limit is not None:
        cases = cases[: args.limit]

    args.output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.spec, args.output / args.spec.name)
    geometry = spec.get("geometry", "arch_branched")

    print(f"Spec       : {args.spec}")
    print(f"Mode       : {spec['mode']}")
    print(f"Cases      : {len(cases)}")
    print(f"Output     : {args.output}")
    print(f"Generator  : {args.generator}")
    print(f"Blender    : {args.blender}")

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
    raise SystemExit(main())
