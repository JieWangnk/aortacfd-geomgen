#!/usr/bin/env python3
"""Compare two STL outputs of the same vessel geometry from different tools.

Designed for head-to-head benchmarking of our generators against VMTK,
SimVascular, or any other tool that emits an STL.  Reports a fixed metric set
that maps cleanly to peer-reviewer questions:

    - lumen volume                  (mm^3)
    - surface area                  (mm^2)
    - bounding box                  (mm)
    - n_components / n_holes        (topological health)
    - watertight / manifold flags
    - triangle count + aspect ratios (mesh quality)
    - Hausdorff distance ours <-> reference (worst-case surface deviation)
    - mean nearest-point distance   (average surface deviation)
    - centreline RMS distance       (after symmetric centreline extraction)

Outputs a single JSON report and a compact human-readable table.

Usage
-----
    python compare.py ours/bifurcation_murray.stl vmtk/bifurcation_murray.stl
    python compare.py --json out.json A.stl B.stl
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import trimesh


# ---------------------------------------------------------------------------
# Mesh loading and topology
# ---------------------------------------------------------------------------

def load_mesh(stl_path: Path) -> trimesh.Trimesh:
    """Load a (possibly multi-solid) STL and return one concatenated mesh."""
    obj = trimesh.load(str(stl_path))
    if isinstance(obj, trimesh.Scene):
        meshes = list(obj.geometry.values())
        if not meshes:
            raise ValueError(f"{stl_path} has no geometry")
        m = trimesh.util.concatenate(meshes)
    else:
        m = obj
    m.merge_vertices()
    return m


def topology_metrics(m: trimesh.Trimesh) -> dict:
    edges = m.edges_sorted
    _, c = np.unique(edges, axis=0, return_counts=True)
    bbox = m.bounds[1] - m.bounds[0]
    return {
        "n_vertices": int(len(m.vertices)),
        "n_faces": int(len(m.faces)),
        "n_edges_unique": int(len(c)),
        "n_boundary_edges": int((c == 1).sum()),
        "n_nonmanifold_edges": int((c >= 3).sum()),
        "n_components": int(m.body_count),
        "is_watertight": bool(m.is_watertight),
        "is_winding_consistent": bool(m.is_winding_consistent),
        "lumen_volume_mm3": float(m.volume),
        "surface_area_mm2": float(m.area),
        "bbox_mm": [float(x) for x in bbox],
        "centroid_mm": [float(x) for x in m.centroid],
    }


# ---------------------------------------------------------------------------
# Mesh quality
# ---------------------------------------------------------------------------

def triangle_aspect_ratios(m: trimesh.Trimesh) -> np.ndarray:
    """Return per-triangle aspect ratio = longest_edge / shortest_edge."""
    v = m.vertices[m.faces]
    e0 = np.linalg.norm(v[:, 1] - v[:, 0], axis=1)
    e1 = np.linalg.norm(v[:, 2] - v[:, 1], axis=1)
    e2 = np.linalg.norm(v[:, 0] - v[:, 2], axis=1)
    edges = np.stack([e0, e1, e2], axis=1)
    edges = np.maximum(edges, 1e-12)
    return edges.max(axis=1) / edges.min(axis=1)


def triangle_areas(m: trimesh.Trimesh) -> np.ndarray:
    return m.area_faces


def quality_metrics(m: trimesh.Trimesh) -> dict:
    ar = triangle_aspect_ratios(m)
    areas = triangle_areas(m)
    return {
        "aspect_ratio_min": float(ar.min()),
        "aspect_ratio_mean": float(ar.mean()),
        "aspect_ratio_p95": float(np.percentile(ar, 95)),
        "aspect_ratio_max": float(ar.max()),
        "n_slivers_AR_gt_5": int((ar > 5.0).sum()),
        "n_slivers_AR_gt_20": int((ar > 20.0).sum()),
        "area_min_mm2": float(areas.min()),
        "area_mean_mm2": float(areas.mean()),
        "area_max_mm2": float(areas.max()),
        "area_total_mm2": float(areas.sum()),
    }


# ---------------------------------------------------------------------------
# Surface comparison (Hausdorff and mean nearest distance)
# ---------------------------------------------------------------------------

def surface_distances(a: trimesh.Trimesh, b: trimesh.Trimesh,
                      n_samples: int = 50_000) -> dict:
    """Sample-based Hausdorff and mean nearest-point distances.

    Distance is symmetric: max(d(A->B), d(B->A)) for Hausdorff,
    0.5*(mean(d(A->B)) + mean(d(B->A))) for the mean.
    """
    rng = np.random.default_rng(seed=42)
    pts_a = trimesh.sample.sample_surface(a, n_samples)[0]
    pts_b = trimesh.sample.sample_surface(b, n_samples)[0]
    # Closest point on each mesh
    _, dist_a_to_b, _ = trimesh.proximity.closest_point(b, pts_a)
    _, dist_b_to_a, _ = trimesh.proximity.closest_point(a, pts_b)
    return {
        "hausdorff_mm": float(max(dist_a_to_b.max(), dist_b_to_a.max())),
        "mean_nearest_mm": float(0.5 * (dist_a_to_b.mean() + dist_b_to_a.mean())),
        "p95_a_to_b_mm": float(np.percentile(dist_a_to_b, 95)),
        "p95_b_to_a_mm": float(np.percentile(dist_b_to_a, 95)),
        "n_samples_per_mesh": int(n_samples),
    }


# ---------------------------------------------------------------------------
# Centreline extraction (lightweight: principal-axis sweep + per-slice centroid)
# ---------------------------------------------------------------------------

def extract_centreline(m: trimesh.Trimesh, axis: np.ndarray = None,
                       n_slices: int = 50) -> np.ndarray:
    """Extract a centreline by slicing along a given axis and taking centroids.

    Returns (k, 3) array.  ``axis`` defaults to the mesh's longest principal
    axis.  Handles single-tube geometry well; for branched geometry (Y or
    aorta) the slices through bifurcation points lump multiple lumens
    together — the user should pass an axis aligned with the main vessel.
    """
    if axis is None:
        # Use the longest bounding-box dimension as a coarse axis
        bbox = m.bounds[1] - m.bounds[0]
        axis = np.zeros(3); axis[int(np.argmax(bbox))] = 1.0
    axis = axis / np.linalg.norm(axis)
    proj = m.vertices @ axis
    s_min, s_max = proj.min(), proj.max()
    centres = []
    for u in np.linspace(s_min, s_max, n_slices):
        # Soft slice: vertices within +/- step of plane u
        step = (s_max - s_min) / n_slices
        mask = np.abs(proj - u) < step
        if mask.sum() < 6:
            continue
        centres.append(m.vertices[mask].mean(axis=0))
    return np.asarray(centres)


def centreline_rms(a: trimesh.Trimesh, b: trimesh.Trimesh,
                   axis: np.ndarray = None, n_slices: int = 50) -> dict:
    cl_a = extract_centreline(a, axis=axis, n_slices=n_slices)
    cl_b = extract_centreline(b, axis=axis, n_slices=n_slices)
    if len(cl_a) < 3 or len(cl_b) < 3:
        return {"centreline_rms_mm": None, "note": "insufficient slices"}
    # Match by arclength fraction along the axis
    n = min(len(cl_a), len(cl_b))
    cl_a = cl_a[np.linspace(0, len(cl_a) - 1, n).astype(int)]
    cl_b = cl_b[np.linspace(0, len(cl_b) - 1, n).astype(int)]
    diffs = np.linalg.norm(cl_a - cl_b, axis=1)
    return {
        "centreline_rms_mm": float(np.sqrt(np.mean(diffs ** 2))),
        "centreline_max_mm": float(diffs.max()),
        "n_slices_used": int(n),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def compare(stl_a: Path, stl_b: Path, label_a: str, label_b: str) -> dict:
    a = load_mesh(stl_a)
    b = load_mesh(stl_b)

    report = {
        "inputs": {label_a: str(stl_a), label_b: str(stl_b)},
        f"{label_a}_topology": topology_metrics(a),
        f"{label_b}_topology": topology_metrics(b),
        f"{label_a}_quality": quality_metrics(a),
        f"{label_b}_quality": quality_metrics(b),
        "surface_distance": surface_distances(a, b),
        "centreline": centreline_rms(a, b),
    }

    # Volume / area deltas
    va, vb = report[f"{label_a}_topology"]["lumen_volume_mm3"], report[f"{label_b}_topology"]["lumen_volume_mm3"]
    aa, ab = report[f"{label_a}_topology"]["surface_area_mm2"], report[f"{label_b}_topology"]["surface_area_mm2"]
    report["deltas"] = {
        "volume_diff_mm3":     float(vb - va),
        "volume_diff_pct":     float(100 * (vb - va) / va) if va > 0 else None,
        "area_diff_mm2":       float(ab - aa),
        "area_diff_pct":       float(100 * (ab - aa) / aa) if aa > 0 else None,
    }

    return report


def print_table(report: dict, label_a: str, label_b: str) -> None:
    print()
    print(f"  {'metric':<32s}  {label_a:>16s}  {label_b:>16s}  {'delta':>12s}")
    print("  " + "-" * 80)
    ta = report[f"{label_a}_topology"]
    tb = report[f"{label_b}_topology"]
    qa = report[f"{label_a}_quality"]
    qb = report[f"{label_b}_quality"]
    rows = [
        ("lumen volume (mm^3)",      ta["lumen_volume_mm3"], tb["lumen_volume_mm3"]),
        ("surface area (mm^2)",      ta["surface_area_mm2"], tb["surface_area_mm2"]),
        ("vertex count",             ta["n_vertices"],       tb["n_vertices"]),
        ("face count",               ta["n_faces"],          tb["n_faces"]),
        ("boundary edges",           ta["n_boundary_edges"], tb["n_boundary_edges"]),
        ("non-manifold edges",       ta["n_nonmanifold_edges"], tb["n_nonmanifold_edges"]),
        ("aspect ratio mean",        qa["aspect_ratio_mean"], qb["aspect_ratio_mean"]),
        ("aspect ratio p95",         qa["aspect_ratio_p95"], qb["aspect_ratio_p95"]),
        ("slivers (AR > 5)",         qa["n_slivers_AR_gt_5"], qb["n_slivers_AR_gt_5"]),
        ("min triangle area (mm^2)", qa["area_min_mm2"],     qb["area_min_mm2"]),
    ]
    for name, va, vb in rows:
        delta = vb - va if isinstance(va, (int, float)) else None
        if isinstance(va, float):
            print(f"  {name:<32s}  {va:>16.4g}  {vb:>16.4g}  {delta:>12.4g}")
        else:
            print(f"  {name:<32s}  {va:>16d}  {vb:>16d}  {delta:>12d}")
    print()
    sd = report["surface_distance"]
    print(f"  Surface comparison ({sd['n_samples_per_mesh']} samples / mesh)")
    print(f"    Hausdorff distance:   {sd['hausdorff_mm']:.4g} mm")
    print(f"    mean nearest:         {sd['mean_nearest_mm']:.4g} mm")
    print(f"    p95(A->B):            {sd['p95_a_to_b_mm']:.4g} mm")
    print(f"    p95(B->A):            {sd['p95_b_to_a_mm']:.4g} mm")
    cl = report["centreline"]
    if cl.get("centreline_rms_mm") is not None:
        print(f"    centreline RMS:       {cl['centreline_rms_mm']:.4g} mm "
              f"({cl['n_slices_used']} slices)")
    print()
    d = report["deltas"]
    print(f"  Volume delta: {d['volume_diff_mm3']:+.2f} mm^3 ({d['volume_diff_pct']:+.2f}%)")
    print(f"  Area delta:   {d['area_diff_mm2']:+.2f} mm^2 ({d['area_diff_pct']:+.2f}%)")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("stl_a", type=Path)
    parser.add_argument("stl_b", type=Path)
    parser.add_argument("--label-a", default="A")
    parser.add_argument("--label-b", default="B")
    parser.add_argument("--json", type=Path, default=None,
                        help="optional output JSON path")
    args = parser.parse_args()

    if not args.stl_a.exists() or not args.stl_b.exists():
        print(f"Missing input(s): {args.stl_a} or {args.stl_b}")
        sys.exit(2)

    report = compare(args.stl_a, args.stl_b, args.label_a, args.label_b)
    print_table(report, args.label_a, args.label_b)
    if args.json is not None:
        args.json.write_text(json.dumps(report, indent=2))
        print(f"  Wrote JSON report: {args.json}")


if __name__ == "__main__":
    main()
