"""Split a watertight aortic arch STL into wall + cap patches (Block A helper).

Ported from ``~/GitHub/AortaCFD-X/scripts/split_patches.py`` with two
changes for AortaCFD-app compatibility:

  1. Patch names follow the AortaCFD-app convention:
       wall_aorta.stl, inlet.stl, outlet1.stl, outlet2.stl, ...
     (instead of wall / outlet_desc / outlet_branch1 / outlet_branch2 / ...)
  2. ``split_stl(stl_path, output_dir=None)`` is importable as a Python
     function, so the orchestrator (``cli.py``) does not need to spawn
     a subprocess.

Requires: numpy, numpy-stl, scipy.

Standalone CLI is preserved::

    python split_patches.py /path/to/geometry.stl
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from stl import mesh as stl_mesh


def build_adjacency(vectors):
    edge_to_tri = defaultdict(list)
    for i, tri in enumerate(vectors):
        for j in range(3):
            p1 = tuple(tri[j])
            p2 = tuple(tri[(j + 1) % 3])
            edge = tuple(sorted([p1, p2]))
            edge_to_tri[edge].append(i)
    adj = defaultdict(set)
    for tris in edge_to_tri.values():
        if len(tris) == 2:
            adj[tris[0]].add(tris[1])
            adj[tris[1]].add(tris[0])
    return adj


def unit_normals(normals):
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    return normals / np.clip(lengths, 1e-12, None)


def flood_fill_cap(seed_idx, adj, unit_norms, centroids, ref_normal, center, max_radius, normal_tol=0.85):
    cap: set[int] = set()
    queue = [seed_idx]
    while queue:
        tri = queue.pop()
        if tri in cap:
            continue
        dot = abs(np.dot(unit_norms[tri], ref_normal))
        if dot < normal_tol:
            continue
        if np.linalg.norm(centroids[tri] - center) > max_radius:
            continue
        cap.add(tri)
        for nb in adj[tri]:
            if nb not in cap:
                queue.append(nb)
    return sorted(cap)


def find_cap(centroids, unit_norms, adj, center, expected_normal, search_radius, normal_tol=0.85):
    dists = np.linalg.norm(centroids - center, axis=1)
    candidates = np.where(dists < search_radius)[0]
    if len(candidates) == 0:
        raise RuntimeError(f"No triangles within {search_radius} mm of {center}")
    dots = np.abs(np.sum(unit_norms[candidates] * expected_normal, axis=1))
    good = candidates[dots > normal_tol]
    if len(good) == 0:
        good = candidates[dots > 0.7]
    if len(good) == 0:
        raise RuntimeError(
            f"No triangle near {center} with normal aligned to {expected_normal} "
            f"(best dot = {dots.max():.3f})"
        )
    seed = good[np.argmin(dists[good])]
    cap_tris = flood_fill_cap(
        seed,
        adj,
        unit_norms,
        centroids,
        ref_normal=unit_norms[seed],
        center=center,
        max_radius=search_radius * 1.5,
        normal_tol=normal_tol,
    )
    return cap_tris


def extract_submesh(full_mesh, indices):
    sub = stl_mesh.Mesh(np.zeros(len(indices), dtype=stl_mesh.Mesh.dtype))
    for i, idx in enumerate(indices):
        sub.vectors[i] = full_mesh.vectors[idx]
        sub.normals[i] = full_mesh.normals[idx]
    return sub


def split_stl(stl_path: Path | str, output_dir: Path | str | None = None) -> dict[str, Path]:
    """Split a watertight aortic STL into named AortaCFD-app patches.

    Reads ``stl_path`` and its sidecar ``<stl_path>.json`` (or
    ``stl_path.with_suffix('.json')``), identifies the cap regions via
    flood-fill, and writes wall + cap STLs into ``output_dir``
    (defaults to the STL's parent directory).

    Returns
    -------
    Dict mapping logical patch name (``"inlet"``, ``"outlet1"``,
    ``"wall_aorta"``, …) to written STL path.
    """
    stl_path = Path(stl_path)
    json_path = stl_path.with_suffix(".json")
    if not stl_path.exists():
        raise FileNotFoundError(f"STL not found: {stl_path}")
    if not json_path.exists():
        raise FileNotFoundError(f"Metadata JSON not found: {json_path}")

    out_dir = Path(output_dir) if output_dir is not None else stl_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    m = stl_mesh.Mesh.from_file(str(stl_path))
    meta: dict[str, Any] = json.loads(json_path.read_text())
    n_tri = len(m.vectors)
    centroids = m.vectors.mean(axis=1)
    u_normals = unit_normals(m.normals)
    adj = build_adjacency(m.vectors)

    main_diam = meta["main_diameter_mm"]
    main_radius = main_diam / 2.0
    branch_diam = meta["branch_diameter_mm"]
    branch_radius = branch_diam / 2.0
    asc_length = meta["ascending_length_mm"]  # noqa: F841 (kept for future use)
    desc_length = meta["descending_length_mm"]  # noqa: F841
    arch_span = meta["arch_span_mm"]

    inlet_center = np.array([0.0, 0.0, 0.0])
    inlet_normal = np.array([0.0, 0.0, -1.0])

    desc_x_mask = np.abs(centroids[:, 0] - arch_span) < main_diam
    outlet_z = centroids[desc_x_mask, 2].min() if desc_x_mask.sum() > 0 else centroids[:, 2].min()
    outlet_center = np.array([arch_span, 0.0, outlet_z])
    outlet_normal = np.array([0.0, 0.0, -1.0])

    cap_specs = [
        # (aortacfd_name, center, normal, radius)
        ("inlet", inlet_center, inlet_normal, main_radius + 2),
        ("outlet1", outlet_center, outlet_normal, main_radius + 2),  # descending = outlet1
    ]
    branch_specs = list(meta.get("branches") or [])
    for k, br in enumerate(branch_specs, start=2):
        origin = np.array(br["origin_xyz_mm"])
        end = np.array(br["end_xyz_mm"])
        direction = end - origin
        direction /= np.linalg.norm(direction)
        cap_specs.append((f"outlet{k}", end, direction, branch_radius + 3))

    all_cap_indices: set[int] = set()
    cap_results: dict[str, list[int]] = {}
    for name, center, normal, radius in cap_specs:
        indices = find_cap(centroids, u_normals, adj, center=center, expected_normal=normal, search_radius=radius)
        cap_results[name] = indices
        all_cap_indices.update(indices)

    wall_indices = sorted(set(range(n_tri)) - all_cap_indices)

    paths: dict[str, Path] = {}
    wall_path = out_dir / "wall_aorta.stl"
    extract_submesh(m, wall_indices).save(str(wall_path))
    paths["wall_aorta"] = wall_path

    for name, indices in cap_results.items():
        p = out_dir / f"{name}.stl"
        extract_submesh(m, indices).save(str(p))
        paths[name] = p

    return paths


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 1
    stl_path = Path(argv[0])
    out_dir = Path(argv[1]) if len(argv) > 1 else None
    paths = split_stl(stl_path, out_dir)
    for name, p in paths.items():
        print(f"{name:14s} -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
