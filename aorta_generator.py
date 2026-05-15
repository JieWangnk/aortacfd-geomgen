#!/usr/bin/env python3
"""Pure-Python watertight aortic arch generator.

Mirrors the bifurcation-union approach from ``vessel_generator.py``: build
the arch + three supra-aortic branches as capped tubes, run a boolean
union via ``trimesh`` + ``manifold3d``, then relabel patches by cap-plane
proximity.  Produces a single watertight STL with named solids
``wall`` / ``inlet`` / ``descending_outlet`` /
``BCA_outlet`` / ``LCCA_outlet`` / ``LSA_outlet`` -- directly consumable by
snappyHexMesh.

Requires ``pip install trimesh manifold3d``.

Compared with ``blender_aorta_like_generator.py``, this version:
- Has no Blender dependency (runs anywhere Python does).
- Uses RMF tube construction from ``vessel_generator.py``.
- Does a genuine boolean union (not a surface overlap).

Run
---
    python aorta_generator.py --demo --out /tmp/aorta/
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, field
from typing import List

import numpy as np

from vessel_generator import (
    build_tube,
    write_multisolid_stl,
    _trimesh_from_arrays,
    _label_cap_faces,
)

EPS = 1e-9


# =============================================================================
# Parameter dataclass (matches the 16-parameter browser model)
# =============================================================================

@dataclass
class Aorta:
    # Centreline (4)
    ascending_length: float = 45.0         # mm
    arch_span: float = 70.0                # mm
    arch_height: float = 35.0              # mm
    descending_length: float = 80.0        # mm

    # Diameter (3)
    d_ascending: float = 30.0              # mm
    taper_ratio: float = 0.73              # d_descending / d_ascending
    taper_exponent: float = 1.0

    # Branches (4)
    bca_ratio: float = 0.45                # d_BCA / d_ascending
    lcca_ratio: float = 0.23
    lsa_ratio: float = 0.33
    branch_spacing: float = 13.0           # mm (arc-length spacing between origins)

    # Coarctation (4)
    coarctation_area_reduction: float = 0.0
    coarctation_offset_from_lsa: float = 10.0   # mm distal from LSA origin
    coarctation_shape: float = 0.5               # 0 = shelf, 1 = hourglass
    coarctation_length: float = 20.0             # mm

    # Hypoplasia (1)
    proximal_hypoplasia: float = 0.0

    # Branch geometry (fixed shape parameters, not part of the 16)
    branch_length: float = 40.0            # mm
    branch_tilt_deg: float = 18.0          # tilt from arch tangent toward +z
    branch_splay_deg: float = 15.0         # lateral splay amplitude (y-axis)


# =============================================================================
# Centreline construction (XZ plane, y = 0)
# =============================================================================

def _cubic_bezier(p0, p1, p2, p3, t):
    s = 1.0 - t
    return (s ** 3) * p0 + 3 * (s ** 2) * t * p1 + 3 * s * (t ** 2) * p2 + (t ** 3) * p3


def _build_centreline(a: Aorta, n: int = 200) -> np.ndarray:
    """Piecewise centreline: ascending line -> arch crown bezier ->
    arch-to-descending bezier -> descending line.
    """
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([0.0, 0.0, a.ascending_length])

    p2 = np.array([a.arch_span, 0.0, a.ascending_length + 0.20 * a.arch_height])
    # G^1 at p1: h1 lies on the line through p0 -> p1 (directly above p1).
    h1 = np.array([0.0, 0.0, a.ascending_length + a.arch_height])
    # G^1 at p2: h2 directly above p2, so the arch crown arrives at p2 along -z.
    h2 = np.array([a.arch_span, 0.0, a.ascending_length + a.arch_height])

    transition_drop = min(0.30 * a.descending_length, 40.0)
    p3 = np.array([a.arch_span, 0.0, a.ascending_length - transition_drop])
    # G^1 at p2: h3 directly below p2 along -z, parallel to (p2 - h2).
    h3 = np.array([a.arch_span, 0.0, a.ascending_length + 0.10 * a.arch_height])
    # G^1 at p3: h4 directly above p3 along -z, parallel to (p4 - p3).
    h4 = np.array([a.arch_span, 0.0, a.ascending_length - 0.30 * transition_drop])

    p4 = np.array([a.arch_span, 0.0, a.ascending_length - a.descending_length])

    n1 = max(12, int(n * 0.15))
    n2 = max(36, int(n * 0.50))
    n3 = max(18, int(n * 0.20))
    n4 = max(8, n - n1 - n2 - n3)

    seg1 = np.array([p0 + (p1 - p0) * t for t in np.linspace(0, 1, n1, endpoint=False)])
    seg2 = np.array([_cubic_bezier(p1, h1, h2, p2, t) for t in np.linspace(0, 1, n2, endpoint=False)])
    seg3 = np.array([_cubic_bezier(p2, h3, h4, p3, t) for t in np.linspace(0, 1, n3, endpoint=False)])
    seg4 = np.array([p3 + (p4 - p3) * t for t in np.linspace(0, 1, n4)])
    return np.concatenate([seg1, seg2, seg3, seg4])


def _arc_lengths(centreline: np.ndarray) -> np.ndarray:
    d = np.diff(centreline, axis=0)
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(d, axis=1))])


# =============================================================================
# Radius profile (taper + hypoplasia + coarctation)
# =============================================================================

def _radius_profile(a: Aorta, centreline: np.ndarray,
                    s_lsa: float) -> np.ndarray:
    """Return per-point radius along the centreline.

    Applied in order: taper -> proximal hypoplasia -> coarctation.
    ``s_lsa`` is the arc-length coordinate of the LSA branch origin, used
    to anchor the coarctation location distal to it.
    """
    s = _arc_lengths(centreline)
    S = max(s[-1], EPS)
    xi = s / S
    r_base = 0.5 * a.d_ascending
    r = r_base * (1.0 - (1.0 - a.taper_ratio) * np.power(xi, a.taper_exponent))

    # Proximal hypoplasia (transverse arch narrowing, raised cosine envelope)
    eta = np.clip(a.proximal_hypoplasia, 0.0, 0.6)
    if eta > 1e-6:
        lo, hi = 0.52 * S, 0.80 * S
        mask = (s > lo) & (s < hi)
        if mask.any():
            x = (s[mask] - lo) / (hi - lo)
            smooth = 0.5 * (1.0 - np.cos(math.pi * x))
            r[mask] *= (1.0 - eta * smooth)
        mask_after = s >= hi
        r[mask_after] *= (1.0 - eta)

    # Coarctation: placed `coarctation_offset_from_lsa` mm distal to LSA
    alpha = np.clip(a.coarctation_area_reduction, 0.0, 0.95)
    if alpha > 1e-6:
        s_coa_centre = s_lsa + a.coarctation_offset_from_lsa
        half = 0.5 * max(a.coarctation_length, 1.0)
        mask = np.abs(s - s_coa_centre) <= half
        if mask.any():
            xi_c = np.clip((s[mask] - (s_coa_centre - half)) / (2 * half), 0.0, 1.0)
            # Asymmetric power-law shape (same as single-vessel power_law)
            sig = np.clip(a.coarctation_shape, 0.0, 1.0)
            a_p = 1.5 + 0.5 * sig
            a_d = 4.0 - 2.0 * sig
            f = np.where(xi_c <= 0.5,
                         (2.0 * xi_c) ** a_p,
                         (2.0 * (1.0 - xi_c)) ** a_d)
            r_throat = r[mask] * math.sqrt(max(1.0 - alpha, 1e-6))
            r[mask] = r[mask] - (r[mask] - r_throat) * f

    return r


# =============================================================================
# Branch placement
# =============================================================================

def _interp_point_by_arc_length(centreline: np.ndarray, s_arr: np.ndarray,
                                s_target: float) -> tuple:
    """Linear interpolation of a point and tangent at arc-length s_target."""
    s_target = float(np.clip(s_target, 0.0, s_arr[-1]))
    i = int(np.searchsorted(s_arr, s_target))
    i = max(1, min(i, len(centreline) - 1))
    denom = max(s_arr[i] - s_arr[i - 1], EPS)
    t = (s_target - s_arr[i - 1]) / denom
    pt = centreline[i - 1] * (1 - t) + centreline[i] * t
    tangent = centreline[i] - centreline[i - 1]
    tangent = tangent / (np.linalg.norm(tangent) + EPS)
    return pt, tangent


def _branch_direction(tangent: np.ndarray,
                      tilt_deg: float, splay_deg: float) -> np.ndarray:
    """Branch direction: tilt in arch-plane toward +z, splay laterally in y."""
    t = tangent / (np.linalg.norm(tangent) + EPS)
    up = np.array([0.0, 0.0, 1.0])
    up_proj = up - t * np.dot(t, up)
    if np.linalg.norm(up_proj) < EPS:
        up_proj = np.array([0.0, 1.0, 0.0])
    up_proj /= np.linalg.norm(up_proj)

    side = np.cross(t, up)
    if np.linalg.norm(side) < EPS:
        side = np.array([0.0, 1.0, 0.0])
    side /= np.linalg.norm(side)

    tilt = math.radians(tilt_deg)
    splay = math.radians(splay_deg)
    d = math.cos(tilt) * t + math.sin(tilt) * up_proj
    d = math.cos(splay) * d + math.sin(splay) * side
    return d / (np.linalg.norm(d) + EPS)


def _compute_branch_specs(a: Aorta, centreline: np.ndarray
                          ) -> List[dict]:
    """Three supra-aortic branches anchored at arch-crown fractions (BCA, LCCA, LSA)."""
    s_arr = _arc_lengths(centreline)
    S = s_arr[-1]
    # Anchor fractions match the 16-parameter web model.
    # Crown is roughly 0.38 of total arc length; space branches by branch_spacing mm.
    s_crown = 0.38 * S
    delta = a.branch_spacing
    specs = [
        ("BCA",  s_crown - 1.5 * delta, a.bca_ratio,  -1),
        ("LCCA", s_crown - 0.5 * delta, a.lcca_ratio,  0),
        ("LSA",  s_crown + 0.5 * delta, a.lsa_ratio,  +1),
    ]
    out = []
    for name, s_target, d_ratio, splay_sign in specs:
        origin, tangent = _interp_point_by_arc_length(centreline, s_arr, s_target)
        direction = _branch_direction(
            tangent,
            tilt_deg=a.branch_tilt_deg,
            splay_deg=splay_sign * a.branch_splay_deg,
        )
        d_branch = a.d_ascending * d_ratio
        # Embed origin slightly inside the arch wall so the union has a clean cut.
        embed = 0.6 * 0.5 * a.d_ascending
        start = origin - direction * embed
        end = origin + direction * a.branch_length
        out.append({
            "name": name,
            "s_origin": s_target,
            "origin": origin,
            "tangent": tangent,
            "direction": direction,
            "diameter": d_branch,
            "start": start,
            "end": end,
            "length": a.branch_length + embed,
        })
    return out


# =============================================================================
# Main generator
# =============================================================================

def generate_aorta_union(a: Aorta,
                         n_sectors: int = 48,
                         n_rings_arch: int = 200,
                         n_rings_branch: int = 60,
                         out: str = "aorta_union.stl") -> dict:
    """Build a watertight aortic arch STL via boolean union.

    Parameters
    ----------
    a : Aorta
    n_sectors : int
    n_rings_arch, n_rings_branch : int
    out : str

    Returns
    -------
    info : dict with ``is_watertight``, ``n_faces``, ``groups`` (dict of
        patch name -> face-count), ``branch_specs``.
    """
    import trimesh
    from trimesh.boolean import union as tm_union

    # --- Main arch (capped at both ends) --------------------------------------
    arch_line = _build_centreline(a, n=n_rings_arch)
    branches = _compute_branch_specs(a, arch_line)
    s_lsa = next(b["s_origin"] for b in branches if b["name"] == "LSA")
    arch_r = _radius_profile(a, arch_line, s_lsa=s_lsa)

    a_v, a_wall, a_inlet, a_outlet, *_ = build_tube(arch_line, arch_r, n_sectors=n_sectors)
    arch_mesh = _trimesh_from_arrays(a_v, [a_wall, a_inlet, a_outlet])

    # --- Three branches (straight capped cylinders) ---------------------------
    branch_meshes = []
    branch_planes = []
    for b in branches:
        line = b["start"][None, :] + np.linspace(0.0, b["length"], n_rings_branch)[:, None] * b["direction"][None, :]
        r = np.full(n_rings_branch, 0.5 * b["diameter"])
        bv, bw, bi, bo, *_ = build_tube(line, r, n_sectors=n_sectors)
        branch_meshes.append(_trimesh_from_arrays(bv, [bw, bi, bo]))
        branch_planes.append((f"{b['name']}_outlet", b["end"], b["direction"],
                              0.7 * b["diameter"]))

    # --- Boolean union --------------------------------------------------------
    merged = tm_union([arch_mesh, *branch_meshes], engine="manifold")

    # --- Relabel patches by plane proximity -----------------------------------
    # Inlet cap: at arch_line[0], normal points out (-z at the ascending start)
    inlet_normal = arch_line[0] - arch_line[1]
    inlet_normal /= np.linalg.norm(inlet_normal) + EPS
    # Descending outlet cap: at arch_line[-1], normal points out (-z at descending end)
    desc_normal = arch_line[-1] - arch_line[-2]
    desc_normal /= np.linalg.norm(desc_normal) + EPS

    planes = [
        ("inlet",              arch_line[0],  inlet_normal, 0.7 * a.d_ascending),
        ("descending_outlet",  arch_line[-1], desc_normal,  0.7 * a.d_ascending),
        *branch_planes,
    ]
    groups = _label_cap_faces(merged, planes, tol_dist=5e-2, tol_normal=0.80)
    groups_ordered = {"wall": groups["wall"]}
    for name, *_rest in planes:
        groups_ordered[name] = groups[name]

    write_multisolid_stl(out, merged.vertices, groups_ordered)

    return {
        "vertices": merged.vertices,
        "groups": groups_ordered,
        "branch_specs": branches,
        "is_watertight": bool(merged.is_watertight),
        "n_faces": len(merged.faces),
    }


# =============================================================================
# Demo
# =============================================================================

def demo_aorta(out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    # Healthy default
    a_healthy = Aorta()
    info_h = generate_aorta_union(a_healthy,
                                  out=os.path.join(out_dir, "aorta_healthy.stl"))
    g = info_h["groups"]
    print(f"  Healthy:   watertight={info_h['is_watertight']}, {info_h['n_faces']} faces")
    print(f"             wall={len(g['wall'])}  inlet={len(g['inlet'])}  "
          f"desc={len(g['descending_outlet'])}  "
          f"BCA={len(g['BCA_outlet'])} LCCA={len(g['LCCA_outlet'])} LSA={len(g['LSA_outlet'])}")

    # Coarctation case
    a_coa = Aorta(
        coarctation_area_reduction=0.70,
        coarctation_shape=0.2,
        coarctation_length=18.0,
    )
    info_c = generate_aorta_union(a_coa,
                                  out=os.path.join(out_dir, "aorta_coarctation.stl"))
    print(f"  CoA 70%:   watertight={info_c['is_watertight']}, {info_c['n_faces']} faces")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="aorta_out")
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args()
    if args.demo:
        demo_aorta(args.out)
    else:
        demo_aorta(args.out)


if __name__ == "__main__":
    main()
