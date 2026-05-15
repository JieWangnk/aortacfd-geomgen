#!/usr/bin/env python3
"""Parametric single-vessel geometry generator (1D -> 3D).

Builds a watertight, CFD-ready ASCII STL from a centreline and a radius
profile, with an optional stenosis. No Blender dependency; the STL can be
imported into Blender, ParaView, or snappyHexMesh as-is.

Features
--------
- Rotation-minimising frames via Wang et al. (2008) double reflection.
- Shared-vertex indexed mesh (watertight; no ring-seam gaps).
- Named patches: ``inlet``, ``outlet``, ``wall`` via multi-solid ASCII STL.
- Pluggable stenosis shapes:
    * ``cosine``    -- Young & Tsai (1973), the CFD literature default
    * ``power_law`` -- asymmetric, shelf/hourglass via sigma
    * ``gaussian``  -- smooth bump
- Built-in validation presets:
    * FDA sudden-expansion nozzle (piecewise conical)
    * Young & Tsai 50% severity straight stenosis
    * Curved vessel (torus arc) to exercise RMF

Usage
-----
Run the built-in demos::

    python vessel_generator.py --demo all --out out/

Call as a library::

    from vessel_generator import (
        generate_vessel, Stenosis, straight_centreline, arc_centreline,
    )
    centreline = arc_centreline(R=50.0, arc_deg=120.0, n=200)
    generate_vessel(
        centreline=centreline,
        radius=5.0,
        stenosis=Stenosis(position=0.5, severity=0.5, length=20.0, shape="cosine"),
        n_sectors=48,
        out="curved_50pct.stl",
    )
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np


# =============================================================================
# Centreline helpers
# =============================================================================

def straight_centreline(length: float, n: int = 200) -> np.ndarray:
    """Straight vessel along +x of given length, n points."""
    x = np.linspace(0.0, length, n)
    return np.column_stack([x, np.zeros(n), np.zeros(n)])


def arc_centreline(R: float, arc_deg: float, n: int = 200) -> np.ndarray:
    """Planar circular arc in the xy-plane, starting at the origin tangent to +x."""
    theta = np.linspace(0.0, math.radians(arc_deg), n)
    x = R * np.sin(theta)
    y = R * (1.0 - np.cos(theta))
    z = np.zeros_like(theta)
    return np.column_stack([x, y, z])


def catmull_rom(ctrl: np.ndarray, n: int = 200) -> np.ndarray:
    """Uniform Catmull-Rom spline through control points (N>=4)."""
    ctrl = np.asarray(ctrl, dtype=float)
    if len(ctrl) < 4:
        raise ValueError("Catmull-Rom needs >=4 control points")
    segs = len(ctrl) - 3
    out = []
    for s in range(segs):
        p0, p1, p2, p3 = ctrl[s:s + 4]
        ts = np.linspace(0.0, 1.0, max(2, n // segs), endpoint=(s == segs - 1))
        for t in ts:
            t2 = t * t
            t3 = t2 * t
            out.append(
                0.5 * (
                    (2.0 * p1) +
                    (-p0 + p2) * t +
                    (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 +
                    (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3
                )
            )
    return np.asarray(out)


def arc_lengths(centreline: np.ndarray) -> np.ndarray:
    d = np.diff(centreline, axis=0)
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(d, axis=1))])


# =============================================================================
# Rotation-minimising frames (double reflection, Wang et al. 2008)
# =============================================================================

def _tangents(centreline: np.ndarray) -> np.ndarray:
    t = np.gradient(centreline, axis=0)
    t /= np.linalg.norm(t, axis=1, keepdims=True) + 1e-15
    return t


def _initial_normal(t0: np.ndarray) -> np.ndarray:
    """Pick a unit vector perpendicular to t0 that is numerically stable."""
    for ref in (np.array([0.0, 0.0, 1.0]),
                np.array([0.0, 1.0, 0.0]),
                np.array([1.0, 0.0, 0.0])):
        n = ref - t0 * np.dot(t0, ref)
        if np.linalg.norm(n) > 1e-6:
            return n / np.linalg.norm(n)
    return np.array([0.0, 1.0, 0.0])


def rmf_double_reflection(centreline: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rotation-minimising frames by double reflection.

    See Wang, Juettler, Zheng, Liu (2008), "Computation of rotation minimizing
    frames", ACM TOG 27(1), Algorithm "double reflection".

    Returns tangents, normals, binormals of shape (N, 3) each.
    """
    c = np.asarray(centreline, dtype=float)
    N = len(c)
    t = _tangents(c)
    normals = np.zeros_like(c)
    binormals = np.zeros_like(c)

    normals[0] = _initial_normal(t[0])
    binormals[0] = np.cross(t[0], normals[0])
    binormals[0] /= np.linalg.norm(binormals[0]) + 1e-15

    for i in range(N - 1):
        v1 = c[i + 1] - c[i]
        c1 = np.dot(v1, v1)
        if c1 < 1e-20:
            normals[i + 1] = normals[i]
            binormals[i + 1] = binormals[i]
            continue
        # First reflection: across plane with normal v1
        rL_n = normals[i] - (2.0 / c1) * np.dot(v1, normals[i]) * v1
        tL_i = t[i] - (2.0 / c1) * np.dot(v1, t[i]) * v1
        # Second reflection: across plane with normal (t[i+1] - tL_i)
        v2 = t[i + 1] - tL_i
        c2 = np.dot(v2, v2)
        if c2 < 1e-20:
            n_next = rL_n
        else:
            n_next = rL_n - (2.0 / c2) * np.dot(v2, rL_n) * v2
        n_next /= np.linalg.norm(n_next) + 1e-15
        normals[i + 1] = n_next
        binormals[i + 1] = np.cross(t[i + 1], n_next)
        binormals[i + 1] /= np.linalg.norm(binormals[i + 1]) + 1e-15

    return t, normals, binormals


# =============================================================================
# Stenosis profiles
# =============================================================================

@dataclass
class Stenosis:
    """Stenosis descriptor.

    Parameters
    ----------
    position : float
        Arc-length fraction (0..1) of the stenosis centre.
    severity : float
        Fractional area reduction at the throat, 0..1.  Throat radius is
        ``r(s) * sqrt(1 - severity)``.
    length : float
        Axial extent in millimetres.
    shape : str
        One of ``cosine``, ``power_law``, ``gaussian``.
    sigma : float
        Shape parameter in [0,1] for ``power_law`` only: 0 = shelf
        (steep proximal, gentle distal), 1 = symmetric hourglass.
    """
    position: float
    severity: float
    length: float
    shape: str = "cosine"
    sigma: float = 1.0


def _stenosis_factor(xi: np.ndarray, stenosis: Stenosis) -> np.ndarray:
    """Return f(xi) in [0,1] with f(0)=f(1)=0, f(0.5)=1 (or peak near centre)."""
    shape = stenosis.shape
    if shape == "cosine":
        # Young & Tsai 1973: r(z) = R0 [1 - (delta/R0)*0.5*(1 + cos(pi*z/L))]
        # Expressed on xi in [0,1], centred at xi=0.5:
        return 0.5 * (1.0 + np.cos(2.0 * math.pi * (xi - 0.5)))
    if shape == "gaussian":
        # Half-width at half-max set to 1/3 of the stenosis length
        return np.exp(-((xi - 0.5) ** 2) / (2.0 * (0.15 ** 2)))
    if shape == "power_law":
        a_p = 1.5 + 0.5 * stenosis.sigma
        a_d = 4.0 - 2.0 * stenosis.sigma
        return np.where(xi <= 0.5,
                        (2.0 * xi) ** a_p,
                        (2.0 * (1.0 - xi)) ** a_d)
    raise ValueError(f"Unknown stenosis shape: {shape!r}")


def apply_stenosis(s: np.ndarray, r: np.ndarray, stenosis: Stenosis) -> np.ndarray:
    """Multiply baseline radius r(s) by a stenosis envelope."""
    s_total = s[-1]
    s_centre = stenosis.position * s_total
    s0 = s_centre - 0.5 * stenosis.length
    s1 = s_centre + 0.5 * stenosis.length
    in_coa = (s >= s0) & (s <= s1)
    r_out = r.copy()
    if not np.any(in_coa):
        return r_out
    xi = np.clip((s[in_coa] - s0) / stenosis.length, 0.0, 1.0)
    f = _stenosis_factor(xi, stenosis)
    r_throat = r[in_coa] * math.sqrt(max(1e-6, 1.0 - stenosis.severity))
    r_out[in_coa] = r[in_coa] - (r[in_coa] - r_throat) * f
    return r_out


# =============================================================================
# Mesh construction (shared-vertex, watertight)
# =============================================================================

def build_tube(centreline: np.ndarray,
               radii: np.ndarray,
               n_sectors: int = 48,
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Build a tube mesh with caps.

    Returns
    -------
    vertices : (V, 3) array
    wall_faces : (Fw, 3) int array
    inlet_faces : (Fi, 3) int array
    outlet_faces : (Fo, 3) int array
    inlet_centre_index : int
    outlet_centre_index : int
    """
    c = np.asarray(centreline, dtype=float)
    r = np.asarray(radii, dtype=float)
    if len(c) != len(r):
        raise ValueError("centreline and radii must have the same length")
    if len(c) < 2:
        raise ValueError("need >=2 centreline points")

    _, normals, binormals = rmf_double_reflection(c)
    N = len(c)

    theta = 2.0 * math.pi * np.arange(n_sectors) / n_sectors
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    # Ring vertices: shape (N, n_sectors, 3) flattened to (N*n_sectors, 3)
    ring_verts = (
        c[:, None, :]
        + r[:, None, None] * (cos_t[None, :, None] * normals[:, None, :]
                              + sin_t[None, :, None] * binormals[:, None, :])
    ).reshape(-1, 3)

    inlet_centre = c[0]
    outlet_centre = c[-1]
    vertices = np.vstack([ring_verts, inlet_centre[None, :], outlet_centre[None, :]])
    inlet_centre_idx = N * n_sectors
    outlet_centre_idx = N * n_sectors + 1

    # Wall triangles
    wall_faces = []
    for i in range(N - 1):
        for j in range(n_sectors):
            jn = (j + 1) % n_sectors
            a = i * n_sectors + j
            b = i * n_sectors + jn
            cc = (i + 1) * n_sectors + jn
            d = (i + 1) * n_sectors + j
            wall_faces.append((a, b, cc))
            wall_faces.append((a, cc, d))
    wall_faces = np.array(wall_faces, dtype=np.int64)

    # Inlet cap: fan from inlet_centre, winding chosen so normal points outward
    # (opposite to t[0], i.e. out of the vessel). We reverse the ring order.
    inlet_faces = []
    for j in range(n_sectors):
        jn = (j + 1) % n_sectors
        inlet_faces.append((inlet_centre_idx, jn, j))
    inlet_faces = np.array(inlet_faces, dtype=np.int64)

    outlet_faces = []
    base = (N - 1) * n_sectors
    for j in range(n_sectors):
        jn = (j + 1) % n_sectors
        outlet_faces.append((outlet_centre_idx, base + j, base + jn))
    outlet_faces = np.array(outlet_faces, dtype=np.int64)

    return vertices, wall_faces, inlet_faces, outlet_faces, inlet_centre_idx, outlet_centre_idx


# =============================================================================
# STL writer (multi-solid ASCII, named patches for snappyHexMesh)
# =============================================================================

def _triangle_normal(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    n = np.cross(p1 - p0, p2 - p0)
    norm = np.linalg.norm(n)
    return n / norm if norm > 1e-15 else np.array([0.0, 0.0, 1.0])


def write_multisolid_stl(path: str,
                         vertices: np.ndarray,
                         groups: dict) -> None:
    """Write an ASCII STL with one ``solid NAME`` block per patch group.

    ``groups`` maps patch name (str) -> (F, 3) int array of triangle indices.
    """
    lines = []
    for name, faces in groups.items():
        lines.append(f"solid {name}")
        for (i0, i1, i2) in faces:
            p0, p1, p2 = vertices[i0], vertices[i1], vertices[i2]
            n = _triangle_normal(p0, p1, p2)
            lines.append(f"  facet normal {n[0]:.9e} {n[1]:.9e} {n[2]:.9e}")
            lines.append("    outer loop")
            lines.append(f"      vertex {p0[0]:.9e} {p0[1]:.9e} {p0[2]:.9e}")
            lines.append(f"      vertex {p1[0]:.9e} {p1[1]:.9e} {p1[2]:.9e}")
            lines.append(f"      vertex {p2[0]:.9e} {p2[1]:.9e} {p2[2]:.9e}")
            lines.append("    endloop")
            lines.append("  endfacet")
        lines.append(f"endsolid {name}")
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# =============================================================================
# Main entry point
# =============================================================================

def generate_vessel(centreline: np.ndarray,
                    radius,
                    stenosis: Stenosis | None = None,
                    n_sectors: int = 48,
                    out: str = "vessel.stl") -> dict:
    """Generate a named-patch STL from a centreline and radius description.

    Parameters
    ----------
    centreline : (N, 3) array
        Sampled centreline points.
    radius : float, (N,) array, or callable s -> r
        Baseline radius description.
    stenosis : Stenosis, optional
        If provided, applied on top of the baseline radius.
    n_sectors : int
        Circumferential resolution.
    out : str
        Output STL path.

    Returns
    -------
    info : dict with keys ``vertices``, ``wall_faces``, ``inlet_faces``,
        ``outlet_faces``, ``radii``, ``arc_length``.
    """
    c = np.asarray(centreline, dtype=float)
    s = arc_lengths(c)

    if callable(radius):
        r = np.array([radius(si) for si in s], dtype=float)
    elif np.isscalar(radius):
        r = np.full(len(c), float(radius))
    else:
        r = np.asarray(radius, dtype=float)
        if len(r) != len(c):
            raise ValueError("radius array must match centreline length")

    if stenosis is not None:
        r = apply_stenosis(s, r, stenosis)

    verts, wall, inlet, outlet, *_ = build_tube(c, r, n_sectors=n_sectors)
    write_multisolid_stl(out, verts, {"wall": wall, "inlet": inlet, "outlet": outlet})

    return {
        "vertices": verts,
        "wall_faces": wall,
        "inlet_faces": inlet,
        "outlet_faces": outlet,
        "radii": r,
        "arc_length": s,
    }


# =============================================================================
# Bifurcation primitive (Y-junction of parent + two daughters)
# =============================================================================

@dataclass
class Bifurcation:
    """Y-junction descriptor for ``generate_bifurcation``.

    The parent vessel ends at a junction point J; two daughters start at J
    and radiate outward at given angles.  Murray's law relates parent and
    daughter diameters via ``d_parent^3 = d_1^3 + d_2^3``; set
    ``use_murray=True`` to derive ``d_1`` and ``d_2`` from ``d_parent`` and
    the split ratio.

    Parameters
    ----------
    d_parent : float
        Parent diameter (mm).
    L_parent : float
        Parent length (mm).
    d_1, d_2 : float
        Daughter diameters (mm).  Ignored when ``use_murray`` is True.
    L_1, L_2 : float
        Daughter lengths (mm).
    angle_total : float
        Full angle between daughters (degrees).
    asymmetry : float
        Angular asymmetry in [-1, 1]; 0 is symmetric.
    use_murray : bool
        If True, enforce d_parent^3 = d_1^3 + d_2^3 with split ``murray_split``.
    murray_split : float
        Fraction of cube-sum carried by daughter 1, in [0.1, 0.9].
    """
    d_parent: float = 20.0
    L_parent: float = 80.0
    d_1: float = 15.87
    d_2: float = 15.87
    L_1: float = 60.0
    L_2: float = 60.0
    angle_total: float = 60.0
    asymmetry: float = 0.0
    use_murray: bool = False
    murray_split: float = 0.5

    def resolve(self) -> "Bifurcation":
        """Return a copy with daughter diameters set via Murray if enabled."""
        if not self.use_murray:
            return self
        total_cube = self.d_parent ** 3
        d1 = (self.murray_split * total_cube) ** (1.0 / 3.0)
        d2 = ((1.0 - self.murray_split) * total_cube) ** (1.0 / 3.0)
        return Bifurcation(
            d_parent=self.d_parent, L_parent=self.L_parent,
            d_1=d1, d_2=d2,
            L_1=self.L_1, L_2=self.L_2,
            angle_total=self.angle_total, asymmetry=self.asymmetry,
            use_murray=False, murray_split=self.murray_split,
        )


def _trimesh_from_arrays(verts: np.ndarray, faces: list) -> "trimesh.Trimesh":
    """Stack face arrays and wrap in a trimesh.Trimesh (closed, oriented)."""
    import trimesh
    all_faces = np.vstack(faces)
    mesh = trimesh.Trimesh(vertices=verts.copy(), faces=all_faces.copy(),
                           process=True)
    # process=True already merges duplicate vertices and removes unreferenced ones
    mesh.fix_normals()
    return mesh


def _label_cap_faces(mesh: "trimesh.Trimesh",
                     planes: list,
                     tol_dist: float = 1e-3,
                     tol_normal: float = 0.9) -> dict:
    """Split triangles of a merged mesh into named patches by plane proximity.

    ``planes`` is a list of ``(name, point, normal, radius)``.  A face is
    assigned to a patch if all three vertices lie within ``tol_dist`` of the
    plane and the face normal is (anti-)aligned with the plane normal.
    Everything else is labelled ``wall``.
    """
    V = mesh.vertices
    F = mesh.faces
    Fn = mesh.face_normals
    centroids = V[F].mean(axis=1)
    assigned = np.full(len(F), -1, dtype=int)

    for i, (_name, pt, nrm, _rad) in enumerate(planes):
        nrm = np.asarray(nrm, dtype=float)
        nrm /= np.linalg.norm(nrm)
        pt = np.asarray(pt, dtype=float)

        d_vert = np.abs((V[F] - pt) @ nrm).max(axis=1)  # max |dist| over triangle
        cos_n = np.abs(Fn @ nrm)
        mask = (d_vert < tol_dist) & (cos_n > tol_normal) & (assigned < 0)
        assigned[mask] = i

    wall_mask = assigned < 0
    groups = {"wall": F[wall_mask]}
    for i, (name, *_rest) in enumerate(planes):
        groups[name] = F[assigned == i]
    return groups


def generate_bifurcation_union(bif: Bifurcation,
                               n_sectors: int = 48,
                               n_rings: int = 120,
                               out: str = "bifurcation_union.stl") -> dict:
    """Build a truly watertight Y-junction via trimesh boolean union.

    Requires ``trimesh`` and ``manifold3d`` (``pip install trimesh manifold3d``).
    Unions three capped tubes, then relabels faces into
    ``inlet`` / ``outlet_1`` / ``outlet_2`` / ``wall`` patches by matching
    against the known cap planes.  Produces a single watertight surface with
    no overlapping walls at the junction -- CFD-ready.
    """
    import trimesh
    from trimesh.boolean import union as tm_union

    b = bif.resolve()

    # ---- Parent tube (capped) ------------------------------------------------
    parent_line = straight_centreline(length=b.L_parent, n=n_rings)
    r_parent = np.full(n_rings, 0.5 * b.d_parent)
    p_v, p_wall, p_inlet, p_outlet, *_ = build_tube(parent_line, r_parent, n_sectors=n_sectors)
    parent_mesh = _trimesh_from_arrays(p_v, [p_wall, p_inlet, p_outlet])

    # ---- Daughter directions -------------------------------------------------
    junction = parent_line[-1]
    t_p = np.array([1.0, 0.0, 0.0])
    n_axis = np.array([0.0, 0.0, 1.0])
    b_axis = np.cross(t_p, n_axis)
    half = 0.5 * math.radians(b.angle_total)
    theta1 = half + half * b.asymmetry
    theta2 = -half + half * b.asymmetry

    def _rotate(v, axis, theta):
        c_, s_ = math.cos(theta), math.sin(theta)
        return v * c_ + np.cross(axis, v) * s_ + axis * np.dot(axis, v) * (1 - c_)

    d1_dir = _rotate(t_p, b_axis, theta1)
    d2_dir = _rotate(t_p, b_axis, theta2)

    # ---- Daughters (capped, extended back through J so CSG crotch is clean) --
    # Extending daughters backward by overlap_back avoids knife-edges where the
    # daughter cross-section meets the parent wall obliquely.
    overlap_back = 0.6 * b.d_parent
    d1_start = junction - d1_dir * overlap_back
    d2_start = junction - d2_dir * overlap_back
    d1_end = junction + d1_dir * b.L_1
    d2_end = junction + d2_dir * b.L_2
    L1_ext = b.L_1 + overlap_back
    L2_ext = b.L_2 + overlap_back

    d1_line = d1_start[None, :] + np.linspace(0.0, L1_ext, n_rings)[:, None] * d1_dir[None, :]
    d2_line = d2_start[None, :] + np.linspace(0.0, L2_ext, n_rings)[:, None] * d2_dir[None, :]
    r_d1 = np.full(n_rings, 0.5 * b.d_1)
    r_d2 = np.full(n_rings, 0.5 * b.d_2)

    d1v, d1w, d1i, d1o, *_ = build_tube(d1_line, r_d1, n_sectors=n_sectors)
    d2v, d2w, d2i, d2o, *_ = build_tube(d2_line, r_d2, n_sectors=n_sectors)
    d1_mesh = _trimesh_from_arrays(d1v, [d1w, d1i, d1o])
    d2_mesh = _trimesh_from_arrays(d2v, [d2w, d2i, d2o])

    # ---- Boolean union -------------------------------------------------------
    merged = tm_union([parent_mesh, d1_mesh, d2_mesh], engine="manifold")
    if not merged.is_watertight:
        print(f"  warning: union produced non-watertight mesh "
              f"({len(merged.faces)} faces, {len(merged.fill_holes())} holes filled)")

    # ---- Relabel patches by plane proximity ----------------------------------
    # Each cap is a flat disk in a known plane.
    planes = [
        ("inlet",    parent_line[0], -t_p,      0.6 * b.d_parent),
        ("outlet_1", d1_end,          d1_dir,   0.6 * b.d_1),
        ("outlet_2", d2_end,          d2_dir,   0.6 * b.d_2),
    ]
    groups = _label_cap_faces(merged, planes, tol_dist=1e-2, tol_normal=0.85)
    groups_ordered = {
        "wall":     groups["wall"],
        "inlet":    groups["inlet"],
        "outlet_1": groups["outlet_1"],
        "outlet_2": groups["outlet_2"],
    }

    write_multisolid_stl(out, merged.vertices, groups_ordered)

    return {
        "vertices": merged.vertices,
        "groups": groups_ordered,
        "junction": junction,
        "d1_direction": d1_dir,
        "d2_direction": d2_dir,
        "resolved_diameters": (b.d_parent, b.d_1, b.d_2),
        "is_watertight": bool(merged.is_watertight),
        "n_faces": len(merged.faces),
    }


def generate_bifurcation_smooth(bif: Bifurcation,
                                smoothness: float = 3.0,
                                grid_res: float = 0.8,
                                out: str = "bifurcation_smooth.stl") -> dict:
    """Build a smooth watertight Y-junction via SDF + marching cubes.

    Unlike ``generate_bifurcation_union`` (which boolean-unions three triangle
    meshes and leaves T-junctions where all three walls meet at the crotch),
    this builds a signed distance field for each tube, blends them via a
    polynomial smooth-min, and extracts the zero-level set with marching
    cubes.  The resulting surface is genuinely manifold (no T-junctions) with
    a smooth fillet at the crotch controlled by ``smoothness``.

    Parameters
    ----------
    bif : Bifurcation
    smoothness : float
        Fillet radius at the crotch in mm.  0.1 -> sharp union, 5 -> very
        rounded.  Default 3 mm is a visually pleasant aortic-style blend.
    grid_res : float
        Voxel size for marching cubes in mm.  0.8 mm gives ~50k triangles for
        a default geometry.  Halving it doubles resolution at 8x cost.
    out : str
        Output STL path.

    Requires ``pip install scipy scikit-image trimesh``.
    """
    from scipy.spatial import cKDTree
    from skimage.measure import marching_cubes
    import trimesh

    b = bif.resolve()

    # --- Build three centrelines (densified for accurate nearest-point SDF) --
    n_dense = 200
    parent_line = straight_centreline(length=b.L_parent, n=n_dense)
    junction = parent_line[-1]

    t_p = np.array([1.0, 0.0, 0.0])
    n_axis = np.array([0.0, 0.0, 1.0])
    b_axis = np.cross(t_p, n_axis)
    half = 0.5 * math.radians(b.angle_total)
    theta1 = half + half * b.asymmetry
    theta2 = -half + half * b.asymmetry

    def _rotate(v, axis, theta):
        c_, s_ = math.cos(theta), math.sin(theta)
        return v * c_ + np.cross(axis, v) * s_ + axis * np.dot(axis, v) * (1 - c_)

    d1_dir = _rotate(t_p, b_axis, theta1)
    d2_dir = _rotate(t_p, b_axis, theta2)

    d1_end = junction + d1_dir * b.L_1
    d2_end = junction + d2_dir * b.L_2
    d1_line = junction[None, :] + np.linspace(0.0, b.L_1, n_dense)[:, None] * d1_dir[None, :]
    d2_line = junction[None, :] + np.linspace(0.0, b.L_2, n_dense)[:, None] * d2_dir[None, :]

    r_parent = 0.5 * b.d_parent
    r_d1 = 0.5 * b.d_1
    r_d2 = 0.5 * b.d_2

    # --- Build voxel grid bounded by the union of the three capsule bboxes --
    # Pad generously: the grid must contain the rounded capsule extensions
    # with room to spare, or the isosurface will touch the grid boundary and
    # MC will leave open edges.  Rule of thumb: 3x max radius + smoothness.
    pad = 15.0 + 2.0 * smoothness + 2.5 * max(r_parent, r_d1, r_d2)
    all_pts = np.vstack([parent_line, d1_line, d2_line])
    mins = all_pts.min(axis=0) - pad
    maxs = all_pts.max(axis=0) + pad
    xs = np.arange(mins[0], maxs[0] + grid_res, grid_res)
    ys = np.arange(mins[1], maxs[1] + grid_res, grid_res)
    zs = np.arange(mins[2], maxs[2] + grid_res, grid_res)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    query = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    def _sdf_capsule(pts, p0, p1, r):
        """SDF to a constant-radius capsule (rounded-ends cylinder) from p0 to p1."""
        d = p1 - p0
        L2 = float(d @ d)
        t = np.clip(((pts - p0) @ d) / L2, 0.0, 1.0)
        proj = p0 + t[:, None] * d
        return np.linalg.norm(pts - proj, axis=1) - r

    def _sdf_taper_end(pts, p0, p1, r0, r_end, taper_frac=0.15):
        """Capsule with constant radius r0 for the first (1-taper_frac) of its
        length, then smoothstep transition to r_end at p1."""
        d = p1 - p0
        L2 = float(d @ d)
        t = np.clip(((pts - p0) @ d) / L2, 0.0, 1.0)
        proj = p0 + t[:, None] * d
        u = np.clip((t - (1.0 - taper_frac)) / taper_frac, 0.0, 1.0)
        u = u * u * (3.0 - 2.0 * u)
        r_local = r0 + (r_end - r0) * u
        return np.linalg.norm(pts - proj, axis=1) - r_local

    def _sdf_taper_start(pts, p0, p1, r_start, r1, taper_frac=0.25):
        """Capsule flared at the base: radius = r_start at p0 (t=0), smoothstep
        transition to r1 within the first taper_frac, constant after."""
        d = p1 - p0
        L2 = float(d @ d)
        t = np.clip(((pts - p0) @ d) / L2, 0.0, 1.0)
        proj = p0 + t[:, None] * d
        u = np.clip(t / taper_frac, 0.0, 1.0)
        u = u * u * (3.0 - 2.0 * u)
        r_local = r_start + (r1 - r_start) * u
        return np.linalg.norm(pts - proj, axis=1) - r_local

    def _smin_poly(a, b, k):
        """Inigo Quilez polynomial smooth-min (C1 continuous)."""
        h = np.maximum(k - np.abs(a - b), 0.0) / max(k, 1e-6)
        return np.minimum(a, b) - 0.25 * h * h * k

    # Constant-radius capsules.  The polynomial smin gradually widens the
    # junction to ~1.15*r_parent at J -- this is physiologically correct
    # (Murray's law: equal-daughter bifurcation has 1.26x cross-section at
    # the split) and gives a smooth monotonic transition with no neck.
    sdf_p = _sdf_capsule(query, parent_line[0], junction, r_parent)
    sdf_1 = _sdf_capsule(query, junction,       d1_end,   r_d1)
    sdf_2 = _sdf_capsule(query, junction,       d2_end,   r_d2)
    sdf = _smin_poly(_smin_poly(sdf_p, sdf_1, smoothness), sdf_2, smoothness)
    sdf = sdf.reshape(X.shape)

    # --- Marching cubes ------------------------------------------------------
    verts_mc, faces_mc, _normals_mc, _vals = marching_cubes(
        sdf, level=0.0, spacing=(grid_res, grid_res, grid_res))
    verts_mc += mins  # shift grid coords back to world

    # CRITICAL: process=False preserves marching cubes' manifold guarantee.
    # With process=True trimesh merges near-identical vertices with its default
    # tolerance, which glues adjacent MC triangles together and creates false
    # non-manifold edges where the surface was previously clean.
    mesh = trimesh.Trimesh(vertices=verts_mc, faces=faces_mc, process=False)
    mesh.fix_normals()

    # Weld near-identical vertices (MC can produce sub-epsilon duplicates
    # that collapse to degenerate triangles on STL reload).  Merging here
    # with 1e-5 mm tolerance lets us drop the resulting zero-area faces
    # explicitly, so the written STL round-trips cleanly.
    mesh.merge_vertices(digits_vertex=5)
    f = mesh.faces
    non_deg = (f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])
    if (~non_deg).any():
        mesh.update_faces(non_deg)
        mesh.remove_unreferenced_vertices()

    # Replace the rounded capsule ends with flat caps by plane slicing.  Each
    # slice trims the mesh at the cap plane and closes the cut with shapely-
    # triangulated cap polygons -- so the final STL has flat inlet / outlet
    # patches suitable for CFD boundary conditions, while keeping the smooth
    # filleted crotch.
    from trimesh.intersections import slice_mesh_plane
    mesh = slice_mesh_plane(mesh, plane_normal=t_p,
                            plane_origin=parent_line[0], cap=True)
    mesh = slice_mesh_plane(mesh, plane_normal=-d1_dir,
                            plane_origin=d1_end, cap=True)
    mesh = slice_mesh_plane(mesh, plane_normal=-d2_dir,
                            plane_origin=d2_end, cap=True)
    mesh.fix_normals()

    # --- Label patches by plane proximity (caps are now flat) ----------------
    plane_specs = [
        ("inlet",    parent_line[0], -t_p,    0.7 * b.d_parent),
        ("outlet_1", d1_end,          d1_dir, 0.7 * b.d_1),
        ("outlet_2", d2_end,          d2_dir, 0.7 * b.d_2),
    ]
    groups = _label_cap_faces(mesh, plane_specs,
                              tol_dist=0.05, tol_normal=0.80)
    groups_ordered = {
        "wall":     groups["wall"],
        "inlet":    groups["inlet"],
        "outlet_1": groups["outlet_1"],
        "outlet_2": groups["outlet_2"],
    }

    write_multisolid_stl(out, mesh.vertices, groups_ordered)

    # Topology diagnostic
    from collections import Counter
    edges = mesh.edges_sorted
    _u, counts = np.unique(edges, axis=0, return_counts=True)

    return {
        "vertices": mesh.vertices,
        "groups": groups_ordered,
        "junction": junction,
        "d1_direction": d1_dir,
        "d2_direction": d2_dir,
        "resolved_diameters": (b.d_parent, b.d_1, b.d_2),
        "is_watertight": bool(mesh.is_watertight),
        "n_faces": len(mesh.faces),
        "boundary_edges": int((counts == 1).sum()),
        "nonmanifold_edges": int((counts >= 3).sum()),
    }


def generate_bifurcation(bif: Bifurcation,
                         n_sectors: int = 48,
                         n_rings: int = 120,
                         out: str = "bifurcation.stl") -> dict:
    """Build a Y-junction STL with named patches.

    Writes three tubes (parent + two daughters) into a single ASCII STL with
    six named solids: ``wall_parent``, ``wall_d1``, ``wall_d2``, ``inlet``,
    ``outlet_1``, ``outlet_2``.  The walls overlap near the junction point;
    downstream meshers (snappyHexMesh) union them implicitly when building
    the volume mesh.  For CFD-grade watertight stitching use a CSG backend
    (trimesh + manifold3d) as a post-process.

    Parameters
    ----------
    bif : Bifurcation
        Junction descriptor.  Murray's law is applied via ``bif.resolve()``.
    n_sectors, n_rings : int
        Circumferential / axial resolution per tube.
    out : str
        Output STL path.

    Returns
    -------
    info : dict with per-patch vertex and face arrays, plus ``junction`` point
        and daughter direction vectors.
    """
    b = bif.resolve()

    # --- Parent: straight along +x from origin to junction at x = L_parent ---
    parent_line = straight_centreline(length=b.L_parent, n=n_rings)
    r_parent = np.full(n_rings, 0.5 * b.d_parent)
    p_verts, p_wall, p_inlet, _p_outlet, _i_c, _o_c = build_tube(
        parent_line, r_parent, n_sectors=n_sectors)
    # We do NOT cap the parent's downstream end (it opens into the junction).
    # Keep only inlet cap; discard outlet cap faces.
    p_outlet_discarded = _p_outlet  # unused

    junction = parent_line[-1]
    t_p = np.array([1.0, 0.0, 0.0])                 # parent tangent at J
    n_axis = np.array([0.0, 0.0, 1.0])              # reference "up"
    b_axis = np.cross(t_p, n_axis)                  # rotation axis for daughters

    half = 0.5 * math.radians(b.angle_total)
    theta1 = half + half * b.asymmetry              # daughter 1 angle from t_p
    theta2 = -half + half * b.asymmetry             # daughter 2 angle from t_p

    def _rotate(v: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
        """Rodrigues rotation of v about axis by theta (axis assumed unit)."""
        c_, s_ = math.cos(theta), math.sin(theta)
        return v * c_ + np.cross(axis, v) * s_ + axis * np.dot(axis, v) * (1 - c_)

    d1_dir = _rotate(t_p, b_axis, theta1)
    d2_dir = _rotate(t_p, b_axis, theta2)

    d1_line = junction[None, :] + np.linspace(0.0, b.L_1, n_rings)[:, None] * d1_dir[None, :]
    d2_line = junction[None, :] + np.linspace(0.0, b.L_2, n_rings)[:, None] * d2_dir[None, :]
    r_d1 = np.full(n_rings, 0.5 * b.d_1)
    r_d2 = np.full(n_rings, 0.5 * b.d_2)

    d1_verts, d1_wall, _d1_inlet, d1_outlet, *_ = build_tube(
        d1_line, r_d1, n_sectors=n_sectors)
    d2_verts, d2_wall, _d2_inlet, d2_outlet, *_ = build_tube(
        d2_line, r_d2, n_sectors=n_sectors)

    # Concatenate with vertex offsets
    n_p = len(p_verts)
    n_d1 = len(d1_verts)
    verts = np.vstack([p_verts, d1_verts, d2_verts])
    d1_wall_o = d1_wall + n_p
    d1_outlet_o = d1_outlet + n_p
    d2_wall_o = d2_wall + n_p + n_d1
    d2_outlet_o = d2_outlet + n_p + n_d1

    groups = {
        "wall_parent": p_wall,
        "wall_d1":     d1_wall_o,
        "wall_d2":     d2_wall_o,
        "inlet":       p_inlet,
        "outlet_1":    d1_outlet_o,
        "outlet_2":    d2_outlet_o,
    }
    write_multisolid_stl(out, verts, groups)

    return {
        "vertices": verts,
        "groups": groups,
        "junction": junction,
        "d1_direction": d1_dir,
        "d2_direction": d2_dir,
        "resolved_diameters": (b.d_parent, b.d_1, b.d_2),
    }


# =============================================================================
# Validation demos
# =============================================================================

def demo_young_tsai(out_dir: str) -> None:
    """50% diameter reduction cosine stenosis on a straight vessel (Young & Tsai 1973)."""
    L = 100.0
    c = straight_centreline(length=L, n=400)
    # 50% diameter reduction == 75% area reduction
    stenosis = Stenosis(position=0.5, severity=0.75, length=20.0, shape="cosine")
    info = generate_vessel(c, radius=5.0, stenosis=stenosis, n_sectors=64,
                           out=os.path.join(out_dir, "young_tsai_50pct.stl"))
    print(f"  Young & Tsai 50%: r_min={info['radii'].min():.3f} mm "
          f"(expected {5.0 * math.sqrt(0.25):.3f} mm)")


def demo_fda_nozzle(out_dir: str) -> None:
    """FDA sudden-expansion nozzle (piecewise conical).

    Geometry from the FDA CFD benchmark: D_upstream=12 mm, D_throat=4 mm,
    throat length 40 mm, upstream/downstream extensions and conical transitions.
    """
    L_up = 80.0
    L_cone_in = 22.3  # length of 18.1 deg half-angle contraction for r: 6 -> 2
    L_throat = 40.0
    L_cone_out = 38.0  # sudden expansion would be ~0; FDA uses gradual version
    L_down = 120.0
    total = L_up + L_cone_in + L_throat + L_cone_out + L_down

    n = 1200
    s = np.linspace(0.0, total, n)
    r = np.empty_like(s)
    R_up, R_th = 6.0, 2.0
    for i, si in enumerate(s):
        if si <= L_up:
            r[i] = R_up
        elif si <= L_up + L_cone_in:
            u = (si - L_up) / L_cone_in
            r[i] = R_up + (R_th - R_up) * u
        elif si <= L_up + L_cone_in + L_throat:
            r[i] = R_th
        elif si <= L_up + L_cone_in + L_throat + L_cone_out:
            u = (si - L_up - L_cone_in - L_throat) / L_cone_out
            r[i] = R_th + (R_up - R_th) * u
        else:
            r[i] = R_up

    c = straight_centreline(length=total, n=n)
    info = generate_vessel(c, radius=r, n_sectors=64,
                           out=os.path.join(out_dir, "fda_nozzle.stl"))
    print(f"  FDA nozzle: throat radius {info['radii'].min():.3f} mm (expected 2.000 mm)")


def demo_curved_vessel(out_dir: str) -> None:
    """120 degree curved vessel to exercise RMF."""
    c = arc_centreline(R=50.0, arc_deg=120.0, n=300)
    stenosis = Stenosis(position=0.5, severity=0.5, length=20.0,
                        shape="power_law", sigma=0.0)
    info = generate_vessel(c, radius=5.0, stenosis=stenosis, n_sectors=64,
                           out=os.path.join(out_dir, "curved_shelf_50pct.stl"))
    print(f"  Curved: r_min={info['radii'].min():.3f} mm "
          f"(expected {5.0 * math.sqrt(0.5):.3f} mm)")


def demo_bifurcation(out_dir: str) -> None:
    """Symmetric Y-junction with Murray's law (aortic-like).

    Writes two STLs: the simple overlap-tubes version (for visualisation) and
    the boolean-unioned watertight version (for CFD meshing).
    """
    bif = Bifurcation(
        d_parent=20.0, L_parent=80.0,
        L_1=60.0, L_2=60.0,
        angle_total=60.0, asymmetry=0.0,
        use_murray=True, murray_split=0.5,
    )
    info = generate_bifurcation(bif, n_sectors=48, n_rings=100,
                                out=os.path.join(out_dir, "bifurcation_overlap.stl"))
    d_p, d1, d2 = info["resolved_diameters"]
    murray_lhs = d_p ** 3
    murray_rhs = d1 ** 3 + d2 ** 3
    print(f"  Bifurcation (overlap): d_p={d_p:.2f}, d_1={d1:.2f}, d_2={d2:.2f} mm")
    print(f"                         Murray residual: {abs(murray_lhs - murray_rhs):.3e}")

    try:
        info_u = generate_bifurcation_union(
            bif, n_sectors=48, n_rings=100,
            out=os.path.join(out_dir, "bifurcation_union.stl"),
        )
        g = info_u["groups"]
        print(f"  Bifurcation (union):   watertight={info_u['is_watertight']}, "
              f"{info_u['n_faces']} faces")
        print(f"                         wall={len(g['wall'])} "
              f"inlet={len(g['inlet'])} "
              f"outlet_1={len(g['outlet_1'])} "
              f"outlet_2={len(g['outlet_2'])}")
    except ImportError as e:
        print(f"  Bifurcation (union):   skipped ({e})."
              f"  Install: pip install trimesh manifold3d")

    try:
        info_s = generate_bifurcation_smooth(
            bif, smoothness=3.0, grid_res=0.8,
            out=os.path.join(out_dir, "bifurcation_smooth.stl"),
        )
        g = info_s["groups"]
        print(f"  Bifurcation (smooth):  watertight={info_s['is_watertight']}, "
              f"{info_s['n_faces']} faces, "
              f"boundary={info_s['boundary_edges']}, "
              f"nonmanifold={info_s['nonmanifold_edges']}")
        print(f"                         wall={len(g['wall'])} "
              f"inlet={len(g['inlet'])} "
              f"outlet_1={len(g['outlet_1'])} "
              f"outlet_2={len(g['outlet_2'])}")
    except ImportError as e:
        print(f"  Bifurcation (smooth):  skipped ({e})."
              f"  Install: pip install scipy scikit-image")


def demo_all(out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    print(f"Writing demos to {out_dir}/")
    demo_young_tsai(out_dir)
    demo_fda_nozzle(out_dir)
    demo_curved_vessel(out_dir)
    demo_bifurcation(out_dir)
    print("Done.")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--demo", choices=["all", "young_tsai", "fda", "curved", "bifurcation"],
                        default="all", help="which built-in demo to run")
    parser.add_argument("--out", default="out_vessels",
                        help="output directory")
    args = parser.parse_args()

    if args.demo == "all":
        demo_all(args.out)
    elif args.demo == "young_tsai":
        os.makedirs(args.out, exist_ok=True)
        demo_young_tsai(args.out)
    elif args.demo == "fda":
        os.makedirs(args.out, exist_ok=True)
        demo_fda_nozzle(args.out)
    elif args.demo == "curved":
        os.makedirs(args.out, exist_ok=True)
        demo_curved_vessel(args.out)
    elif args.demo == "bifurcation":
        os.makedirs(args.out, exist_ok=True)
        demo_bifurcation(args.out)


if __name__ == "__main__":
    main()
