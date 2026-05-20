#!/usr/bin/env python3
"""Healthy aorta geometry generator (v2) — three segments, no branches, no coarctation.

Designed to be the geometrical companion of the SynthAorta dataset
(Bošnjak et al. 2025, arXiv:2409.08635) inside our parameter-sweep
framework. Produces a watertight, single-solid aortic arch tube:

    inlet (s=0)  ──ascending──  arch  ──descending── outlet (s=1)

with three per-segment radii (r_ascending, r_arch, r_descending) joined
by a smoothstep / linear / piecewise taper, and the arch defined by a
circular arc with clinical radius of curvature ``arch_R_c`` and
subtended angle ``arch_angle_deg``.

Run (inside Blender, headless)::

  blender -b -P blender_aorta_v2.py -- \\
      --r_ascending 13.7 --r_arch 13.0 --r_descending 12.2 \\
      --ascending_length 50 --arch_R_c 40.4 --arch_angle_deg 180 \\
      --descending_length 200 --taper_mode smoothstep \\
      --segments_radial 64 --curve_samples 220 \\
      --metadata --triangulate \\
      --output /tmp/v2_baseline.stl

Compared with ``blender_aorta_like_generator.py`` (v1), this script:
  - has no supra-aortic branches (single inlet + single outlet)
  - has no coarctation / hypoplasia / roughness
  - varies radius along the centreline (v1's main lumen is one diameter)
  - parameterises curvature with a clinical R_c + subtended angle
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from math import pi
from typing import Sequence

import bpy
from mathutils import Vector

EPS = 1.0e-9


# -----------------------------------------------------------------------------
# Args
# -----------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="blender_aorta_v2", description=__doc__)

    # Per-segment radii (replaces v1's single --diameter)
    parser.add_argument("--r_ascending", type=float, default=13.7,
                        help="Ascending aorta radius [mm]. Default from SynthAorta Table I "
                             "(Schäfer 2018; Wolak 2008).")
    parser.add_argument("--r_arch", type=float, default=13.0,
                        help="Arch radius [mm]. Default from SynthAorta Table I.")
    parser.add_argument("--r_descending", type=float, default=12.2,
                        help="Descending aorta radius [mm]. Default from SynthAorta Table I.")
    parser.add_argument("--taper_mode", choices=["piecewise", "linear", "smoothstep"],
                        default="smoothstep",
                        help="How to blend r_ascending → r_arch → r_descending across "
                             "segment boundaries. smoothstep = cosine-Hermite (C¹).")

    # Lengths
    parser.add_argument("--ascending_length", type=float, default=50.0,
                        help="Ascending aorta length [mm].")
    parser.add_argument("--descending_length", type=float, default=200.0,
                        help="Descending aorta length [mm].")

    # Arch geometry — clinical R_c + subtended angle
    parser.add_argument("--arch_R_c", type=float, default=40.4,
                        help="Arch radius of curvature [mm]. Default from SynthAorta Table I "
                             "Gumbel(40.4, 2.4). Clinical: Choi 2017, Saitta 2022.")
    parser.add_argument("--arch_angle_deg", type=float, default=180.0,
                        help="Subtended angle of the arch arc [deg]. 180=U-arch, "
                             "<180=shallow, >180=over-arched. Engineering range [120, 200].")
    parser.add_argument("--arch_tilt_deg", type=float, default=0.0,
                        help="RIGID rotation of the arch+descending segments around the "
                             "inlet z-axis [deg]. 0=arch lies in xz-plane (default). "
                             ">0=arch tilts toward +y. Anatomically the real aortic "
                             "arch tilts ~5-15° to the patient's left. Arch stays "
                             "PLANAR (just in a rotated plane).")
    parser.add_argument("--arch_twist_deg", type=float, default=0.0,
                        help="GRADUAL twist around the z-axis along the arch [deg]. "
                             "0=no twist. Ramps linearly from 0 at the ascending-arch "
                             "boundary to twist_deg at the arch-descending boundary, "
                             "then held constant through descending. Result: arch "
                             "becomes a NON-PLANAR 3D curve (unlike arch_tilt_deg "
                             "which keeps the arch planar in a rotated plane).")

    # Non-planar Fourier multipliers (SynthAorta Eq 13)
    parser.add_argument("--delta_3", type=float, default=0.0,
                        help="SynthAorta non-planar Fourier δ_3: scales the "
                             "cos(2w·||x||) second-harmonic out-of-plane displacement. "
                             "Default 0.0 = strictly planar (backwards-compat). "
                             "1.0 = SynthAorta nominal centreline shape. "
                             "Sample around 1.0 with std≈0.09 for SynthAorta variability "
                             "(Bošnjak et al. 2025, Table I).")
    parser.add_argument("--delta_4", type=float, default=0.0,
                        help="SynthAorta non-planar Fourier δ_4: scales the "
                             "sin(2w·||x||) second-harmonic. See --delta_3.")

    # Mesh resolution
    parser.add_argument("--segments_radial", type=int, default=96,
                        help="Circumferential ring vertices. Higher = smoother tube wall.")
    parser.add_argument("--curve_samples", type=int, default=300,
                        help="Total centreline sample count (split across the 5 sub-segments "
                             "proportionally to their arc length). Higher = smoother arch.")
    parser.add_argument("--junction_blend_mm", type=float, default=12.0,
                        help="Length of the cubic-Bezier blend zone inserted at each "
                             "ascending→arch and arch→descending junction. 0 = no blend "
                             "(sharp circular-arc corners — visible C² fold). Default "
                             "12 mm gives a visually smooth transition while keeping the "
                             "middle of the arch on the exact circular arc.")

    # Output
    parser.add_argument("--output", required=True, help="Output STL path.")
    parser.add_argument("--triangulate", action="store_true",
                        help="Triangulate before export (CFD-friendly).")
    parser.add_argument("--metadata", action="store_true",
                        help="Write <output>.json sidecar with parameters + derived geometry.")
    parser.add_argument("--save_blend", action="store_true",
                        help="Save .blend file alongside the STL (debugging).")

    argv = sys.argv
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    else:
        argv = argv[1:]
    return parser.parse_args(argv)


# -----------------------------------------------------------------------------
# Geometry math
# -----------------------------------------------------------------------------


def derive_arch_geometry(R_c: float, angle_deg: float) -> dict:
    """Compute arch span, height, peak position from R_c and subtended angle.

    The arch is a circular arc lying in the xz-plane, with the inlet at
    the origin (0,0,0), ascending up along +z to (0,0,asc_length), then
    arcing in the +x direction. Circle centre is at (R_c, 0, asc_length).

    Arc point at parameter φ ∈ [0, θ]::

        P(φ) = (R_c*(1 - cos φ), 0, asc_length + R_c*sin φ)

    Returns a dict with span, height (above asc top), end_dz (signed), peak_dz.
    """
    theta = math.radians(angle_deg)
    span = R_c * (1.0 - math.cos(theta))
    end_dz = R_c * math.sin(theta)
    # Peak occurs at φ=π/2 if θ ≥ π/2, else at φ=θ
    peak_dz = R_c * 1.0 if theta >= math.pi / 2.0 else R_c * math.sin(theta)
    return {"span": span, "end_dz": end_dz, "peak_dz": peak_dz, "theta_rad": theta}


def smoothstep(t: float) -> float:
    """C¹ Hermite smoothstep: 3t² - 2t³. Clamped to [0,1]."""
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


def radius_along_arc(s_arc: float, total_arc: float, asc_len: float, arch_len: float,
                     r_asc: float, r_arch: float, r_desc: float,
                     taper_mode: str, blend_window: float = 30.0) -> float:
    """Return the lumen radius at arc-length position ``s_arc`` along the full centreline.

    Three segments: ascending [0, asc_len], arch [asc_len, asc_len+arch_len],
    descending [asc_len+arch_len, total_arc]. Two boundaries to blend across.
    """
    s1 = asc_len                    # ascending → arch boundary
    s2 = asc_len + arch_len         # arch → descending boundary

    if taper_mode == "piecewise":
        if s_arc < s1:
            return r_asc
        elif s_arc < s2:
            return r_arch
        else:
            return r_desc

    # Identify segment + blend interpolant
    def _segment_radius(s):
        if s < s1:
            return r_asc
        elif s < s2:
            return r_arch
        else:
            return r_desc

    if taper_mode == "linear":
        blend = lambda t: t
    elif taper_mode == "smoothstep":
        blend = smoothstep
    else:
        raise ValueError(f"Unknown taper_mode: {taper_mode!r}")

    half_w = 0.5 * blend_window

    # Near the ascending→arch boundary
    if abs(s_arc - s1) < half_w:
        t = (s_arc - (s1 - half_w)) / blend_window
        return r_asc + (r_arch - r_asc) * blend(t)
    # Near the arch→descending boundary
    if abs(s_arc - s2) < half_w:
        t = (s_arc - (s2 - half_w)) / blend_window
        return r_arch + (r_desc - r_arch) * blend(t)
    return _segment_radius(s_arc)


# ── SynthAorta Eq 13 non-planar Fourier displacement ─────────────────────
# Constants from Bošnjak et al. 2025, Table II (fitted to the base patient
# centreline they share publicly). We apply the full Eq 13 form, with γ_1,
# γ_2 fixed at 1.0 (paper sensitivity analysis identifies them as
# non-influential). δ_3, δ_4 are exposed as user-tunable.
_SYNTHAORTA_A1 = -0.798
_SYNTHAORTA_B1 = -0.453
_SYNTHAORTA_A2 = 1.517
_SYNTHAORTA_B2 = 2.699
_SYNTHAORTA_W = 0.027  # series frequency [1/mm]


def apply_arch_tilt(points: list[Vector], tilt_deg: float,
                     pivot_index: int) -> list[Vector]:
    """Rotate arch + descending segments around the inlet z-axis by tilt_deg.

    The ascending segment (points[:pivot_index]) stays unchanged on the
    z-axis. Every point from pivot_index onward (arch + descending) is
    rotated around the z-axis through the ascending-top pivot point so
    the tube's continuity at the join is preserved.

    tilt_deg = 0 → no rotation (backwards-compat).
    tilt_deg > 0 → arch tilts toward +y.
    """
    if tilt_deg == 0.0:
        return points

    angle = math.radians(tilt_deg)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    pivot = points[pivot_index]  # last ascending point (on the z-axis)

    out: list[Vector] = []
    for i, p in enumerate(points):
        if i < pivot_index:
            out.append(p)
            continue
        dx = p.x - pivot.x
        dy = p.y - pivot.y
        rx = dx * cos_a - dy * sin_a
        ry = dx * sin_a + dy * cos_a
        out.append(Vector((pivot.x + rx, pivot.y + ry, p.z)))
    return out


def apply_arch_twist(points: list[Vector], twist_deg: float,
                      n_asc: int, n_desc: int) -> list[Vector]:
    """Gradual z-axis rotation that ramps linearly from 0 across the
    arch region, then is held at ``twist_deg`` through descending.

    Difference from ``apply_arch_tilt``:
      - ``apply_arch_tilt`` applies a single CONSTANT angle to every
        arch+descending point — arch stays planar (in a rotated plane).
      - ``apply_arch_twist`` applies a VARYING angle that increases
        along the arch — arch becomes a non-planar 3D curve. The
        descending tube ends up offset in the same direction as
        arch_tilt would, but the arch path between is helical.

    Both transforms can be applied together; they compose
    multiplicatively in the rotation group.
    """
    if twist_deg == 0.0:
        return points

    n_total = len(points)
    arch_first = n_asc                  # first non-ascending centreline point
    arch_last = n_total - n_desc - 1   # last non-descending centreline point
    if arch_last < arch_first:
        return points

    twist_rad = math.radians(twist_deg)
    denom = max(1, arch_last - arch_first)

    out: list[Vector] = []
    for i, p in enumerate(points):
        if i < arch_first:
            out.append(p)            # ascending — unchanged
            continue
        if i > arch_last:
            angle = twist_rad        # descending — full twist
        else:
            progress = (i - arch_first) / denom
            angle = twist_rad * progress
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        # Rotate around global z-axis (which passes through the inlet at origin)
        rx = p.x * cos_a - p.y * sin_a
        ry = p.x * sin_a + p.y * cos_a
        out.append(Vector((rx, ry, p.z)))
    return out


def apply_nonplanar_displacement(points: list[Vector],
                                  delta_3: float, delta_4: float) -> list[Vector]:
    """Apply SynthAorta Eq 13 non-planar Fourier displacement along the y-axis.

    The planar centreline (built by ``build_centreline``) lies entirely in
    the xz-plane. This adds a y-component:

        d_y(x) = A_1 · cos(w·||x||) + B_1 · sin(w·||x||)
               + A_2 · δ_3 · cos(2w·||x||) + B_2 · δ_4 · sin(2w·||x||)

    where ||x|| is the Euclidean norm of the point's position vector
    (per the paper). When ``delta_3 == 0 and delta_4 == 0``, no
    displacement is applied (backwards-compatible planar behaviour).

    Note: with δ_3=δ_4=1, this reproduces the SynthAorta nominal
    non-planar shape (a few mm of y-wobble along the centreline).
    """
    if delta_3 == 0.0 and delta_4 == 0.0:
        return points

    A1, B1 = _SYNTHAORTA_A1, _SYNTHAORTA_B1
    A2, B2 = _SYNTHAORTA_A2, _SYNTHAORTA_B2
    w = _SYNTHAORTA_W

    out: list[Vector] = []
    for p in points:
        s = p.length  # ||x||
        ws = w * s
        d_y = (A1 * math.cos(ws) + B1 * math.sin(ws)
               + A2 * delta_3 * math.cos(2.0 * ws)
               + B2 * delta_4 * math.sin(2.0 * ws))
        out.append(Vector((p.x, p.y + d_y, p.z)))
    return out


def _cubic_bezier(p0: Vector, p1: Vector, p2: Vector, p3: Vector, t: float) -> Vector:
    u = 1.0 - t
    return (u * u * u) * p0 + (3.0 * u * u * t) * p1 + (3.0 * u * t * t) * p2 + (t * t * t) * p3


def build_centreline(asc_len: float, R_c: float, angle_deg: float, desc_len: float,
                     curve_samples: int, junction_blend_mm: float = 12.0,
                     ) -> tuple[list[Vector], list[float], dict]:
    """Build the centreline as 5 sub-segments with cubic-Bezier corner blends.

    Composition (in arc-length order)::

        [Ascending straight] [Blend1] [Circular arc] [Blend2] [Descending straight]
                              ↑ width = blend/2 each side of the original sharp junction

    The straight ascending portion is shortened by blend/2; the arc starts
    at φ = blend/(2·R_c) instead of φ = 0. A cubic Bezier with G¹-continuous
    tangents at both endpoints fills the blend zone. Same construction at
    the arch→descending join. The middle of the arc remains *exactly* on
    the circular path (radius R_c), so the geometric arch is still the
    designed circular arc — only the small junction regions are smoothed.

    junction_blend_mm = 0 reverts to the original sharp ascending/arc/desc
    composition (backwards-compat for users who want exact corners).

    Returns (points, arc_lengths, derived_geometry).
    """
    arch_geom = derive_arch_geometry(R_c, angle_deg)
    theta = arch_geom["theta_rad"]
    arch_len_arc = R_c * theta  # exact circular arc length

    # Clamp the blend width to fit inside every segment.
    blend = max(0.0, float(junction_blend_mm))
    blend = min(blend, asc_len * 0.5, desc_len * 0.5, arch_len_arc * 0.5)
    half = blend * 0.5

    # Convenience: phi range of the blend at the arc end (half-blend in arc length)
    phi_blend = (half / R_c) if R_c > 0 else 0.0

    # Endpoints + tangents
    p_inlet = Vector((0.0, 0.0, 0.0))
    p_asc_top_full = Vector((0.0, 0.0, asc_len))          # original (sharp-corner) top
    asc_tangent = Vector((0.0, 0.0, 1.0))                 # ascending heads +z

    # Where the *straight* ascending now ends (start of Blend1):
    p_asc_blend_start = Vector((0.0, 0.0, asc_len - half))
    # Where the *circular arc* now starts (end of Blend1):
    p_arc_start = Vector(
        (R_c * (1.0 - math.cos(phi_blend)), 0.0,
         asc_len + R_c * math.sin(phi_blend))
    )
    arc_start_tangent = Vector(
        (math.sin(phi_blend), 0.0, math.cos(phi_blend))
    )

    # Where the *circular arc* now ends (start of Blend2):
    phi_arc_end = theta - phi_blend
    p_arc_end = Vector(
        (R_c * (1.0 - math.cos(phi_arc_end)), 0.0,
         asc_len + R_c * math.sin(phi_arc_end))
    )
    arc_end_tangent = Vector(
        (math.sin(phi_arc_end), 0.0, math.cos(phi_arc_end))
    )

    # Full arch endpoint (where Blend2 lands on the descending tangent line):
    arch_end_full = Vector(
        (R_c * (1.0 - math.cos(theta)), 0.0,
         asc_len + R_c * math.sin(theta))
    )
    desc_natural_tangent = Vector(
        (math.sin(theta), 0.0, math.cos(theta))
    )
    # Descending direction (continuation of arch tangent; flip if pointing up)
    desc_dir = Vector(
        (desc_natural_tangent.x, desc_natural_tangent.y, desc_natural_tangent.z)
    )
    if desc_dir.z > 0:
        desc_dir = -desc_dir
    p_desc_blend_end = arch_end_full + desc_dir * half  # start of descending straight
    p_outlet = arch_end_full + desc_dir * desc_len

    # Sample-count allocation across the 5 sub-segments by approximate length.
    seg_lens = {
        "asc": max(0.0, asc_len - half),       # straight ascending (shortened)
        "b1":  blend,                          # Blend1 chord ≈ blend
        "arc": max(0.0, arch_len_arc - blend), # circular arc (shortened on both ends)
        "b2":  blend,                          # Blend2 chord ≈ blend
        "desc": max(0.0, desc_len - half),     # straight descending (shortened)
    }
    total_eff = sum(seg_lens.values())
    def alloc(name: str, floor: int) -> int:
        return max(floor, int(round(curve_samples * seg_lens[name] / total_eff))) \
               if total_eff > 0 else floor
    n_asc = alloc("asc", 2)
    n_b1  = alloc("b1", 6 if blend > 0 else 0)
    n_arc = alloc("arc", 4)
    n_b2  = alloc("b2", 6 if blend > 0 else 0)
    n_desc = alloc("desc", 2)

    points: list[Vector] = []
    arc_s: list[float] = []
    s = 0.0  # cumulative arc-length so far

    # 1) Ascending straight: p_inlet → p_asc_blend_start
    for i in range(n_asc):
        t = i / (n_asc - 1) if n_asc > 1 else 0.0
        points.append(p_inlet.lerp(p_asc_blend_start, t))
        arc_s.append(seg_lens["asc"] * t)
    s = seg_lens["asc"]

    # 2) Blend1: cubic Bezier p_asc_blend_start → p_arc_start
    if blend > 0 and n_b1 > 0:
        chord = (p_arc_start - p_asc_blend_start).length
        L = chord / 3.0
        h0 = p_asc_blend_start + L * asc_tangent
        h1 = p_arc_start - L * arc_start_tangent
        for i in range(1, n_b1):  # skip i=0 (duplicate of last ascending point)
            t = i / (n_b1 - 1)
            points.append(_cubic_bezier(p_asc_blend_start, h0, h1, p_arc_start, t))
            arc_s.append(s + blend * t)
        s += blend
    else:
        # No blend zone: p_asc_blend_start == p_asc_top_full and p_arc_start
        # is also at the ascending top — the points coincide.
        pass

    # 3) Circular arc: p_arc_start → p_arc_end (φ from phi_blend to phi_arc_end)
    arc_extent = phi_arc_end - phi_blend
    for i in range(1, n_arc):  # skip i=0 (duplicate of Blend1 end)
        t = i / (n_arc - 1)
        phi = phi_blend + arc_extent * t
        pt = Vector(
            (R_c * (1.0 - math.cos(phi)), 0.0,
             asc_len + R_c * math.sin(phi))
        )
        points.append(pt)
        arc_s.append(s + (arc_extent * R_c) * t)
    s += arc_extent * R_c  # = arch_len_arc - blend

    # 4) Blend2: cubic Bezier p_arc_end → p_desc_blend_end
    if blend > 0 and n_b2 > 0:
        chord = (p_desc_blend_end - p_arc_end).length
        L = chord / 3.0
        h0 = p_arc_end + L * arc_end_tangent
        # Use the descending direction for the second handle (G¹ at the
        # descending side). For θ in [90°, 270°] this matches the natural
        # arc-end tangent; for shallow θ the defensive flip in desc_dir
        # would create a discontinuity, but our working range is fine.
        h1 = p_desc_blend_end - L * desc_dir
        for i in range(1, n_b2):
            t = i / (n_b2 - 1)
            points.append(_cubic_bezier(p_arc_end, h0, h1, p_desc_blend_end, t))
            arc_s.append(s + blend * t)
        s += blend
    else:
        pass

    # 5) Descending straight: p_desc_blend_end → p_outlet
    for i in range(1, n_desc + 1):
        t = i / n_desc
        points.append(p_desc_blend_end + desc_dir * (seg_lens["desc"] * t))
        arc_s.append(s + seg_lens["desc"] * t)

    # Recompute arc-length from actual point-to-point distances
    arc_s_exact = [0.0]
    for i in range(1, len(points)):
        arc_s_exact.append(arc_s_exact[-1] + (points[i] - points[i - 1]).length)
    total_arc_exact = arc_s_exact[-1]

    arch_geom["arch_len"] = arch_len_arc
    arch_geom["total_arc"] = total_arc_exact
    arch_geom["n_asc"] = n_asc
    arch_geom["n_arc"] = n_arc
    arch_geom["n_desc"] = n_desc
    arch_geom["n_blend1"] = n_b1
    arch_geom["n_blend2"] = n_b2
    arch_geom["junction_blend_mm"] = blend
    arch_geom["desc_dir"] = (desc_dir.x, desc_dir.y, desc_dir.z)
    arch_geom["outlet_xyz"] = (points[-1].x, points[-1].y, points[-1].z)
    return points, arc_s_exact, arch_geom


# -----------------------------------------------------------------------------
# RMF tube construction (adapted from blender_aorta_like_generator.py:235-323)
# -----------------------------------------------------------------------------


def choose_perpendicular(tangent: Vector) -> Vector:
    tangent = tangent.normalized()
    for ref in (Vector((0, 0, 1)), Vector((0, 1, 0)), Vector((1, 0, 0))):
        if abs(tangent.dot(ref)) < 0.95:
            n = ref - tangent * tangent.dot(ref)
            if n.length > EPS:
                return n.normalized()
    return tangent.orthogonal().normalized()


def compute_frames(centres: Sequence[Vector]) -> tuple[list[Vector], list[Vector], list[Vector]]:
    """Rotation-minimising frames via parallel transport (Wang 2008 simplification)."""
    n = len(centres)
    tangents: list[Vector] = []
    for i in range(n):
        if i == 0:
            t = centres[1] - centres[0]
        elif i == n - 1:
            t = centres[-1] - centres[-2]
        else:
            t = centres[i + 1] - centres[i - 1]
        if t.length < EPS:
            t = Vector((1, 0, 0))
        tangents.append(t.normalized())

    normals: list[Vector] = [choose_perpendicular(tangents[0])]
    binormals: list[Vector] = [tangents[0].cross(normals[0]).normalized()]
    for i in range(1, n):
        t = tangents[i]
        prev_n = normals[-1]
        nv = prev_n - t * prev_n.dot(t)
        if nv.length < EPS:
            nv = choose_perpendicular(t)
        else:
            nv.normalize()
        bv = t.cross(nv)
        if bv.length < EPS:
            bv = tangents[i - 1].cross(nv)
        bv.normalize()
        normals.append(nv)
        binormals.append(bv)
    return tangents, normals, binormals


def build_tube_mesh(name: str, centres: Sequence[Vector], radii: Sequence[float],
                    segments_radial: int = 64) -> bpy.types.Object:
    if len(centres) != len(radii):
        raise ValueError("centres and radii must have same length")
    if len(centres) < 2:
        raise ValueError("need at least two centre points")

    _, normals, binormals = compute_frames(centres)

    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)

    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []
    rings: list[list[int]] = []

    for i, centre in enumerate(centres):
        ring = []
        r = radii[i]
        nv = normals[i]
        bv = binormals[i]
        for j in range(segments_radial):
            theta = 2.0 * pi * j / segments_radial
            p = centre + r * (math.cos(theta) * nv + math.sin(theta) * bv)
            ring.append(len(vertices))
            vertices.append((p.x, p.y, p.z))
        rings.append(ring)

    for i in range(len(centres) - 1):
        a = rings[i]
        b = rings[i + 1]
        for j in range(segments_radial):
            jn = (j + 1) % segments_radial
            faces.append((a[j], a[jn], b[jn], b[j]))

    # Cap fans
    start_idx = len(vertices)
    vertices.append((centres[0].x, centres[0].y, centres[0].z))
    end_idx = len(vertices)
    vertices.append((centres[-1].x, centres[-1].y, centres[-1].z))
    ring0 = rings[0]
    ring1 = rings[-1]
    for j in range(segments_radial):
        jn = (j + 1) % segments_radial
        # Inlet (s=0): wind so normals face outward (-tangent direction at s=0)
        faces.append((start_idx, ring0[jn], ring0[j]))
        # Outlet (s=1): wind so normals face outward (+tangent direction at s=1)
        faces.append((end_idx, ring1[j], ring1[jn]))

    mesh.from_pydata(vertices, [], faces)
    mesh.update(calc_edges=True)
    return obj


# -----------------------------------------------------------------------------
# Blender scene helpers
# -----------------------------------------------------------------------------


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    for block_collection in (bpy.data.meshes, bpy.data.materials, bpy.data.curves):
        for block in list(block_collection):
            block_collection.remove(block)


def recalc_normals(obj: bpy.types.Object) -> None:
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.mesh.select_all(action="SELECT")
    bpy.ops.mesh.normals_make_consistent(inside=False)
    bpy.ops.object.mode_set(mode="OBJECT")


def triangulate(obj: bpy.types.Object) -> None:
    bpy.context.view_layer.objects.active = obj
    mod = obj.modifiers.new(name="Triangulate", type="TRIANGULATE")
    bpy.ops.object.modifier_apply(modifier=mod.name)


def export_stl(obj: bpy.types.Object, output_path: str, do_triangulate: bool,
               save_blend: bool) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    bpy.context.view_layer.objects.active = obj
    if do_triangulate:
        triangulate(obj)
    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    try:
        bpy.ops.wm.stl_export(filepath=output_path, export_selected_objects=True)
    except AttributeError:
        bpy.ops.export_mesh.stl(filepath=output_path, use_selection=True, ascii=False)

    if save_blend:
        blend_path = os.path.splitext(output_path)[0] + ".blend"
        bpy.ops.wm.save_as_mainfile(filepath=blend_path)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> int:
    args = parse_args()

    # Validate
    if args.arch_angle_deg < 30 or args.arch_angle_deg > 270:
        raise SystemExit(f"arch_angle_deg {args.arch_angle_deg} out of geometric "
                         f"sanity range [30, 270] deg")
    if args.arch_R_c <= 0:
        raise SystemExit(f"arch_R_c must be > 0, got {args.arch_R_c}")
    if min(args.r_ascending, args.r_arch, args.r_descending) <= 0:
        raise SystemExit("All segment radii must be > 0")

    # Arch sanity vs ascending length: peak should not exceed asc top by too much
    # for the descending segment to reach below the inlet plane.
    arch_geom_pre = derive_arch_geometry(args.arch_R_c, args.arch_angle_deg)
    if args.descending_length < 1.0:
        raise SystemExit("descending_length must be ≥ 1 mm")

    clean_scene()

    centres, arc_s, arch_geom = build_centreline(
        asc_len=args.ascending_length,
        R_c=args.arch_R_c,
        angle_deg=args.arch_angle_deg,
        desc_len=args.descending_length,
        curve_samples=args.curve_samples,
        junction_blend_mm=args.junction_blend_mm,
    )

    # Rotate arch+descending around inlet z-axis by arch_tilt_deg (rigid).
    centres = apply_arch_tilt(centres, args.arch_tilt_deg,
                              pivot_index=arch_geom["n_asc"] - 1)

    # Gradually twist along the arch by arch_twist_deg (non-planar 3D curve).
    centres = apply_arch_twist(centres, args.arch_twist_deg,
                                n_asc=arch_geom["n_asc"],
                                n_desc=arch_geom["n_desc"])

    # SynthAorta non-planar Fourier displacement (skipped when both δ_3=δ_4=0)
    centres = apply_nonplanar_displacement(centres, args.delta_3, args.delta_4)

    total_arc = arch_geom["total_arc"]
    arch_len = arch_geom["arch_len"]
    radii = [
        radius_along_arc(
            s_arc=s,
            total_arc=total_arc,
            asc_len=args.ascending_length,
            arch_len=arch_len,
            r_asc=args.r_ascending,
            r_arch=args.r_arch,
            r_desc=args.r_descending,
            taper_mode=args.taper_mode,
        )
        for s in arc_s
    ]

    obj = build_tube_mesh(
        name="wall_aorta",
        centres=centres,
        radii=radii,
        segments_radial=args.segments_radial,
    )
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    recalc_normals(obj)

    export_stl(
        obj=obj,
        output_path=args.output,
        do_triangulate=args.triangulate,
        save_blend=args.save_blend,
    )

    if args.metadata:
        # Sidecar JSON — keys compatible with split_patches.py expectations
        # so the v1 splitter can be reused with branches=[]. We also include
        # explicit inlet/outlet positions + normals so the splitter can use
        # them directly (necessary for non-180° arches where the descending
        # column is angled and the heuristic outlet-by-x fails).
        #
        # IMPORTANT: when non-planar Fourier (δ_3/δ_4 != 0) is applied, the
        # actual cap positions in the STL differ from the planar predictions
        # in arch_geom. Use the final post-displacement centreline endpoints
        # for the cap positions.
        nonplanar_on = (args.delta_3 != 0.0 or args.delta_4 != 0.0)
        inlet_pt = centres[0]
        outlet_pt = centres[-1]
        outlet_xyz = (outlet_pt.x, outlet_pt.y, outlet_pt.z)
        inlet_xyz = (inlet_pt.x, inlet_pt.y, inlet_pt.z)
        # Apply the same z-rotation to the planar desc_dir that we applied to
        # the points: arch_tilt_deg (rigid) plus full arch_twist_deg (gradual
        # — by the time we reach descending, the twist has accumulated to
        # the full value).
        planar_desc_dir = arch_geom["desc_dir"]
        total_z_rot_deg = args.arch_tilt_deg + args.arch_twist_deg
        if total_z_rot_deg != 0.0:
            t = math.radians(total_z_rot_deg)
            cos_t, sin_t = math.cos(t), math.sin(t)
            dx, dy, dz = planar_desc_dir
            desc_dir = (dx * cos_t - dy * sin_t, dx * sin_t + dy * cos_t, dz)
        else:
            desc_dir = planar_desc_dir
        meta = {
            "schema_version": "2.0",
            "generator": "blender_aorta_v2",
            "geometry": "healthy_arch_v2",
            "args": {
                "r_ascending": args.r_ascending,
                "r_arch": args.r_arch,
                "r_descending": args.r_descending,
                "taper_mode": args.taper_mode,
                "ascending_length": args.ascending_length,
                "arch_R_c": args.arch_R_c,
                "arch_angle_deg": args.arch_angle_deg,
                "arch_tilt_deg": args.arch_tilt_deg,
                "arch_twist_deg": args.arch_twist_deg,
                "junction_blend_mm": args.junction_blend_mm,
                "descending_length": args.descending_length,
                "delta_3": args.delta_3,
                "delta_4": args.delta_4,
                "segments_radial": args.segments_radial,
                "curve_samples": args.curve_samples,
            },
            "derived": {
                "arch_span_mm": arch_geom["span"],
                "arch_peak_height_above_asc_top_mm": arch_geom["peak_dz"],
                "arch_end_height_above_asc_top_mm": arch_geom["end_dz"],
                "arch_arc_length_mm": arch_geom["arch_len"],
                "total_centreline_length_mm": arch_geom["total_arc"],
                "outlet_xyz_mm": list(outlet_xyz),
                "nonplanar_active": nonplanar_on,
            },
            # ─── Keys consumed by split_patches.split_stl ────────────────────
            # The splitter sizes its cap-search radii from main_diameter_mm and
            # uses arch_span_mm to position the outlet centre.
            "main_diameter_mm": 2.0 * max(args.r_ascending, args.r_arch, args.r_descending),
            "branch_diameter_mm": 0.0,
            "ascending_length_mm": args.ascending_length,
            "descending_length_mm": args.descending_length,
            "arch_span_mm": arch_geom["span"],
            "branches": [],  # explicit empty list → no outletN beyond outlet1
            # Explicit cap positions + normals (v2 addition; used by splitter
            # when present). Outlet normal points OUT of the geometry along
            # the descending tangent.
            "inlet_xyz_mm": list(inlet_xyz),
            "inlet_normal_xyz": [0.0, 0.0, -1.0],
            "outlet_xyz_mm": list(outlet_xyz),
            "outlet_normal_xyz": list(desc_dir),
        }
        sidecar = os.path.splitext(args.output)[0] + ".json"
        with open(sidecar, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
