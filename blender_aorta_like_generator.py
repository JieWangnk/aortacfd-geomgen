#!/usr/bin/env python3
"""
Blender aorta-like benchmark geometry generator for snappyHexMesh stability studies.

Designed for controlled cardiovascular CFD benchmark generation, especially
arch / coarctation / branch-takeoff problems where refinement transitions may
interact with curvature and junctions.

Run inside Blender, for example:

  blender -b -P blender_aorta_like_generator.py -- \
    --geometry aorta_suite \
    --output_dir /tmp/aorta_suite \
    --metadata --triangulate

All dimensions are in millimetres.
"""

import os
import sys
import json
import math
import argparse
from math import pi

import bpy
from mathutils import Vector, noise

EPS = 1.0e-9


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args():
    argv = sys.argv
    if "--" in argv:
        argv = argv[argv.index("--") + 1 :]
    else:
        argv = []

    parser = argparse.ArgumentParser(description="Generate aorta-like CFD benchmark geometries in Blender.")
    parser.add_argument(
        "--geometry",
        default="arch_branched_coarctation",
        choices=[
            "arch_simple",
            "arch_branched",
            "arch_coarctation",
            "arch_branched_coarctation",
            "rough_arch_branched_coarctation",
            "aorta_suite",
        ],
    )

    # Global main-tube dimensions
    parser.add_argument("--diameter", type=float, default=24.0, help="Nominal main lumen diameter [mm].")
    parser.add_argument("--ascending_length", type=float, default=45.0, help="Ascending aorta length [mm].")
    parser.add_argument("--arch_span", type=float, default=70.0, help="Arch span from ascending to descending [mm].")
    parser.add_argument("--arch_height", type=float, default=35.0, help="Arch rise above the ascending top [mm].")
    parser.add_argument("--descending_length", type=float, default=80.0, help="Descending aorta length [mm].")

    # Branch configuration
    parser.add_argument("--branch_count", type=int, default=3, choices=[1, 2, 3], help="Number of supra-aortic branches.")
    parser.add_argument("--branch_diameter_ratio", type=float, default=0.42, help="Branch diameter / main diameter.")
    parser.add_argument("--branch_length", type=float, default=35.0, help="Branch outlet extension length [mm].")
    parser.add_argument("--branch_spacing", type=float, default=14.0, help="Approximate spacing of branch origins along arch centreline [mm].")
    parser.add_argument("--branch_tilt_deg", type=float, default=60.0, help="Branch take-off tilt from +x direction toward +z [deg].")
    parser.add_argument("--branch_splay_deg", type=float, default=22.0, help="Out-of-plane branch splay magnitude [deg].")

    # Coarctation controls
    parser.add_argument("--coarctation", action="store_true", help="Enable isthmus coarctation for single-geometry runs.")
    parser.add_argument("--coarctation_area_reduction", type=float, default=0.65,
                        help="Area reduction fraction at coarctation throat, between 0 and 0.95.")
    parser.add_argument("--coarctation_length", type=float, default=30.0, help="Axial length of the smooth coarctation [mm].")
    parser.add_argument("--coarctation_centre_fraction", type=float, default=0.72,
                        help="Centre location of coarctation along main centreline arc length [0-1].")
    parser.add_argument("--proximal_hypoplasia", type=float, default=0.0,
                        help="Optional smooth proximal arch diameter reduction fraction over distal arch segment [0-0.5].")

    # Surface roughness (to mimic imperfect segmentation / rough STL)
    parser.add_argument("--roughness", action="store_true", help="Enable mild radial surface roughness.")
    parser.add_argument("--noise_amplitude", type=float, default=0.35, help="Roughness amplitude [mm].")
    parser.add_argument("--noise_scale", type=float, default=0.08, help="Noise spatial scale.")

    # Mesh/sampling controls
    parser.add_argument("--segments_radial", type=int, default=64, help="Circumferential ring resolution.")
    parser.add_argument("--curve_samples", type=int, default=220, help="Main centreline sample count.")
    parser.add_argument("--cylinder_vertices", type=int, default=64, help="Cylinder vertices for branches.")

    # Cleanup/export
    parser.add_argument("--voxel_remesh", type=float, default=0.0, help="Optional voxel remesh size [mm].")
    parser.add_argument("--smooth_iterations", type=int, default=0, help="Optional smoothing iterations.")
    parser.add_argument("--triangulate", action="store_true", help="Triangulate before export.")
    parser.add_argument("--output", default="", help="Output path (.stl or .obj) for single geometry.")
    parser.add_argument("--output_dir", default="", help="Output directory for suite mode.")
    parser.add_argument("--metadata", action="store_true", help="Write JSON metadata sidecars.")
    parser.add_argument("--save_blend", action="store_true", help="Save .blend next to exports.")
    return parser.parse_args(argv)


# -----------------------------------------------------------------------------
# Scene / export utilities
# -----------------------------------------------------------------------------

def ensure_dir(path):
    if path:
        os.makedirs(path, exist_ok=True)


def clean_scene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    for data_block in (bpy.data.meshes, bpy.data.curves, bpy.data.materials,
                       bpy.data.textures, bpy.data.images, bpy.data.objects):
        for block in list(data_block):
            if getattr(block, "users", 0) == 0:
                data_block.remove(block)


def set_active(obj):
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj


def apply_object_transforms(obj):
    set_active(obj)
    bpy.ops.object.transform_apply(location=False, rotation=True, scale=True)


def recalc_normals(obj):
    set_active(obj)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.normals_make_consistent(inside=False)
    bpy.ops.object.mode_set(mode='OBJECT')


def shade_smooth(obj):
    set_active(obj)
    bpy.ops.object.shade_smooth()


def optional_triangulate(obj):
    mod = obj.modifiers.new(name="Triangulate", type='TRIANGULATE')
    set_active(obj)
    bpy.ops.object.modifier_apply(modifier=mod.name)


def optional_voxel_remesh(obj, voxel_size):
    if voxel_size <= 0:
        return obj
    mod = obj.modifiers.new(name="VoxelRemesh", type='REMESH')
    mod.mode = 'VOXEL'
    mod.voxel_size = voxel_size
    mod.use_smooth_shade = True
    mod.adaptivity = 0.0
    set_active(obj)
    bpy.ops.object.modifier_apply(modifier=mod.name)
    return obj


def optional_smooth(obj, iterations):
    if iterations <= 0:
        return obj
    mod = obj.modifiers.new(name="Smooth", type='SMOOTH')
    mod.iterations = iterations
    mod.factor = 0.5
    set_active(obj)
    bpy.ops.object.modifier_apply(modifier=mod.name)
    return obj


def export_object(obj, output_path, triangulate=False, save_blend=False):
    output_path = os.path.abspath(output_path)
    ensure_dir(os.path.dirname(output_path))

    if triangulate:
        optional_triangulate(obj)

    set_active(obj)
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)

    ext = os.path.splitext(output_path)[1].lower()
    if ext == ".stl":
        try:
            bpy.ops.export_mesh.stl(filepath=output_path, use_selection=True, ascii=False)
        except Exception:
            bpy.ops.wm.stl_export(filepath=output_path, export_selected_objects=True)
    elif ext == ".obj":
        try:
            bpy.ops.export_scene.obj(filepath=output_path, use_selection=True, use_mesh_modifiers=True)
        except Exception:
            bpy.ops.wm.obj_export(filepath=output_path, export_selected_objects=True)
    else:
        raise ValueError("Only .stl and .obj export are currently supported")

    if save_blend:
        blend_path = os.path.splitext(output_path)[0] + ".blend"
        bpy.ops.wm.save_as_mainfile(filepath=blend_path)


def maybe_write_metadata(output_path, metadata):
    sidecar = os.path.splitext(output_path)[0] + ".json"
    with open(sidecar, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


def finalize_object(obj, args):
    optional_voxel_remesh(obj, args.voxel_remesh)
    optional_smooth(obj, args.smooth_iterations)
    recalc_normals(obj)
    shade_smooth(obj)
    return obj


# -----------------------------------------------------------------------------
# Geometry helpers
# -----------------------------------------------------------------------------

def choose_perpendicular(tangent):
    tangent = tangent.normalized()
    refs = [Vector((0, 0, 1)), Vector((0, 1, 0)), Vector((1, 0, 0))]
    for ref in refs:
        if abs(tangent.dot(ref)) < 0.95:
            n = ref - tangent * tangent.dot(ref)
            if n.length > EPS:
                return n.normalized()
    return tangent.orthogonal().normalized()


def compute_tangents(centres):
    tangents = []
    n = len(centres)
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
    return tangents


def compute_frames(centres):
    tangents = compute_tangents(centres)
    normals = [choose_perpendicular(tangents[0])]
    binormals = [tangents[0].cross(normals[0]).normalized()]
    for i in range(1, len(centres)):
        t = tangents[i]
        prev_n = normals[-1]
        n = prev_n - t * prev_n.dot(t)
        if n.length < EPS:
            n = choose_perpendicular(t)
        else:
            n.normalize()
        b = t.cross(n)
        if b.length < EPS:
            b = tangents[i - 1].cross(n)
        b.normalize()
        normals.append(n)
        binormals.append(b)
    return tangents, normals, binormals


def build_tube_mesh(name, centres, radii, segments_radial=64, cap_ends=True):
    if len(centres) != len(radii):
        raise ValueError("centres and radii must have same length")
    if len(centres) < 2:
        raise ValueError("need at least two centre points")

    _, normals, binormals = compute_frames(centres)

    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)

    vertices = []
    faces = []
    ring_vertex_indices = []

    for i, centre in enumerate(centres):
        ring = []
        r = radii[i]
        n = normals[i]
        b = binormals[i]
        for j in range(segments_radial):
            theta = 2.0 * pi * j / segments_radial
            p = centre + r * (math.cos(theta) * n + math.sin(theta) * b)
            ring.append(len(vertices))
            vertices.append((p.x, p.y, p.z))
        ring_vertex_indices.append(ring)

    for i in range(len(centres) - 1):
        ring_a = ring_vertex_indices[i]
        ring_b = ring_vertex_indices[i + 1]
        for j in range(segments_radial):
            jn = (j + 1) % segments_radial
            faces.append((ring_a[j], ring_a[jn], ring_b[jn], ring_b[j]))

    if cap_ends:
        start_idx = len(vertices)
        vertices.append(tuple(centres[0]))
        end_idx = len(vertices)
        vertices.append(tuple(centres[-1]))
        ring0 = ring_vertex_indices[0]
        ring1 = ring_vertex_indices[-1]
        for j in range(segments_radial):
            jn = (j + 1) % segments_radial
            faces.append((start_idx, ring0[jn], ring0[j]))
            faces.append((end_idx, ring1[j], ring1[jn]))

    mesh.from_pydata(vertices, [], faces)
    mesh.update(calc_edges=True)
    obj["ring_count"] = len(ring_vertex_indices)
    obj["segments_radial"] = segments_radial
    return obj, ring_vertex_indices


def add_cylinder_between(p0, p1, radius, name, vertices=64):
    p0 = Vector(p0)
    p1 = Vector(p1)
    direction = p1 - p0
    length = direction.length
    if length < EPS:
        raise ValueError(f"Cylinder '{name}' has zero length")

    bpy.ops.mesh.primitive_cylinder_add(
        vertices=vertices,
        radius=radius,
        depth=length,
        end_fill_type='NGON',
        location=(0, 0, 0),
    )
    obj = bpy.context.object
    obj.name = name
    obj.location = (p0 + p1) * 0.5
    obj.rotation_mode = 'QUATERNION'
    obj.rotation_quaternion = direction.to_track_quat('Z', 'Y')
    apply_object_transforms(obj)
    return obj


def boolean_union(objs, name="union"):
    if not objs:
        raise ValueError("No objects for boolean union")
    base = objs[0]
    set_active(base)
    for other in objs[1:]:
        mod = base.modifiers.new(name=f"Boolean_{other.name}", type='BOOLEAN')
        mod.operation = 'UNION'
        mod.solver = 'EXACT'
        mod.object = other
        bpy.ops.object.modifier_apply(modifier=mod.name)
        if other.name in bpy.data.objects:
            bpy.data.objects.remove(other, do_unlink=True)
    base.name = name
    recalc_normals(base)
    return base


def cumulative_lengths(points):
    s = [0.0]
    for i in range(1, len(points)):
        s.append(s[-1] + (points[i] - points[i - 1]).length)
    return s


def cubic_bezier(p0, p1, p2, p3, t):
    u = 1.0 - t
    return (u**3) * p0 + 3.0 * (u**2) * t * p1 + 3.0 * u * (t**2) * p2 + (t**3) * p3


def line_points(p0, p1, n, include_end=True):
    p0 = Vector(p0)
    p1 = Vector(p1)
    if n < 2:
        pts = [p0, p1]
    else:
        pts = [p0.lerp(p1, i / (n - 1)) for i in range(n)]
    return pts if include_end else pts[:-1]


def bezier_points(p0, p1, p2, p3, n, include_end=True):
    pts = [cubic_bezier(Vector(p0), Vector(p1), Vector(p2), Vector(p3), i / (n - 1)) for i in range(n)]
    return pts if include_end else pts[:-1]


# -----------------------------------------------------------------------------
# Aorta-like centreline and radius models
# -----------------------------------------------------------------------------

def build_arch_centres(ascending_length, arch_span, arch_height, descending_length, curve_samples):
    # Piecewise path: ascending straight -> smooth arch crown -> smooth
    # arch-to-descending transition -> descending straight.
    # Layout in the XZ plane; y = 0 for the main lumen.
    #
    # The arch bezier now includes a tangent-continuous transition into
    # the descending segment via a short second bezier that smoothly
    # turns the curve from horizontal to vertical, eliminating the
    # sharp crease at the old p2 junction.
    n1 = max(12, int(curve_samples * 0.15))
    n2 = max(36, int(curve_samples * 0.50))
    n3 = max(18, int(curve_samples * 0.20))
    n4 = max(8, curve_samples - n1 - n2 - n3 + 3)

    p0 = Vector((0.0, 0.0, 0.0))
    p1 = Vector((0.0, 0.0, ascending_length))

    # Arch crown (bezier from ascending top to apex)
    p2 = Vector((arch_span, 0.0, ascending_length + 0.20 * arch_height))
    # G^1 at p1: h1 directly above p1, so the arch leaves p1 along +z (matches the ascending line).
    h1 = Vector((0.0, 0.0, ascending_length + arch_height))
    # G^1 at p2: h2 directly above p2, so the arch arrives at p2 along -z.
    h2 = Vector((arch_span, 0.0, ascending_length + arch_height))

    # Smooth transition from arch into descending (second bezier).
    # All four control points stay at x=arch_span: the transition descends purely along -z.
    transition_drop = min(0.30 * descending_length, 40.0)
    p3 = Vector((arch_span, 0.0, ascending_length - transition_drop))
    # h3: directly below p2, same x — G^1 at p2 (parallel to p2-h2).
    h3 = Vector((arch_span, 0.0, ascending_length + 0.10 * arch_height))
    # h4: directly above p3, same x — G^1 at p3 (parallel to p4-p3).
    h4 = Vector((arch_span, 0.0, ascending_length - 0.30 * transition_drop))

    # Descending straight
    p4 = Vector((arch_span, 0.0, ascending_length - descending_length))

    seg1 = line_points(p0, p1, n1, include_end=False)
    seg2 = bezier_points(p1, h1, h2, p2, n2, include_end=False)
    seg3 = bezier_points(p2, h3, h4, p3, n3, include_end=False)
    seg4 = line_points(p3, p4, n4, include_end=True)
    centres = seg1 + seg2 + seg3 + seg4
    return centres


def arch_radii(centres, diameter, coarctation=False, coarctation_area_reduction=0.65,
               coarctation_length=30.0, coarctation_centre_fraction=0.72,
               proximal_hypoplasia=0.0):
    base_r = 0.5 * diameter
    s = cumulative_lengths(centres)
    total = max(s[-1], EPS)
    centre_s = min(max(coarctation_centre_fraction, 0.0), 1.0) * total

    # Proximal arch hypoplasia applied smoothly over distal arch / isthmus region.
    hyp_frac = max(0.0, min(proximal_hypoplasia, 0.5))

    throat_r = base_r * math.sqrt(max(1.0 - min(max(coarctation_area_reduction, 0.0), 0.95), 1.0e-6))
    half_len = 0.5 * max(coarctation_length, 1.0)

    radii = []
    for si in s:
        r = base_r

        # Optional smooth distal arch hypoplasia before the discrete coarctation.
        if hyp_frac > 0.0:
            a = 0.52 * total
            b = 0.80 * total
            if a < si < b:
                x = (si - a) / (b - a)
                smooth = 0.5 * (1.0 - math.cos(pi * x))
                r *= (1.0 - hyp_frac * smooth)
            elif si >= b:
                r *= (1.0 - hyp_frac)

        # Smooth coarctation taper using a raised-cosine-squared profile.
        # cos^2 gives a gentler shoulder than the old cos(2πx) which had
        # abrupt inflection points at the start and end of the narrowing.
        if coarctation and abs(si - centre_s) <= half_len:
            x = (si - (centre_s - half_len)) / max(2.0 * half_len, EPS)
            constriction = math.sin(pi * x) ** 2
            local_base = r
            r = local_base - (local_base - throat_r) * constriction

        radii.append(r)
    return radii


def interpolate_point_by_fraction(centres, fraction):
    fraction = min(max(fraction, 0.0), 1.0)
    s = cumulative_lengths(centres)
    total = max(s[-1], EPS)
    target = fraction * total
    for i in range(1, len(centres)):
        if s[i] >= target:
            denom = max(s[i] - s[i - 1], EPS)
            t = (target - s[i - 1]) / denom
            return centres[i - 1].lerp(centres[i], t), (centres[i] - centres[i - 1]).normalized()
    return centres[-1].copy(), (centres[-1] - centres[-2]).normalized()


def branch_direction_from_tangent(tangent, tilt_deg, splay_deg):
    # Make the branch primarily upward (+z) and slightly out-of-plane in y.
    tangent = tangent.normalized()
    up = Vector((0.0, 0.0, 1.0))
    side = tangent.cross(up)
    if side.length < EPS:
        side = Vector((0.0, 1.0, 0.0))
    side.normalize()

    tilt = math.radians(tilt_deg)
    splay = math.radians(splay_deg)
    # Remove tangent component from up to avoid shooting back into the arch.
    up_proj = up - tangent * up.dot(tangent)
    if up_proj.length < EPS:
        up_proj = Vector((0.0, 1.0, 0.0))
    up_proj.normalize()

    d = math.cos(tilt) * tangent + math.sin(tilt) * up_proj
    d = (math.cos(splay) * d + math.sin(splay) * side).normalized()
    return d


def apply_ring_roughness(obj, centres, rings, amplitude, scale):
    for i in range(1, len(rings) - 1):
        centre = centres[i]
        for vidx in rings[i]:
            v = obj.data.vertices[vidx]
            radial = (v.co - centre)
            if radial.length < EPS:
                continue
            scalar = noise.noise(v.co * scale)
            v.co += radial.normalized() * (amplitude * scalar)
    obj.data.update()
    recalc_normals(obj)


# -----------------------------------------------------------------------------
# Generators
# -----------------------------------------------------------------------------

def generate_arch_main(args, with_coarctation=False, rough=False):
    centres = build_arch_centres(
        ascending_length=args.ascending_length,
        arch_span=args.arch_span,
        arch_height=args.arch_height,
        descending_length=args.descending_length,
        curve_samples=max(args.curve_samples, 80),
    )
    radii = arch_radii(
        centres,
        diameter=args.diameter,
        coarctation=with_coarctation,
        coarctation_area_reduction=args.coarctation_area_reduction,
        coarctation_length=args.coarctation_length,
        coarctation_centre_fraction=args.coarctation_centre_fraction,
        proximal_hypoplasia=args.proximal_hypoplasia,
    )
    name = "arch_main"
    obj, rings = build_tube_mesh(name, centres, radii, segments_radial=args.segments_radial)
    if rough:
        apply_ring_roughness(obj, centres, rings, args.noise_amplitude, args.noise_scale)
    return obj, centres, radii


def add_branches_to_arch(main_obj, centres, args):
    branch_radius = 0.5 * args.diameter * args.branch_diameter_ratio
    branch_specs = []

    # Branch origins on the arch crown (fractions of total arc length).
    # All three sit on the ascending-to-apex portion of the arch bezier,
    # before the arch-to-descending transition begins.
    fractions = {1: [0.38], 2: [0.34, 0.44], 3: [0.30, 0.38, 0.46]}[args.branch_count]
    splay_signs = {1: [0], 2: [-1, 1], 3: [-1, 0, 1]}[args.branch_count]

    objs = [main_obj]
    for idx, (f, sign) in enumerate(zip(fractions, splay_signs), start=1):
        origin, tangent = interpolate_point_by_fraction(centres, f)
        direction = branch_direction_from_tangent(
            tangent,
            tilt_deg=args.branch_tilt_deg,
            splay_deg=sign * args.branch_splay_deg,
        )
        end = origin + args.branch_length * direction
        cyl = add_cylinder_between(origin, end, branch_radius, f"branch_{idx}", vertices=args.cylinder_vertices)
        objs.append(cyl)
        branch_specs.append({
            "index": idx,
            "origin_fraction": f,
            "origin_xyz_mm": [origin.x, origin.y, origin.z],
            "end_xyz_mm": [end.x, end.y, end.z],
            "diameter_mm": 2.0 * branch_radius,
        })

    union = boolean_union(objs, name="arch_branched")
    return union, branch_specs


def generate_one(args):
    geom = args.geometry

    if geom == "arch_simple":
        obj, centres, radii = generate_arch_main(args, with_coarctation=False, rough=False)
        meta = {
            "type": geom,
            "main_diameter_mm": args.diameter,
            "ascending_length_mm": args.ascending_length,
            "arch_span_mm": args.arch_span,
            "arch_height_mm": args.arch_height,
            "descending_length_mm": args.descending_length,
            "branch_count": 0,
            "coarctation": False,
        }
        return obj, meta

    if geom == "arch_branched":
        main_obj, centres, radii = generate_arch_main(args, with_coarctation=False, rough=False)
        obj, branch_specs = add_branches_to_arch(main_obj, centres, args)
        meta = {
            "type": geom,
            "main_diameter_mm": args.diameter,
            "ascending_length_mm": args.ascending_length,
            "arch_span_mm": args.arch_span,
            "arch_height_mm": args.arch_height,
            "descending_length_mm": args.descending_length,
            "branch_count": args.branch_count,
            "branch_diameter_mm": args.diameter * args.branch_diameter_ratio,
            "branch_length_mm": args.branch_length,
            "branches": branch_specs,
            "coarctation": False,
        }
        return obj, meta

    if geom == "arch_coarctation":
        obj, centres, radii = generate_arch_main(args, with_coarctation=True, rough=False)
        meta = {
            "type": geom,
            "main_diameter_mm": args.diameter,
            "ascending_length_mm": args.ascending_length,
            "arch_span_mm": args.arch_span,
            "arch_height_mm": args.arch_height,
            "descending_length_mm": args.descending_length,
            "branch_count": 0,
            "coarctation": True,
            "coarctation_area_reduction": args.coarctation_area_reduction,
            "coarctation_length_mm": args.coarctation_length,
            "coarctation_centre_fraction": args.coarctation_centre_fraction,
            "proximal_hypoplasia": args.proximal_hypoplasia,
        }
        return obj, meta

    if geom in ["arch_branched_coarctation", "rough_arch_branched_coarctation"]:
        rough = geom.startswith("rough") or args.roughness
        main_obj, centres, radii = generate_arch_main(args, with_coarctation=True, rough=rough)
        obj, branch_specs = add_branches_to_arch(main_obj, centres, args)
        meta = {
            "type": geom,
            "main_diameter_mm": args.diameter,
            "ascending_length_mm": args.ascending_length,
            "arch_span_mm": args.arch_span,
            "arch_height_mm": args.arch_height,
            "descending_length_mm": args.descending_length,
            "branch_count": args.branch_count,
            "branch_diameter_mm": args.diameter * args.branch_diameter_ratio,
            "branch_length_mm": args.branch_length,
            "branches": branch_specs,
            "coarctation": True,
            "coarctation_area_reduction": args.coarctation_area_reduction,
            "coarctation_length_mm": args.coarctation_length,
            "coarctation_centre_fraction": args.coarctation_centre_fraction,
            "proximal_hypoplasia": args.proximal_hypoplasia,
            "roughness": rough,
            "noise_amplitude_mm": args.noise_amplitude if rough else 0.0,
            "noise_scale": args.noise_scale if rough else 0.0,
        }
        return obj, meta

    raise ValueError(f"Unsupported geometry: {geom}")


# -----------------------------------------------------------------------------
# Suite presets tailored to aorta / CoA study
# -----------------------------------------------------------------------------

def override_args(base_args, spec):
    namespace = argparse.Namespace(**vars(base_args))
    for k, v in spec.items():
        if k != "name":
            setattr(namespace, k, v)
    return namespace


def aorta_suite_specs(args):
    d = args.diameter
    return [
        {
            "name": "A01_arch_simple",
            "geometry": "arch_simple",
            "diameter": d,
            "branch_count": 0,
        },
        {
            "name": "A02_arch_branched",
            "geometry": "arch_branched",
            "diameter": d,
            "branch_count": 3,
        },
        {
            "name": "A03_arch_coarctation",
            "geometry": "arch_coarctation",
            "diameter": d,
            "coarctation": True,
            "coarctation_area_reduction": 0.50,
            "coarctation_length": 14.0,
            "coarctation_centre_fraction": 0.72,
        },
        {
            "name": "A04_arch_branched_coarctation_mild",
            "geometry": "arch_branched_coarctation",
            "diameter": d,
            "branch_count": 3,
            "coarctation": True,
            "coarctation_area_reduction": 0.50,
            "coarctation_length": 14.0,
            "coarctation_centre_fraction": 0.72,
        },
        {
            "name": "A05_arch_branched_coarctation_severe",
            "geometry": "arch_branched_coarctation",
            "diameter": d,
            "branch_count": 3,
            "coarctation": True,
            "coarctation_area_reduction": 0.70,
            "coarctation_length": 16.0,
            "coarctation_centre_fraction": 0.74,
        },
        {
            "name": "A06_arch_branched_hypoplastic_coarctation",
            "geometry": "arch_branched_coarctation",
            "diameter": d,
            "branch_count": 3,
            "coarctation": True,
            "coarctation_area_reduction": 0.65,
            "coarctation_length": 18.0,
            "coarctation_centre_fraction": 0.74,
            "proximal_hypoplasia": 0.15,
        },
        {
            "name": "A07_rough_arch_branched_coarctation",
            "geometry": "rough_arch_branched_coarctation",
            "diameter": d,
            "branch_count": 3,
            "coarctation": True,
            "coarctation_area_reduction": 0.65,
            "coarctation_length": 16.0,
            "coarctation_centre_fraction": 0.74,
            "noise_amplitude": max(args.noise_amplitude, 0.25),
            "noise_scale": args.noise_scale,
        },
        {
            "name": "A08_arch_one_branch_coarctation",
            "geometry": "arch_branched_coarctation",
            "diameter": d,
            "branch_count": 1,
            "coarctation": True,
            "coarctation_area_reduction": 0.65,
            "coarctation_length": 16.0,
            "coarctation_centre_fraction": 0.74,
        },
    ]


def run_suite(args):
    out_dir = args.output_dir or os.path.join(os.getcwd(), "aorta_suite")
    ensure_dir(out_dir)
    manifest = []

    for spec in aorta_suite_specs(args):
        clean_scene()
        local_args = override_args(args, spec)
        obj, meta = generate_one(local_args)
        finalize_object(obj, local_args)
        output_path = os.path.join(out_dir, f"{spec['name']}.stl")
        export_object(obj, output_path, triangulate=local_args.triangulate, save_blend=local_args.save_blend)
        meta["file"] = os.path.basename(output_path)
        if local_args.metadata:
            maybe_write_metadata(output_path, meta)
        manifest.append(meta)

    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"Wrote {len(manifest)} geometries to: {out_dir}")


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    if args.geometry == "aorta_suite":
        run_suite(args)
        return

    clean_scene()
    obj, meta = generate_one(args)
    finalize_object(obj, args)

    if args.output:
        export_object(obj, args.output, triangulate=args.triangulate, save_blend=args.save_blend)
        if args.metadata:
            maybe_write_metadata(args.output, meta)
        print(f"Exported: {os.path.abspath(args.output)}")
    else:
        print("Geometry generated in scene. No file exported because --output was not provided.")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
