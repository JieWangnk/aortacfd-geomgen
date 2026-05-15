#!/usr/bin/env python3
"""
Blender benchmark geometry generator for CFD meshing studies.

Generates capped, watertight lumen-style geometries suitable for exporting to STL
for snappyHexMesh benchmark work. The script is designed to run inside Blender:

    blender -b -P blender_geometry_generator.py -- --geometry u_bend --output /tmp/u_bend.stl

Supported geometry types
------------------------
- straight_tube
- u_bend
- rough_u_bend
- t_junction
- y_bifurcation
- stenosis
- suite   (writes a small benchmark set into --output_dir)

All dimensions are interpreted in millimetres.
"""

import os
import sys
import json
import math
import argparse
from math import pi, cos, sin

import bpy
import bmesh
from mathutils import Vector, noise

EPS = 1.0e-9


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------

def parse_args():
    argv = sys.argv
    if "--" in argv:
        argv = argv[argv.index("--") + 1 :]
    else:
        argv = []

    parser = argparse.ArgumentParser(description="Generate benchmark CFD geometries in Blender.")
    parser.add_argument("--geometry", default="u_bend",
                        choices=["straight_tube", "u_bend", "rough_u_bend", "t_junction",
                                 "y_bifurcation", "stenosis", "suite"])

    # Core dimensions
    parser.add_argument("--diameter", type=float, default=5.0, help="Main tube diameter [mm].")
    parser.add_argument("--length", type=float, default=60.0, help="Straight tube total length [mm].")
    parser.add_argument("--inlet_length", type=float, default=30.0, help="Straight inlet extension [mm].")
    parser.add_argument("--outlet_length", type=float, default=30.0, help="Straight outlet extension [mm].")
    parser.add_argument("--bend_radius", type=float, default=25.0, help="U-bend centreline radius [mm].")
    parser.add_argument("--bend_angle", type=float, default=180.0, help="Bend angle [deg].")

    # Branching / stenosis controls
    parser.add_argument("--branch_length", type=float, default=30.0, help="Branch extension [mm].")
    parser.add_argument("--branch_angle", type=float, default=30.0, help="Y-branch angle between daughters [deg].")
    parser.add_argument("--branch_diameter_ratio", type=float, default=1.0,
                        help="Branch diameter / main diameter.")
    parser.add_argument("--stenosis_length", type=float, default=15.0,
                        help="Length of smooth stenotic segment [mm].")
    parser.add_argument("--stenosis_area_reduction", type=float, default=0.5,
                        help="Area reduction fraction at throat, between 0 and 0.95.")

    # Mesh / sampling controls
    parser.add_argument("--segments_radial", type=int, default=48, help="Circumferential resolution.")
    parser.add_argument("--segments_axial", type=int, default=120, help="Axial sampling for swept shapes.")
    parser.add_argument("--arc_segments", type=int, default=64, help="Arc sampling for bends.")
    parser.add_argument("--cylinder_vertices", type=int, default=64, help="Cylinder vertices for junctions.")

    # Roughness / cleanup
    parser.add_argument("--noise_amplitude", type=float, default=0.10,
                        help="Radial roughness amplitude [mm] for rough_u_bend.")
    parser.add_argument("--noise_scale", type=float, default=0.15,
                        help="Noise spatial scale for rough_u_bend.")
    parser.add_argument("--voxel_remesh", type=float, default=0.0,
                        help="Optional voxel size [mm] for remesh cleanup, 0 disables.")
    parser.add_argument("--smooth_iterations", type=int, default=0,
                        help="Optional smoothing iterations after generation.")
    parser.add_argument("--triangulate", action="store_true", help="Triangulate before export.")

    # Output / metadata
    parser.add_argument("--output", default="", help="Output file path, usually .stl.")
    parser.add_argument("--output_dir", default="", help="Directory used by --geometry suite.")
    parser.add_argument("--metadata", action="store_true", help="Write JSON metadata next to export.")
    parser.add_argument("--save_blend", action="store_true", help="Also save a .blend file next to the export.")
    return parser.parse_args(argv)


def ensure_dir(path):
    if path:
        os.makedirs(path, exist_ok=True)


def clean_scene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)

    # Remove orphan data blocks for cleaner repeated runs.
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


def choose_perpendicular(tangent):
    tangent = tangent.normalized()
    refs = [Vector((0, 0, 1)), Vector((0, 1, 0)), Vector((1, 0, 0))]
    for ref in refs:
        if abs(tangent.dot(ref)) < 0.95:
            n = ref - tangent * tangent.dot(ref)
            if n.length > EPS:
                return n.normalized()
    # Fallback
    return tangent.orthogonal().normalized()


# -----------------------------------------------------------------------------
# Swept-tube mesh generator (good for straight, bends, stenosis)
# -----------------------------------------------------------------------------

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


def build_tube_mesh(name, centres, radii, segments_radial=48, cap_ends=True):
    if len(centres) != len(radii):
        raise ValueError("centres and radii must have the same length")
    if len(centres) < 2:
        raise ValueError("at least two centre points are required")

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

    n_rings = len(centres)
    for i in range(n_rings - 1):
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

    # Store a lightweight manifest for later processing.
    obj["ring_count"] = len(ring_vertex_indices)
    obj["segments_radial"] = segments_radial
    return obj, ring_vertex_indices


# -----------------------------------------------------------------------------
# Primitive / boolean helpers (good for T/Y junctions)
# -----------------------------------------------------------------------------

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
        raise ValueError("No objects provided for boolean union")
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


# -----------------------------------------------------------------------------
# Centreline builders
# -----------------------------------------------------------------------------

def line_points(p0, p1, n):
    p0 = Vector(p0)
    p1 = Vector(p1)
    if n < 2:
        return [p0, p1]
    return [p0.lerp(p1, i / (n - 1)) for i in range(n)]


def build_straight_centres(length, n):
    return line_points((0.0, 0.0, 0.0), (length, 0.0, 0.0), n)


def build_u_bend_centres(bend_radius, inlet_length, outlet_length, bend_angle_deg, axial_samples, arc_samples):
    angle = math.radians(bend_angle_deg)
    if bend_radius <= 0:
        raise ValueError("bend_radius must be positive")
    if angle <= 0 or angle > 2 * pi:
        raise ValueError("bend_angle must be in (0, 360]")

    # Use a 180-degree U-bend layout in the XZ plane by default, with the circular arc
    # centred at the origin. General angle support preserves the same construction logic.
    start_angle = pi
    end_angle = max(pi - angle, -pi)

    arc_pts = []
    for i in range(max(arc_samples, 3)):
        t = i / (max(arc_samples, 3) - 1)
        a = start_angle + t * (end_angle - start_angle)
        arc_pts.append(Vector((bend_radius * cos(a), 0.0, bend_radius * sin(a))))

    inlet_end = arc_pts[0]
    outlet_start = arc_pts[-1]

    # Tangent directions at each end of the arc.
    inlet_tangent = (arc_pts[1] - arc_pts[0]).normalized()
    outlet_tangent = (arc_pts[-1] - arc_pts[-2]).normalized()

    inlet_start = inlet_end - inlet_tangent * inlet_length
    outlet_end = outlet_start + outlet_tangent * outlet_length

    inlet_pts = line_points(inlet_start, inlet_end, max(2, axial_samples // 4))[:-1]
    outlet_pts = line_points(outlet_start, outlet_end, max(2, axial_samples // 4))[1:]
    centres = inlet_pts + arc_pts + outlet_pts
    return centres


def build_stenosis_centres_and_radii(length, diameter, stenosis_length, area_reduction, axial_samples):
    area_reduction = min(max(area_reduction, 0.0), 0.95)
    r0 = 0.5 * diameter
    r_throat = r0 * math.sqrt(max(1.0 - area_reduction, 1.0e-6))
    xs = [length * i / (axial_samples - 1) for i in range(axial_samples)]
    x0 = 0.5 * (length - stenosis_length)
    x1 = 0.5 * (length + stenosis_length)
    centres = []
    radii = []
    for x in xs:
        centres.append(Vector((x, 0.0, 0.0)))
        if x <= x0 or x >= x1:
            radii.append(r0)
        else:
            s = (x - x0) / max(stenosis_length, EPS)
            constriction = 0.5 * (1.0 - math.cos(2.0 * pi * s))
            radii.append(r0 - (r0 - r_throat) * constriction)
    return centres, radii


# -----------------------------------------------------------------------------
# Geometry generators
# -----------------------------------------------------------------------------

def generate_straight_tube(args):
    centres = build_straight_centres(args.length, max(args.segments_axial, 12))
    radii = [0.5 * args.diameter] * len(centres)
    obj, _ = build_tube_mesh("straight_tube", centres, radii, segments_radial=args.segments_radial)
    return obj, {"type": "straight_tube", "diameter_mm": args.diameter, "length_mm": args.length}


def generate_u_bend(args, rough=False):
    centres = build_u_bend_centres(
        bend_radius=args.bend_radius,
        inlet_length=args.inlet_length,
        outlet_length=args.outlet_length,
        bend_angle_deg=args.bend_angle,
        axial_samples=max(args.segments_axial, 24),
        arc_samples=max(args.arc_segments, 12),
    )
    radii = [0.5 * args.diameter] * len(centres)
    name = "rough_u_bend" if rough else "u_bend"
    obj, rings = build_tube_mesh(name, centres, radii, segments_radial=args.segments_radial)

    if rough:
        for i in range(1, len(rings) - 1):
            centre = centres[i]
            for vidx in rings[i]:
                v = obj.data.vertices[vidx]
                radial = (v.co - centre)
                if radial.length < EPS:
                    continue
                scalar = noise.noise(v.co * args.noise_scale)
                v.co += radial.normalized() * (args.noise_amplitude * scalar)
        obj.data.update()
        recalc_normals(obj)

    meta = {
        "type": name,
        "diameter_mm": args.diameter,
        "bend_radius_mm": args.bend_radius,
        "bend_angle_deg": args.bend_angle,
        "inlet_length_mm": args.inlet_length,
        "outlet_length_mm": args.outlet_length,
    }
    if rough:
        meta.update({
            "noise_amplitude_mm": args.noise_amplitude,
            "noise_scale": args.noise_scale,
        })
    return obj, meta


def generate_t_junction(args):
    r_main = 0.5 * args.diameter
    r_branch = r_main * args.branch_diameter_ratio

    main_start = Vector((-args.inlet_length, 0.0, 0.0))
    main_end = Vector((args.outlet_length, 0.0, 0.0))
    branch_start = Vector((0.0, 0.0, 0.0))
    branch_end = Vector((0.0, 0.0, args.branch_length))

    main = add_cylinder_between(main_start, main_end, r_main, "main", vertices=args.cylinder_vertices)
    branch = add_cylinder_between(branch_start, branch_end, r_branch, "branch", vertices=args.cylinder_vertices)
    obj = boolean_union([main, branch], name="t_junction")

    meta = {
        "type": "t_junction",
        "main_diameter_mm": args.diameter,
        "branch_diameter_mm": args.diameter * args.branch_diameter_ratio,
        "inlet_length_mm": args.inlet_length,
        "outlet_length_mm": args.outlet_length,
        "branch_length_mm": args.branch_length,
        "junction_angle_deg": 90.0,
    }
    return obj, meta


def generate_y_bifurcation(args):
    r_main = 0.5 * args.diameter
    r_branch = r_main * args.branch_diameter_ratio
    half_angle = math.radians(args.branch_angle * 0.5)

    parent_start = Vector((-args.inlet_length, 0.0, 0.0))
    junction = Vector((0.0, 0.0, 0.0))
    branch_dir_1 = Vector((math.cos(half_angle), 0.0, math.sin(half_angle))).normalized()
    branch_dir_2 = Vector((math.cos(half_angle), 0.0, -math.sin(half_angle))).normalized()
    daughter_end_1 = junction + args.branch_length * branch_dir_1
    daughter_end_2 = junction + args.branch_length * branch_dir_2

    parent = add_cylinder_between(parent_start, junction, r_main, "parent", vertices=args.cylinder_vertices)
    d1 = add_cylinder_between(junction, daughter_end_1, r_branch, "daughter_1", vertices=args.cylinder_vertices)
    d2 = add_cylinder_between(junction, daughter_end_2, r_branch, "daughter_2", vertices=args.cylinder_vertices)
    obj = boolean_union([parent, d1, d2], name="y_bifurcation")

    meta = {
        "type": "y_bifurcation",
        "main_diameter_mm": args.diameter,
        "branch_diameter_mm": args.diameter * args.branch_diameter_ratio,
        "inlet_length_mm": args.inlet_length,
        "branch_length_mm": args.branch_length,
        "branch_angle_deg": args.branch_angle,
    }
    return obj, meta


def generate_stenosis(args):
    centres, radii = build_stenosis_centres_and_radii(
        length=args.length,
        diameter=args.diameter,
        stenosis_length=args.stenosis_length,
        area_reduction=args.stenosis_area_reduction,
        axial_samples=max(args.segments_axial, 24),
    )
    obj, _ = build_tube_mesh("stenosis", centres, radii, segments_radial=args.segments_radial)
    meta = {
        "type": "stenosis",
        "diameter_mm": args.diameter,
        "length_mm": args.length,
        "stenosis_length_mm": args.stenosis_length,
        "stenosis_area_reduction": args.stenosis_area_reduction,
    }
    return obj, meta


# -----------------------------------------------------------------------------
# Export / reporting
# -----------------------------------------------------------------------------

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
# Suite generation
# -----------------------------------------------------------------------------

def suite_specs(args):
    d = args.diameter
    return [
        {
            "name": "straight_tube",
            "geometry": "straight_tube",
            "diameter": d,
            "length": max(args.length, 12.0 * d),
        },
        {
            "name": "u_bend_gentle",
            "geometry": "u_bend",
            "diameter": d,
            "bend_radius": max(10.0 * d, args.bend_radius),
            "inlet_length": max(5.0 * d, args.inlet_length),
            "outlet_length": max(10.0 * d, args.outlet_length),
        },
        {
            "name": "u_bend_medium",
            "geometry": "u_bend",
            "diameter": d,
            "bend_radius": max(5.0 * d, 0.5 * args.bend_radius),
            "inlet_length": max(5.0 * d, args.inlet_length),
            "outlet_length": max(10.0 * d, args.outlet_length),
        },
        {
            "name": "u_bend_tight",
            "geometry": "u_bend",
            "diameter": d,
            "bend_radius": max(2.0 * d, 0.25 * args.bend_radius),
            "inlet_length": max(5.0 * d, args.inlet_length),
            "outlet_length": max(10.0 * d, args.outlet_length),
        },
        {
            "name": "rough_u_bend",
            "geometry": "rough_u_bend",
            "diameter": d,
            "bend_radius": max(5.0 * d, 0.5 * args.bend_radius),
            "inlet_length": max(5.0 * d, args.inlet_length),
            "outlet_length": max(10.0 * d, args.outlet_length),
            "noise_amplitude": max(args.noise_amplitude, 0.02 * d),
            "noise_scale": args.noise_scale,
        },
        {
            "name": "t_junction_equal",
            "geometry": "t_junction",
            "diameter": d,
            "branch_diameter_ratio": 1.0,
            "inlet_length": max(5.0 * d, args.inlet_length),
            "outlet_length": max(5.0 * d, args.outlet_length),
            "branch_length": max(5.0 * d, args.branch_length),
        },
        {
            "name": "t_junction_small_branch",
            "geometry": "t_junction",
            "diameter": d,
            "branch_diameter_ratio": 0.4,
            "inlet_length": max(5.0 * d, args.inlet_length),
            "outlet_length": max(5.0 * d, args.outlet_length),
            "branch_length": max(5.0 * d, args.branch_length),
        },
        {
            "name": "y_bifurcation",
            "geometry": "y_bifurcation",
            "diameter": d,
            "branch_diameter_ratio": 0.8,
            "branch_angle": args.branch_angle,
            "inlet_length": max(5.0 * d, args.inlet_length),
            "branch_length": max(5.0 * d, args.branch_length),
        },
        {
            "name": "stenosis",
            "geometry": "stenosis",
            "diameter": d,
            "length": max(args.length, 14.0 * d),
            "stenosis_length": max(args.stenosis_length, 2.0 * d),
            "stenosis_area_reduction": args.stenosis_area_reduction,
        },
    ]


def override_args(base_args, spec):
    namespace = argparse.Namespace(**vars(base_args))
    for k, v in spec.items():
        if k != "name":
            setattr(namespace, k, v)
    return namespace


def generate_one(args):
    geom = args.geometry
    if geom == "straight_tube":
        return generate_straight_tube(args)
    if geom == "u_bend":
        return generate_u_bend(args, rough=False)
    if geom == "rough_u_bend":
        return generate_u_bend(args, rough=True)
    if geom == "t_junction":
        return generate_t_junction(args)
    if geom == "y_bifurcation":
        return generate_y_bifurcation(args)
    if geom == "stenosis":
        return generate_stenosis(args)
    raise ValueError(f"Unsupported geometry: {geom}")


def run_suite(args):
    out_dir = args.output_dir or os.path.join(os.getcwd(), "geometry_suite")
    ensure_dir(out_dir)
    manifest = []

    for spec in suite_specs(args):
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
    if args.geometry == "suite":
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
