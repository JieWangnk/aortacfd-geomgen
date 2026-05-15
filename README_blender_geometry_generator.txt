Blender geometry generator for benchmark CFD shapes
==================================================

Main script
-----------
/mnt/data/blender_geometry_generator.py

Run inside Blender
------------------
blender -b -P /mnt/data/blender_geometry_generator.py -- [options]

Examples
--------

1) Straight tube
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry straight_tube \
  --diameter 5 \
  --length 80 \
  --output /tmp/straight_tube.stl \
  --metadata

2) Clean U-bend
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry u_bend \
  --diameter 5 \
  --bend_radius 25 \
  --inlet_length 30 \
  --outlet_length 50 \
  --output /tmp/u_bend.stl \
  --metadata

3) Rough U-bend
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry rough_u_bend \
  --diameter 5 \
  --bend_radius 25 \
  --inlet_length 30 \
  --outlet_length 50 \
  --noise_amplitude 0.08 \
  --noise_scale 0.15 \
  --output /tmp/rough_u_bend.stl \
  --metadata

4) T-junction (equal branch)
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry t_junction \
  --diameter 5 \
  --branch_diameter_ratio 1.0 \
  --inlet_length 30 \
  --outlet_length 30 \
  --branch_length 30 \
  --output /tmp/t_junction.stl \
  --triangulate

5) T-junction (small branch)
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry t_junction \
  --diameter 5 \
  --branch_diameter_ratio 0.4 \
  --inlet_length 30 \
  --outlet_length 30 \
  --branch_length 30 \
  --output /tmp/t_junction_small_branch.stl \
  --triangulate

6) Y-bifurcation
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry y_bifurcation \
  --diameter 5 \
  --branch_diameter_ratio 0.8 \
  --branch_angle 30 \
  --inlet_length 30 \
  --branch_length 35 \
  --output /tmp/y_bifurcation.stl \
  --triangulate

7) Smooth stenosis
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry stenosis \
  --diameter 5 \
  --length 80 \
  --stenosis_length 15 \
  --stenosis_area_reduction 0.5 \
  --output /tmp/stenosis.stl \
  --metadata

8) Full benchmark suite
blender -b -P /mnt/data/blender_geometry_generator.py -- \
  --geometry suite \
  --diameter 5 \
  --output_dir /tmp/benchmark_suite \
  --metadata \
  --triangulate

Notes
-----
- All dimensions are in millimetres.
- The script exports capped lumen-style solids suitable for STL export.
- For T/Y junctions, triangulation is recommended.
- Optional cleanup is available via:
    --voxel_remesh <mm>
    --smooth_iterations <N>
- To save a .blend file next to the export, add:
    --save_blend
