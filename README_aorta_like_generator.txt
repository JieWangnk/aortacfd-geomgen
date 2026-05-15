Aorta-like Blender generator for snappyHexMesh / CoA benchmark study
===================================================================

Main script
-----------
/mnt/data/blender_aorta_like_generator.py

Purpose
-------
Generate controlled aorta-like benchmark geometries for mesh-stability studies:
- arch_simple
- arch_branched
- arch_coarctation
- arch_branched_coarctation
- rough_arch_branched_coarctation
- aorta_suite

Recommended first full export
-----------------------------
blender -b -P /mnt/data/blender_aorta_like_generator.py -- \
  --geometry aorta_suite \
  --output_dir /tmp/aorta_suite \
  --metadata --triangulate

Recommended main benchmark case
-------------------------------
blender -b -P /mnt/data/blender_aorta_like_generator.py -- \
  --geometry arch_branched_coarctation \
  --diameter 24 \
  --ascending_length 45 \
  --arch_span 70 \
  --arch_height 35 \
  --descending_length 80 \
  --branch_count 3 \
  --branch_diameter_ratio 0.42 \
  --branch_length 35 \
  --coarctation_area_reduction 0.65 \
  --coarctation_length 16 \
  --coarctation_centre_fraction 0.74 \
  --output /tmp/arch_branched_coa.stl \
  --metadata --triangulate

Roughened realism case
----------------------
blender -b -P /mnt/data/blender_aorta_like_generator.py -- \
  --geometry rough_arch_branched_coarctation \
  --diameter 24 \
  --noise_amplitude 0.35 \
  --noise_scale 0.08 \
  --output /tmp/rough_arch_branched_coa.stl \
  --metadata --triangulate

Hypoplastic arch variant
------------------------
blender -b -P /mnt/data/blender_aorta_like_generator.py -- \
  --geometry arch_branched_coarctation \
  --diameter 24 \
  --proximal_hypoplasia 0.15 \
  --coarctation_area_reduction 0.65 \
  --coarctation_length 18 \
  --coarctation_centre_fraction 0.74 \
  --output /tmp/hypoplastic_arch_coa.stl \
  --metadata --triangulate

Suggested study roles
---------------------
A01_arch_simple
  Control curvature with no branch take-off and no narrowing.

A02_arch_branched
  Tests whether supra-aortic branches alone create vulnerable refinement zones.

A03_arch_coarctation
  Isolates narrowing without branch junctions.

A04/A05 arch_branched_coarctation
  Main mechanistic cases: curvature + branches + isthmus narrowing.

A06 hypoplastic_coarctation
  Mimics congenital arch underdevelopment plus focal narrowing.

A07 rough_arch_branched_coarctation
  Realism check for segmentation-like surface irregularity.

A08 arch_one_branch_coarctation
  Reduced-complexity arch for transition tracing and debugging.

Suggested study order
---------------------
1. A01 vs A02
2. A03 vs A04
3. A04 mild vs A05 severe
4. A04 vs A06
5. A04 vs A07
6. Sweep cpd / level on A04 as the main aorta-like benchmark

Notes
-----
- Dimensions are in mm.
- Outputs are closed lumen solids suitable for STL export.
- For snappyHexMesh testing, keep a consistent cpd definition based on main lumen diameter.
- The suite is intended for controlled benchmark work, not patient-anatomical fidelity.
