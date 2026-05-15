# Vessel Generator Benchmark Harness

Head-to-head geometry comparison between our pure-Python parametric
generators and **VMTK**, the de-facto open-source standard for
centreline-based vascular surface modelling.

The harness fixes the *input geometry* (centreline + per-point radius) so both
tools are fed identical data; differences in output therefore reflect the
algorithms, not the inputs.

## Why VMTK as the reference

VMTK is a serious, mature tool: it operates on patient-derived centrelines
extracted from medical images, builds its surface via maximum-inscribed-sphere
implicit fields, has a full C++/VTK back-end, an active community, and 15+
years of clinical CFD use behind it.  It is the right tool when your starting
point is a 3D scan.

We are **not trying to replace VMTK**.  Our generators serve a different need:
*parametric* synthetic geometry for ML training datasets, where you want
hundreds or thousands of geometries from a closed parameter space, scriptable
end-to-end without a GUI, runnable on any Python install with no compiled
dependencies beyond `pip install trimesh manifold3d`.

The benchmark answers a fair, narrow question: *given the same idealised
centreline and radius profile, how does our generator's surface compare with
VMTK's reference output*?  It does **not** answer "which tool should you use",
because that depends on whether you start from imaging or from parameters.

## Layout

```
benchmarks/
  test_cases/         Canonical input cases (JSON)
    bifurcation_murray.json     Symmetric Y, Murray's law
    aorta_normal.json            16-parameter healthy aortic arch
  ours/
    run_ours.py        Generates STLs via vessel_generator + aorta_generator
                       (bifurcation: SDF + MC; aorta: union + Taubin)
  vmtk/
    run_vmtk.sh        VMTK pipeline (centreline -> level-set -> MC -> STL)
  compare.py           Pairwise STL comparison + JSON report
```

## Quick start

```bash
# 1. Generate ours (already done; outputs in ours/)
cd benchmarks/ours && python3 run_ours.py

# 2. Generate VMTK reference (requires `pip install vmtk` or conda env)
cd ../vmtk && bash run_vmtk.sh

# 3. Compare
cd ..
python3 compare.py ours/bifurcation_murray.stl vmtk/bifurcation_murray.stl \
        --label-a ours --label-b vmtk \
        --json reports/ours_vs_vmtk_bifurcation.json
```

## Metrics reported

For each pair (A, B), `compare.py` reports:

| Category | Metric | What it shows |
|---|---|---|
| Topology | lumen volume, surface area, bbox | Macro-scale agreement |
|  | n boundary edges, n non-manifold | Mesh "health" (CFD-readiness) |
|  | n connected components | Whether the model is one body |
| Quality | aspect ratio min / mean / p95 | Triangle slivers (high = bad) |
|  | n triangles with AR > 5, > 20 | Pathological triangles |
|  | min / mean triangle area | Mesh refinement |
| Geometric agreement | Hausdorff distance | Worst-case surface deviation |
|  | mean nearest distance | Average surface deviation |
|  | p95(A->B), p95(B->A) | Robust deviation |
|  | centreline RMS | Cross-check via re-extracted centrelines |
| Deltas | volume_diff (mm^3, %) | Headline number for reviewers |
|  | area_diff (mm^2, %) | Surface-area discrepancy |

## Honest caveats

1. **Identical inputs are not always realistic.** VMTK is built for
   per-patient centrelines extracted from images; our generator is
   parametric.  A "fair" comparison forces both onto the same idealised
   inputs, which removes VMTK's main strength (working from real patient
   data).  If you want to benchmark patient-realism, you need a different
   protocol that starts from images.

2. **No CFD-output comparison.** The harness compares *geometry*.  A peer
   reviewer will ask for *hemodynamic* comparison: same mesh-density, same
   solver, same boundary conditions, compare WSS / pressure drop / flow
   split.  That's a separate benchmark; this one is the prerequisite.

3. **Mesh quality is sensitive to discretisation.**  The `aspect_ratio_p95`
   and `n_slivers` numbers vary with the chosen MC resolution / NURBS
   sampling / VMTK voxel grid.  Document the parameters you used; do not
   compare across different tool settings.

4. **Centreline RMS uses a slice-based extraction** that lumps multiple
   lumens together at bifurcation points.  Use it as a sanity check, not as
   ground truth.

## What the harness is for

This is the scaffolding that turns "we have a different method" into a
reviewable claim.  The metrics, test cases, and protocol are fixed; running
the comparison once VMTK is installed (or by a collaborator) produces a JSON
report you can drop straight into a paper's supplementary materials.
