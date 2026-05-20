"""Tests for the v2 geometry math: R_c↔(span,height,arc_length) and taper.

These tests run without Blender. The math functions are imported directly
from blender_aorta_v2 — the file contains an `import bpy` at the top that
would normally fail outside Blender, so we test via a try/except shim.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType

import pytest

HERE = Path(__file__).resolve().parent.parent


def _load_v2_math_without_bpy() -> ModuleType:
    """Load blender_aorta_v2.py with a fake `bpy` so the math helpers import.

    The math helpers (derive_arch_geometry, smoothstep, radius_along_arc,
    build_centreline) don't use Blender — they're pure numpy + math. The
    file-level `import bpy` is the only Blender dependency at module load.
    """
    # Stub bpy + mathutils so the import doesn't fail outside Blender
    if "bpy" not in sys.modules:
        bpy_stub = ModuleType("bpy")
        bpy_stub.data = type("D", (), {"meshes": None, "objects": None, "materials": None, "curves": None})()
        bpy_stub.context = type("C", (), {"collection": None, "view_layer": None})()
        bpy_stub.ops = type("O", (), {})()
        bpy_stub.types = type("T", (), {"Object": object})()
        sys.modules["bpy"] = bpy_stub

    if "mathutils" not in sys.modules:
        # Provide a minimal Vector with the operations we use
        import numpy as np

        class _Vec:
            def __init__(self, xyz):
                self._v = np.array(xyz, dtype=float)

            @property
            def x(self):
                return float(self._v[0])

            @property
            def y(self):
                return float(self._v[1])

            @property
            def z(self):
                return float(self._v[2])

            def __add__(self, o):
                return _Vec(self._v + (o._v if isinstance(o, _Vec) else o))

            def __sub__(self, o):
                return _Vec(self._v - (o._v if isinstance(o, _Vec) else o))

            def __mul__(self, s):
                return _Vec(self._v * float(s))

            __rmul__ = __mul__

            def __neg__(self):
                return _Vec(-self._v)

            @property
            def length(self):
                return float(np.linalg.norm(self._v))

            def normalize(self):
                n = self.length
                if n > 0:
                    self._v = self._v / n

            def normalized(self):
                n = self.length
                return _Vec(self._v / n if n > 0 else self._v)

            def dot(self, o):
                return float(np.dot(self._v, o._v))

            def cross(self, o):
                return _Vec(np.cross(self._v, o._v))

            def orthogonal(self):
                # Pick any vector not parallel
                if abs(self._v[0]) < 0.9:
                    other = np.array([1.0, 0.0, 0.0])
                else:
                    other = np.array([0.0, 1.0, 0.0])
                return _Vec(np.cross(self._v, other))

        m = ModuleType("mathutils")
        m.Vector = _Vec
        sys.modules["mathutils"] = m

    spec = importlib.util.spec_from_file_location(
        "blender_aorta_v2", HERE / "blender_aorta_v2.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def v2():
    return _load_v2_math_without_bpy()


# ── derive_arch_geometry ───────────────────────────────────────────────────


def test_arch_geom_180_degrees_canonical_U(v2) -> None:
    out = v2.derive_arch_geometry(R_c=40.0, angle_deg=180.0)
    assert math.isclose(out["span"], 80.0, abs_tol=1e-9)
    assert math.isclose(out["peak_dz"], 40.0, abs_tol=1e-9)
    assert math.isclose(out["end_dz"], 0.0, abs_tol=1e-9)


def test_arch_geom_120_degrees_shallow(v2) -> None:
    out = v2.derive_arch_geometry(R_c=40.0, angle_deg=120.0)
    # span = R_c * (1 - cos 120°) = 40 * 1.5 = 60
    assert math.isclose(out["span"], 60.0, abs_tol=1e-9)
    # peak occurs at φ=90° (since 120° > 90°), peak_dz = R_c
    assert math.isclose(out["peak_dz"], 40.0, abs_tol=1e-9)
    # end_dz = R_c * sin 120° ≈ 34.64
    assert math.isclose(out["end_dz"], 40.0 * math.sin(math.radians(120)), abs_tol=1e-9)


def test_arch_geom_90_degrees_quarter_arc(v2) -> None:
    out = v2.derive_arch_geometry(R_c=40.0, angle_deg=90.0)
    assert math.isclose(out["span"], 40.0, abs_tol=1e-9)
    # peak at the end of arc when θ ≤ 90°
    assert math.isclose(out["peak_dz"], 40.0, abs_tol=1e-9)
    assert math.isclose(out["end_dz"], 40.0, abs_tol=1e-9)


def test_arch_geom_200_degrees_over_arched(v2) -> None:
    out = v2.derive_arch_geometry(R_c=40.0, angle_deg=200.0)
    # span = R_c * (1 - cos 200°) = 40 * 1.94 ≈ 77.6
    assert out["span"] > 75.0 and out["span"] < 80.0
    # peak_dz = R_c (since θ > 90°)
    assert math.isclose(out["peak_dz"], 40.0, abs_tol=1e-9)
    # end_dz = R_c * sin 200° < 0 (descending start is BELOW asc top)
    assert out["end_dz"] < 0


# ── smoothstep ─────────────────────────────────────────────────────────────


def test_smoothstep_endpoints(v2) -> None:
    assert v2.smoothstep(0.0) == 0.0
    assert v2.smoothstep(1.0) == 1.0


def test_smoothstep_midpoint_is_half(v2) -> None:
    assert math.isclose(v2.smoothstep(0.5), 0.5, abs_tol=1e-9)


def test_smoothstep_clamps_outside_unit_interval(v2) -> None:
    assert v2.smoothstep(-1.0) == 0.0
    assert v2.smoothstep(2.0) == 1.0


def test_smoothstep_monotonic(v2) -> None:
    import numpy as np
    xs = np.linspace(0, 1, 50)
    ys = [v2.smoothstep(x) for x in xs]
    for a, b in zip(ys, ys[1:]):
        assert b >= a - 1e-12


# ── radius_along_arc ───────────────────────────────────────────────────────


def test_radius_piecewise_returns_segment_values(v2) -> None:
    # asc 0-50, arch 50-150, desc 150-200; r = 14, 13, 12
    r = lambda s: v2.radius_along_arc(s, total_arc=200, asc_len=50, arch_len=100,
                                       r_asc=14.0, r_arch=13.0, r_desc=12.0,
                                       taper_mode="piecewise")
    assert r(10.0) == 14.0  # ascending
    assert r(100.0) == 13.0  # arch
    assert r(180.0) == 12.0  # descending


def test_radius_smoothstep_blends_across_boundaries(v2) -> None:
    r_asc, r_arch, r_desc = 14.0, 13.0, 12.0
    r = lambda s: v2.radius_along_arc(s, total_arc=200, asc_len=50, arch_len=100,
                                       r_asc=r_asc, r_arch=r_arch, r_desc=r_desc,
                                       taper_mode="smoothstep", blend_window=10.0)
    # At s=50 (asc→arch boundary midpoint): should be midway between r_asc and r_arch
    assert math.isclose(r(50.0), 0.5 * (r_asc + r_arch), abs_tol=1e-6)
    # At s=150 (arch→desc boundary midpoint)
    assert math.isclose(r(150.0), 0.5 * (r_arch + r_desc), abs_tol=1e-6)
    # Inside ascending segment (outside blend window)
    assert r(20.0) == r_asc
    # Inside descending segment
    assert r(180.0) == r_desc


def test_radius_smoothstep_monotonic_across_decreasing_segments(v2) -> None:
    """For r_asc > r_arch > r_desc, the smoothstep transition should never increase."""
    r_asc, r_arch, r_desc = 16.0, 13.0, 10.0
    rs = [
        v2.radius_along_arc(s, total_arc=200, asc_len=50, arch_len=100,
                             r_asc=r_asc, r_arch=r_arch, r_desc=r_desc,
                             taper_mode="smoothstep", blend_window=20.0)
        for s in [10, 30, 45, 50, 55, 70, 90, 130, 145, 150, 155, 180]
    ]
    for a, b in zip(rs, rs[1:]):
        assert b <= a + 1e-9, f"non-monotonic: {rs}"


def test_radius_linear_taper_midpoint(v2) -> None:
    r_asc, r_arch = 14.0, 12.0
    out = v2.radius_along_arc(50.0, total_arc=200, asc_len=50, arch_len=100,
                              r_asc=r_asc, r_arch=r_arch, r_desc=11.0,
                              taper_mode="linear", blend_window=10.0)
    assert math.isclose(out, 0.5 * (r_asc + r_arch), abs_tol=1e-6)


# ── build_centreline ───────────────────────────────────────────────────────


def test_build_centreline_180_endpoints(v2) -> None:
    points, arc_s, geom = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=200,
    )
    # First point at origin
    p0 = points[0]
    assert math.isclose(p0.x, 0.0, abs_tol=1e-9)
    assert math.isclose(p0.z, 0.0, abs_tol=1e-9)
    # Last point: arch ends at (2*R_c, 0, asc_len) → descending goes straight down
    # to (2*R_c, 0, asc_len - desc_len) = (80, 0, -150)
    pN = points[-1]
    assert math.isclose(pN.x, 80.0, abs_tol=1e-3)
    assert math.isclose(pN.z, 50.0 - 200.0, abs_tol=1e-3)
    # Total arc length
    expected_total = 50.0 + 40.0 * math.pi + 200.0
    assert math.isclose(geom["total_arc"], expected_total, rel_tol=1e-3)


def test_build_centreline_120_outlet_position(v2) -> None:
    """For 120° arch, descending column is angled — outlet position is non-trivial."""
    points, _, geom = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=120.0, desc_len=200.0, curve_samples=200,
    )
    # Arch end at (R_c*(1-cos 120°), 0, asc_len + R_c sin 120°) = (60, 0, 50+34.64)
    # Tangent direction: (sin 120°, 0, cos 120°) = (0.866, 0, -0.5), going down+forward.
    # Outlet ≈ arch_end + 200 * (0.866, 0, -0.5) = (60+173.2, 0, 84.64-100) ≈ (233, 0, -15.36)
    pN = points[-1]
    assert pN.x > 200.0
    assert pN.z < 0.0


def test_build_centreline_arc_length_matches_R_theta(v2) -> None:
    """Arch arc length should equal R_c * θ (radians)."""
    _, _, geom = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=300,
    )
    expected_arch_len = 40.0 * math.pi
    assert math.isclose(geom["arch_len"], expected_arch_len, rel_tol=1e-6)


# ── apply_nonplanar_displacement ───────────────────────────────────────────


def test_nonplanar_displacement_zero_means_planar(v2) -> None:
    """δ_3 = δ_4 = 0 must produce identical points (backwards-compat invariant)."""
    points, _, _ = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=80,
    )
    out = v2.apply_nonplanar_displacement(points, delta_3=0.0, delta_4=0.0)
    # Same objects returned (early-return path)
    assert out is points or all(p is q for p, q in zip(points, out))
    # All y-coords stay zero
    assert all(abs(p.y) < 1e-12 for p in out)


def test_nonplanar_displacement_nominal_produces_y_wobble(v2) -> None:
    """δ_3 = δ_4 = 1 (SynthAorta nominal) must give non-trivial y-offsets."""
    import numpy as np

    points, _, _ = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=120,
    )
    out = v2.apply_nonplanar_displacement(points, delta_3=1.0, delta_4=1.0)
    ys = np.array([p.y for p in out])
    # Not all zero
    assert (np.abs(ys) > 1e-6).any(), "Non-planar displacement produced zero y everywhere"
    # Should be bounded by sum of paper amplitudes (~|A_1|+|B_1|+|A_2|+|B_2| ≈ 5.5 mm)
    assert ys.max() - ys.min() < 12.0, f"y-range too large: {ys.max() - ys.min()}"
    assert ys.max() - ys.min() > 0.5, f"y-range too small: {ys.max() - ys.min()}"


def test_nonplanar_displacement_preserves_xz(v2) -> None:
    """x and z coordinates must be unchanged by the Fourier displacement."""
    points, _, _ = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=60,
    )
    out = v2.apply_nonplanar_displacement(points, delta_3=1.0, delta_4=1.0)
    for p, q in zip(points, out):
        assert math.isclose(p.x, q.x, abs_tol=1e-12)
        assert math.isclose(p.z, q.z, abs_tol=1e-12)


def test_arch_tilt_zero_is_identity(v2) -> None:
    """arch_tilt_deg=0 must return the input points unchanged."""
    points, _, _ = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=60,
    )
    out = v2.apply_arch_tilt(points, tilt_deg=0.0, pivot_index=10)
    for p, q in zip(points, out):
        assert math.isclose(p.x, q.x, abs_tol=1e-12)
        assert math.isclose(p.y, q.y, abs_tol=1e-12)
        assert math.isclose(p.z, q.z, abs_tol=1e-12)


def test_arch_tilt_rotates_arch_and_descending_only(v2) -> None:
    """Ascending stays on z-axis; arch+descending rotate around z through pivot."""
    points, _, geom = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=80,
    )
    pivot_i = geom["n_asc"] - 1
    out = v2.apply_arch_tilt(points, tilt_deg=15.0, pivot_index=pivot_i)

    # Ascending segment (i < pivot_i) untouched
    for i in range(pivot_i):
        assert math.isclose(points[i].x, out[i].x, abs_tol=1e-12)
        assert math.isclose(points[i].y, out[i].y, abs_tol=1e-12)
        assert math.isclose(points[i].z, out[i].z, abs_tol=1e-12)
    # The pivot point itself is unchanged (rotation centre)
    assert math.isclose(points[pivot_i].x, out[pivot_i].x, abs_tol=1e-12)
    assert math.isclose(points[pivot_i].y, out[pivot_i].y, abs_tol=1e-12)
    # Arch peak (some point after pivot) now has nonzero y (was strictly 0 before)
    post = [out[i] for i in range(pivot_i + 1, len(out))]
    assert any(abs(p.y) > 1e-6 for p in post), "tilt 15° produced zero y everywhere"
    # z unchanged everywhere (rotation is around z-axis)
    for p, q in zip(points, out):
        assert math.isclose(p.z, q.z, abs_tol=1e-12)


def test_arch_tilt_preserves_distances(v2) -> None:
    """A rotation must preserve point-to-point distances within the rotated set."""
    points, _, geom = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=60,
    )
    pivot_i = geom["n_asc"] - 1
    out = v2.apply_arch_tilt(points, tilt_deg=20.0, pivot_index=pivot_i)
    # Pairwise distances in the rotated subset (arch + descending) are unchanged
    for i in range(pivot_i, len(points) - 1):
        d_before = (points[i + 1] - points[i]).length
        d_after = (out[i + 1] - out[i]).length
        assert math.isclose(d_before, d_after, rel_tol=1e-9), \
            f"distance changed at i={i}: {d_before} → {d_after}"


def test_nonplanar_displacement_scales_with_delta(v2) -> None:
    """Doubling delta_4 should ~roughly double its contribution (linear in δ).

    Specifically: y_displacement(δ_3=0, δ_4=2) - y_displacement(δ_3=0, δ_4=1)
                = (A_2·0 + B_2·1·sin) - 0 — only the B_2·δ_4·sin part is left,
    which doubles with δ_4. Also the A_1, B_1 terms are constant (γ_1=γ_2=1),
    so the difference is purely the B_2 second-harmonic contribution.
    """
    import numpy as np

    points, _, _ = v2.build_centreline(
        asc_len=50.0, R_c=40.0, angle_deg=180.0, desc_len=200.0, curve_samples=80,
    )
    y_d1 = np.array([p.y for p in v2.apply_nonplanar_displacement(points, 0.0, 1.0)])
    y_d2 = np.array([p.y for p in v2.apply_nonplanar_displacement(points, 0.0, 2.0)])
    # The second-harmonic difference (d_4=2 minus d_4=1) is exactly B_2·sin(2w·||x||)
    # which has nonzero magnitude — i.e., the two y arrays should NOT be equal.
    assert not np.allclose(y_d1, y_d2)
