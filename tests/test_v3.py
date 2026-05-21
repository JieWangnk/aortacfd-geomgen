"""Tests for cli_v3 — 5-knob minimal interface."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from cli_v3 import (  # noqa: E402
    PARAMETERS,
    V2_FIXED,
    V3_TO_V2,
    _parse_param_override,
    expand_cases,
    translate_v3_to_v2,
    validate_spec,
)


# ── PARAMETERS schema ───────────────────────────────────────────────────────


def test_v3_has_10_parameters() -> None:
    # 5 primary (radii + arch) + twist_deg + arch_shape + arch_radius_mm
    # + 2 optional (lengths)
    assert len(PARAMETERS) == 10


def test_v3_primary_knobs_present() -> None:
    for name in ("r_inlet", "r_outlet", "arch_width_mm", "arch_height_mm", "torsion_deg"):
        assert name in PARAMETERS


def test_v3_no_internal_v2_params_exposed() -> None:
    """v3 deliberately hides taper_mode, junction_blend_mm, delta_*, mesh res, etc."""
    for hidden in ("taper_mode", "junction_blend_mm", "delta_3", "delta_4",
                   "segments_radial", "curve_samples", "arch_R_c", "arch_angle_deg"):
        assert hidden not in PARAMETERS


# ── v3 → v2 translation ─────────────────────────────────────────────────────


def test_translate_baseline_matches_user_reference() -> None:
    """The v3 baseline should translate to v2 params consistent with the
    user-supplied reference (outputs/v2_dim/baseline_v2)."""
    v3 = {
        "r_inlet": 14.0,
        "r_outlet": 10.0,
        "arch_width_mm": 90.0,
        "arch_height_mm": 45.0,
        "torsion_deg": 0.0,
    }
    v2 = translate_v3_to_v2(v3)
    assert v2["r_ascending"] == 14.0
    assert v2["r_descending"] == 10.0
    assert v2["arch_tilt_deg"] == 0.0
    # r_arch is the midpoint of inlet and outlet
    assert v2["r_arch"] == 12.0
    # The user passes (W, H) and translate_v3_to_v2 calls _resolve_arch_params,
    # which converts to (R_c, angle). With H=45, W=90 → R_c=45, θ=180°.
    assert v2["arch_R_c"] == pytest.approx(45.0)
    assert v2["arch_angle_deg"] == pytest.approx(180.0, abs=1e-6)
    # Fixed defaults injected
    assert v2["taper_mode"] == "smoothstep"
    assert v2["delta_3"] == 0.0
    assert v2["delta_4"] == 0.0
    assert v2["junction_blend_mm"] == V2_FIXED["junction_blend_mm"]
    # v3 keys are stripped
    for v3_only in ("arch_width_mm", "arch_height_mm"):
        assert v3_only not in v2


def test_translate_handles_partial_overrides() -> None:
    """Setting just torsion still works — other defaults pulled from V2_FIXED."""
    v3 = {"r_inlet": 13.7, "r_outlet": 12.0, "arch_width_mm": 80.0,
          "arch_height_mm": 40.0, "torsion_deg": 15.0}
    v2 = translate_v3_to_v2(v3)
    assert v2["arch_tilt_deg"] == 15.0
    assert v2["r_arch"] == pytest.approx(0.5 * (13.7 + 12.0))


def test_translate_rejects_unknown_v3_key() -> None:
    with pytest.raises(ValueError, match="Unknown v3 parameter"):
        translate_v3_to_v2({"r_inletx": 14.0})


# ── arch_radius_mm convenience shortcut ─────────────────────────────────────


def test_arch_radius_expands_to_U_arch() -> None:
    """arch_radius_mm=30 → arch_width_mm=60, arch_height_mm=30 (canonical U)."""
    v2 = translate_v3_to_v2({"r_inlet": 14.0, "r_outlet": 10.0,
                              "arch_radius_mm": 30.0})
    # After translation + _resolve_arch_params, the v2 params should reflect
    # arch_R_c=30 and arch_angle_deg=180 (U-arch)
    assert v2["arch_R_c"] == pytest.approx(30.0)
    assert v2["arch_angle_deg"] == pytest.approx(180.0, abs=1e-6)
    # arch_radius_mm itself doesn't end up in v2
    assert "arch_radius_mm" not in v2


def test_arch_radius_zero_ignored() -> None:
    """Default arch_radius_mm=0 means 'not used' — width/height take effect."""
    v2 = translate_v3_to_v2({"r_inlet": 14.0, "r_outlet": 10.0,
                              "arch_radius_mm": 0.0,
                              "arch_width_mm": 90.0, "arch_height_mm": 45.0})
    assert v2["arch_R_c"] == pytest.approx(45.0)
    assert v2["arch_angle_deg"] == pytest.approx(180.0, abs=1e-6)


def test_arch_radius_conflicts_with_width_height() -> None:
    """Setting arch_radius_mm AND arch_width_mm should raise."""
    with pytest.raises(ValueError, match="cannot be combined"):
        translate_v3_to_v2({"arch_radius_mm": 30.0, "arch_width_mm": 60.0})
    with pytest.raises(ValueError, match="cannot be combined"):
        translate_v3_to_v2({"arch_radius_mm": 30.0, "arch_height_mm": 30.0})


def test_arch_radius_in_ellipse_alone_degenerates_to_circle() -> None:
    """arch_radius_mm alone in ellipse mode produces a degenerate ellipse
    (a = b = R) which is geometrically a circle U-arch."""
    v2 = translate_v3_to_v2({"r_inlet": 14.0, "r_outlet": 10.0,
                              "arch_radius_mm": 30.0, "arch_shape": "ellipse"})
    # Ellipse passes through with span+height untouched (no R_c conversion)
    assert v2["arch_shape"] == "ellipse"
    assert v2["arch_span_mm"] == pytest.approx(60.0)   # = 2R
    assert v2["arch_height_mm"] == pytest.approx(30.0)  # = R


def test_arch_radius_in_ellipse_with_width_derives_height() -> None:
    """In ellipse mode: R_peak + W → H = W² / (4·R_peak)."""
    v2 = translate_v3_to_v2({"arch_radius_mm": 30.0, "arch_shape": "ellipse",
                              "arch_width_mm": 60.0})  # R=30, W=60 → H = 3600/120 = 30
    assert v2["arch_span_mm"] == pytest.approx(60.0)
    assert v2["arch_height_mm"] == pytest.approx(30.0)


def test_arch_radius_in_ellipse_with_height_derives_width() -> None:
    """In ellipse mode: R_peak + H → W = 2·√(R_peak · H)."""
    v2 = translate_v3_to_v2({"arch_radius_mm": 20.0, "arch_shape": "ellipse",
                              "arch_height_mm": 80.0})  # R=20, H=80 → W = 2·√(1600) = 80
    assert v2["arch_span_mm"] == pytest.approx(80.0)
    assert v2["arch_height_mm"] == pytest.approx(80.0)


def test_arch_radius_in_ellipse_rejects_overdetermined() -> None:
    """R + W + H all set in ellipse mode is over-determined."""
    with pytest.raises(ValueError, match="over-determined"):
        translate_v3_to_v2({"arch_radius_mm": 30.0, "arch_shape": "ellipse",
                            "arch_width_mm": 60.0, "arch_height_mm": 40.0})


# ── validate_spec ───────────────────────────────────────────────────────────


def _v3_spec(**kw):
    base = {"schema_version": "3.0", "mode": "single", "geometry": "healthy_arch_v3"}
    base.update(kw)
    return base


def test_validator_accepts_minimal_single() -> None:
    validate_spec(_v3_spec(params={"r_inlet": 14.0}))


def test_validator_rejects_v2_geometry_string() -> None:
    with pytest.raises(ValueError, match="healthy_arch_v3"):
        validate_spec(_v3_spec(geometry="healthy_arch_v2", params={"r_inlet": 14.0}))


def test_validator_rejects_sample_mode() -> None:
    """v3 only supports single + sweep — sample/grid users go to v2."""
    with pytest.raises(ValueError, match="single.*sweep"):
        validate_spec(_v3_spec(mode="sample", params={"r_inlet": {}}))


def test_validator_typo_hint() -> None:
    with pytest.raises(ValueError, match="r_intle"):
        validate_spec(_v3_spec(params={"r_intle": 14.0}))


# ── --param parsing ─────────────────────────────────────────────────────────


def test_param_parse_float() -> None:
    k, v = _parse_param_override("arch_height_mm=42.5")
    assert k == "arch_height_mm" and v == 42.5


def test_param_parse_unknown_key_with_suggestion() -> None:
    with pytest.raises(ValueError, match="torson_deg"):
        _parse_param_override("torson_deg=10")


# ── expand_cases ────────────────────────────────────────────────────────────


def test_expand_single() -> None:
    spec = _v3_spec(case_id="x",
                    params={"r_inlet": 13.5, "r_outlet": 11.0,
                            "arch_width_mm": 85.0, "arch_height_mm": 42.5,
                            "torsion_deg": 0.0})
    cases = expand_cases(spec)
    assert len(cases) == 1
    assert cases[0]["case_id"] == "x"
    assert cases[0]["params"]["r_inlet"] == 13.5


def test_expand_sweep_torsion() -> None:
    spec = _v3_spec(mode="sweep", case_prefix="t",
                    sweep={"param": "torsion_deg", "low": -20.0, "high": 20.0, "n": 5})
    cases = expand_cases(spec)
    assert len(cases) == 5
    assert cases[0]["params"]["torsion_deg"] == -20.0
    assert cases[-1]["params"]["torsion_deg"] == 20.0
    assert cases[0]["case_id"] == "t_001"


# ── Shipped specs validate ──────────────────────────────────────────────────


def test_shipped_v3_specs_validate() -> None:
    import json
    specs_dir = HERE / "specs_v3"
    files = sorted(specs_dir.glob("*.json"))
    assert len(files) >= 2
    for f in files:
        validate_spec(json.loads(f.read_text()), source=f.name)


def test_v3_round_trip_to_v2() -> None:
    """Translating the v3 baseline produces all the v2 keys blender_aorta_v2.py needs."""
    v2 = translate_v3_to_v2({"r_inlet": 14.0, "r_outlet": 10.0,
                              "arch_width_mm": 90.0, "arch_height_mm": 45.0,
                              "torsion_deg": 0.0})
    required_v2_keys = {
        "r_ascending", "r_arch", "r_descending", "taper_mode",
        "ascending_length", "descending_length",
        "arch_R_c", "arch_angle_deg", "arch_tilt_deg", "junction_blend_mm",
        "delta_3", "delta_4", "segments_radial", "curve_samples",
    }
    # The v3 → v2 translation should cover everything blender_aorta_v2.py exposes
    # as a flag, EXCEPT v3 doesn't surface ascending_length/descending_length as
    # primary knobs in the baseline spec — they come from V2_FIXED-style defaults.
    # But translate_v3_to_v2 does inject ALL the fixed defaults, so the result
    # may not include ascending_length/descending_length if v3 spec didn't supply them.
    # That's fine — blender_aorta_v2.py has its own defaults for those.
    must_have = required_v2_keys - {"ascending_length", "descending_length"}
    assert must_have.issubset(set(v2.keys()))
