"""Tests for cli_v2.py spec validator, --param override, and PARAMETERS schema."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from cli_v2 import (  # noqa: E402
    PARAMETERS,
    _format_params_table,
    _parse_param_override,
    _resolve_arch_params,
    _validate_distribution,
    apply_param_overrides,
    validate_spec,
)


# ── PARAMETERS dict invariants ─────────────────────────────────────────────


def test_parameters_dict_has_expected_count() -> None:
    # 18 parameters total: 3 radii + 1 taper_mode + 2 lengths + 6 curvature
    # (R_c, angle, tilt, twist, junction_blend, arch_shape) + 2 direct curvature
    # alternatives (arch_span_mm, arch_height_mm) + 2 non-planar Fourier
    # (δ_3, δ_4) + 2 mesh
    assert len(PARAMETERS) == 18


def test_every_parameter_has_required_keys() -> None:
    for name, info in PARAMETERS.items():
        for k in ("type", "default", "group", "description"):
            assert k in info, f"{name} missing {k}"
        assert info["type"] in {"float", "int", "str", "bool"}


def test_default_distribution_consistent_with_workshop_range() -> None:
    """Default distributions should not put substantial mass outside the workshop range."""
    for name, info in PARAMETERS.items():
        dist = info.get("default_dist")
        if dist is None:
            continue
        if dist["type"] in {"normal", "uniform", "gumbel"}:
            # Truncation bounds (when present) should be within or equal to the
            # workshop range, never wider.
            if "low" in dist and info.get("min") is not None:
                assert dist["low"] >= info["min"] - 1e-9, (
                    f"{name} default_dist.low {dist['low']} < workshop min {info['min']}"
                )
            if "high" in dist and info.get("max") is not None:
                assert dist["high"] <= info["max"] + 1e-9, (
                    f"{name} default_dist.high {dist['high']} > workshop max {info['max']}"
                )


def test_format_params_table_smoke() -> None:
    text = _format_params_table(markdown=False)
    assert "r_ascending" in text
    assert "arch_R_c" in text
    md = _format_params_table(markdown=True)
    assert "| `r_ascending`" in md
    assert "Schäfer" in md  # citations included in markdown


# ── validate_spec ──────────────────────────────────────────────────────────


def _spec(**kw):
    base = {"schema_version": "2.0", "mode": "single", "geometry": "healthy_arch_v2"}
    base.update(kw)
    return base


def test_validator_accepts_minimal_single_spec() -> None:
    validate_spec(_spec(params={"r_ascending": 13.7}))


def test_validator_rejects_unknown_geometry() -> None:
    with pytest.raises(ValueError, match="geometry"):
        validate_spec(_spec(geometry="arch_branched_coarctation",
                            params={"r_ascending": 13.7}))


def test_validator_rejects_unknown_parameter_with_suggestion() -> None:
    with pytest.raises(ValueError, match="r_ascendng"):
        validate_spec(_spec(params={"r_ascendng": 13.7}))  # typo
    try:
        validate_spec(_spec(params={"r_ascendng": 13.7}))
    except ValueError as e:
        assert "r_ascending" in str(e)


def test_validator_rejects_bad_mode() -> None:
    with pytest.raises(ValueError, match="mode"):
        validate_spec(_spec(mode="nope", params={"r_ascending": 13.7}))


def test_validator_sweep_low_lt_high() -> None:
    with pytest.raises(ValueError, match="low.*high"):
        validate_spec(_spec(mode="sweep",
                            sweep={"param": "r_ascending", "low": 18, "high": 10, "n": 5}))


def test_validator_sweep_n_min() -> None:
    with pytest.raises(ValueError, match="n must be at least 2"):
        validate_spec(_spec(mode="sweep",
                            sweep={"param": "r_ascending", "low": 10, "high": 18, "n": 1}))


def test_validator_sample_n_cases_min() -> None:
    with pytest.raises(ValueError, match="n_cases must be at least 4"):
        validate_spec(_spec(mode="sample",
                            params={"r_ascending": {}},
                            n_cases=2))


def test_validator_distribution_override_normal_requires_std() -> None:
    spec = _spec(mode="sample", params={"r_ascending": {}}, n_cases=8,
                 distribution_overrides={"r_ascending": {"type": "normal", "mean": 13.7}})
    with pytest.raises(ValueError, match="std required"):
        validate_spec(spec)


def test_validator_distribution_override_unknown_type() -> None:
    spec = _spec(mode="sample", params={"r_ascending": {}}, n_cases=8,
                 distribution_overrides={"r_ascending": {"type": "exponential", "rate": 1.0}})
    with pytest.raises(ValueError, match="must be one of"):
        validate_spec(spec)


def test_validator_distribution_override_unknown_param() -> None:
    spec = _spec(mode="sample", params={"r_ascending": {}}, n_cases=8,
                 distribution_overrides={"foo": {"type": "normal", "mean": 1, "std": 1}})
    with pytest.raises(ValueError, match="foo"):
        validate_spec(spec)


def test_validator_distribution_helper_uniform_low_lt_high() -> None:
    with pytest.raises(ValueError, match="low must be < high"):
        _validate_distribution({"type": "uniform", "low": 1, "high": 0}, source="x")


def test_validator_grid_requires_two_values() -> None:
    with pytest.raises(ValueError, match=">= 2"):
        validate_spec(_spec(mode="grid",
                            grid={"params": {"r_ascending": [13.7]}}))


# ── --param parsing ────────────────────────────────────────────────────────


def test_param_parse_float() -> None:
    k, v = _parse_param_override("r_arch=14.5")
    assert k == "r_arch" and v == 14.5


def test_param_parse_int() -> None:
    k, v = _parse_param_override("segments_radial=96")
    assert k == "segments_radial" and v == 96 and isinstance(v, int)


def test_param_parse_str() -> None:
    k, v = _parse_param_override("taper_mode=linear")
    assert k == "taper_mode" and v == "linear"


def test_param_parse_unknown_key() -> None:
    with pytest.raises(ValueError, match="r_asc"):
        _parse_param_override("r_asc=13")


def test_param_parse_missing_equals() -> None:
    with pytest.raises(ValueError, match="key=value"):
        _parse_param_override("r_arch")


# ── apply_param_overrides ──────────────────────────────────────────────────


def test_overrides_go_to_params_in_single_mode() -> None:
    spec = _spec(params={"r_ascending": 13.7})
    out = apply_param_overrides(spec, {"r_arch": 14.0})
    assert out["params"]["r_arch"] == 14.0
    assert spec["params"] == {"r_ascending": 13.7}  # input not mutated


def test_overrides_go_to_fixed_in_sweep_mode() -> None:
    spec = _spec(mode="sweep",
                 sweep={"param": "r_ascending", "low": 10, "high": 18, "n": 5},
                 fixed={"r_arch": 13.0})
    out = apply_param_overrides(spec, {"r_descending": 12.0})
    assert out["fixed"]["r_descending"] == 12.0
    assert out["fixed"]["r_arch"] == 13.0


# ── Example specs from specs_v2/ ───────────────────────────────────────────


@pytest.fixture
def specs_v2_dir() -> Path:
    return HERE / "specs_v2"


def test_all_shipped_specs_validate(specs_v2_dir: Path) -> None:
    import json

    files = sorted(specs_v2_dir.glob("*.json"))
    assert len(files) >= 2, f"Expected ≥2 example specs in specs_v2/, found {len(files)}"
    for f in files:
        payload = json.loads(f.read_text())
        validate_spec(payload, source=f.name)


def test_delta_params_default_to_zero_for_backwards_compat() -> None:
    """The scalar defaults of δ_3, δ_4 must be 0.0 so old specs that don't
    set them produce identical planar geometry to pre-2026-05-20 v2."""
    assert PARAMETERS["delta_3"]["default"] == 0.0
    assert PARAMETERS["delta_4"]["default"] == 0.0


def test_delta_default_sampling_distribution_matches_synthaorta() -> None:
    """When sampled, δ_3 / δ_4 should default to Normal(1, 0.09) (SynthAorta Table I)."""
    for name in ("delta_3", "delta_4"):
        d = PARAMETERS[name]["default_dist"]
        assert d["type"] == "normal"
        assert d["mean"] == 1.0
        assert d["std"] == 0.09


# ── _resolve_arch_params (direct span+height → R_c+angle inverse) ──────────


def test_resolve_arch_passthrough_when_no_alt_params() -> None:
    """If neither arch_span_mm nor arch_height_mm is in params, return as-is."""
    p_in = {"r_ascending": 13.7, "arch_R_c": 40.4, "arch_angle_deg": 180.0}
    p_out = _resolve_arch_params(p_in)
    assert p_out == p_in


def test_resolve_arch_u_arch_matches_canonical() -> None:
    """U-arch (S=2H): R_c=H, angle=180°."""
    p = _resolve_arch_params({"arch_span_mm": 80.8, "arch_height_mm": 40.4})
    assert p["arch_R_c"] == pytest.approx(40.4, abs=1e-9)
    assert p["arch_angle_deg"] == pytest.approx(180.0, abs=1e-6)
    # Direct keys removed
    assert "arch_span_mm" not in p
    assert "arch_height_mm" not in p


def test_resolve_arch_shallow_120deg() -> None:
    """S=1.5H → arccos(-0.5) = 120°."""
    p = _resolve_arch_params({"arch_span_mm": 60.0, "arch_height_mm": 40.0})
    assert p["arch_R_c"] == pytest.approx(40.0, abs=1e-9)
    assert p["arch_angle_deg"] == pytest.approx(120.0, abs=1e-6)


def test_resolve_arch_quarter_circle_90deg() -> None:
    """S=H → arccos(0) = 90°."""
    p = _resolve_arch_params({"arch_span_mm": 30.0, "arch_height_mm": 30.0})
    assert p["arch_R_c"] == pytest.approx(30.0, abs=1e-9)
    assert p["arch_angle_deg"] == pytest.approx(90.0, abs=1e-6)


def test_resolve_arch_rejects_unpaired() -> None:
    """Setting only span without height (or vice versa) is an error."""
    with pytest.raises(ValueError, match="must be set together"):
        _resolve_arch_params({"arch_span_mm": 80.0})
    with pytest.raises(ValueError, match="must be set together"):
        _resolve_arch_params({"arch_height_mm": 40.0})


def test_resolve_arch_rejects_span_lt_height() -> None:
    """Span must be at least equal to height (otherwise θ < 90° which we don't solve)."""
    with pytest.raises(ValueError, match="must satisfy"):
        _resolve_arch_params({"arch_span_mm": 30.0, "arch_height_mm": 40.0})


def test_resolve_arch_rejects_span_gt_2height() -> None:
    """Span > 2·height implies over-arched (θ > 180°) — unsolvable from peak height."""
    with pytest.raises(ValueError, match="must satisfy"):
        _resolve_arch_params({"arch_span_mm": 90.0, "arch_height_mm": 40.0})


def test_resolve_arch_ellipse_mode_skips_inverse() -> None:
    """When arch_shape='ellipse', span+height are passed through unchanged
    (no closed-form inverse to R_c/angle)."""
    p = _resolve_arch_params({
        "arch_shape": "ellipse",
        "arch_span_mm": 30.0,
        "arch_height_mm": 80.0,  # W < H — would be rejected in circle mode
    })
    assert p["arch_shape"] == "ellipse"
    assert p["arch_span_mm"] == 30.0
    assert p["arch_height_mm"] == 80.0
    # R_c / angle keys NOT added (Blender uses span+height in ellipse mode)
    assert "arch_R_c" not in p
    assert "arch_angle_deg" not in p


def test_resolve_arch_ellipse_requires_both_dims() -> None:
    with pytest.raises(ValueError, match="requires both"):
        _resolve_arch_params({"arch_shape": "ellipse", "arch_span_mm": 30.0})
    with pytest.raises(ValueError, match="requires both"):
        _resolve_arch_params({"arch_shape": "ellipse", "arch_height_mm": 80.0})


def test_resolve_arch_ellipse_rejects_nonpositive() -> None:
    with pytest.raises(ValueError, match="> 0"):
        _resolve_arch_params({"arch_shape": "ellipse",
                              "arch_span_mm": 30.0, "arch_height_mm": 0.0})


def test_resolve_arch_round_trip_through_expand_cases() -> None:
    """expand_cases should auto-convert span+height to R_c+angle per case."""
    from cli_v2 import expand_cases
    spec = {
        "schema_version": "2.0", "mode": "single", "geometry": "healthy_arch_v2",
        "case_id": "u", "params": {"arch_span_mm": 80.8, "arch_height_mm": 40.4},
    }
    cases = expand_cases(spec)
    assert len(cases) == 1
    p = cases[0]["params"]
    assert "arch_span_mm" not in p
    assert "arch_height_mm" not in p
    assert p["arch_R_c"] == pytest.approx(40.4)
    assert p["arch_angle_deg"] == pytest.approx(180.0, abs=1e-6)
