"""Tests for the spec validator + --param override + grid mode + cost warning.

Covers the workshop-UX additions in cli.py: PARAMETERS dict, validate_spec,
typo detection, --param key=value parsing, grid mode expansion, cost warning.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from cli import (  # noqa: E402
    PARAMETERS,
    _format_params_table,
    _parse_param_override,
    apply_param_overrides,
    estimate_case_count,
    expand_cases,
    validate_spec,
    warn_if_expensive,
)


# ── PARAMETERS dict invariants ─────────────────────────────────────────────


def test_parameters_dict_has_all_required_keys() -> None:
    for name, info in PARAMETERS.items():
        assert "type" in info, f"{name} missing 'type'"
        assert "default" in info, f"{name} missing 'default'"
        assert "min" in info, f"{name} missing 'min'"
        assert "max" in info, f"{name} missing 'max'"
        assert "group" in info, f"{name} missing 'group'"
        assert "description" in info, f"{name} missing 'description'"
        assert info["type"] in {"float", "int", "bool"}, f"{name} has bad type {info['type']!r}"


def test_parameters_dict_covers_blender_flags() -> None:
    """Every key in the DIRECT_FLAGS / BOOL_FLAGS used by build_blender_cmd
    should appear in PARAMETERS (so --list-params is the canonical reference)."""
    from cli import BOOL_FLAGS, DIRECT_FLAGS
    for k in DIRECT_FLAGS:
        assert k in PARAMETERS, f"DIRECT_FLAGS entry {k!r} missing from PARAMETERS"
    # coarctation is auto-derived from coarctation_area_reduction; roughness is explicit
    assert "roughness" in PARAMETERS


# ── --list-params formatting ───────────────────────────────────────────────


def test_list_params_table_terminal() -> None:
    s = _format_params_table(markdown=False)
    assert "Available parameters" in s
    assert "diameter" in s
    assert "coarctation_area_reduction" in s
    assert "Main tube" in s


def test_list_params_table_markdown() -> None:
    s = _format_params_table(markdown=True)
    assert s.startswith("# Aorta geometry generator")
    assert "| `diameter` |" in s
    assert "## Main tube" in s
    assert "## Coarctation" in s


# ── Validator: typo detection ──────────────────────────────────────────────


def test_validate_typo_in_params_suggests_correction() -> None:
    with pytest.raises(ValueError, match="Did you mean 'diameter'"):
        validate_spec({"mode": "single", "params": {"diametr": 24.0}})


def test_validate_typo_in_fixed_suggests_correction() -> None:
    with pytest.raises(ValueError, match="Did you mean 'arch_height'"):
        validate_spec({"mode": "sweep", "fixed": {"arch_hieght": 30.0},
                       "sweep": {"param": "diameter", "low": 22, "high": 36, "n": 5}})


def test_validate_typo_in_sweep_param_suggests_correction() -> None:
    with pytest.raises(ValueError, match="Did you mean 'coarctation_area_reduction'"):
        validate_spec({"mode": "sweep",
                       "sweep": {"param": "coarc_area_reduction", "low": 0, "high": 0.9, "n": 5}})


def test_validate_unknown_param_no_close_match() -> None:
    with pytest.raises(ValueError, match="Unknown parameter 'xyzzy42'"):
        validate_spec({"mode": "single", "params": {"xyzzy42": 0.0}})


# ── Validator: mode + range checks ─────────────────────────────────────────


def test_validate_invalid_mode() -> None:
    with pytest.raises(ValueError, match="mode must be one of single/sweep/sample/grid"):
        validate_spec({"mode": "bogus"})


def test_validate_sweep_low_must_be_less_than_high() -> None:
    with pytest.raises(ValueError, match="low.*must be < high"):
        validate_spec({"mode": "sweep",
                       "sweep": {"param": "diameter", "low": 30, "high": 20, "n": 5}})


def test_validate_sweep_n_too_small() -> None:
    with pytest.raises(ValueError, match="sweep.n must be at least 2"):
        validate_spec({"mode": "sweep",
                       "sweep": {"param": "diameter", "low": 22, "high": 36, "n": 1}})


def test_validate_sample_requires_low_high() -> None:
    with pytest.raises(ValueError, match=r"params\['diameter'\] must be"):
        validate_spec({"mode": "sample", "n_cases": 10,
                       "params": {"diameter": [22, 36]}})


def test_validate_sample_n_cases_too_small() -> None:
    with pytest.raises(ValueError, match="n_cases must be at least 4"):
        validate_spec({"mode": "sample", "n_cases": 2,
                       "params": {"diameter": {"low": 22, "high": 36}}})


# ── Validator: grid mode ───────────────────────────────────────────────────


def test_validate_grid_requires_at_least_two_values() -> None:
    with pytest.raises(ValueError, match=r"params\['diameter'\] must be a list of >= 2"):
        validate_spec({"mode": "grid", "grid": {"params": {"diameter": [24]}}})


def test_validate_grid_typo_in_param_name() -> None:
    with pytest.raises(ValueError, match="Did you mean"):
        validate_spec({"mode": "grid",
                       "grid": {"params": {"diametr": [22, 36], "arch_height": [30, 40]}}})


# ── --param key=value parsing ──────────────────────────────────────────────


def test_param_override_float() -> None:
    k, v = _parse_param_override("diameter=28")
    assert k == "diameter"
    assert v == 28.0
    assert isinstance(v, float)


def test_param_override_int() -> None:
    k, v = _parse_param_override("branch_count=2")
    assert k == "branch_count"
    assert v == 2
    assert isinstance(v, int)


def test_param_override_bool_true() -> None:
    k, v = _parse_param_override("roughness=true")
    assert v is True


def test_param_override_bool_false() -> None:
    k, v = _parse_param_override("roughness=no")
    assert v is False


def test_param_override_unknown_key_typo_hint() -> None:
    with pytest.raises(ValueError, match="Did you mean 'diameter'"):
        _parse_param_override("diametr=28")


def test_param_override_no_equals() -> None:
    with pytest.raises(ValueError, match="expects key=value"):
        _parse_param_override("diameter28")


# ── apply_param_overrides ──────────────────────────────────────────────────


def test_apply_overrides_single_mode_goes_to_params() -> None:
    spec = {"mode": "single", "params": {"diameter": 24}}
    out = apply_param_overrides(spec, {"arch_height": 40})
    assert out["params"]["diameter"] == 24
    assert out["params"]["arch_height"] == 40


def test_apply_overrides_sweep_mode_goes_to_fixed() -> None:
    spec = {"mode": "sweep", "fixed": {"diameter": 24},
            "sweep": {"param": "coarctation_area_reduction", "low": 0, "high": 0.9, "n": 5}}
    out = apply_param_overrides(spec, {"arch_height": 40})
    assert out["fixed"]["diameter"] == 24
    assert out["fixed"]["arch_height"] == 40
    assert "arch_height" not in out["sweep"]


def test_apply_overrides_does_not_mutate_input() -> None:
    spec = {"mode": "single", "params": {"diameter": 24}}
    apply_param_overrides(spec, {"arch_height": 40})
    assert "arch_height" not in spec["params"]


# ── Cost estimation & warning ──────────────────────────────────────────────


def test_estimate_case_count_for_each_mode() -> None:
    assert estimate_case_count({"mode": "single"}) == 1
    assert estimate_case_count({"mode": "sweep", "sweep": {"n": 10}}) == 10
    assert estimate_case_count({"mode": "sample", "n_cases": 50}) == 50
    assert estimate_case_count({"mode": "grid",
                                "grid": {"params": {"a": [1, 2, 3], "b": [10, 20]}}}) == 6


def test_warn_if_expensive_silent_below_threshold() -> None:
    buf = io.StringIO()
    warn_if_expensive(10, threshold=30, stream=buf)
    assert buf.getvalue() == ""


def test_warn_if_expensive_emits_above_threshold() -> None:
    buf = io.StringIO()
    warn_if_expensive(60, threshold=30, stream=buf)
    out = buf.getvalue()
    assert "60 cases" in out
    assert "minutes" in out


# ── Grid mode end-to-end via expand_cases ──────────────────────────────────


def test_grid_mode_expands_to_cartesian_product() -> None:
    spec = {
        "mode": "grid",
        "case_prefix": "g",
        "fixed": {"diameter": 24},
        "grid": {"params": {"coarctation_area_reduction": [0.0, 0.5],
                            "arch_height": [30, 40]}},
    }
    cases = expand_cases(spec)
    assert len(cases) == 4
    params_set = {(c["params"]["coarctation_area_reduction"],
                   c["params"]["arch_height"]) for c in cases}
    assert params_set == {(0.0, 30), (0.0, 40), (0.5, 30), (0.5, 40)}
    # fixed propagates
    for c in cases:
        assert c["params"]["diameter"] == 24
        assert c["mode"] == "grid"
