"""Tests for cli_v2.expand_cases — single / sweep / sample / grid modes."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from cli_v2 import expand_cases  # noqa: E402


def test_single_mode_one_case() -> None:
    spec = {
        "schema_version": "2.0", "mode": "single", "geometry": "healthy_arch_v2",
        "case_id": "base", "params": {"r_ascending": 13.7, "r_arch": 13.0},
    }
    cases = expand_cases(spec)
    assert len(cases) == 1
    assert cases[0]["case_id"] == "base"
    assert cases[0]["params"]["r_ascending"] == 13.7
    assert cases[0]["mode"] == "single"


def test_single_mode_merges_fixed_and_params() -> None:
    spec = {
        "schema_version": "2.0", "mode": "single", "geometry": "healthy_arch_v2",
        "case_id": "x",
        "fixed": {"r_arch": 13.0, "r_descending": 12.0},
        "params": {"r_ascending": 14.0},
    }
    cases = expand_cases(spec)
    p = cases[0]["params"]
    assert p["r_ascending"] == 14.0
    assert p["r_arch"] == 13.0
    assert p["r_descending"] == 12.0


def test_sweep_mode_n_cases() -> None:
    spec = {
        "schema_version": "2.0", "mode": "sweep", "geometry": "healthy_arch_v2",
        "case_prefix": "r",
        "sweep": {"param": "r_ascending", "low": 10.0, "high": 18.0, "n": 9},
        "fixed": {"r_arch": 13.0},
    }
    cases = expand_cases(spec)
    assert len(cases) == 9
    assert cases[0]["params"]["r_ascending"] == 10.0
    assert cases[-1]["params"]["r_ascending"] == 18.0
    # IDs
    assert cases[0]["case_id"] == "r_001"
    assert cases[-1]["case_id"] == "r_009"
    # Fixed parameter present in every case
    for c in cases:
        assert c["params"]["r_arch"] == 13.0


def test_sweep_mode_linear_spacing() -> None:
    spec = {
        "schema_version": "2.0", "mode": "sweep", "geometry": "healthy_arch_v2",
        "sweep": {"param": "arch_R_c", "low": 25.0, "high": 60.0, "n": 8},
    }
    cases = expand_cases(spec)
    vals = [c["params"]["arch_R_c"] for c in cases]
    # Constant step
    deltas = [b - a for a, b in zip(vals, vals[1:])]
    expected = (60.0 - 25.0) / 7
    for d in deltas:
        assert abs(d - expected) < 1e-9


def test_sample_mode_returns_n_cases() -> None:
    spec = {
        "schema_version": "2.0", "mode": "sample", "geometry": "healthy_arch_v2",
        "sampler": "sobol", "n_cases": 16, "seed": 0,
        "params": {"r_ascending": {}, "r_arch": {}, "r_descending": {}},
        "fixed": {"taper_mode": "smoothstep"},
    }
    cases = expand_cases(spec)
    assert len(cases) == 16
    for c in cases:
        assert "r_ascending" in c["params"]
        assert "r_arch" in c["params"]
        assert "r_descending" in c["params"]
        # Fixed param merged
        assert c["params"]["taper_mode"] == "smoothstep"
    assert cases[0]["case_id"] == "case_001"
    assert cases[-1]["case_id"] == "case_016"


def test_grid_mode_cartesian() -> None:
    spec = {
        "schema_version": "2.0", "mode": "grid", "geometry": "healthy_arch_v2",
        "grid": {"params": {"r_ascending": [12.0, 14.0],
                            "arch_R_c": [30.0, 40.0, 50.0]}},
        "fixed": {"r_arch": 13.0},
    }
    cases = expand_cases(spec)
    assert len(cases) == 2 * 3
    pairs = {(c["params"]["r_ascending"], c["params"]["arch_R_c"]) for c in cases}
    assert pairs == {(12.0, 30.0), (12.0, 40.0), (12.0, 50.0),
                     (14.0, 30.0), (14.0, 40.0), (14.0, 50.0)}
    for c in cases:
        assert c["params"]["r_arch"] == 13.0


def test_unknown_mode_raises() -> None:
    with pytest.raises(ValueError, match="mode"):
        expand_cases({"schema_version": "2.0", "mode": "bogus",
                      "geometry": "healthy_arch_v2"})


# ── End-to-end with shipped specs ──────────────────────────────────────────


@pytest.mark.parametrize(
    ("spec_name", "expected_n"),
    [
        ("single_baseline_v2.json", 1),
        ("single_span_height_v2.json", 1),
        ("sweep_r_ascending_v2.json", 10),
        ("sweep_R_c_v2.json", 10),
        ("sweep_arch_angle_v2.json", 10),
        ("sweep_arch_tilt_v2.json", 10),
        ("sample_sobol_synthaorta_v2.json", 256),
        ("sample_sobol_synthaorta_nonplanar_v2.json", 256),
    ],
)
def test_shipped_spec_produces_expected_case_count(spec_name: str, expected_n: int) -> None:
    import json

    spec_path = HERE / "specs_v2" / spec_name
    spec = json.loads(spec_path.read_text())
    cases = expand_cases(spec)
    assert len(cases) == expected_n, (
        f"{spec_name} produced {len(cases)} cases, expected {expected_n}"
    )
