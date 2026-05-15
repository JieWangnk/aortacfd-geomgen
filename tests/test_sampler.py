"""Tests for sampler.py and cli.expand_cases."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from sampler import SobolSampler, LHSSampler, RandomSampler, get_sampler, map_unit_to_physical, linear_sweep  # noqa: E402
from cli import expand_cases, load_spec  # noqa: E402


def test_sobol_returns_unit_cube_samples() -> None:
    s = SobolSampler()
    out = s.sample(8, 4, seed=0)
    assert out.shape == (8, 4)
    assert np.all(out >= 0.0)
    assert np.all(out <= 1.0)


def test_sobol_handles_non_power_of_two() -> None:
    s = SobolSampler()
    out = s.sample(7, 3, seed=0)
    assert out.shape == (7, 3)


def test_lhs_returns_unit_cube_samples() -> None:
    s = LHSSampler()
    out = s.sample(20, 3, seed=0)
    assert out.shape == (20, 3)
    assert np.all((out >= 0) & (out <= 1))


def test_random_is_deterministic_with_seed() -> None:
    s = RandomSampler()
    a = s.sample(10, 5, seed=123)
    b = s.sample(10, 5, seed=123)
    c = s.sample(10, 5, seed=124)
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)


def test_get_sampler_unknown_raises() -> None:
    with pytest.raises(ValueError, match="Unknown sampler"):
        get_sampler("definitely_not_a_sampler")


def test_map_unit_to_physical_linear() -> None:
    unit = np.array([[0.0, 0.5, 1.0], [1.0, 0.0, 0.5]])
    params = {
        "diameter": {"low": 20.0, "high": 40.0},
        "arch_height": {"low": 25.0, "high": 45.0},
        "severity": {"low": 0.0, "high": 0.9},
    }
    physical = map_unit_to_physical(unit, params)
    assert len(physical) == 2
    assert physical[0]["diameter"] == 20.0
    assert physical[0]["arch_height"] == 35.0
    assert physical[0]["severity"] == pytest.approx(0.9)
    assert physical[1]["diameter"] == 40.0
    assert physical[1]["arch_height"] == 25.0
    assert physical[1]["severity"] == pytest.approx(0.45)


def test_map_unit_to_physical_log_scale() -> None:
    unit = np.array([[0.0], [0.5], [1.0]])
    params = {"x": {"low": 1.0, "high": 100.0, "scale": "log"}}
    physical = map_unit_to_physical(unit, params)
    assert physical[0]["x"] == pytest.approx(1.0)
    assert physical[1]["x"] == pytest.approx(10.0)
    assert physical[2]["x"] == pytest.approx(100.0)


def test_linear_sweep_endpoints() -> None:
    values = linear_sweep(0.0, 0.9, 10)
    assert len(values) == 10
    assert values[0] == 0.0
    assert values[-1] == 0.9


def test_expand_cases_single_mode(tmp_path: Path) -> None:
    spec = {
        "schema_version": "1.0",
        "name": "test",
        "mode": "single",
        "geometry": "arch_branched",
        "params": {"diameter": 24.0, "arch_height": 35.0},
        "fixed": {"branch_count": 3},
    }
    cases = expand_cases(spec)
    assert len(cases) == 1
    assert cases[0]["params"]["diameter"] == 24.0
    assert cases[0]["params"]["branch_count"] == 3
    assert cases[0]["mode"] == "single"


def test_expand_cases_sweep_mode() -> None:
    spec = {
        "schema_version": "1.0",
        "name": "test",
        "mode": "sweep",
        "geometry": "arch_branched_coarctation",
        "sweep": {"param": "coarctation_area_reduction", "low": 0.0, "high": 0.9, "n": 10},
        "fixed": {"diameter": 24.0, "branch_count": 3},
    }
    cases = expand_cases(spec)
    assert len(cases) == 10
    assert cases[0]["params"]["coarctation_area_reduction"] == 0.0
    assert cases[-1]["params"]["coarctation_area_reduction"] == pytest.approx(0.9)
    # fixed propagates to every case
    for c in cases:
        assert c["params"]["diameter"] == 24.0
        assert c["params"]["branch_count"] == 3


def test_expand_cases_sample_mode() -> None:
    spec = {
        "schema_version": "1.0",
        "name": "test",
        "mode": "sample",
        "geometry": "arch_branched_coarctation",
        "sampler": "sobol",
        "n_cases": 8,
        "seed": 42,
        "params": {
            "diameter": {"low": 22.0, "high": 36.0},
            "area_reduction": {"low": 0.0, "high": 0.9},
        },
        "fixed": {"branch_count": 3},
    }
    cases = expand_cases(spec)
    assert len(cases) == 8
    for c in cases:
        assert 22.0 <= c["params"]["diameter"] <= 36.0
        assert 0.0 <= c["params"]["area_reduction"] <= 0.9
        assert c["params"]["branch_count"] == 3
        assert c["sampler"] == "sobol"
        assert c["seed"] == 42


def test_expand_cases_sample_mode_reproducible() -> None:
    spec = {
        "mode": "sample",
        "sampler": "sobol",
        "n_cases": 4,
        "seed": 7,
        "params": {"a": {"low": 0.0, "high": 1.0}, "b": {"low": 10.0, "high": 20.0}},
    }
    cases_1 = expand_cases(spec)
    cases_2 = expand_cases(spec)
    for c1, c2 in zip(cases_1, cases_2):
        assert c1["params"] == c2["params"]


def test_invalid_mode_raises() -> None:
    with pytest.raises(ValueError, match="spec.mode must be"):
        expand_cases({"mode": "bogus"})


def test_load_spec_round_trip(tmp_path: Path) -> None:
    spec_path = tmp_path / "spec.json"
    spec_path.write_text(json.dumps({
        "mode": "single",
        "geometry": "arch_branched",
        "params": {"diameter": 24.0},
    }))
    loaded = load_spec(spec_path)
    assert loaded["mode"] == "single"
    assert loaded["params"]["diameter"] == 24.0
