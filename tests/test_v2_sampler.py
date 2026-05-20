"""Tests for cli_v2.py distribution sampling (_map_unit_to_distribution + sample_cases)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from cli_v2 import (  # noqa: E402
    PARAMETERS,
    _map_unit_to_distribution,
    sample_cases,
)


# ── _map_unit_to_distribution ──────────────────────────────────────────────


def test_uniform_endpoints_correct() -> None:
    dist = {"type": "uniform", "low": 10.0, "high": 20.0}
    u = np.array([0.0, 0.5, 1.0])
    out = _map_unit_to_distribution(u, dist)
    np.testing.assert_allclose(out, [10.0, 15.0, 20.0])


def test_normal_truncation_respects_bounds() -> None:
    dist = {"type": "normal", "mean": 13.7, "std": 2.3, "low": 8.0, "high": 22.0}
    rng = np.random.default_rng(0)
    u = rng.random(10000)
    samples = _map_unit_to_distribution(u, dist)
    assert (samples >= 8.0 - 1e-9).all()
    assert (samples <= 22.0 + 1e-9).all()
    # Mean should be near 13.7 (truncation slightly shifts it, but tolerance is loose)
    assert abs(samples.mean() - 13.7) < 0.3


def test_normal_mean_std_match_target_without_truncation() -> None:
    dist = {"type": "normal", "mean": 100.0, "std": 5.0,
            "low": float("-inf"), "high": float("inf")}
    rng = np.random.default_rng(0)
    u = rng.random(20000)
    samples = _map_unit_to_distribution(u, dist)
    assert abs(samples.mean() - 100.0) < 0.2
    assert abs(samples.std() - 5.0) < 0.2


def test_gumbel_truncation_respects_bounds() -> None:
    dist = {"type": "gumbel", "loc": 40.4, "scale": 2.4, "low": 25.0, "high": 60.0}
    rng = np.random.default_rng(0)
    u = rng.random(5000)
    samples = _map_unit_to_distribution(u, dist)
    assert (samples >= 25.0).all()
    assert (samples <= 60.0).all()
    # Gumbel mode is at `loc`; mean ≈ loc + scale*euler ≈ 40.4 + 2.4*0.5772 ≈ 41.8
    assert 40.0 < samples.mean() < 43.0


def test_unknown_distribution_type_raises() -> None:
    with pytest.raises(ValueError, match="Unknown distribution"):
        _map_unit_to_distribution(np.array([0.5]),
                                  {"type": "exponential", "rate": 1.0})


# ── sample_cases (end-to-end) ──────────────────────────────────────────────


def test_sample_cases_uses_default_distributions_when_params_empty() -> None:
    spec = {
        "schema_version": "2.0",
        "mode": "sample",
        "geometry": "healthy_arch_v2",
        "sampler": "sobol",
        "n_cases": 64,
        "seed": 42,
        "params": {
            "r_ascending": {},
            "r_arch": {},
            "r_descending": {},
            "arch_R_c": {},
        },
    }
    cases = sample_cases(spec)
    assert len(cases) == 64
    # All values within workshop bounds
    for c in cases:
        assert 8.0 <= c["r_ascending"] <= 22.0
        assert 8.0 <= c["r_arch"] <= 20.0
        assert 8.0 <= c["r_descending"] <= 20.0
        assert 25.0 <= c["arch_R_c"] <= 60.0

    # Sample means should be near the default distribution means (loose tol)
    ras = np.array([c["r_ascending"] for c in cases])
    assert abs(ras.mean() - 13.7) < 1.0  # 64 cases, sobol — small but biased


def test_sample_cases_respects_distribution_overrides() -> None:
    spec = {
        "schema_version": "2.0",
        "mode": "sample",
        "geometry": "healthy_arch_v2",
        "sampler": "sobol",
        "n_cases": 32,
        "seed": 0,
        "params": {"r_ascending": {}},
        "distribution_overrides": {
            "r_ascending": {"type": "uniform", "low": 15.0, "high": 17.0},
        },
    }
    cases = sample_cases(spec)
    for c in cases:
        assert 15.0 <= c["r_ascending"] <= 17.0


def test_sample_cases_respects_explicit_low_high_in_params() -> None:
    spec = {
        "schema_version": "2.0",
        "mode": "sample",
        "geometry": "healthy_arch_v2",
        "sampler": "sobol",
        "n_cases": 16,
        "seed": 0,
        # Explicit low/high in params overrides default distribution
        "params": {"arch_R_c": {"low": 30.0, "high": 35.0}},
    }
    cases = sample_cases(spec)
    for c in cases:
        assert 30.0 <= c["arch_R_c"] <= 35.0


def test_sample_cases_deterministic_for_same_seed() -> None:
    spec = {
        "schema_version": "2.0",
        "mode": "sample",
        "geometry": "healthy_arch_v2",
        "sampler": "sobol",
        "n_cases": 8,
        "seed": 123,
        "params": {"r_ascending": {}, "r_arch": {}, "r_descending": {}},
    }
    a = sample_cases(spec)
    b = sample_cases(spec)
    for ca, cb in zip(a, b):
        for k in ca:
            assert ca[k] == cb[k]
