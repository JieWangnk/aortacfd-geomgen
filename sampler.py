"""Parameter samplers for the aorta geometry generator (Block A).

Pure Python — scipy + numpy only, no Blender, no AortaCFD-app dependency.

Three samplers:

  - ``SobolSampler``  : low-discrepancy quasi-random (recommended for ML
                        training data and global sensitivity studies)
  - ``LHSSampler``    : Latin Hypercube — guaranteed marginal coverage,
                        less uniform joint coverage
  - ``RandomSampler`` : plain uniform random (fastest, baseline)

All three return a (n_cases, n_params) array of points in [0,1]^M. The
caller maps to physical units.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Sampler(Protocol):
    def sample(self, n_cases: int, n_params: int, seed: int = 0) -> np.ndarray:
        """Return an (n_cases, n_params) array of samples in [0,1]."""
        ...


@dataclass
class SobolSampler:
    """Sobol quasi-random sequence (scrambled).

    Works best when ``n_cases`` is a power of two; otherwise we round
    up to the next power, generate, and slice. The low-discrepancy
    property is preserved for the prefix.
    """

    def sample(self, n_cases: int, n_params: int, seed: int = 0) -> np.ndarray:
        from scipy.stats.qmc import Sobol

        sampler = Sobol(d=n_params, scramble=True, seed=seed)
        m_pow = int(np.ceil(np.log2(max(n_cases, 2))))
        n_sobol = 2 ** m_pow
        raw = sampler.random(n_sobol)
        return np.asarray(raw[:n_cases], dtype=float)


@dataclass
class LHSSampler:
    """Latin Hypercube sampling — every marginal interval is hit once."""

    def sample(self, n_cases: int, n_params: int, seed: int = 0) -> np.ndarray:
        from scipy.stats.qmc import LatinHypercube

        sampler = LatinHypercube(d=n_params, seed=seed)
        return np.asarray(sampler.random(n=n_cases), dtype=float)


@dataclass
class RandomSampler:
    """Plain uniform random — baseline, not recommended for production sweeps."""

    def sample(self, n_cases: int, n_params: int, seed: int = 0) -> np.ndarray:
        rng = np.random.default_rng(seed)
        return rng.random((n_cases, n_params))


SAMPLERS: dict[str, type[Sampler]] = {
    "sobol": SobolSampler,
    "lhs": LHSSampler,
    "random": RandomSampler,
}


def get_sampler(name: str) -> Sampler:
    name = (name or "sobol").lower()
    if name not in SAMPLERS:
        raise ValueError(f"Unknown sampler {name!r}. Valid: {sorted(SAMPLERS)}.")
    return SAMPLERS[name]()


def map_unit_to_physical(unit_samples: np.ndarray, params: dict[str, dict]) -> list[dict]:
    """Map [0,1]^M unit samples to physical units using per-parameter ranges.

    Parameters
    ----------
    unit_samples
        Shape (n_cases, n_params). Order of columns matches insertion
        order of ``params``.
    params
        Dict mapping param name -> {"low": ..., "high": ..., "scale": "linear"|"log"}.
        ``scale`` defaults to ``"linear"``.

    Returns
    -------
    List of dicts, one per case, mapping name -> physical value.
    """
    names = list(params.keys())
    if unit_samples.shape[1] != len(names):
        raise ValueError(
            f"unit_samples has {unit_samples.shape[1]} cols but {len(names)} params given"
        )
    out: list[dict] = []
    for row in unit_samples:
        case: dict = {}
        for j, name in enumerate(names):
            p = params[name]
            lo, hi = float(p["low"]), float(p["high"])
            scale = (p.get("scale") or "linear").lower()
            u = float(row[j])
            if scale == "log":
                if lo <= 0 or hi <= 0:
                    raise ValueError(f"log-scale param {name!r} needs strictly positive bounds")
                v = float(np.exp(np.log(lo) + (np.log(hi) - np.log(lo)) * u))
            else:
                v = lo + (hi - lo) * u
            case[name] = v
        out.append(case)
    return out


def linear_sweep(low: float, high: float, n: int) -> np.ndarray:
    """Return n equally-spaced values from low to high (inclusive)."""
    return np.linspace(low, high, n)


__all__ = [
    "Sampler",
    "SobolSampler",
    "LHSSampler",
    "RandomSampler",
    "SAMPLERS",
    "get_sampler",
    "map_unit_to_physical",
    "linear_sweep",
]
