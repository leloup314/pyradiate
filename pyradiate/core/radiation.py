from dataclasses import dataclass

from pyradiate.core.nuclide import NuclideIdentifier


@dataclass(frozen=True)
class Radiation:
    energy_kev: float


@dataclass(frozen=True)
class Alpha(Radiation):
    """Alpha particle with a parsed energy (ENSDF A record, field E)."""


@dataclass(frozen=True)
class Gamma(Radiation):
    """Gamma transition with energy and intensity (ENSDF G record, E + RI/TI)."""

    intensity: float


@dataclass(frozen=True)
class Beta(Radiation):
    """Beta transition with energy and branch intensity (ENSDF B record, E + IB)."""

    intensity: float


@dataclass(frozen=True)
class Xray(Radiation):
    """Atomic X-ray from an ENSDF tG table (energy, intensity, shell line label)."""

    intensity: float
    label: str


@dataclass(frozen=True)
class Decay:
    """One evaluated decay pathway (e.g. EC decay of 89Zr)."""

    mode: str
    dataset_id: str
    daughter_nuclide: NuclideIdentifier
    half_life_s: float
    branches: tuple[tuple[str, float], ...]
    alphas: tuple[Alpha, ...]
    betas: tuple[Beta, ...]
    gammas: tuple[Gamma, ...]
    xrays: tuple[Xray, ...]
