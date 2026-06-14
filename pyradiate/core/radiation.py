from dataclasses import dataclass
from enum import StrEnum

from pyradiate.core.nuclide_id import NuclideIdentifier
from pyradiate.tools.ensdf.helpers import BranchIdentifier


@dataclass(frozen=True)
class DecayBranch:
    """One competing decay mode within a parent-state decay."""

    identifier: BranchIdentifier
    fraction: float
    daughter_nuclide: NuclideIdentifier


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


class BetaKind(StrEnum):
    """Charged beta particle emitted in the transition."""

    ELECTRON = "electron"
    POSITRON = "positron"


@dataclass(frozen=True)
class Beta(Radiation):
    """Beta transition with energy and branch intensity (ENSDF B or E record)."""

    intensity: float
    kind: BetaKind


@dataclass(frozen=True)
class Xray(Radiation):
    """Atomic X-ray from an ENSDF tG table (energy, intensity, shell line label)."""

    intensity: float
    label: str


@dataclass(frozen=True)
class Decay:
    """One parent-state decay with one or more competing branch modes."""

    dataset_id: str
    half_life_s: float
    branches: tuple[DecayBranch, ...]
    alphas: tuple[Alpha, ...]
    betas: tuple[Beta, ...]
    gammas: tuple[Gamma, ...]
    xrays: tuple[Xray, ...]

    @property
    def mode(self) -> str:
        """ENSDF-style summary label derived from branch identifiers."""
        labels = sorted(branch.identifier.value for branch in self.branches)
        return "+".join(labels) + " DECAY" if labels else ""
