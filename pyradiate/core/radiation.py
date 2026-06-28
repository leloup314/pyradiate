from dataclasses import dataclass

from pyradiate.core.nuclide_id import NuclideIdentifier
from pyradiate.tools.ensdf.helpers import BranchIdentifier


@dataclass(frozen=True)
class Radiation:
    energy_kev: float


@dataclass(frozen=True)
class Alpha(Radiation):
    """Alpha particle with a parsed energy (ENSDF A record, field E)."""

    def __repr__(self):
        return super().__repr__().replace("Alpha", "\u03b1")


@dataclass(frozen=True)
class BetaPlus(Radiation):
    """Beta (Positron) transition with energy and branch intensity (ENSDF B or E record)."""

    intensity: float

    def __repr__(self):
        return super().__repr__().replace("BetaPlus", "\u03b2+")


@dataclass(frozen=True)
class BetaMinus(Radiation):
    """Beta (Electron) transition with energy and branch intensity (ENSDF B or E record)."""

    intensity: float

    def __repr__(self):
        return super().__repr__().replace("BetaMinus", "\u03b2-")


@dataclass(frozen=True)
class Gamma(Radiation):
    """Gamma transition with energy and intensity (ENSDF G record, E + RI/TI)."""

    intensity: float

    def __repr__(self):
        return super().__repr__().replace("Gamma", "\u03b3")


@dataclass(frozen=True)
class Xray(Radiation):
    """Atomic X-ray from an ENSDF tG table (energy, intensity, shell line label)."""

    intensity: float
    label: str

    def __repr__(self):
        return super().__repr__().replace("Xray", "x")


@dataclass(frozen=True)
class DecayBranch:
    """One competing decay mode within a parent-state decay."""

    identifier: BranchIdentifier
    fraction: float
    daughter_nuclide: NuclideIdentifier


@dataclass(frozen=True)
class Decay:
    """One parent-state decay with one or more competing branch modes."""

    dataset_id: str
    half_life_s: float
    branches: tuple[DecayBranch, ...]
    alphas: tuple[Alpha, ...]
    betas: tuple[BetaPlus | BetaMinus, ...]
    gammas: tuple[Gamma, ...]
    xrays: tuple[Xray, ...]

    @property
    def mode(self) -> str:
        """ENSDF-style summary label derived from branch identifiers."""
        labels = sorted(branch.identifier.value for branch in self.branches)
        return "+".join(labels) + " DECAY" if labels else ""
