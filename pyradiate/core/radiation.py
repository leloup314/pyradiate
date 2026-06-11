import re
from dataclasses import dataclass
from enum import StrEnum

from pyradiate.core.nuclide import NuclideIdentifier


class DecayBranch(StrEnum):
    """ENSDF decay branch identifier (e.g. %B-=, %IT=, %EC+%B+=)."""

    C14 = "14C"
    B_MINUS_2 = "2B-"
    EC_2 = "2EC"
    N_2 = "2N"
    P_2 = "2P"
    ALPHA = "A"
    B_PLUS = "B+"
    B_PLUS_2 = "2B+"
    B_PLUS_2P = "B+2P"
    B_PLUS_3P = "B+3P"
    B_PLUS_ALPHA = "B+A"
    B_PLUS_PROTON = "B+P"
    B_MINUS = "B-"
    B_MINUS_2N = "B-2N"
    B_MINUS_ALPHA = "B-A"
    B_MINUS_N = "B-N"
    B_MINUS_PROTON = "B-P"
    EC = "EC"
    EC_B_PLUS = "EC+B+"
    EC_2P = "EC2P"
    EC_3P = "EC3P"
    EC_ALPHA = "ECA"
    EC_PROTON = "ECP"
    IT = "IT"
    NEUTRON = "N"
    PROTON = "P"
    SF = "SF"

    @classmethod
    def from_string(cls, ensdf_decay_string: str):
        normalized = re.sub(r"\s+DECAY.*", "", ensdf_decay_string.strip(), flags=re.IGNORECASE).strip().upper()
        for branch in cls:
            if branch.value == normalized:
                return branch
        raise ValueError(f"Unknown decay branch identifier: {ensdf_decay_string!r}")


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
    branch: DecayBranch
    branch_fraction: float
    alphas: tuple[Alpha, ...]
    betas: tuple[Beta, ...]
    gammas: tuple[Gamma, ...]
    xrays: tuple[Xray, ...]
