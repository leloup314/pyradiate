"""Minimal parser for ENSDF decay datasets."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

# Half-life unit multipliers to seconds (ENSDF field T)
_HALFLIFE_TO_SECONDS = {
    "YS": 1e-24,
    "ZS": 1e-21,
    "AS": 1e-18,
    "FS": 1e-15,
    "PS": 1e-12,
    "NS": 1e-9,
    "US": 1e-6,
    "MS": 1e-3,
    "S": 1.0,
    "M": 60.0,
    "MIN": 60.0,
    "H": 3600.0,
    "D": 86400.0,
    "Y": 31557600.0,  # Julian year approximation
    "KY": 31557600e3,
    "MY": 31557600e6,
}

_DECAY_DATASET_RE = re.compile(
    r"^\s*(\S{1,5})\s{4,}(.+?\b([A-Z][A-Z0-9+\-]*)\s+DECAY\b.*?)\s{2,}\S",
    re.IGNORECASE,
)

_DECAY_MODE_RE = re.compile(
    r"^(\d+\w*)\s+((?:B[+-]|EC|A|IT|SF|EP|B-N|B\+N|B-2N|B-3N|B\+EC|B\+A|"
    r"B-EC|B\+B-|B-DECAY|B\+B\+|B\+A|B\+N|B\+P|B\+EC\+A|"
    r"EC\+B\+|EC\+A|A\+EC|A\+B\+|A\+B-|A\+EC\+B\+)\s+DECAY)",
    re.IGNORECASE,
)

_BRANCH_RE = re.compile(r"%([A-Z][A-Z0-9+\-]*)=([^\s$%]+)")

_RADIATION_RECORDS = frozenset({"G", "B", "A"})


@dataclass(frozen=True)
class AlphaRay:
    """Alpha particle with a parsed energy (ENSDF A record, field E)."""

    energy_kev: float


@dataclass(frozen=True)
class GammaRay:
    """Gamma transition with energy and intensity (ENSDF G record, E + RI/TI)."""

    energy_kev: float
    intensity: float


@dataclass(frozen=True)
class BetaRay:
    """Beta transition with energy and branch intensity (ENSDF B record, E + IB)."""

    energy_kev: float
    intensity: float


@dataclass(frozen=True)
class DecayMode:
    """One evaluated decay pathway (e.g. EC decay of 89Zr)."""

    mode: str
    dataset_id: str
    daughter_nuclide: str
    half_life_s: float
    branches: tuple[tuple[str, float], ...]
    alphas: tuple[AlphaRay, ...]
    betas: tuple[BetaRay, ...]
    gammas: tuple[GammaRay, ...]

    def __iter__(self) -> Iterator[AlphaRay | BetaRay | GammaRay]:
        yield from self.alphas
        yield from self.betas
        yield from self.gammas


@dataclass(frozen=True)
class Radionuclide:
    """Radionuclide with one or more evaluated decay modes."""

    nuclide_id: str
    mass_number: int
    symbol: str
    half_life_s: float
    decay_modes: tuple[DecayMode, ...]

    def __iter__(self) -> Iterator[DecayMode]:
        yield from self.decay_modes


def _pad80(line: str) -> str:
    return line if len(line) >= 80 else line.ljust(80)


def _slice(line: str, start: int, end: int) -> str:
    """1-based inclusive ENSDF column slice."""
    return line[start - 1 : end].strip()


def _parse_value(text: str) -> float | None:
    """Parse a single ENSDF numeric value; return None if not available."""
    if not text:
        return None
    raw = text.strip()
    upper = raw.upper()
    if upper in ("", "?", "UNKNOWN", "AP", "SY") or "STABLE" in upper:
        return None

    token = raw.split()[0].lstrip("<>")
    if token.upper() in ("?", "AP", "SY"):
        return None

    m = re.match(r"^(\d+\.?\d*(?:[Ee][+-]?\d+)?)\+(\d+)-(\d+)$", token.replace(" ", ""))
    if m:
        return float(m.group(1))

    parts = raw.split()
    value_s = parts[0].lstrip("<>")
    try:
        value = float(value_s)
    except ValueError:
        return None

    if len(parts) >= 2 and parts[1].upper() in ("E2", "E3", "E4", "E5", "E6"):
        value *= 10 ** int(parts[1][1])

    return value


def _parse_energy(line: str) -> float | None:
    """Cols 10-19 (E), including E2/E3 multipliers."""
    e_text = _slice(line, 10, 19)
    value = _parse_value(e_text)
    if value is not None and e_text:
        m = re.search(r"E(\d)\s*$", e_text, re.IGNORECASE)
        if m:
            value *= 10 ** int(m.group(1))
    return value


def _parse_intensity(line: str, record: str) -> float | None:
    """IB (beta/alpha) cols 22-29; RI or TI for gamma."""
    if record == "G":
        text = _slice(line, 22, 29) or _slice(line, 65, 74)
    else:
        text = _slice(line, 22, 29)
    return _parse_value(text)


def _parse_halflife(text: str) -> float | None:
    """Convert ENSDF half-life field (cols 40-49) to seconds."""
    if not text or "STABLE" in text.upper():
        return None

    parts = text.split()
    unit = None
    value_parts: list[str] = []
    for p in parts:
        pu = p.upper()
        if pu in _HALFLIFE_TO_SECONDS:
            unit = pu
        else:
            value_parts.append(p)

    if not value_parts:
        return None

    value = _parse_value(" ".join(value_parts[:2]))
    if value is None:
        return None
    if unit is None:
        return value
    return value * _HALFLIFE_TO_SECONDS[unit]


def _parse_alpha(line: str) -> AlphaRay | None:
    energy = _parse_energy(line)
    if energy is None:
        return None
    return AlphaRay(energy_kev=energy)


def _parse_gamma(line: str) -> GammaRay | None:
    energy = _parse_energy(line)
    intensity = _parse_intensity(line, "G")
    if energy is None or intensity is None:
        return None
    return GammaRay(energy_kev=energy, intensity=intensity)


def _parse_beta(line: str) -> BetaRay | None:
    energy = _parse_energy(line)
    intensity = _parse_intensity(line, "B")
    if energy is None or intensity is None:
        return None
    return BetaRay(energy_kev=energy, intensity=intensity)


def _record_type(line: str) -> str | None:
    """Return ENSDF record letter in column 8, or None for comments / supplementals."""
    line = _pad80(line)
    if len(line) < 8:
        return None
    if line[6] == "c":
        return None
    if line[5] in "SHNQDCRdtx":
        return None
    rtype = line[7]
    if not rtype.isalpha() or rtype in "SHNQDCRdtx":
        return None
    return rtype


def _is_dataset_header(line: str) -> bool:
    line = _pad80(line)
    return line[5:7] == "  " and line[7] == " "


def _extract_decay_mode_from_title(title: str) -> str | None:
    m = _DECAY_MODE_RE.search(title.strip())
    if m:
        return m.group(2).upper().replace(" ", " ")
    m2 = re.search(r"(\S+)\s+DECAY", title, re.IGNORECASE)
    return m2.group(1).upper() if m2 else None


def _parse_nuclide_id(nucid: str) -> tuple[int, str]:
    nucid = nucid.strip()
    m = re.match(r"(\d+)([A-Za-z]*)", nucid)
    if not m:
        raise ValueError(f"Invalid NUCID: {nucid!r}")
    return int(m.group(1)), m.group(2).capitalize() or "?"


def _split_datasets(lines: list[str]) -> list[tuple[str, str, list[str]]]:
    datasets: list[tuple[str, str, list[str]]] = []
    current_title: str | None = None
    current_daughter: str | None = None
    current_body: list[str] = []

    def flush() -> None:
        nonlocal current_title, current_daughter, current_body
        if current_title and current_body:
            datasets.append((current_daughter or "", current_title, current_body))
        current_title = None
        current_daughter = None
        current_body = []

    for line in lines:
        if _is_dataset_header(line):
            m = _DECAY_DATASET_RE.match(line)
            if m:
                flush()
                current_daughter = m.group(1).strip()
                current_title = m.group(2).strip()
                current_body = []
            else:
                flush()
            continue
        if current_title is not None:
            current_body.append(line)

    flush()
    return datasets


def _parse_decay_dataset(
    daughter: str, title: str, lines: list[str]
) -> tuple[str, DecayMode] | None:
    if " DECAY" not in title.upper():
        return None

    mode_label = _extract_decay_mode_from_title(title) or "DECAY"
    parent_id: str | None = None
    half_life_s: float | None = None
    branch_items: list[tuple[str, float]] = []
    alphas: list[AlphaRay] = []
    betas: list[BetaRay] = []
    gammas: list[GammaRay] = []

    for line in lines:
        line = _pad80(line.rstrip("\n"))
        rtype = _record_type(line)
        if rtype is None:
            continue

        if rtype == "P":
            if line[5] != " " or line[6] != " ":
                continue
            parent_id = _slice(line, 1, 5)
            hl = _parse_halflife(_slice(line, 40, 49))
            if hl is not None:
                half_life_s = hl
        elif rtype == "L":
            if line[5] not in (" ", ""):
                text = line[9:].strip() if len(line) > 9 else ""
                for bm, val in _BRANCH_RE.findall(text):
                    v = _parse_value(val.replace("$", ""))
                    if v is not None:
                        branch_items.append((bm.upper(), v))
        elif rtype in _RADIATION_RECORDS:
            if line[5] != " " or line[6] != " ":
                continue
            if rtype == "A":
                alpha = _parse_alpha(line)
                if alpha is not None:
                    alphas.append(alpha)
            elif rtype == "G":
                gamma = _parse_gamma(line)
                if gamma is not None:
                    gammas.append(gamma)
            elif rtype == "B":
                beta = _parse_beta(line)
                if beta is not None:
                    betas.append(beta)

    if parent_id is None:
        m = _DECAY_MODE_RE.match(title.strip())
        if m:
            parent_id = m.group(1).upper()
        else:
            return None

    if half_life_s is None:
        m = re.search(r"\(([^)]+)\)", title)
        if m:
            half_life_s = _parse_halflife(m.group(1))

    if half_life_s is None:
        return None

    return parent_id.strip(), DecayMode(
        mode=mode_label,
        dataset_id=title,
        daughter_nuclide=daughter.strip(),
        half_life_s=half_life_s,
        branches=tuple(branch_items),
        alphas=tuple(alphas),
        betas=tuple(betas),
        gammas=tuple(gammas),
    )


def parse_ensdf_file(path: str | Path) -> dict[str, Radionuclide]:
    """
    Parse one ENSDF mass file and return radionuclides keyed by NUCID.

    Only decay datasets with a parseable half-life are included.
    Radiation records are kept only when all required fields for their type parse.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    pending: dict[str, list[DecayMode]] = {}

    for daughter, title, body in _split_datasets(lines):
        if " DECAY" not in title.upper():
            continue
        parsed = _parse_decay_dataset(daughter, title, body)
        if parsed is None:
            continue
        parent_id, mode = parsed
        pending.setdefault(parent_id, []).append(mode)

    nuclides: dict[str, Radionuclide] = {}
    for parent_id, modes in pending.items():
        a, sym = _parse_nuclide_id(parent_id)
        nuclides[parent_id] = Radionuclide(
            nuclide_id=parent_id,
            mass_number=a,
            symbol=sym,
            half_life_s=modes[0].half_life_s,
            decay_modes=tuple(modes),
        )
    return nuclides


def parse_ensdf_directory(path: str | Path) -> dict[str, Radionuclide]:
    """Parse all ensdf.* files in a directory."""
    path = Path(path)
    merged_modes: dict[str, list[DecayMode]] = {}
    for fp in sorted(path.glob("ensdf.*")):
        for nid, nuc in parse_ensdf_file(fp).items():
            merged_modes.setdefault(nid, []).extend(nuc.decay_modes)

    return {
        nid: Radionuclide(
            nuclide_id=nid,
            mass_number=_parse_nuclide_id(nid)[0],
            symbol=_parse_nuclide_id(nid)[1],
            half_life_s=modes[0].half_life_s,
            decay_modes=tuple(modes),
        )
        for nid, modes in merged_modes.items()
        if modes
    }


def iter_radionuclides(path: str | Path) -> Iterator[Radionuclide]:
    """Iterate all radionuclides from a file or directory."""
    path = Path(path)
    data = parse_ensdf_directory(path) if path.is_dir() else parse_ensdf_file(path)
    yield from data.values()
