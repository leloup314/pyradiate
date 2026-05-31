"""Minimal parser for ENSDF decay datasets."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

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

# Dataset titles like "89SE B- DECAY (0.43 S)" or "288FL A DECAY (0.64 S)"
_DECAY_DATASET_RE = re.compile(
    r"^\s*(\S{1,5})\s{4,}(.+?\b([A-Z][A-Z0-9+\-]*)\s+DECAY\b.*?)\s{2,}\S",
    re.IGNORECASE,
)

# Parent / mode from title: "89ZR EC DECAY"
_DECAY_MODE_RE = re.compile(
    r"^(\d+\w*)\s+((?:B[+-]|EC|A|IT|SF|EP|B-N|B\+N|B-2N|B-3N|B\+EC|B\+A|"
    r"B-EC|B\+B-|B-DECAY|B\+B\+|B\+A|B\+N|B\+P|B\+EC\+A|"
    r"EC\+B\+|EC\+A|A\+EC|A\+B\+|A\+B-|A\+EC\+B\+)\s+DECAY)",
    re.IGNORECASE,
)

# Continuation lines with branching: %B-=100$ %EC=50
_BRANCH_RE = re.compile(r"%([A-Z][A-Z0-9+\-]*)=([^\s$%]+)")

# Numeric token: value, optional uncertainty, optional unit / qualifier
_NUM_RE = re.compile(
    r"^([<>]?\s*\d+\.?\d*(?:[Ee][+-]?\d+)?|AP|SY|STABLE|UNSTABLE|"
    r"\d+\.\d+\+\d+-\d+|\d+\+\d+-\d+)"
    r"(?:\s+([A-Z][A-Z0-9+\-]*|\d+))?"
    r"(?:\s+([A-Z]{1,3}|SY|AP))?",
    re.IGNORECASE,
)

_RADIATION_RECORDS = frozenset({"G", "B", "A", "E"})


@dataclass(frozen=True)
class Radiation:
    """Alpha, beta, or gamma line from an ENSDF decay dataset."""

    kind: str  # "alpha", "beta", "gamma"
    energy_kev: Optional[float]
    intensity: Optional[float]
    energy_uncertainty_kev: Optional[float] = None
    intensity_uncertainty: Optional[float] = None
    raw_energy: Optional[str] = None
    raw_intensity: Optional[str] = None


@dataclass
class DecayMode:
    """One decay pathway (e.g. EC decay of 89Zr)."""

    mode: str
    dataset_id: str
    daughter_nuclide: str
    half_life_s: Optional[float]
    branches: dict[str, float] = field(default_factory=dict)
    _alphas: list[Radiation] = field(default_factory=list)
    _betas: list[Radiation] = field(default_factory=list)
    _gammas: list[Radiation] = field(default_factory=list)

    @property
    def alphas(self) -> tuple[Radiation, ...]:
        return tuple(self._alphas)

    @property
    def betas(self) -> tuple[Radiation, ...]:
        return tuple(self._betas)

    @property
    def gammas(self) -> tuple[Radiation, ...]:
        return tuple(self._gammas)

    def __iter__(self) -> Iterator[Radiation]:
        """Iterate all radiation (alphas, then betas, then gammas)."""
        yield from self._alphas
        yield from self._betas
        yield from self._gammas


@dataclass
class Radionuclide:
    """Radionuclide with one or more evaluated decay modes."""

    nuclide_id: str
    mass_number: int
    symbol: str
    half_life_s: Optional[float] = None
    _modes: list[DecayMode] = field(default_factory=list)

    @property
    def decay_modes(self) -> tuple[DecayMode, ...]:
        return tuple(self._modes)

    def __iter__(self) -> Iterator[DecayMode]:
        yield from self._modes


def _pad80(line: str) -> str:
    return line if len(line) >= 80 else line.ljust(80)


def _slice(line: str, start: int, end: int) -> str:
    """1-based inclusive ENSDF column slice."""
    return line[start - 1 : end].strip()


def _parse_numeric_field(text: str) -> tuple[Optional[float], Optional[float], Optional[str]]:
    """Parse an ENSDF numeric field (energy, intensity, half-life fragment)."""
    if not text:
        return None, None, None
    raw = text.strip()
    upper = raw.upper()
    if upper in ("", "?", "UNKNOWN") or "STABLE" in upper:
        return None, None, raw

    # Strip trailing qualifiers on the line (coincidence flags etc.)
    token = raw.split()[0] if raw.split() else ""
    if token.upper() in ("?", "AP", "SY"):
        return None, None, raw

    # Asymmetric uncertainty: 108+20-14
    m = re.match(r"^(\d+\.?\d*(?:[Ee][+-]?\d+)?)\+(\d+)-(\d+)$", token.replace(" ", ""))
    if m:
        return float(m.group(1)), float(m.group(2)), raw

    parts = raw.split()
    value_s = parts[0].lstrip("<>")
    try:
        value = float(value_s)
    except ValueError:
        return None, None, raw

    uncertainty = None
    if len(parts) >= 2:
        try:
            uncertainty = float(parts[1])
        except ValueError:
            pass

    # Energy multiplier suffix glued to number (e.g. 9.21E3)
    if len(parts) >= 2 and parts[1].upper() in ("E2", "E3", "E4", "E5", "E6"):
        exp = int(parts[1][1])
        value *= 10**exp

    return value, uncertainty, raw


def _parse_halflife(text: str) -> Optional[float]:
    """Convert ENSDF half-life field (cols 40-49) to seconds."""
    if not text:
        return None
    raw = text.strip()
    if "STABLE" in raw.upper():
        return None

    parts = raw.split()
    if not parts:
        return None

    # Find unit token (last alphabetic token)
    unit = None
    value_parts: list[str] = []
    for p in parts:
        pu = p.upper().rstrip("0123456789")
        if pu in _HALFLIFE_TO_SECONDS:
            unit = pu
        elif re.match(r"^[A-Z]{1,3}$", p.upper()) and p.upper() in _HALFLIFE_TO_SECONDS:
            unit = p.upper()
        else:
            value_parts.append(p)

    if not value_parts:
        return None

    value, _, _ = _parse_numeric_field(" ".join(value_parts[:2]))
    if value is None:
        return None

    if unit is None:
        # Bare number — assume seconds if no unit (unusual)
        return value

    return value * _HALFLIFE_TO_SECONDS[unit]


def _parse_energy(line: str) -> tuple[Optional[float], Optional[float], Optional[str]]:
    """Cols 10-19 (E), 20-21 (DE)."""
    e_text = _slice(line, 10, 19)
    de_text = _slice(line, 20, 21)
    value, unc, raw = _parse_numeric_field(e_text)

    # E3-style multiplier in column 19
    if e_text and value is not None:
        m = re.search(r"E(\d)\s*$", e_text, re.IGNORECASE)
        if m:
            value *= 10 ** int(m.group(1))

    if value is not None and de_text:
        try:
            de = float(de_text)
            # ENSDF uncertainty applies to last digits
            if de < value:
                unc = de
        except ValueError:
            pass

    return value, unc, raw or e_text or None


def _parse_intensity(line: str, record: str) -> tuple[Optional[float], Optional[float], Optional[str]]:
    """IB (beta/alpha) cols 22-29; RI/TI for gamma."""
    if record == "G":
        ri = _slice(line, 22, 29)
        ti = _slice(line, 65, 74)
        text = ri or ti
    else:
        text = _slice(line, 22, 29)
    return _parse_numeric_field(text)


def _record_type(line: str) -> Optional[str]:
    """Return ENSDF record letter in column 8, or None for comments / supplementals."""
    line = _pad80(line)
    if len(line) < 8:
        return None
    # Comment record: col 7 (index 6) is 'c'
    if line[6] == "c":
        return None
    # Supplemental / non-body records (col 6)
    if line[5] in "SHNQDCRdtx":
        return None
    rtype = line[7]
    if not rtype.isalpha():
        return None
    if rtype in "SHNQDCRdtx":
        return None
    return rtype


def _is_dataset_header(line: str) -> bool:
    """Identification record: col 8 blank, cols 6-7 blank."""
    line = _pad80(line)
    return line[5:7] == "  " and line[7] == " "


def _extract_decay_mode_from_title(title: str) -> Optional[str]:
    m = _DECAY_MODE_RE.search(title.strip())
    if m:
        return m.group(2).upper().replace(" ", " ")
    # Fallback: token before DECAY
    m2 = re.search(r"(\S+)\s+DECAY", title, re.IGNORECASE)
    return m2.group(1).upper() if m2 else None


def _parse_nuclide_id(nucid: str) -> tuple[int, str]:
    nucid = nucid.strip()
    m = re.match(r"(\d+)([A-Za-z]*)", nucid)
    if not m:
        raise ValueError(f"Invalid NUCID: {nucid!r}")
    return int(m.group(1)), m.group(2).capitalize() or "?"


def _split_datasets(lines: list[str]) -> list[tuple[str, str, list[str]]]:
    """Return (daughter_nucid, dataset_title, body_lines)."""
    datasets: list[tuple[str, str, list[str]]] = []
    current_title: Optional[str] = None
    current_daughter: Optional[str] = None
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
                continue
            # Non-decay dataset header ends previous decay dataset
            if current_title and " DECAY" in current_title.upper():
                flush()
            else:
                flush()
            continue

        if current_title is not None:
            current_body.append(line)

    flush()
    return datasets


def _parse_decay_dataset(
    daughter: str, title: str, lines: list[str]
) -> Optional[tuple[str, DecayMode]]:
    if " DECAY" not in title.upper():
        return None

    mode_label = _extract_decay_mode_from_title(title) or "DECAY"
    parent_id: Optional[str] = None
    half_life_s: Optional[float] = None
    branches: dict[str, float] = {}
    alphas: list[Radiation] = []
    betas: list[Radiation] = []
    gammas: list[Radiation] = []

    for line in lines:
        line = _pad80(line.rstrip("\n"))
        rtype = _record_type(line)
        if rtype is None:
            continue

        nucid = _slice(line, 1, 5)

        if rtype == "P":
            if line[5] != " " or line[6] != " ":
                continue
            parent_id = nucid
            t_field = _slice(line, 40, 49)
            hl = _parse_halflife(t_field)
            if hl is not None:
                half_life_s = hl
        elif rtype == "L":
            # Branching fractions on continuation L records
            cont = line[5] if len(line) > 5 else " "
            if cont not in (" ", ""):
                text = line[9:].strip() if len(line) > 9 else ""
                for bm, val in _BRANCH_RE.findall(text):
                    v, _, _ = _parse_numeric_field(val.replace("$", ""))
                    if v is not None:
                        branches[bm.upper()] = v
        elif rtype in _RADIATION_RECORDS:
            # Only standard data records: cols 6-7 blank
            if line[5] != " " or line[6] != " ":
                continue
            energy, e_unc, raw_e = _parse_energy(line)
            intensity, i_unc, raw_i = _parse_intensity(line, rtype)

            if energy is None and intensity is None:
                continue

            if rtype == "G":
                kind = "gamma"
            elif rtype == "A":
                kind = "alpha"
            elif rtype == "B":
                kind = "beta"
            else:
                continue

            rad = Radiation(
                kind=kind,
                energy_kev=energy,
                intensity=intensity,
                energy_uncertainty_kev=e_unc,
                intensity_uncertainty=i_unc,
                raw_energy=raw_e,
                raw_intensity=raw_i,
            )
            if kind == "alpha":
                alphas.append(rad)
            elif kind == "beta":
                betas.append(rad)
            else:
                gammas.append(rad)

    if parent_id is None:
        # Infer parent from title (first token before mode)
        m = _DECAY_MODE_RE.match(title.strip())
        if m:
            parent_id = m.group(1).upper()
        else:
            return None

    if half_life_s is None:
        m = re.search(r"\(([^)]+)\)\s*$", title)
        if m:
            half_life_s = _parse_halflife(m.group(1))

    decay = DecayMode(
        mode=mode_label,
        dataset_id=title,
        daughter_nuclide=daughter.strip(),
        half_life_s=half_life_s,
        branches=branches,
    )
    decay._alphas = alphas
    decay._betas = betas
    decay._gammas = gammas
    return parent_id.strip(), decay


def parse_ensdf_file(path: str | Path) -> dict[str, Radionuclide]:
    """
    Parse one ENSDF mass file and return radionuclides keyed by NUCID.

    Only decay datasets (with a parent record) are included.
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    nuclides: dict[str, Radionuclide] = {}

    for daughter, title, body in _split_datasets(lines):
        if " DECAY" not in title.upper():
            continue
        parsed = _parse_decay_dataset(daughter, title, body)
        if parsed is None:
            continue
        parent_id, mode = parsed

        if parent_id not in nuclides:
            a, sym = _parse_nuclide_id(parent_id)
            nuclides[parent_id] = Radionuclide(
                nuclide_id=parent_id,
                mass_number=a,
                symbol=sym,
            )

        nuc = nuclides[parent_id]
        nuc._modes.append(mode)
        if mode.half_life_s is not None and nuc.half_life_s is None:
            nuc.half_life_s = mode.half_life_s

    return nuclides


def parse_ensdf_directory(path: str | Path) -> dict[str, Radionuclide]:
    """Parse all ensdf.* files in a directory."""
    path = Path(path)
    merged: dict[str, Radionuclide] = {}
    for fp in sorted(path.glob("ensdf.*")):
        for nid, nuc in parse_ensdf_file(fp).items():
            if nid in merged:
                merged[nid]._modes.extend(nuc._modes)
                if merged[nid].half_life_s is None:
                    merged[nid].half_life_s = nuc.half_life_s
            else:
                merged[nid] = nuc
    return merged


def iter_radionuclides(path: str | Path) -> Iterator[Radionuclide]:
    """Iterate all radionuclides from a file or directory."""
    path = Path(path)
    if path.is_dir():
        data = parse_ensdf_directory(path)
    else:
        data = parse_ensdf_file(path)
    yield from data.values()
