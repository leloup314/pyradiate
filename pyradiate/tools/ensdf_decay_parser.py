"""Minimal parser for ENSDF decay datasets."""

import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from pyradiate import ensdf_path
from pyradiate.core import radiation
from pyradiate.core.nuclide import NuclideIdentifier, NuclideIdentifierError

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

_DATASET_HEADER_RE = re.compile(
    r"^\s*(\S{1,5})\s{4,}(.+?)\s{2,}\S",
    re.IGNORECASE,
)

_DECAY_DATASET_RE = re.compile(
    r"^\s*(\S{1,5})\s{4,}(.+?\b([A-Z][A-Z0-9+\-]*)\s+DECAY\b.*?)\s{2,}\S",
    re.IGNORECASE,
)

_DECAY_MODE_RE = re.compile(
    r"^(\d+\w*)\s+("
    + r"(?:"
    + r"B[+-]|EC|A|IT|SF|EP|"
    + r"B-N|B\+N|B-2N|B-3N|"
    + r"B\+EC|B\+A|B-EC|"
    + r"B\+B-|B\+B\+|"
    + r"B\+P|"
    + r"B\+EC\+A|"
    + r"EC\+B\+|EC\+A|"
    + r"A\+EC|A\+B\+|A\+B-|A\+EC\+B\+"
    + r")\s+DECAY"
    + r")",
    re.IGNORECASE,
)

_BRANCH_RE = re.compile(r"%([0-9]*[A-Z][A-Z0-9+\-]*)=([0-9.]+)")

_RADIATION_RECORDS = frozenset({"G", "B", "A"})

# tG table rows: "102.024 {I20}    23.2 {I14}  K|a{-2}| x ray"
_ENSDF_TABLE_NUMBER = re.compile(r"(\d+\.?\d*(?:[Ee][+-]?\d+)?)\s*(?:\{I[^}]+\})?")

_XRAY_BLOCK_START = re.compile(
    r"x\s*ray|x-ray",
    re.IGNORECASE,
)
_XRAY_BLOCK_END = re.compile(
    r"^\$[-=]{3,}|^\$Other |\|g-ray intensities|^\|a\|g|\|a\(|coincidence data",
    re.IGNORECASE,
)
_XRAY_SKIP_LINE = re.compile(
    r"^E\(x[- ]?ray\)|^I\(x[- ]?ray\)|^-{3,}|^\(\% per|intensity\s*$|^\s*energy\s*$",
    re.IGNORECASE,
)


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


def _parse_alpha(line: str) -> radiation.Alpha | None:
    energy = _parse_energy(line)
    if energy is None:
        return None
    return radiation.Alpha(energy_kev=energy)


def _parse_gamma(line: str) -> radiation.Gamma | None:
    energy = _parse_energy(line)
    intensity = _parse_intensity(line, "G")
    if energy is None or intensity is None:
        return None
    return radiation.Gamma(energy_kev=energy, intensity=intensity)


def _parse_beta(line: str) -> radiation.Beta | None:
    energy = _parse_energy(line)
    intensity = _parse_intensity(line, "B")
    if energy is None or intensity is None:
        return None
    return radiation.Beta(energy_kev=energy, intensity=intensity)


def _is_tg_line(line: str) -> bool:
    line = _pad80(line)
    return len(line) >= 8 and line[6] == "t" and line[7] == "G"


def _tg_body(line: str) -> str:
    return line[8:].strip()


def _table_numeric_region(body: str) -> str:
    """Limit number parsing to the table columns, not shell labels like K|b{-3}."""
    m = re.search(r"\s[KLMNO]\|", body)
    if m:
        return body[: m.start()]
    m = re.search(r"x\s*ray", body, re.I)
    if m:
        return body[: m.start()]
    return body[:50]


def _parse_table_number_spans(body: str) -> list[tuple[float, int, int]]:
    """Return (value, start, end) for each ENSDF number token in a tG row."""
    region = _table_numeric_region(body)
    spans: list[tuple[float, int, int]] = []
    for m in _ENSDF_TABLE_NUMBER.finditer(region):
        v = _parse_value(m.group(1))
        if v is not None:
            spans.append((v, m.start(), m.end()))
    return spans


def _extract_xray_label(body: str) -> str:
    cleaned = re.sub(r"\{I[^}]+\}", "", body).replace("$", "")
    shell = re.search(
        r"(?:K\|[^+\s]+|L\|[^+\s]+|M\|[^+\s]+|N\|[^+\s]+|K-O\{[^}]+\})"
        r"(?:\s*\+\s*(?:K\|[^+\s]+|L\|[^+\s]+|K-O\{[^}]+\}))*"
        r"(?:\s+x\s*ray)?",
        cleaned,
        re.I,
    )
    if shell:
        return shell.group(0).strip()
    if re.search(r"x\s*ray", cleaned, re.I):
        return "x ray"
    return ""


def _tg_starts_xray_block(body: str) -> bool:
    if _XRAY_BLOCK_END.search(body):
        return False
    if body.startswith("$") and _XRAY_BLOCK_START.search(body):
        return True
    if re.search(r"E\(x[- ]?ray\)\s+I\(x[- ]?ray\)", body, re.I):
        return True
    if re.search(r"x\s*ray.*measured|measured.*x\s*ray", body, re.I):
        return True
    return False


def _tg_ends_xray_block(body: str) -> bool:
    return bool(_XRAY_BLOCK_END.search(body))


def _parse_tg_xray_rows(lines: list[str]) -> list[radiation.Xray]:
    """Parse X-ray lines from tG table blocks inside a decay dataset."""
    xrays: list[radiation.Xray] = []
    in_block = False
    pending_energy: float | None = None
    pending_label = ""

    for line in lines:
        if not _is_tg_line(line):
            if in_block:
                in_block = False
                pending_energy = None
            continue

        body = _tg_body(line)
        if not body:
            continue

        if _tg_ends_xray_block(body):
            in_block = False
            pending_energy = None
            continue

        if _tg_starts_xray_block(body):
            in_block = True
            pending_energy = None
            continue

        if not in_block:
            continue

        if _XRAY_SKIP_LINE.search(body) or body.startswith("("):
            continue
        if set(body) <= {"-", " "}:
            continue

        spans = _parse_table_number_spans(body)
        values = [s[0] for s in spans]

        if len(spans) >= 2:
            label = _extract_xray_label(body)
            if not label:
                label = "x ray"
            xrays.append(
                radiation.Xray(
                    energy_kev=spans[0][0],
                    intensity=spans[1][0],
                    label=label,
                )
            )
            pending_energy = None
            continue

        if len(spans) == 1 and pending_energy is not None:
            value = spans[0][0]
            label = _extract_xray_label(body)
            # Continuation row: intensity only (typically much smaller than X-ray energy)
            if value < pending_energy:
                xrays.append(
                    radiation.Xray(
                        energy_kev=pending_energy,
                        intensity=value,
                        label=label or pending_label or "x ray",
                    )
                )
                pending_energy = None
                pending_label = ""
            elif label:
                pending_energy = value
                pending_label = label
            continue

        if len(spans) == 1 and pending_energy is None:
            label = _extract_xray_label(body)
            if label:
                pending_energy = spans[0][0]
                pending_label = label
            continue

    return xrays


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


@dataclass(frozen=True)
class _AdoptedLevelData:
    half_life_s: float | None
    branches: tuple[tuple[str, float], ...]


@dataclass(frozen=True)
class _AdoptedNuclideData:
    parent_id: str
    half_life_s: float | None
    branches: tuple[tuple[str, float], ...]
    levels: tuple[_AdoptedLevelData, ...]
    decay_xref_titles: frozenset[str]


@dataclass(frozen=True)
class ParsedRadioNuclide:
    decays: tuple[radiation.Decay, ...]
    recommended_half_life_s: float


@dataclass(frozen=True)
class _EnsdfFileData:
    decays_by_parent: dict[str, list[radiation.Decay]]
    adopted_map: dict[str, _AdoptedNuclideData]
    adopted_xrefs: frozenset[str]


def _extract_branches_from_text(text: str) -> list[tuple[str, float]]:
    branches: list[tuple[str, float]] = []
    for mode, value_s in _BRANCH_RE.findall(text):
        if re.search(r"[\$\{\|\(;<]", mode):
            continue
        value = _parse_value(value_s)
        if value is None:
            continue
        try:
            branch_id = radiation.DecayBranch.from_string(mode)
        except ValueError:
            continue
        branches.append((branch_id.value, value))
    return branches


def _is_xref_line(line: str) -> bool:
    line = _pad80(line)
    return len(line) >= 9 and line[7] == "X" and line[8].isalpha()


def _parse_xref_body(line: str) -> str:
    return _pad80(line)[9:].strip()


def _normalize_decay_title(title: str) -> str:
    return re.sub(r"\s+", " ", title.upper().strip())


def _is_ground_state_level(line: str) -> bool:
    energy = _parse_value(_slice(_pad80(line), 10, 19))
    return energy is None or energy < 1.0


def _parse_adopted_dataset(nucid: str, title: str, lines: list[str]) -> _AdoptedNuclideData | None:
    if "ADOPTED" not in title.upper():
        return None

    try:
        parent_id = NuclideIdentifier.from_string(nucid.strip()).identifier
    except NuclideIdentifierError:
        return None

    half_life_s: float | None = None
    ground_state_branches: list[tuple[str, float]] = []
    decay_xrefs: set[str] = set()
    levels: list[dict[str, object]] = []
    current_level_idx: int | None = None

    for line in lines:
        line = _pad80(line.rstrip("\n"))

        if _is_xref_line(line):
            xref = _parse_xref_body(line)
            if " DECAY" in xref.upper():
                decay_xrefs.add(_normalize_decay_title(xref))
            continue

        rtype = _record_type(line)
        if rtype is None:
            continue

        if rtype == "P" and line[5] == " " and line[6] == " ":
            hl = _parse_halflife(_slice(line, 40, 49))
            if hl is not None:
                half_life_s = hl
        elif rtype == "L":
            if line[5] in (" ", ""):
                hl = _parse_halflife(_slice(line, 40, 49))
                levels.append({"half_life_s": hl, "branches": [], "is_ground_state": _is_ground_state_level(line)})
                current_level_idx = len(levels) - 1
                if levels[-1]["is_ground_state"] and hl is not None:
                    half_life_s = hl
            elif line[5] not in (" ", "") and current_level_idx is not None:
                text = line[9:].strip() if len(line) > 9 else ""
                level_branches = _extract_branches_from_text(text)
                if level_branches:
                    level = levels[current_level_idx]
                    level["branches"] = [*level["branches"], *level_branches]
                    if level["is_ground_state"]:
                        ground_state_branches.extend(level_branches)

    adopted_levels = tuple(
        _AdoptedLevelData(
            half_life_s=level["half_life_s"],
            branches=tuple(level["branches"]),
        )
        for level in levels
        if level["branches"] or level["half_life_s"] is not None
    )

    return _AdoptedNuclideData(
        parent_id=parent_id,
        half_life_s=half_life_s,
        branches=tuple(ground_state_branches),
        levels=adopted_levels,
        decay_xref_titles=frozenset(decay_xrefs),
    )


def _split_datasets(lines: list[str]) -> list[tuple[str, str, list[str]]]:
    datasets: list[tuple[str, str, list[str]]] = []
    current_title: str | None = None
    current_nucid: str | None = None
    current_body: list[str] = []

    def flush() -> None:
        nonlocal current_title, current_nucid, current_body
        if current_title and current_body:
            datasets.append((current_nucid or "", current_title, current_body))
        current_title = None
        current_nucid = None
        current_body = []

    for line in lines:
        if _is_dataset_header(line):
            m = _DATASET_HEADER_RE.match(line)
            if m:
                flush()
                current_nucid = m.group(1).strip()
                current_title = m.group(2).strip()
                current_body = []
            else:
                flush()
            continue
        if current_title is not None:
            current_body.append(line)

    flush()
    return datasets


def _decay_identity_key(parent_id: str, decay: radiation.Decay) -> tuple:
    mode = _extract_decay_mode_from_title(decay.mode) or decay.mode
    return (parent_id, mode, decay.daughter_nuclide.identifier, round(decay.half_life_s, 6))


def _radiation_count(decay: radiation.Decay) -> int:
    return len(decay.alphas) + len(decay.betas) + len(decay.gammas) + len(decay.xrays)


def _select_preferred_decay(
    candidates: list[tuple[str, radiation.Decay, str, int]],
    adopted_xrefs: frozenset[str],
) -> tuple[str, radiation.Decay, str, int]:
    def score(item: tuple[str, radiation.Decay, str, int]) -> int:
        _, _, title, rad_count = item
        value = rad_count
        if _normalize_decay_title(title) in adopted_xrefs:
            value += 10000
        if re.search(r"DATA SET\s*#\d", title, re.IGNORECASE):
            value -= 50
        return value

    return max(candidates, key=score)


def _hl_close(a: float, b: float, rtol: float = 0.05) -> bool:
    if a <= 0 or b <= 0:
        return False
    ratio = a / b
    return (1 - rtol) <= ratio <= (1 + rtol)


def _apply_adopted_data(decay: radiation.Decay, adopted: _AdoptedNuclideData) -> radiation.Decay:
    if adopted.half_life_s is None or not _hl_close(decay.half_life_s, adopted.half_life_s):
        return decay

    return radiation.Decay(
        mode=decay.mode,
        dataset_id=decay.dataset_id,
        daughter_nuclide=decay.daughter_nuclide,
        half_life_s=adopted.half_life_s,
        branch=decay.branch,
        branch_fraction=decay.branch_fraction,
        alphas=decay.alphas,
        betas=decay.betas,
        gammas=decay.gammas,
        xrays=decay.xrays,
    )


def _decay_mode_key(mode: str) -> str:
    m = _DECAY_MODE_RE.search(mode.strip())
    if m:
        return m.group(1).upper()
    return re.sub(r"\s+DECAY.*", "", mode, flags=re.IGNORECASE).strip().upper()


_BRANCH_MODE_ALIASES: dict[str, tuple[str, ...]] = {
    "EC+B+": ("EC+%B+", "EC+B+", "B+"),
    "EC+B-": ("EC+%B-", "EC+B-"),
    "EC": ("EC", "EC+%B+", "EC+%B-"),
}


def _lookup_branch_fraction(mode_key: str, level_branches: tuple[tuple[str, float], ...]) -> float | None:
    branch_map = {bm.upper(): val for bm, val in level_branches}
    if mode_key in branch_map:
        return branch_map[mode_key]
    for alias_key, aliases in _BRANCH_MODE_ALIASES.items():
        if mode_key == alias_key:
            for name in aliases:
                if name in branch_map:
                    return branch_map[name]
    return None


def _level_branches_for_decay(
    decay: radiation.Decay,
    adopted: _AdoptedNuclideData | None,
) -> tuple[tuple[str, float], ...]:
    if adopted is None:
        return ()

    for level in adopted.levels:
        if level.half_life_s is not None and _hl_close(decay.half_life_s, level.half_life_s):
            if level.branches:
                return level.branches

    if adopted.half_life_s is not None and _hl_close(decay.half_life_s, adopted.half_life_s):
        return adopted.branches

    return ()


def _finalize_decay_branches(
    decay: radiation.Decay,
    adopted: _AdoptedNuclideData | None,
    *,
    sibling_count: int,
) -> radiation.Decay:
    branch = radiation.DecayBranch.from_string(_decay_mode_key(decay.mode))
    level_branches = _level_branches_for_decay(decay, adopted)
    branch_fraction = _lookup_branch_fraction(branch.value, level_branches)
    if branch_fraction is None and sibling_count == 1:
        branch_fraction = 100.0
    if branch_fraction is None:
        branch_fraction = decay.branch_fraction

    return radiation.Decay(
        mode=decay.mode,
        dataset_id=decay.dataset_id,
        daughter_nuclide=decay.daughter_nuclide,
        half_life_s=decay.half_life_s,
        branch=branch,
        branch_fraction=branch_fraction,
        alphas=decay.alphas,
        betas=decay.betas,
        gammas=decay.gammas,
        xrays=decay.xrays,
    )


def _finalize_parent_decays(
    decays: list[radiation.Decay],
    adopted: _AdoptedNuclideData | None,
) -> list[radiation.Decay]:
    siblings_by_hl: dict[float, int] = {}
    for decay in decays:
        hl_key = round(decay.half_life_s, 6)
        siblings_by_hl[hl_key] = siblings_by_hl.get(hl_key, 0) + 1

    return [
        _finalize_decay_branches(
            decay,
            adopted,
            sibling_count=siblings_by_hl[round(decay.half_life_s, 6)],
        )
        for decay in decays
    ]


def _recommended_half_life_s(
    decays: list[radiation.Decay],
    adopted: _AdoptedNuclideData | None,
) -> float:
    if adopted is not None and adopted.half_life_s is not None:
        return adopted.half_life_s
    return decays[0].half_life_s


def _dedupe_decays_for_parent(
    parent_id: str,
    decays: list[radiation.Decay],
    adopted_xrefs: frozenset[str],
) -> list[radiation.Decay]:
    grouped: dict[tuple, list[radiation.Decay]] = {}
    for decay in decays:
        key = _decay_identity_key(parent_id, decay)
        grouped.setdefault(key, []).append(decay)

    deduped: list[radiation.Decay] = []
    for group in grouped.values():
        if len(group) == 1:
            deduped.append(group[0])
            continue
        candidates = [(parent_id, decay, decay.dataset_id, _radiation_count(decay)) for decay in group]
        _, best, _, _ = _select_preferred_decay(candidates, adopted_xrefs)
        deduped.append(best)
    return deduped


def _parse_decay_dataset(daughter: str, title: str, lines: list[str]) -> tuple[str, radiation.Decay] | None:
    if " DECAY" not in title.upper():
        return None

    mode_label = _extract_decay_mode_from_title(title) or "DECAY"
    branch = radiation.DecayBranch.from_string(_decay_mode_key(mode_label))
    parent_id: str | None = None
    half_life_s: float | None = None
    alphas: list[radiation.Alpha] = []
    betas: list[radiation.Beta] = []
    gammas: list[radiation.Gamma] = []
    xrays = _parse_tg_xray_rows(lines)

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
    try:
        nuclid_identifier = NuclideIdentifier.from_string(parent_id.strip()).identifier
        return nuclid_identifier, radiation.Decay(
            mode=mode_label,
            dataset_id=title,
            daughter_nuclide=NuclideIdentifier.from_string(daughter.strip()),
            half_life_s=half_life_s,
            branch=branch,
            branch_fraction=100.0,
            alphas=tuple(alphas),
            betas=tuple(betas),
            gammas=tuple(gammas),
            xrays=tuple(xrays),
        )
    except NuclideIdentifierError:
        return None


def parse_radio_nuclide(nuclide_identifier: NuclideIdentifier) -> ParsedRadioNuclide:
    mass_data = parse_ensdf_directory(ensdf_path)
    meta = mass_data.get(nuclide_identifier.identifier)
    if meta is None:
        raise NuclideIdentifierError(f"No decay data for {nuclide_identifier.identifier}")
    return meta


@cache
def parse_ensdf_file(path: str | Path) -> _EnsdfFileData:
    """
    Parse one ENSDF mass file.

    Returns per-parent decay lists (not yet branch-finalized) and any adopted-level
    data found in this file.  Decay datasets for a parent may also appear in other
    mass files (keyed by the daughter nucleus); use parse_ensdf_directory for the
    complete merged result.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    adopted_map: dict[str, _AdoptedNuclideData] = {}
    adopted_xrefs: set[str] = set()
    raw_decays: list[tuple[str, radiation.Decay, str, int]] = []

    for nucid, title, body in _split_datasets(lines):
        title_upper = title.upper()
        if "ADOPTED" in title_upper:
            adopted = _parse_adopted_dataset(nucid, title, body)
            if adopted is not None:
                adopted_map[adopted.parent_id] = adopted
                adopted_xrefs.update(adopted.decay_xref_titles)
        elif " DECAY" in title_upper:
            parsed = _parse_decay_dataset(nucid, title, body)
            if parsed is None:
                continue
            parent_id, decay = parsed
            raw_decays.append((parent_id, decay, title, _radiation_count(decay)))

    xref_set = frozenset(adopted_xrefs)
    normalized_decays: list[tuple[str, radiation.Decay, str, int]] = []
    for parent_id, decay, title, rad_count in raw_decays:
        if parent_id in adopted_map:
            decay = _apply_adopted_data(decay, adopted_map[parent_id])
        normalized_decays.append((parent_id, decay, title, rad_count))

    grouped: dict[tuple, list[tuple[str, radiation.Decay, str, int]]] = {}
    for item in normalized_decays:
        parent_id, decay, title, rad_count = item
        key = _decay_identity_key(parent_id, decay)
        grouped.setdefault(key, []).append(item)

    decays_by_parent: dict[str, list[radiation.Decay]] = {}
    for candidates in grouped.values():
        parent_id, decay, _, _ = _select_preferred_decay(candidates, xref_set)
        decays_by_parent.setdefault(parent_id, []).append(decay)

    return _EnsdfFileData(
        decays_by_parent=decays_by_parent,
        adopted_map=adopted_map,
        adopted_xrefs=xref_set,
    )


def parse_ensdf_directory(path: str | Path) -> dict[str, ParsedRadioNuclide]:
    """Parse all ensdf.* files in a directory."""
    path = Path(path)
    all_decays: dict[str, list[radiation.Decay]] = {}
    adopted_map: dict[str, _AdoptedNuclideData] = {}
    adopted_xrefs: set[str] = set()

    for fp in sorted(path.glob("ensdf.*")):
        file_data = parse_ensdf_file(fp)
        adopted_map.update(file_data.adopted_map)
        adopted_xrefs.update(file_data.adopted_xrefs)
        for parent_id, decays in file_data.decays_by_parent.items():
            all_decays.setdefault(parent_id, []).extend(decays)

    xref_set = frozenset(adopted_xrefs)
    result: dict[str, ParsedRadioNuclide] = {}
    for parent_id, parent_decays in all_decays.items():
        adopted = adopted_map.get(parent_id)
        normalized = [_apply_adopted_data(decay, adopted) if adopted is not None else decay for decay in parent_decays]
        deduped = _dedupe_decays_for_parent(parent_id, normalized, xref_set)
        finalized = _finalize_parent_decays(deduped, adopted)
        result[parent_id] = ParsedRadioNuclide(
            decays=tuple(finalized),
            recommended_half_life_s=_recommended_half_life_s(finalized, adopted),
        )

    return result
