"""Build the NIST X-ray mass-coefficient archive and Compounds enum."""

import json
import re
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from tqdm import tqdm

from pyradiate import logger
from pyradiate.core.elements import Elements
from pyradiate.tools.nist.coefficients import NIST_ARCHIVE_FILE, TABLE_DTYPE

NIST_BASE_URL = "https://physics.nist.gov/PhysRefData/XrayMassCoef"
NIST_ELEMENTS_INDEX_URL = f"{NIST_BASE_URL}/tab3.html"
NIST_COMPOUNDS_INDEX_URL = f"{NIST_BASE_URL}/tab4.html"

NIST_COMPOUNDS_MODULE = Path(__file__).parent / "compounds.py"

_SCI_RE = re.compile(r"\d+\.\d+E[+-]\d+", re.IGNORECASE)
_ELEM_HREF_RE = re.compile(r"ElemTab/z(\d+)\.html$", re.IGNORECASE)
_COMP_HREF_RE = re.compile(r"ComTab/([^/]+)\.html$", re.IGNORECASE)
_ENUM_UNSAFE_RE = re.compile(r"[^A-Za-z0-9]+")
_LEADING_DIGITS_RE = re.compile(r"^\d+")
_LETTER_UNDERSCORE_DIGITS_RE = re.compile(r"([A-Z])_(\d+)")


class _LinkParser(HTMLParser):
    """Collect (href, text) pairs from an HTML page."""

    def __init__(self):
        super().__init__()
        self.links: list[tuple[str, str]] = []
        self._href: str | None = None
        self._text: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self._href = href
                self._text = []

    def handle_data(self, data):
        if self._href is not None:
            self._text.append(data)

    def handle_endtag(self, tag):
        if tag == "a" and self._href is not None:
            self.links.append((self._href, "".join(self._text)))
            self._href = None


class _CompoundCellParser(HTMLParser):
    """Collect compound links together with the surrounding table-cell label."""

    def __init__(self):
        super().__init__()
        self.cells: list[tuple[str, str]] = []
        self._in_td = False
        self._href: str | None = None
        self._text: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag == "td":
            self._in_td = True
            self._href = None
            self._text = []
        elif tag == "a" and self._in_td:
            href = dict(attrs).get("href")
            if href and _COMP_HREF_RE.search(href):
                self._href = href

    def handle_data(self, data):
        if self._in_td:
            self._text.append(data)

    def handle_endtag(self, tag):
        if tag == "td" and self._in_td:
            if self._href is not None:
                self.cells.append((self._href, "".join(self._text)))
            self._in_td = False


def fetch_url(url: str) -> str:
    """Download *url* and return the decoded page body."""
    request = Request(url, headers={"User-Agent": "pyradiate (NIST X-ray mass coefficient builder)"})
    with urlopen(request) as response:
        assert response.status == 200, f"Incorrect return code {response.status} for {url}"
        raw = response.read()
    for encoding in ("utf-8", "latin-1"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def parse_element_index(html: str) -> list[tuple[Elements, str]]:
    """Return ``(Elements member, relative href)`` pairs from Table 3 HTML."""
    parser = _LinkParser()
    parser.feed(html)
    parser.close()

    elements: list[tuple[Elements, str]] = []
    seen: set[int] = set()
    for href, _text in parser.links:
        match = _ELEM_HREF_RE.search(href)
        if match is None:
            continue
        atomic_number = int(match.group(1))
        if atomic_number in seen:
            continue
        seen.add(atomic_number)
        elements.append((Elements.from_atomic_number(atomic_number), href))

    elements.sort(key=lambda item: item[0].atomic_number)
    return elements


def parse_compound_index(html: str) -> list[tuple[str, str, str]]:
    """Return ``(enum_name, display_name, nist_id)`` triples from Table 4 HTML."""
    parser = _CompoundCellParser()
    parser.feed(html)
    parser.close()

    compounds: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    for href, text in parser.cells:
        match = _COMP_HREF_RE.search(href)
        if match is None:
            continue
        nist_id = match.group(1)
        if nist_id in seen:
            continue
        seen.add(nist_id)
        display_name = " ".join(text.split())
        compounds.append((compound_enum_name(display_name, nist_id), display_name, nist_id))
    return uniquify_compound_enum_names(compounds)


def compound_enum_name(display_name: str, nist_id: str, *, keep_parens: bool = False) -> str:
    """Build a stable ``Compounds`` member name from a NIST display name."""
    name = display_name if keep_parens else re.sub(r"\([^)]*\)", "", display_name)
    name = _ENUM_UNSAFE_RE.sub("_", name)
    name = re.sub(r"_+", "_", name).strip("_").upper()
    name = _LEADING_DIGITS_RE.sub("", name).strip("_")
    name = _LETTER_UNDERSCORE_DIGITS_RE.sub(r"\1\2", name)
    name = re.sub(r"^(\d+_)*MMOL_L\d+_", "", name)
    if not name or not name[0].isalpha():
        name = nist_id.upper()
    return name


def uniquify_compound_enum_names(compounds: list[tuple[str, str, str]]) -> list[tuple[str, str, str]]:
    """Keep the first short name; later collisions keep parenthetical NIST notes."""
    counts = Counter(name for name, _display, _nist_id in compounds)
    seen_short: set[str] = set()
    unique: list[tuple[str, str, str]] = []
    for name, display_name, nist_id in compounds:
        if counts[name] > 1:
            if name in seen_short:
                name = compound_enum_name(display_name, nist_id, keep_parens=True)
            else:
                seen_short.add(name)
        unique.append((name, display_name, nist_id))

    still_duplicated = {name for name, count in Counter(item[0] for item in unique).items() if count > 1}
    if still_duplicated:
        unique = [
            (f"{name}_{nist_id.upper()}" if name in still_duplicated else name, display_name, nist_id)
            for name, display_name, nist_id in unique
        ]

    names = [name for name, _display, _nist_id in unique]
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate Compounds enum names: {sorted(n for n, c in Counter(names).items() if c > 1)}")
    return unique


def parse_coefficient_table(html: str) -> np.ndarray:
    """Parse a NIST element/compound page into a structured array of raw table rows."""
    pre_match = re.search(r"<PRE>(.*?)</PRE>", html, flags=re.IGNORECASE | re.DOTALL)
    if pre_match is None:
        raise ValueError("NIST coefficient page has no ASCII <PRE> table")

    plain = re.sub(r"<[^>]+>", "", pre_match.group(1))
    rows: list[tuple[float, float, float]] = []
    for line in plain.splitlines():
        numbers = _SCI_RE.findall(line)
        if len(numbers) < 3:
            continue
        energy_mev, mu_rho, mu_en_rho = (float(value) for value in numbers[:3])
        rows.append((energy_mev, mu_rho, mu_en_rho))

    if not rows:
        raise ValueError("No coefficient rows found in NIST ASCII table")
    return np.array(rows, dtype=TABLE_DTYPE)


def _write_compounds_module(compounds: list[tuple[str, str, str]], path: Path = NIST_COMPOUNDS_MODULE) -> None:
    lines = [
        '"""NIST Table 4 compounds. Generated by pyradiate.tools.nist.lib_builder. Do not edit."""',
        "",
        "from dataclasses import dataclass",
        "from enum import Enum",
        "",
        "",
        "@dataclass(frozen=True)",
        "class Compound:",
        "    name: str",
        "    nist_id: str",
        "",
        "",
        "class Compounds(Enum):",
        "    @property",
        "    def name(self):",
        "        return self.value.name",
        "",
        "    @property",
        "    def nist_id(self):",
        "        return self.value.nist_id",
        "",
    ]
    for enum_name, display_name, nist_id in sorted(compounds, key=lambda item: item[0]):
        lines.append(f"    {enum_name} = Compound({json.dumps(display_name)}, {json.dumps(nist_id)})")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_nist_archive(
    archive_file: Path = NIST_ARCHIVE_FILE,
    compounds_module: Path = NIST_COMPOUNDS_MODULE,
) -> Path:
    """
    Scrape NIST Tables 3 and 4 and write the numpy archive plus ``compounds.py``.

    Energies are stored in MeV as published by NIST. Coefficients are cm^2/g.
    """
    logger.info("Fetching NIST X-ray mass coefficient indexes...")
    element_entries = parse_element_index(fetch_url(NIST_ELEMENTS_INDEX_URL))
    compound_entries = parse_compound_index(fetch_url(NIST_COMPOUNDS_INDEX_URL))

    if not element_entries:
        raise RuntimeError(f"No element links found at {NIST_ELEMENTS_INDEX_URL}")
    if not compound_entries:
        raise RuntimeError(f"No compound links found at {NIST_COMPOUNDS_INDEX_URL}")

    payload: dict[str, np.ndarray] = {}
    element_symbols: list[str] = []
    for element, href in tqdm(element_entries, desc="NIST elements", unit="elem"):
        table = parse_coefficient_table(fetch_url(f"{NIST_BASE_URL}/{href}"))
        payload[f"element/{element.symbol}"] = table
        element_symbols.append(element.symbol)
        logger.debug(f"Parsed {element.symbol} ({len(table)} rows)")

    compound_nist_ids: list[str] = []
    compound_names: list[str] = []
    for enum_name, display_name, nist_id in tqdm(compound_entries, desc="NIST compounds", unit="comp"):
        table = parse_coefficient_table(fetch_url(f"{NIST_BASE_URL}/ComTab/{nist_id}.html"))
        payload[f"compound/{nist_id}"] = table
        compound_nist_ids.append(nist_id)
        compound_names.append(display_name)
        logger.debug(f"Parsed {enum_name} ({len(table)} rows)")

    payload["_element_symbols"] = np.asarray(element_symbols)
    payload["_compound_nist_ids"] = np.asarray(compound_nist_ids)
    payload["_compound_names"] = np.asarray(compound_names)
    payload["_energy_unit"] = np.asarray("MeV")
    payload["_coefficient_unit"] = np.asarray("cm2/g")
    payload["_source_elements"] = np.asarray(NIST_ELEMENTS_INDEX_URL)
    payload["_source_compounds"] = np.asarray(NIST_COMPOUNDS_INDEX_URL)

    archive_file.parent.mkdir(parents=True, exist_ok=True)
    np.savez(archive_file, **payload)
    _write_compounds_module(compound_entries, path=compounds_module)

    logger.info(
        f"Wrote NIST archive ({len(element_symbols)} elements, {len(compound_nist_ids)} compounds) to {archive_file}"
    )
    logger.info(f"Wrote Compounds enum to {compounds_module}")
    return archive_file


if __name__ == "__main__":
    build_nist_archive()
