from functools import cache
from pathlib import Path

import numpy as np

from pyradiate.core.elements import Elements
from pyradiate.tools.nist.compounds import Compounds

NIST_ARCHIVE_FILE = Path(__file__).parent / "data" / "xray_mass_coef.npz"

# Raw NIST columns: Energy (MeV), μ/ρ (cm²/g), μ_en/ρ (cm²/g)
TABLE_DTYPE = np.dtype(
    [
        ("energy_mev", np.float64),
        ("mu_rho", np.float64),
        ("mu_en_rho", np.float64),
    ]
)

_KEV_TO_MEV = 1e-3
_MU_RHO = "mu_rho"
_MU_EN_RHO = "mu_en_rho"


@cache
def _archive() -> dict[str, np.ndarray]:
    if not NIST_ARCHIVE_FILE.is_file():
        raise FileNotFoundError(
            f"NIST coefficient archive not found at {NIST_ARCHIVE_FILE}. Run python -m pyradiate.tools.nist.lib_builder"
        )
    with np.load(NIST_ARCHIVE_FILE) as data:
        return {key: data[key] for key in data.files}


def _material_key(material: Elements | Compounds) -> str:
    if isinstance(material, Elements):
        return f"element/{material.symbol}"
    if isinstance(material, Compounds):
        return f"compound/{material.nist_id}"
    raise TypeError(f"material must be Elements or Compounds, got {type(material)!r}")


def _material_label(material: Elements | Compounds) -> str:
    if isinstance(material, Elements):
        return material.symbol
    return material.name


def _table(material: Elements | Compounds) -> np.ndarray:
    key = _material_key(material)
    try:
        return _archive()[key]
    except KeyError:
        raise ValueError(f"No NIST X-ray mass coefficient data for {_material_label(material)}") from None


def load_table(material: Elements | Compounds) -> np.ndarray:
    """
    Raw NIST coefficient table for *material* as a structured array.

    Fields are ``energy_mev``, ``mu_rho`` and ``mu_en_rho`` (cm²/g), as published.
    """
    return _table(material).copy()


def _strictly_increasing(energies: np.ndarray) -> np.ndarray:
    """Nudge duplicate absorption-edge energies so interpolators stay monotonic."""
    out = np.array(energies, dtype=float, copy=True)
    for i in range(1, len(out)):
        if out[i] <= out[i - 1]:
            out[i] = np.nextafter(out[i - 1], np.inf)
    return out


def _interpolate(material: Elements | Compounds, energy_kev: float | np.ndarray, field: str) -> float | np.ndarray:
    table = _table(material)
    energies_kev = np.asarray(energy_kev, dtype=float)
    scalar = energies_kev.ndim == 0
    energies_mev = np.atleast_1d(energies_kev) * _KEV_TO_MEV

    lo, hi = table["energy_mev"][0], table["energy_mev"][-1]
    if np.any(energies_mev < lo) or np.any(energies_mev > hi):
        raise ValueError(
            f"energy_kev must be within [{lo / _KEV_TO_MEV:g}, {hi / _KEV_TO_MEV:g}] keV "
            f"for {_material_label(material)}"
        )

    # Log-log interpolation is the usual scheme for these coefficients.
    # Absorption edges are tabulated twice at the same energy; adjacent floats
    # still collapse under np.log, so the log-energy grid is made strictly
    # increasing. Querying the exact edge energy returns the pre-edge value.
    x = _strictly_increasing(np.log(table["energy_mev"]))
    y = np.log(table[field])
    values = np.exp(np.interp(np.log(energies_mev), x, y))
    return values.item() if scalar else values


def get_mass_attenuation(material: Elements | Compounds, energy_kev: float | np.ndarray) -> float | np.ndarray:
    """
    Mass attenuation coefficient μ/ρ (cm²/g) at *energy_kev*.

    Values are interpolated in log-log space from the NIST X-ray mass
    coefficient tables.
    """
    return _interpolate(material, energy_kev, _MU_RHO)


def get_energy_absorption(material: Elements | Compounds, energy_kev: float | np.ndarray) -> float | np.ndarray:
    """
    Mass energy-absorption coefficient μ_en/ρ (cm²/g) at *energy_kev*.

    Values are interpolated in log-log space from the NIST X-ray mass
    coefficient tables.
    """
    return _interpolate(material, energy_kev, _MU_EN_RHO)
