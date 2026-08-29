import numpy as np
import pytest

from pyradiate.core.elements import Elements
from pyradiate.tools.nist import Compounds, get_energy_absorption, get_mass_attenuation, load_table
from pyradiate.tools.nist.coefficients import TABLE_DTYPE
from pyradiate.tools.nist.lib_builder import (
    NIST_BASE_URL,
    NIST_COMPOUNDS_INDEX_URL,
    NIST_ELEMENTS_INDEX_URL,
    fetch_url,
    parse_coefficient_table,
    parse_compound_index,
    parse_element_index,
)


@pytest.fixture(scope="session")
def nist_elements_index():
    return fetch_url(NIST_ELEMENTS_INDEX_URL)


@pytest.fixture(scope="session")
def nist_compounds_index():
    return fetch_url(NIST_COMPOUNDS_INDEX_URL)


@pytest.fixture(scope="session")
def nist_hydrogen_page():
    return fetch_url(f"{NIST_BASE_URL}/ElemTab/z01.html")


def test_load_table_covers_nist_materials():
    for atomic_number in range(1, 93):
        table = load_table(Elements.from_atomic_number(atomic_number))
        assert table.dtype == TABLE_DTYPE
        assert table["energy_mev"][0] == pytest.approx(1.0e-3)
        assert table["energy_mev"][-1] == pytest.approx(20.0)

    assert len(Compounds) == 48
    for compound in Compounds:
        assert load_table(compound).dtype == TABLE_DTYPE


def test_load_table_returns_an_independent_copy():
    table = load_table(Compounds.WATER_LIQUID)
    table["mu_rho"][0] = -1.0
    assert load_table(Compounds.WATER_LIQUID)["mu_rho"][0] == pytest.approx(4.078e3)


def test_tabulated_hydrogen_and_water_values():
    assert get_mass_attenuation(Elements.H, 1.0) == pytest.approx(7.217)
    assert get_energy_absorption(Elements.H, 1.0) == pytest.approx(6.820)
    assert get_mass_attenuation(Elements.H, 100.0) == pytest.approx(0.2944)
    assert get_mass_attenuation(Compounds.WATER_LIQUID, 100.0) == pytest.approx(0.1707)
    assert get_energy_absorption(Compounds.WATER_LIQUID, 100.0) == pytest.approx(0.02546)


def test_interpolation_is_log_log_between_neighbors():
    table = load_table(Elements.H)
    e0, mu0 = table["energy_mev"][0], table["mu_rho"][0]
    e1, mu1 = table["energy_mev"][1], table["mu_rho"][1]
    mid_mev = np.sqrt(e0 * e1)
    expected = np.exp(np.interp(np.log(mid_mev), np.log([e0, e1]), np.log([mu0, mu1])))
    assert get_mass_attenuation(Elements.H, mid_mev * 1e3) == pytest.approx(expected)


def test_array_and_scalar_energy_inputs():
    scalar = get_mass_attenuation(Elements.Pb, 100.0)
    vector = get_mass_attenuation(Elements.Pb, np.array([10.0, 100.0, 1000.0]))
    assert isinstance(scalar, float)
    assert isinstance(vector, np.ndarray)
    assert vector.shape == (3,)
    assert vector[1] == pytest.approx(scalar)


def test_lead_k_edge_keeps_pre_edge_value():
    table = load_table(Elements.Pb)
    edge_rows = table[np.isclose(table["energy_mev"], 8.80045e-2)]
    assert len(edge_rows) == 2
    pre_energy_kev = edge_rows["energy_mev"][0] * 1e3
    assert get_mass_attenuation(Elements.Pb, pre_energy_kev) == pytest.approx(edge_rows["mu_rho"][0])


def test_unsupported_element_and_out_of_range_energy():
    with pytest.raises(ValueError, match="No NIST"):
        get_mass_attenuation(Elements.Og, 100.0)
    with pytest.raises(ValueError, match="energy_kev"):
        get_mass_attenuation(Elements.H, 0.5)
    with pytest.raises(TypeError):
        get_mass_attenuation("water", 100.0)  # type: ignore[arg-type]


def test_parse_element_index_from_nist(nist_elements_index):
    elements = parse_element_index(nist_elements_index)
    assert [element.symbol for element, _href in elements] == [
        Elements.from_atomic_number(z).symbol for z in range(1, 93)
    ]
    assert elements[0][1] == "ElemTab/z01.html"
    assert elements[-1][1] == "ElemTab/z92.html"


def test_parse_compound_index_from_nist(nist_compounds_index):
    compounds = parse_compound_index(nist_compounds_index)
    by_name = {name: (display, nist_id) for name, display, nist_id in compounds}
    assert len(compounds) == 48
    assert by_name["WATER_LIQUID"] == ("Water, Liquid", "water")
    assert by_name["TISSUE_SOFT"][1] == "tissue"
    assert by_name["TISSUE_SOFT_ICRU_FOUR_COMPONENT"][1] == "tissue4"
    assert by_name["CERIC_AMMONIUM_SULFATE_SOLUTION"][1] == "ceric"
    assert {name for name, _display, _nist_id in compounds} == {member._name_ for member in Compounds}


def test_parse_coefficient_table_from_nist(nist_hydrogen_page):
    table = parse_coefficient_table(nist_hydrogen_page)
    archived = load_table(Elements.H)
    assert table.dtype == TABLE_DTYPE
    assert table["energy_mev"][0] == pytest.approx(1.0e-3)
    assert table["mu_rho"][0] == pytest.approx(7.217)
    assert table["mu_en_rho"][0] == pytest.approx(6.820)
    np.testing.assert_allclose(table["energy_mev"], archived["energy_mev"])
    np.testing.assert_allclose(table["mu_rho"], archived["mu_rho"])
    np.testing.assert_allclose(table["mu_en_rho"], archived["mu_en_rho"])
