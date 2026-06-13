import pytest

from pyradiate.core import Elements, RadioNuclide
from pyradiate.core.nuclide_id import NuclideIdentifierError


@pytest.fixture(params=["Zn65", "Zn_65", "133_ba", "22 Na", "18-f", "1b3a3"])  # Last one weird but allowed
def valid_identifier(request):
    return request.param


@pytest.fixture(params=["Zinc65", "natrium22", "65", "aN22"])
def invalid_identifier(request):
    return request.param


@pytest.fixture(
    params=[
        "3H",
        "11C",
        "14C",
        "18F",
        "22Na",
        "32P",
        "35S",
        "51Cr",
        "54Mn",
        "57Co",
        "60Co",
        "65Zn",
        "68Ga",
        "75Se",
        "85Kr",
        "89Sr",
        "90Sr",
        "90Y",
        "99Tc",
        "99Mo",
        "106Ru",
        "111In",
        "123I",
        "125I",
        "129I",
        "131I",
        "134Cs",
        "137Cs",
        "147Pm",
        "153Sm",
        "177Lu",
        "192Ir",
        "201Tl",
        "223Ra",
        "225Ac",
        "226Ra",
        "232Th",
        "235U",
        "238U",
        "238Pu",
        "239Pu",
        "241Am",
    ]
)
def important_radionuclide(request):
    return request.param


def test_radio_nuclide_creation_from_identifier(valid_identifier):
    rn = RadioNuclide.from_string(valid_identifier)
    assert rn.element in (Elements["Zn"], Elements["Na"], Elements["F"], Elements["Ba"])
    assert rn.half_life == rn.decays[0].half_life_s
    assert sum(branch.fraction for branch in rn.decays[0].branches) == 100.0


def test_radio_nuclide_creation_from_invalid_identifier_fails(invalid_identifier):
    with pytest.raises(NuclideIdentifierError):
        _ = RadioNuclide.from_string(invalid_identifier)


def test_important_radio_nuclides_can_be_created(important_radionuclide):
    _ = RadioNuclide.from_string(important_radionuclide)


def test_radio_nuclide_cashes_instance_per_nuclid_id(valid_identifier):
    rn_one = RadioNuclide.from_string(valid_identifier)
    rn_two = RadioNuclide.from_string(valid_identifier)
    assert rn_one is rn_two
