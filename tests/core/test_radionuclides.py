import pytest

from pyradiate.core.source import RadioNuclide


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


def test_important_radio_nuclides_can_be_created(important_radionuclide):
    _ = RadioNuclide.from_string(important_radionuclide)
