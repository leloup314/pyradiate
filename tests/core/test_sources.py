import pytest

from pyradiate.core.elements import Elements
from pyradiate.core.sources import RadioNuclide
from pyradiate.core.nuclide import NuclideIdentifierError


@pytest.fixture(params=["Zn65", "Zn_65", "65_zn", "65zN", "zn-65"])
def valid_identifier(request):
    return request.param


@pytest.fixture(params=["Zinc65", "natrium22", "65", "22"])
def invalid_identifier(request):
    return request.param


def test_radio_nuclide_from_valid_identifier(valid_identifier):
    rn = RadioNuclide.from_string(valid_identifier)
    assert rn.element is Elements["Zn"]
    assert rn.atomic_number == 30
    assert rn.activity is None


def test_radio_nuclide_from_invalid_identifier_fails(invalid_identifier):
    with pytest.raises(NuclideIdentifierError):
        _ = RadioNuclide.from_string(invalid_identifier)
