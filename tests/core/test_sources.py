import pytest

from pyradiate.core.elements import Elements
from pyradiate.core.sources import RadioNuclide
from pyradiate.core.errors import NuclidIdentifierError


@pytest.fixture(params=["Zn65", "Zn_65", "65_zn", "65zN", "zn-65"])
def valid_identifier(request):
    return request.param


@pytest.fixture(params=["Zinc65", "natrium22", "65", "22"])
def invalid_identifier(request):
    return request.param


def test_radio_nuclide_from_valid_identifier(valid_identifier):
    rn = RadioNuclide.from_identifier(valid_identifier)
    assert rn.element is Elements["Zn"]
    assert rn.atomic_number == 30
    assert rn.activity is None


def test_radio_nuclide_from_invalid_identifier_fails(invalid_identifier):
    with pytest.raises(NuclidIdentifierError):
        _ = RadioNuclide.from_identifier(invalid_identifier)
