import pytest

from pyradiate.core.elements import Elements
from pyradiate.core.sources import RadioNuclide


@pytest.fixture(params=["Zn65", "Zn_65", "65_Zn", "65Zn", "Zn-?/65"])
def identifier(request):
    return request.param


def test_radio_nuclide_from_identifier(identifier):
    rn = RadioNuclide.from_identifier(identifier)
    assert rn.element is Elements["Zn"]
    assert rn.atomic_number == 30
    assert rn.activity is None
