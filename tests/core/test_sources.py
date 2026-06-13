import pytest

from pyradiate.core import RadioNuclide, Source


@pytest.fixture
def radio_nuclides_with_activities():
    return {"133Ba": 1e5, "65Zn": 1e4, "60Co": 1e3}


def test_source_contains_nuclide(radio_nuclides_with_activities):
    source = Source(
        radio_nuclides=[RadioNuclide.from_string(x) for x in radio_nuclides_with_activities],
        activities=list(radio_nuclides_with_activities.values()),
    )

    assert all(RadioNuclide.from_string(x) in source for x in radio_nuclides_with_activities)
    assert source.activity == sum(radio_nuclides_with_activities.values())
