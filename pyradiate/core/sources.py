from dataclasses import dataclass
from pyradiate.core.radio_nuclide import RadioNuclide


@dataclass(frozen=True)
class UnitSource:
    """The most basic verison of a radioactive source: a single radioactive nuclide with given activity"""

    nuclide: RadioNuclide
    activity: float

    def __repr__(self) -> str:
        return f"{self.nuclide.identifier}({self.activity} Bq)"


class Source:
    @property
    def nuclides(self) -> list[RadioNuclide]:
        return [us.nuclide for us in self._unit_sources]

    @property
    def activity(self) -> float:
        return sum(us.activity for us in self._unit_sources)

    def __init__(self, radio_nuclides: list[RadioNuclide], activities: list[float]):
        assert len(radio_nuclides) == len(activities), "Lengths of RadioNuclides and activities do not match"
        self._unit_sources = [
            UnitSource(nuclide=radio_nuclides[i], activity=activities[i]) for i in range(len(activities))
        ]

    def __contains__(self, rn: RadioNuclide):
        return rn in self.nuclides

    def __repr__(self):
        return f"Source[{', '.join(f'{sc!r}' for sc in self._unit_sources)}]"


class CalibratedSource(Source):
    """
    Radioactive source of calibrated activity at given time
    """


class Sample(Source):
    pass
