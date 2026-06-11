from pyradiate.tools.ensdf_decay_parser import parse_radio_nuclide
from pyradiate.core.elements import Elements
from pyradiate.core.nuclide import NuclideIdentifier
from pyradiate.core.radiation import Decay
from pyradiate.physics import decay_constant


class RadioNuclide:
    """
    Definition of a radioactive isotope of an element
    """

    @classmethod
    def from_string(cls, nuclid_string: str, **kwargs):
        nuclid_id = NuclideIdentifier.from_string(nuclid_string)
        return cls.from_identifier(nuclid_identifier=nuclid_id)

    @classmethod
    def from_identifier(cls, nuclid_identifier: NuclideIdentifier, **kwargs):
        return cls(mass_number=nuclid_identifier.mass_number, element=nuclid_identifier.element, **kwargs)

    @property
    def identifier(self) -> NuclideIdentifier:
        return self._nuclide_identifier

    @property
    def atomic_number(self) -> int:
        return self.identifier.element.atomic_number

    @property
    def mass_number(self) -> int:
        return self.identifier.mass_number

    @property
    def element(self) -> Elements:
        return self.identifier.element

    @property
    def half_life(self) -> float:
        return self._recommended_half_life_s

    @property
    def decays(self) -> list[Decay]:
        return self._decays

    @property
    def decay_constant(self) -> float:
        return decay_constant(self.half_life)

    def __init__(self, mass_number: int, element: Elements, activity: float | None = None):
        self._nuclide_identifier = NuclideIdentifier(mass_number=mass_number, element=element)
        parsed = parse_radio_nuclide(self._nuclide_identifier)
        self._decays = list(parsed.decays)
        self._recommended_half_life_s = parsed.recommended_half_life_s
        self.activity = activity

    def __iter__(self):
        yield from self.decays

    def __repr__(self):
        return f"{self.__class__.__name__}[{self.identifier.identifier}, decay_modes={len(self.decays)})"


class Source:
    @property
    def nuclides(self):
        return self._radio_nuclides

    @property
    def activity(self):
        return sum(n.activity for n in self.nuclides)

    def __init__(self, radio_nuclides: list[RadioNuclide]):
        pass


class CalibratedSource(Source):
    """
    Radioactive source of calibrated activity at given time
    """


class Sample(Source):
    pass
