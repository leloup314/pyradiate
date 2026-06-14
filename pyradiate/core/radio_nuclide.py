from pyradiate.core import radiation
from pyradiate.core.elements import Elements
from pyradiate.core.nuclide_id import NuclideIdentifier
from pyradiate.physics import decay_constant
from pyradiate.tools.ensdf.decay_parser import parse_radio_nuclide


class RadioNuclide:
    """
    Definition of a radioactive isotope of an element
    """

    @classmethod
    def from_string(cls, nuclid_string: str):
        nuclid_id = NuclideIdentifier.from_string(nuclid_string)
        return cls.from_identifier(nuclid_identifier=nuclid_id)

    @classmethod
    def from_identifier(cls, nuclid_identifier: NuclideIdentifier):
        return cls(mass_number=nuclid_identifier.mass_number, element=nuclid_identifier.element)

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
    def decays(self) -> list[radiation.Decay]:
        return self._decays

    @property
    def decay_constant(self) -> float:
        return decay_constant(self.half_life)

    # Store instances of already created radio nuclides, i.e. multiton pattern
    _instances = {}

    def __new__(cls, mass_number: int, element: Elements):
        nuclid_id = NuclideIdentifier(mass_number=mass_number, element=element)
        if nuclid_id not in cls._instances:
            cls._instances[nuclid_id] = super().__new__(cls)
        return cls._instances[nuclid_id]

    def __init__(self, mass_number: int, element: Elements):
        self._nuclide_identifier = NuclideIdentifier(mass_number=mass_number, element=element)
        self._parsed_rn = parse_radio_nuclide(self._nuclide_identifier)
        self._decays = list(self._parsed_rn.decays)
        self._recommended_half_life_s = self._parsed_rn.recommended_half_life_s

    def __iter__(self):
        yield from self.decays

    def __repr__(self):
        return f"{self.__class__.__name__}[{self.identifier.identifier}, decay_modes={len(self.decays)})"

    # Convenience access to radiation
    def _radiation(self, kind: str, sort_by: str = "I"):
        rad = [r for d in self.decays for r in getattr(d, kind)]
        rad.sort(key=lambda x: x.intensity if sort_by == "I" else x.energy_kev, reverse=True)
        return rad

    def alphas(self) -> list[radiation.Alpha]:
        return self._radiation(kind="alphas", sort_by="E")

    def betas(self, sort_by: str = "I") -> list[radiation.Beta]:
        return self._radiation(kind="betas", sort_by=sort_by)

    def gammas(self, sort_by: str = "I") -> list[radiation.Gamma]:
        return self._radiation(kind="gammas", sort_by=sort_by)
