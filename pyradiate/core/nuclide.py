from dataclasses import dataclass

from pyradiate.core.elements import Elements
from pyradiate.core import radiation
from pyradiate.physics import decay_constant


class NuclideIdentifierError(Exception):
    pass


@dataclass(frozen=True)
class NuclideIdentifier:
    mass_number: int
    element: Elements

    @classmethod
    def from_string(cls, nuclide_string):
        """
        Generate *NuclideIdentifier* from string *nuclide_string*.
        The string must be a concatenation of <= 2 (adjacent) characters,
        denoting the element (e.g. 'Zn') and 1<= digits <= 3 with their
        integer value <= 300, providing the mass number (e.g. '65').
        The string is insensitive to the order, the letter casing and
        potential separators.

        Examples: '65Zn', 'zn65', '65_zn', 'zn 65'
        """
        # Extract exact element symbol from very permissive identifier string
        symbol = [c for c in nuclide_string if c.isalpha()]
        symbol = "".join([c.upper() if i == 0 else c.lower() for i, c in enumerate(symbol)])
        mass_number = int("".join([c for c in nuclide_string if c.isdigit()]))

        # Construct error strings
        symbol_error = mass_number_error = ""
        error = f"Input '{nuclide_string}' contains invalid {{}} '{{}}':\n"
        if not 1 <= len(symbol) <= 2:
            symbol_error = f"Symbol must contain at least one, at most two characters (has {len(symbol)}."
        elif not hasattr(Elements, symbol):
            symbol_error = f"Symbol '{symbol}' not an element"
        elif not 1 <= mass_number <= 300:
            mass_number_error = f"Mass number must be 1 <= mass_number <= 300, is {mass_number}."
        if symbol_error:
            raise NuclideIdentifierError(error.format("symbol", symbol) + symbol_error)
        elif mass_number_error:
            raise NuclideIdentifierError(error.format("mass number", mass_number) + mass_number_error)

        return cls(mass_number=mass_number, element=Elements[symbol])

    @property
    def identifier(self):
        return f"{self.mass_number}{self.element.symbol}"

    def __repr__(self):
        return f"{self.__class__.__name__}[{self.identifier}]"

    def __str__(self):
        return self.identifier


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
