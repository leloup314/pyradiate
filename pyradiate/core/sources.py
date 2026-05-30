from pyradiate.core.elements import Elements
from pyradiate.core.errors import NuclidIdentifierError
from pyradiate.physics import decay_constant


class RadioNuclide:
    """
    Definition of a radioactive isotope of an element
    """

    @classmethod
    def from_identifier(cls, identifier: str, **kwargs):
        # Extract exact element symbol from very permissive identifier string
        symbol = [c for c in identifier if c.isalpha()]
        symbol = "".join([c.upper() if i == 0 else c.lower() for i, c in enumerate(symbol)])
        mass_number = int("".join([c for c in identifier if c.isdigit()]))

        error_str = ""
        if not 0 < len(symbol) < 3:
            error_str = f"Identifier {identifier} results in invalid element symbol {symbol}"
        elif not hasattr(Elements, symbol):
            error_str = f"Symbol {symbol} not an element"
        if error_str:
            raise NuclidIdentifierError(error_str)

        return cls(Elements[symbol].atomic_number, mass_number, **kwargs)

    @property
    def atomic_number(self) -> int:
        return self._atomic_number

    @property
    def mass_number(self) -> int:
        return self._mass_number

    @property
    def element(self) -> Elements:
        return self._element

    @property
    def half_life(self) -> float:
        return None

    @property
    def decay_constant(self) -> float:
        return decay_constant(self.half_life)

    def __init__(self, atomic_number: int, mass_number: int, activity: float | None = None):
        self._mass_number = mass_number
        self._atomic_number = atomic_number
        self._element = Elements.from_atomic_number(atomic_number)

        # Accesible attributes
        self.activity = activity


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
