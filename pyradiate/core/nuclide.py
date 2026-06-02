from dataclasses import dataclass

from pyradiate.core.elements import Elements


class NuclideIdentifierError(Exception):
    pass


@dataclass(frozen=True)
class NuclideIdentifier:
    mass_number: int
    element: Elements

    @classmethod
    def from_string(cls, radionuclide_string):
        # Extract exact element symbol from very permissive identifier string
        symbol = [c for c in radionuclide_string if c.isalpha()]
        symbol = "".join([c.upper() if i == 0 else c.lower() for i, c in enumerate(symbol)])
        mass_number = int("".join([c for c in radionuclide_string if c.isdigit()]))
        print(radionuclide_string, "#" * 55)
        error = ""
        if not 0 < len(symbol) < 3:
            error = f"radionuclid_string {radionuclide_string} results in invalid element symbol {symbol}"
        elif not hasattr(Elements, symbol):
            error = f"Symbol {symbol} not an element"
        if error:
            raise NuclideIdentifierError(error)
        return cls(mass_number=mass_number, element=Elements[symbol])

    @property
    def identifier(self):
        return f"{self.mass_number}{self.element.symbol}"

    def __repr__(self):
        return f"{self.__class__.__name__}[{self.identifier}]"
