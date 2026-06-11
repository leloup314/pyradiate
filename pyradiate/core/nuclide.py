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
        """
        Generate *NuclideIdentifier* from string *radionuclid_string*.
        The string must be a concatenation of <= 2 (adjacent) characters,
        denoting the element (e.g. 'Zn') and 1<= digits <= 3 with their
        integer value <= 300, providing the mass number (e.g. '65').
        The string is insensitive to the order, the letter casing and
        potential separators.

        Examples: '65Zn', 'zn65', '65_zn', 'zn 65'
        """
        # Extract exact element symbol from very permissive identifier string
        symbol = [c for c in radionuclide_string if c.isalpha()]
        symbol = "".join([c.upper() if i == 0 else c.lower() for i, c in enumerate(symbol)])
        mass_number = int("".join([c for c in radionuclide_string if c.isdigit()]))

        # Construct error strings
        symbol_error = mass_number_error = ""
        error = f"Input '{radionuclide_string}' contains invalid {{}} '{{}}':\n"
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
