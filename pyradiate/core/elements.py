from enum import Enum
from dataclasses import dataclass


@dataclass(frozen=True)
class Element:
    name: str
    symbol: str
    atomic_number: int


class Elements(Enum):

    @staticmethod
    def from_atomic_number(atomic_number: int):
        try:
            return _ELEMENT_BY_ATOMIC_NUMBER[atomic_number]
        except KeyError:
            raise ValueError(f"No element with atomic number {atomic_number}")

    @property
    def symbol(self):
        return self.value.symbol

    @property
    def atomic_number(self):
        return self.value.atomic_number

    @property
    def name(self):
        return self.value.name

    # List of all elements
    H = Element("Hydrogen", "H", 1)
    He = Element("Helium", "He", 2)
    Li = Element("Lithium", "Li", 3)
    Be = Element("Beryllium", "Be", 4)
    B = Element("Boron", "B", 5)
    C = Element("Carbon", "C", 6)
    N = Element("Nitrogen", "N", 7)
    O = Element("Oxygen", "O", 8)
    F = Element("Fluorine", "F", 9)
    Ne = Element("Neon", "Ne", 10)
    Na = Element("Sodium", "Na", 11)
    Mg = Element("Magnesium", "Mg", 12)
    Al = Element("Aluminium", "Al", 13)
    Si = Element("Silicon", "Si", 14)
    P = Element("Phosphorus", "P", 15)
    S = Element("Sulfur", "S", 16)
    Cl = Element("Chlorine", "Cl", 17)
    Ar = Element("Argon", "Ar", 18)
    K = Element("Potassium", "K", 19)
    Ca = Element("Calcium", "Ca", 20)
    Sc = Element("Scandium", "Sc", 21)
    Ti = Element("Titanium", "Ti", 22)
    V = Element("Vanadium", "V", 23)
    Cr = Element("Chromium", "Cr", 24)
    Mn = Element("Manganese", "Mn", 25)
    Fe = Element("Iron", "Fe", 26)
    Co = Element("Cobalt", "Co", 27)
    Ni = Element("Nickel", "Ni", 28)
    Cu = Element("Copper", "Cu", 29)
    Zn = Element("Zinc", "Zn", 30)
    Ga = Element("Gallium", "Ga", 31)
    Ge = Element("Germanium", "Ge", 32)
    As = Element("Arsenic", "As", 33)
    Se = Element("Selenium", "Se", 34)
    Br = Element("Bromine", "Br", 35)
    Kr = Element("Krypton", "Kr", 36)
    Rb = Element("Rubidium", "Rb", 37)
    Sr = Element("Strontium", "Sr", 38)
    Y = Element("Yttrium", "Y", 39)
    Zr = Element("Zirconium", "Zr", 40)
    Nb = Element("Niobium", "Nb", 41)
    Mo = Element("Molybdenum", "Mo", 42)
    Tc = Element("Technetium", "Tc", 43)
    Ru = Element("Ruthenium", "Ru", 44)
    Rh = Element("Rhodium", "Rh", 45)
    Pd = Element("Palladium", "Pd", 46)
    Ag = Element("Silver", "Ag", 47)
    Cd = Element("Cadmium", "Cd", 48)
    In = Element("Indium", "In", 49)
    Sn = Element("Tin", "Sn", 50)
    Sb = Element("Antimony", "Sb", 51)
    Te = Element("Tellurium", "Te", 52)
    I = Element("Iodine", "I", 53)
    Xe = Element("Xenon", "Xe", 54)
    Cs = Element("Caesium", "Cs", 55)
    Ba = Element("Barium", "Ba", 56)
    La = Element("Lanthanum", "La", 57)
    Ce = Element("Cerium", "Ce", 58)
    Pr = Element("Praseodymium", "Pr", 59)
    Nd = Element("Neodymium", "Nd", 60)
    Pm = Element("Promethium", "Pm", 61)
    Sm = Element("Samarium", "Sm", 62)
    Eu = Element("Europium", "Eu", 63)
    Gd = Element("Gadolinium", "Gd", 64)
    Tb = Element("Terbium", "Tb", 65)
    Dy = Element("Dysprosium", "Dy", 66)
    Ho = Element("Holmium", "Ho", 67)
    Er = Element("Erbium", "Er", 68)
    Tm = Element("Thulium", "Tm", 69)
    Yb = Element("Ytterbium", "Yb", 70)
    Lu = Element("Lutetium", "Lu", 71)
    Hf = Element("Hafnium", "Hf", 72)
    Ta = Element("Tantalum", "Ta", 73)
    W = Element("Tungsten", "W", 74)
    Re = Element("Rhenium", "Re", 75)
    Os = Element("Osmium", "Os", 76)
    Ir = Element("Iridium", "Ir", 77)
    Pt = Element("Platinum", "Pt", 78)
    Au = Element("Gold", "Au", 79)
    Hg = Element("Mercury", "Hg", 80)
    Tl = Element("Thallium", "Tl", 81)
    Pb = Element("Lead", "Pb", 82)
    Bi = Element("Bismuth", "Bi", 83)
    Po = Element("Polonium", "Po", 84)
    At = Element("Astatine", "At", 85)
    Rn = Element("Radon", "Rn", 86)
    Fr = Element("Francium", "Fr", 87)
    Ra = Element("Radium", "Ra", 88)
    Ac = Element("Actinium", "Ac", 89)
    Th = Element("Thorium", "Th", 90)
    Pa = Element("Protactinium", "Pa", 91)
    U = Element("Uranium", "U", 92)
    Np = Element("Neptunium", "Np", 93)
    Pu = Element("Plutonium", "Pu", 94)
    Am = Element("Americium", "Am", 95)
    Cm = Element("Curium", "Cm", 96)
    Bk = Element("Berkelium", "Bk", 97)
    Cf = Element("Californium", "Cf", 98)
    Es = Element("Einsteinium", "Es", 99)
    Fm = Element("Fermium", "Fm", 100)
    Md = Element("Mendelevium", "Md", 101)
    No = Element("Nobelium", "No", 102)
    Lr = Element("Lawrencium", "Lr", 103)
    Rf = Element("Rutherfordium", "Rf", 104)
    Db = Element("Dubnium", "Db", 105)
    Sg = Element("Seaborgium", "Sg", 106)
    Bh = Element("Bohrium", "Bh", 107)
    Hs = Element("Hassium", "Hs", 108)
    Mt = Element("Meitnerium", "Mt", 109)
    Ds = Element("Darmstadtium", "Ds", 110)
    Rg = Element("Roentgenium", "Rg", 111)
    Cn = Element("Copernicium", "Cn", 112)
    Nh = Element("Nihonium", "Nh", 113)
    Fl = Element("Flerovium", "Fl", 114)
    Mc = Element("Moscovium", "Mc", 115)
    Lv = Element("Livermorium", "Lv", 116)
    Ts = Element("Tennessine", "Ts", 117)
    Og = Element("Oganesson", "Og", 118)


# Helper map to atomic number
_ELEMENT_BY_ATOMIC_NUMBER = {
    elem.atomic_number: elem for elem in Elements
}
