import re
from enum import StrEnum


class BranchIdentifier(StrEnum):
    """ENSDF decay branch identifier (e.g. %B-=, %IT=, %EC+%B+=)."""

    C14 = "14C"
    B_MINUS_2 = "2B-"
    EC_2 = "2EC"
    N_2 = "2N"
    P_2 = "2P"
    ALPHA = "A"
    B_PLUS = "B+"
    B_PLUS_2 = "2B+"
    B_PLUS_2P = "B+2P"
    B_PLUS_3P = "B+3P"
    B_PLUS_ALPHA = "B+A"
    B_PLUS_PROTON = "B+P"
    B_MINUS = "B-"
    B_MINUS_2N = "B-2N"
    B_MINUS_ALPHA = "B-A"
    B_MINUS_N = "B-N"
    B_MINUS_PROTON = "B-P"
    EC = "EC"
    EC_B_PLUS = "EC+B+"
    EC_2P = "EC2P"
    EC_3P = "EC3P"
    EC_ALPHA = "ECA"
    EC_PROTON = "ECP"
    IT = "IT"
    NEUTRON = "N"
    PROTON = "P"
    SF = "SF"

    @classmethod
    def from_string(cls, ensdf_decay_string: str):
        normalized = re.sub(r"\s+DECAY.*", "", ensdf_decay_string.strip(), flags=re.IGNORECASE).strip().upper()
        normalized = re.sub(r"\+%", "+", normalized)
        for branch in cls:
            if branch.value == normalized:
                return branch
        raise ValueError(f"Unknown decay branch identifier: {ensdf_decay_string!r}")
