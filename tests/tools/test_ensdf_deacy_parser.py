from collections import Counter

import pytest

from pyradiate import ensdf_path
from pyradiate.core.radiation import BranchIdentifier
from pyradiate.tools import ensdf_decay_parser


def _branch_sum(decay) -> float:
    return sum(branch.fraction for branch in decay.branches)


def test_ensdf_parser():
    n_isos = 0
    n_gammas = 0
    n_xrays = 0
    n_alphas = 0
    n_betas = 0
    for iso, parsed in ensdf_decay_parser.parse_ensdf_directory(ensdf_path).items():
        n_isos += 1
        for mode in parsed.decays:
            n_gammas += len(mode.gammas)
            n_betas += len(mode.betas)
            n_alphas += len(mode.alphas)
            n_xrays += len(mode.xrays)

    assert n_isos == 2570
    assert n_gammas == 114641
    assert n_alphas == 2736
    assert n_betas == 439
    assert n_xrays == 546


def test_ensdf_parser_deduplicates_same_decay():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    keys = [
        (parent, round(decay.half_life_s, 6))
        for parent, parsed in data.items()
        for decay in parsed.decays
    ]
    counts = Counter(keys)
    duplicates = [key for key, count in counts.items() if count > 1]
    assert duplicates == []


def test_ensdf_parser_branch_fractions_sum_to_one_hundred():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    for parsed in data.values():
        for decay in parsed.decays:
            assert _branch_sum(decay) == pytest.approx(100.0, rel=0.02, abs=0.5)


def test_ensdf_parser_uses_adopted_data_for_zn65():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    zn65 = data["65Zn"]
    assert len(zn65.decays) == 1
    assert zn65.recommended_half_life_s == 243.93 * 86400
    decay = zn65.decays[0]
    assert decay.half_life_s == zn65.recommended_half_life_s
    assert len(decay.branches) == 1
    assert decay.branches[0].identifier is BranchIdentifier.EC_B_PLUS
    assert decay.branches[0].fraction == 100.0
    assert decay.branches[0].daughter_nuclide.identifier == "65Cu"
    assert decay.mode == "EC+B+ DECAY"


def test_ensdf_parser_parses_composite_ec_b_plus_branch():
    branches = ensdf_decay_parser._extract_branches_from_text("%EC+%B+=100")
    assert branches == [(BranchIdentifier.EC_B_PLUS.value, 100.0)]

    branches = ensdf_decay_parser._extract_branches_from_text("%EC=50 %B+=50")
    assert branches == [(BranchIdentifier.EC.value, 50.0), (BranchIdentifier.B_PLUS.value, 50.0)]


def test_ensdf_parser_branch_fractions_for_multi_decay_nuclide():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    carbon9 = data["9C"]
    assert len(carbon9.decays) == 1
    by_branch = {branch.identifier: branch.fraction for branch in carbon9.decays[0].branches}
    assert by_branch[BranchIdentifier.B_PLUS_ALPHA] == 37.9
    assert by_branch[BranchIdentifier.B_PLUS_PROTON] == 62.0


def test_ensdf_parser_isomer_branch_fractions_for_co60():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    co60 = data["60Co"]
    isomer = next(d for d in co60.decays if d.half_life_s < 1000)
    by_branch = {branch.identifier: branch.fraction for branch in isomer.branches}
    assert co60.recommended_half_life_s > 1e8
    assert by_branch[BranchIdentifier.IT] == 99.75
    assert by_branch[BranchIdentifier.B_MINUS] == 0.25


def test_ensdf_parser_ba133_groups_isomer_branches():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    ba133 = data["133Ba"]
    assert len(ba133.decays) == 2

    ground = next(d for d in ba133.decays if d.half_life_s > 1e8)
    isomer = next(d for d in ba133.decays if d.half_life_s < 1e6)

    assert len(ground.branches) == 1
    assert ground.branches[0].identifier is BranchIdentifier.EC
    assert ground.branches[0].fraction == 100.0

    isomer_branches = {branch.identifier: branch.fraction for branch in isomer.branches}
    assert isomer_branches[BranchIdentifier.IT] == pytest.approx(99.9896, rel=1e-3)
    assert isomer_branches[BranchIdentifier.EC] == pytest.approx(0.0104, rel=1e-2)
    assert _branch_sum(isomer) == pytest.approx(100.0, rel=1e-3)

    assert ground.branches[0].daughter_nuclide.identifier == "133Cs"
    ec_branch = next(b for b in isomer.branches if b.identifier is BranchIdentifier.EC)
    it_branch = next(b for b in isomer.branches if b.identifier is BranchIdentifier.IT)
    assert ec_branch.daughter_nuclide.identifier == "133Cs"
    assert it_branch.daughter_nuclide.identifier == "133Ba"
    assert isomer.mode == "EC+IT DECAY"
