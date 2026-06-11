from collections import Counter

from pyradiate import ensdf_path
from pyradiate.core.radiation import DecayBranch
from pyradiate.tools import ensdf_decay_parser


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
        (parent, decay.mode, decay.daughter_nuclide.identifier, round(decay.half_life_s, 6))
        for parent, parsed in data.items()
        for decay in parsed.decays
    ]
    counts = Counter(keys)
    duplicates = [key for key, count in counts.items() if count > 1]
    assert duplicates == []


def test_ensdf_parser_uses_adopted_data_for_zn65():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    zn65 = data["65Zn"]
    assert len(zn65.decays) == 1
    assert zn65.recommended_half_life_s == 243.93 * 86400
    decay = zn65.decays[0]
    assert decay.half_life_s == zn65.recommended_half_life_s
    assert decay.branch_fraction == 100.0
    assert decay.branch is DecayBranch.EC_B_PLUS


def test_ensdf_parser_branch_fractions_for_multi_decay_nuclide():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    carbon9 = data["9C"]
    by_branch = {d.branch: d for d in carbon9.decays}
    assert by_branch[DecayBranch.B_PLUS_ALPHA].branch_fraction == 37.9
    assert by_branch[DecayBranch.B_PLUS_PROTON].branch_fraction == 62.0


def test_ensdf_parser_isomer_branch_fractions_for_co60():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    co60 = data["60Co"]
    isomer = next(d for d in co60.decays if d.branch is DecayBranch.IT)
    isomer_beta = next(d for d in co60.decays if d.mode == "B- DECAY" and d.half_life_s < 1000)
    assert co60.recommended_half_life_s > 1e8
    assert isomer.branch_fraction == 99.75
    assert isomer_beta.branch_fraction == 0.25
