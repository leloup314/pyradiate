from collections import Counter

from pyradiate import ensdf_path
from pyradiate.tools import ensdf_decay_parser


def test_ensdf_parser():
    n_isos = 0
    n_gammas = 0
    n_xrays = 0
    n_alphas = 0
    n_betas = 0
    for iso, data in ensdf_decay_parser.parse_ensdf_directory(ensdf_path).items():
        n_isos += 1
        for mode in data:
            n_gammas += len(mode.gammas)
            n_betas += len(mode.betas)
            n_alphas += len(mode.alphas)
            n_xrays += len(mode.xrays)

    assert n_isos == 2570
    assert n_gammas == 114642
    assert n_alphas == 2745
    assert n_betas == 439
    assert n_xrays == 546


def test_ensdf_parser_deduplicates_same_decay():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    keys = [
        (parent, decay.mode, decay.daughter_nuclide.identifier, round(decay.half_life_s, 6))
        for parent, decays in data.items()
        for decay in decays
    ]
    counts = Counter(keys)
    duplicates = [key for key, count in counts.items() if count > 1]
    assert duplicates == []


def test_ensdf_parser_uses_adopted_data_for_zn65():
    data = ensdf_decay_parser.parse_ensdf_directory(ensdf_path)
    zn65_decays = data["65Zn"]
    assert len(zn65_decays) == 1
    decay = zn65_decays[0]
    assert decay.half_life_s == 243.93 * 86400
    assert decay.branches
