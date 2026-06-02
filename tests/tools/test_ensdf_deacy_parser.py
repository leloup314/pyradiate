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

    assert n_isos == 2561
    assert n_gammas == 116249
    assert n_alphas == 2760
    assert n_betas == 440
    assert n_xrays == 546
