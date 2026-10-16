"""Pauli decomposition (Sec. 2.3, Eqs. 6-8) and its structural consequences.

Reproduces: Eqs. (6)-(8) equal the directly built matrices (max deviation 0.00),
H2 = H_mass + H_E is diagonal, the Gauss-law ZZ coefficients
c(j,k) = (N-1-k) g^2 a/4 (e.g. 1.50 and 0.25 at N=8), and the 25 distinct
Pauli terms of the N=6 Hamiltonian.
"""
import numpy as np
from _common import Checks, write_json, run_main
from schwinger_core import parts, pauli_terms, matrix_from_terms, A_DEFAULT
from schwinger_circuit_lib import diagonal_terms


def main(quick=False):
    c = Checks('Pauli decomposition')
    out = {}
    for N in (4, 6, 8):
        for (m, g) in [(0.5, 1.0), (1.3, 0.7)]:
            Hh, Hm, He = parts(N, m, g)
            terms = pauli_terms(N, m, g, drop_identity=False)
            Hp = matrix_from_terms(terms)
            dev = np.max(abs(Hp - (Hh + Hm + He)))
            c.less(f'Pauli form = direct matrix, N={N}, m={m}, g={g}', dev, 1e-12, 'Sec. 2.3: 0.00')
            H2 = Hm + He
            off = np.max(abs(H2 - np.diag(np.diag(H2))))
            c.true(f'H2 diagonal, N={N}, m={m}, g={g}', off == 0.0, f'max off-diagonal {off}')
    N = 8
    z, zz = diagonal_terms(N, 0.5, 1.0)
    ok = all(abs(zz[(j, k)] - (N - 1 - k) * 1.0 * A_DEFAULT / 4) < 1e-12 for (j, k) in zz)
    c.true('ZZ coefficients c(j,k) = (N-1-k) g^2 a/4 for all pairs, N=8', ok)
    c.close('c(0,1) at N=8', zz[(0, 1)], 1.50, 1e-12, 'Sec. 5.3')
    c.close('c(5,6) at N=8', zz[(5, 6)], 0.25, 1e-12, 'Sec. 5.3')
    nterms = len(pauli_terms(6, 0.5, 1.0))
    c.true('N=6 Hamiltonian has 25 distinct non-identity Pauli terms', nterms == 25, f'{nterms} terms')
    out['zz_coefficients_N8'] = {f'{j},{k}': v for (j, k), v in zz.items()}
    out['pauli_terms_N6'] = pauli_terms(6, 0.5, 1.0)
    write_json('pauli_decomposition.json', out)
    return c.finish()


if __name__ == '__main__':
    run_main(main)
