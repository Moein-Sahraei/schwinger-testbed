"""Beyond exact diagonalization: matrix-product states (Sec. 8.4, Outlook 1).

The Hamiltonian of Eq. (2) is written as a matrix-product operator from the
same Pauli decomposition used for the circuit (Eqs. 6-8), so the
Gauss-law-induced all-to-all ZZ couplings c(j,k) = (N-1-k) g^2 a/4 are carried
exactly by the MPO.  The ground state is found with two-site DMRG (TeNPy).

Reproduces: the exact-diagonalization ground-state energy for N <= 12 (the
preprint quotes agreement to 1e-14 or better; differences at this level are
floating-point noise, so the check requires < 1e-13), and the extension to
N = 30 at bond dimension chi = 96.  The energy variance <H^2> - <H>^2 is
reported as an independent convergence measure.

Parameters: m = 0.5, g = 1, a = 1, open boundaries.  Only ground states are
computed here; time evolution with this long-range MPO is the open direction
named in the preprint's outlook.
"""
import io, time, logging, warnings, contextlib
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from _common import Checks, write_csv, write_json, run_main
from schwinger_core import pauli_terms, hamiltonian

warnings.filterwarnings('ignore')
logging.getLogger('tenpy').setLevel(logging.ERROR)
from tenpy.networks.site import SpinHalfSite          # noqa: E402
from tenpy.models.model import CouplingMPOModel       # noqa: E402
from tenpy.networks.mps import MPS                    # noqa: E402
from tenpy.algorithms import dmrg                     # noqa: E402

M0, G0 = 0.5, 1.0
TENPY_OP = {'X': 'Sigmax', 'Y': 'Sigmay', 'Z': 'Sigmaz'}   # TeNPy's Pauli matrices; 'up' = Z eigenvalue +1


def split_terms(N, m=M0, g=G0):
    """Pauli terms of Eqs. (6)-(8) without the identity, plus the identity constant."""
    terms = pauli_terms(N, m, g, drop_identity=False)
    const = terms.pop('I' * N, 0.0)
    return terms, const


class SchwingerMPO(CouplingMPOModel):
    """Eq. (2) as a spin-1/2 chain.  Labels are little-endian (rightmost character = site 0)."""

    def init_sites(self, params):
        return SpinHalfSite(conserve=None)

    def init_terms(self, params):
        N = params.get('L', None)
        for label, coef in params.get('terms', None).items():
            ops = sorted((N - 1 - i, ch) for i, ch in enumerate(label) if ch != 'I')
            if len(ops) == 1:
                self.add_onsite_term(coef, ops[0][0], TENPY_OP[ops[0][1]])
            elif len(ops) == 2:
                self.add_coupling_term(coef, ops[0][0], ops[1][0], TENPY_OP[ops[0][1]], TENPY_OP[ops[1][1]])
            else:
                raise ValueError(f'unexpected Pauli string {label}')


def dmrg_ground_state(N, chi, m=M0, g=G0):
    terms, const = split_terms(N, m, g)
    model = SchwingerMPO({'L': N, 'terms': terms, 'bc_MPS': 'finite'})
    vacuum = ['down' if n % 2 == 0 else 'up' for n in range(N)]   # |1010...>, site 0 in |1>
    psi = MPS.from_product_state(model.lat.mps_sites(), vacuum, bc='finite')
    with contextlib.redirect_stdout(io.StringIO()):
        info = dmrg.run(psi, model, {'trunc_params': {'chi_max': chi, 'svd_min': 1e-12},
                                     'max_sweeps': 40, 'mixer': True})
    sz = psi.expectation_value('Sigmaz')
    sigma = sum((-1) ** n * sz[n] for n in range(N)) / (2 * N)
    return {'E0': float(np.real(info['E'])) + const, 'Sigma': float(np.real(sigma)), 'chi': int(max(psi.chi)),
            'variance': float(np.real(model.H_MPO.variance(psi))), 'mpo_chi': int(max(model.H_MPO.chi))}


PAULI_SPARSE = {'I': sp.identity(2, format='csr', dtype=complex),
                'X': sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex)),
                'Y': sp.csr_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex)),
                'Z': sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))}


def sparse_hamiltonian(N, m=M0, g=G0):
    """Eq. (2) as a sparse matrix from its Pauli decomposition, identity included."""
    H = sp.csr_matrix((2 ** N, 2 ** N), dtype=complex)
    for label, coef in pauli_terms(N, m, g, drop_identity=False).items():
        M = sp.identity(1, format='csr', dtype=complex)
        for ch in label:                       # leftmost character = highest qubit
            M = sp.kron(M, PAULI_SPARSE[ch], format='csr')
        H = H + coef * M
    return H


def exact_ground_energy(N):
    return eigsh(sparse_hamiltonian(N), k=1, which='SA', return_eigenvectors=False, tol=0)[0].real


def main(quick=False):
    c = Checks('Matrix-product states beyond exact diagonalization (Sec. 8.4)')
    dev = np.max(abs(sparse_hamiltonian(6).toarray() - hamiltonian(6, M0, G0)))
    c.less('Sparse Pauli-sum Hamiltonian = directly built H (N=6)', dev, 1e-12, 'max elementwise deviation')

    val = []
    for N in ((4, 6, 8) if quick else (4, 6, 8, 10, 12)):
        E_ed = exact_ground_energy(N)
        r = dmrg_ground_state(N, chi=64)
        diff = abs(r['E0'] - E_ed)
        val.append((N, E_ed, r['E0'], diff, r['chi']))
        c.less(f'DMRG ground-state energy = exact diagonalization, N={N}', diff, 1e-13,
               'Sec. 8.4: 1e-14 or better')
    write_csv('mps_validation.csv', ['N', 'E0_exact', 'E0_dmrg', 'abs_difference', 'bond_dimension'], val)

    big = []
    if not quick:
        for N in (16, 20, 24, 30):
            t0 = time.time()
            r = dmrg_ground_state(N, chi=96)
            big.append((N, r['E0'], r['E0'] / N, r['Sigma'], r['chi'], r['variance'], r['mpo_chi'],
                        round(time.time() - t0, 1)))
        write_csv('mps_ground_state.csv', ['N', 'E0', 'E0_per_site', 'Sigma', 'bond_dimension',
                                           'energy_variance', 'mpo_bond_dimension', 'seconds'], big)
        last = big[-1]
        c.true('N = 30 reached at bond dimension chi = 96', last[0] == 30 and last[4] == 96,
               f'max bond dimension {last[4]}')
        c.less('Energy variance at N = 30 (converged eigenstate)', last[5], 1e-9, 'convergence measure')
        c.true('Energy per site decreases monotonically with N', all(b[2] < a[2] for a, b in zip(big, big[1:])),
               ', '.join(f'N={b[0]}: {b[2]:.6f}' for b in big))
    write_json('mps.json', {'validation': [{'N': v[0], 'abs_difference': v[3]} for v in val],
                            'beyond_ed': [{'N': b[0], 'E0': b[1], 'Sigma': b[3], 'chi': b[4],
                                           'variance': b[5]} for b in big]})
    return c.finish()


if __name__ == '__main__':
    run_main(main)
