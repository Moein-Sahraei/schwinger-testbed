"""Reference 1: exact diagonalization (Sec. 2, Table 1, Fig. 1).

Reproduces: Hermiticity, charge conservation, free-field limit, strong-mass
limit, chiral condensate at m=20, Gauss-law consistency, ||[H1,H2]||_2,
and the Nguyen et al. convention correspondence (Eq. 4-5).
Run from the repository root:  python code/schwinger_ed.py
"""
import numpy as np
from scipy.linalg import eigvalsh
from _common import Checks, write_json, run_main
from schwinger_core import (parts, hamiltonian, split_H1_H2, charge_op, condensate_op,
                            site_op, Z, comm, norm2, ground_state, A_DEFAULT)

M0, G0 = 0.5, 1.0


def nguyen_hamiltonian(N, x, mu):
    """Eq. (4): Nguyen et al., PRX Quantum 3, 020324 (2022); sites m=1..N."""
    from schwinger_core import two_site_op, X, Y
    dim = 2 ** N
    H = np.zeros((dim, dim), complex)
    for n in range(N - 1):
        H += x * (two_site_op(X, n, X, n + 1, N) + two_site_op(Y, n, Y, n + 1, N))
    S = np.zeros((dim, dim), complex)
    for n in range(N - 1):
        S = S + site_op(Z, n, N) + (-1) ** (n + 1) * np.eye(dim)   # their m = our n+1
        H += 0.25 * S @ S
    for n in range(N):
        H += mu * (-1) ** (n + 1) * (site_op(Z, n, N) + np.eye(dim)) / 2
    return H


def main(quick=False):
    c = Checks('Reference 1: exact diagonalization')
    out = {}
    N = 6
    Hh, Hm, He = parts(N, M0, G0)
    H = Hh + Hm + He
    herm = max(np.max(abs(M - M.conj().T)) for M in (H, Hh, Hm, He))
    c.less('Hermiticity of H and each part (max |H-H^dag|)', herm, 1e-14, 'Table 1: 0.00')
    c.less('Charge conservation ||[Q,H]||_max', np.max(abs(comm(charge_op(N), H))), 1e-12, 'Table 1: 0.00')

    # Free-field limit m = g = 0 vs sum of negative single-particle energies (1/a) cos(pi j/(N+1))
    free = {}
    for n in (4, 6, 8):
        e0 = eigvalsh(hamiltonian(n, 0.0, 0.0))[0]
        sp = [np.cos(np.pi * j / (n + 1)) / A_DEFAULT for j in range(1, n + 1)]
        ref = sum(e for e in sp if e < 0)
        free[n] = abs(e0 - ref)
        c.less(f'Free-field limit N={n}', free[n], 1e-14, 'Table 1: error < 1e-14')
    out['free_field_error'] = free

    # Strong-mass limit <Z_n> -> -(-1)^n at m=100, chiral condensate -> -1/2
    e, v = ground_state(hamiltonian(N, 100.0, G0))
    dev = max(abs((v.conj() @ site_op(Z, n, N) @ v).real + (-1) ** n) for n in range(N))
    c.less('Strong-mass limit max|<Z_n> + (-1)^n| at m=100', dev, 5e-5, 'Table 1: 0.0000 to 4 d.p.')
    e, v = ground_state(hamiltonian(N, 20.0, G0))
    sig20 = (v.conj() @ condensate_op(N) @ v).real
    c.close('Chiral condensate at m=20', sig20, -0.4997, 1e-4, 'Table 1: -0.4997')
    out['strong_mass_dev_m100'] = dev
    out['condensate_m20'] = sig20

    # Gauss law: <H_E> = g^2 a/2 sum_n <L_n^2> in the ground state
    e, v = ground_state(H)
    dim = 2 ** N
    L = np.zeros((dim, dim), complex)
    acc = 0.0
    for n in range(N - 1):
        L = L + (site_op(Z, n, N) + (-1) ** n * np.eye(dim)) / 2
        acc += (v.conj() @ L @ L @ v).real
    c.less('Gauss law <H_E> = g^2a/2 sum <L_n^2>', abs((v.conj() @ He @ v).real - G0 ** 2 * A_DEFAULT / 2 * acc), 1e-12, 'Table 1: exact match')

    comms = {}
    for n in (4, 6, 8, 10):
        H1, H2 = split_H1_H2(n, M0, G0)
        comms[n] = norm2(comm(H1, H2))
    c.close('||[H1,H2]||_2, N=6', comms[6], 2.5507, 1e-4, 'Table 1 / Table 3')
    c.close('||[H1,H2]||_2, N=4', comms[4], 1.609, 1e-3, 'Table 3')
    c.close('||[H1,H2]||_2, N=8', comms[8], 3.813, 1e-3, 'Table 3')
    c.close('||[H1,H2]||_2, N=10', comms[10], 5.376, 1e-3, 'Table 3')
    out['comm_norm_H1H2'] = comms

    # Convention correspondence with Nguyen et al. (Sec. 2.2)
    for n in (4, 6, 8):
        for (m, g) in [(0.5, 1.0), (1.3, 0.7)]:
            x = 1 / (2 * g ** 2 * A_DEFAULT ** 2)
            mu = 2 * m / (g ** 2 * A_DEFAULT)
            Hpub = nguyen_hamiltonian(n, x, mu)
            Hflip = hamiltonian(n, m, g, stagger=-1)
            diff = np.max(abs(Hflip - Hpub * g ** 2 * A_DEFAULT / 2))
            c.less(f'Nguyen et al. Eq.(4) = flipped-staggering H x 2/(g^2 a), N={n}, m={m}, g={g}', diff, 1e-12, 'Sec. 2.2')
        ea = eigvalsh(hamiltonian(n, 0.5, 1.0))
        eb = eigvalsh(hamiltonian(n, 0.5, 1.0, stagger=-1))
        shift = np.mean(eb - ea)
        c.less(f'Isospectral conventions up to a constant, N={n}', np.max(abs(eb - ea - shift)), 1e-12, 'Sec. 2.2: 1e-13 - 1e-15')

    # Data for Fig. 1
    ms = np.linspace(0, 20, 41)
    fig = {'m': ms.tolist(), 'condensate': [], 'gap': []}
    Q = charge_op(N)
    for m in ms:
        H = hamiltonian(N, m, G0)
        w, V = np.linalg.eigh(H)
        fig['condensate'].append((V[:, 0].conj() @ condensate_op(N) @ V[:, 0]).real)
        q = np.real(np.einsum('ij,ik,kj->j', V.conj(), Q, V))
        e0 = w[np.isclose(q, 0, atol=1e-8)]
        fig['gap'].append(e0[1] - e0[0])
    fig['spectrum_N6'] = eigvalsh(hamiltonian(6, M0, G0)).tolist()
    gs = np.linspace(0, 2, 21)
    fig['g'] = gs.tolist()
    fig['comm_vs_g'] = {str(n): [norm2(comm(*split_H1_H2(n, M0, g))) for g in gs] for n in (4, 6, 8)}
    out['fig1'] = fig
    write_json('ed_reference1.json', out)
    return c.finish()


if __name__ == '__main__':
    run_main(main)
