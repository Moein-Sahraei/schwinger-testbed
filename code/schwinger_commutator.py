"""Testing the commutator-scaling bound (Sec. 4.3, Eq. 10, Result 1).

Reproduces: ||[H1,[H1,H2]]|| = 3.750, ||[H2,[H1,H2]]|| = 6.254, bound coefficient
0.6774, the bound/measured ratio 1.74 (constant in dt), and the numerically
extracted prefactors c1 = 1/24, c2 = 1/12 from H_eff = i log S(dt)/dt over four
parameter sets.
"""
import numpy as np
from scipy.linalg import logm
from _common import Checks, write_json, write_csv, run_main
from schwinger_core import (split_H1_H2, comm, norm2, strang_step, trotter_propagator,
                            exact_propagator)

N, M0, G0, T = 6, 0.5, 1.0, 1.0


def main(quick=False):
    c = Checks('Commutator-scaling bound (Eq. 10, Result 1)')
    H1, H2 = split_H1_H2(N, M0, G0)
    H = H1 + H2
    C = comm(H1, H2)
    n1, n2 = norm2(comm(H1, C)), norm2(comm(H2, C))
    coef = n1 / 24 + n2 / 12
    c.close('||[H1,[H1,H2]]||', n1, 3.750, 5e-4, 'Sec. 4.3')
    c.close('||[H2,[H1,H2]]||', n2, 6.254, 5e-4, 'Sec. 4.3')
    c.close('Bound coefficient', coef, 0.6774, 5e-5, 'Sec. 4.3')

    # bound / measured across dt.  dt from 1/2 down to 1/8192 (3.6 decades); far
    # smaller dt is limited by round-off accumulated over r = t/dt steps.
    ks = range(1, 14) if not quick else range(3, 10)
    Uex = exact_propagator(H, T)
    rows = []
    for k in ks:
        dt = 1 / 2 ** k
        err = norm2(trotter_propagator(H1, H2, T, 2 ** k) - Uex)
        rows.append((dt, err, coef * T * dt ** 2, coef * T * dt ** 2 / err))
    write_csv('bound_vs_measured_N6.csv', ['dt', 'measured_operator_norm_error', 'bound_eq10', 'ratio'], rows)
    fit_rows = [r for r in rows if r[0] <= 1 / 32]
    ratio_fit = np.mean([r[3] for r in fit_rows])
    c.close('Bound / measured ratio (dt <= 1/32)', ratio_fit, 1.74, 0.005, 'Result 1: 1.74')
    spread = max(r[3] for r in fit_rows) - min(r[3] for r in fit_rows)
    c.less('Ratio constant to 3 significant figures for dt <= 1/32', spread, 0.005, 'Result 1')

    # prefactor extraction from the effective Hamiltonian
    pref = []
    for (n, m, g) in [(4, 0.5, 1.0), (6, 0.5, 1.0), (6, 1.0, 1.0), (6, 0.5, 2.0)]:
        A, B = split_H1_H2(n, m, g)
        Hn = A + B
        K1, K2 = comm(A, comm(A, B)), comm(B, comm(A, B))
        dt = 2e-3
        Heff = 1j * logm(strang_step(A, B, dt)) / dt
        D = ((Heff - Hn) / dt ** 2).ravel()
        Mx = np.stack([K1.ravel(), K2.ravel()], 1)
        sol, *_ = np.linalg.lstsq(Mx, D, rcond=None)
        resid = np.linalg.norm(Mx @ sol - D) / np.linalg.norm(D)
        pref.append({'N': n, 'm': m, 'g': g, 'c1': sol[0].real, 'c2': sol[1].real, 'rel_residual': resid})
        c.close(f'c1 (N={n}, m={m}, g={g})', sol[0].real, 1 / 24, 2e-6, 'Sec. 4.3: 0.0416666 = 1/24')
        c.close(f'c2 (N={n}, m={m}, g={g})', sol[1].real, 1 / 12, 4e-6, 'Sec. 4.3: 0.0833336 = 1/12')
    write_json('commutator_bound_N6.json', {'norm_H1_H1H2': n1, 'norm_H2_H1H2': n2,
                                            'bound_coefficient': coef, 'ratio_mean_dt_le_1_32': ratio_fit,
                                            'prefactor_extraction': pref})
    return c.finish()


if __name__ == '__main__':
    run_main(main)
