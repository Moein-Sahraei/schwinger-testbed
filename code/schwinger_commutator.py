"""Testing the commutator-scaling bound (Sec. 4.3, Eq. 10, Result 1).

Reproduces: ||[H1,[H1,H2]]|| = 3.750, ||[H2,[H1,H2]]|| = 6.254, bound coefficient
0.6774; the bound/measured ratio 1.74 at the reference point, constant to three
significant figures over dt in [1/512, 1/16]; the same ratio at all nine
parameter points of Table 3 (1.51 - 6.84, median 2.01; 1.81 -> 2.31 for
N = 4 -> 10; 1.57 -> 6.84 for g = 0.5 -> 2); and the numerically extracted
prefactors c1 = 1/24, c2 = 1/12 from H_eff = i log S(dt)/dt over four parameter sets.
"""
import numpy as np
from scipy.linalg import logm
from _common import Checks, write_json, write_csv, run_main
from schwinger_core import (split_H1_H2, comm, norm2, strang_step, trotter_propagator,
                            exact_propagator)

N, M0, G0, T = 6, 0.5, 1.0, 1.0
# the nine parameter points of Table 3 (N, m, g)
NINE = [(4, 0.5, 1.0), (6, 0.5, 1.0), (8, 0.5, 1.0), (10, 0.5, 1.0), (6, 0.0, 1.0),
        (6, 1.0, 1.0), (6, 2.0, 1.0), (6, 0.5, 0.5), (6, 0.5, 2.0)]
DT_NINE = 1 / 128   # the dt of Table 3; the ratio is dt-independent at leading order


def bound_coefficient(H1, H2):
    C = comm(H1, H2)
    n1, n2 = norm2(comm(H1, C)), norm2(comm(H2, C))
    return n1, n2, n1 / 24 + n2 / 12


def bound_ratio(N, m, g, dt, t=T):
    """(bound of Eq. 10) / (measured operator-norm error) at one parameter point."""
    H1, H2 = split_H1_H2(N, m, g)
    _, _, coef = bound_coefficient(H1, H2)
    r = int(round(t / dt))
    err = norm2(trotter_propagator(H1, H2, t, r) - exact_propagator(H1 + H2, t))
    return coef * t * dt ** 2 / err, coef, err


def main(quick=False):
    c = Checks('Commutator-scaling bound (Eq. 10, Result 1)')
    H1, H2 = split_H1_H2(N, M0, G0)
    H = H1 + H2
    n1, n2, coef = bound_coefficient(H1, H2)
    c.close('||[H1,[H1,H2]]||', n1, 3.750, 5e-4, 'Sec. 4.3')
    c.close('||[H2,[H1,H2]]||', n2, 6.254, 5e-4, 'Sec. 4.3')
    c.close('Bound coefficient', coef, 0.6774, 5e-5, 'Sec. 4.3')

    # bound / measured across dt.  Written for dt from 1/2 down to 1/8192; the
    # preprint's statement concerns dt in [1/512, 1/16].
    ks = range(1, 14) if not quick else range(3, 10)
    Uex = exact_propagator(H, T)
    rows = []
    for k in ks:
        dt = 1 / 2 ** k
        err = norm2(trotter_propagator(H1, H2, T, 2 ** k) - Uex)
        rows.append((dt, err, coef * T * dt ** 2, coef * T * dt ** 2 / err))
    write_csv('bound_vs_measured_N6.csv', ['dt', 'measured_operator_norm_error', 'bound_eq10', 'ratio'], rows)
    sel = [r[3] for r in rows if 1 / 512 <= r[0] <= 1 / 16]
    ratio_ref = float(np.mean(sel))
    c.close('Bound / measured ratio at the reference point, dt in [1/512, 1/16]', ratio_ref, 1.74, 0.005,
            'Result 1: 1.74')
    c.true('Ratio constant to three significant figures over dt in [1/512, 1/16]',
           all(round(x, 2) == 1.74 for x in sel), f'range {min(sel):.4f} - {max(sel):.4f}')

    # Result 1 extended to the nine parameter points of Table 3
    nine = []
    for (n, m, g) in NINE:
        ratio, cf, err = bound_ratio(n, m, g, DT_NINE)
        nine.append((n, m, g, cf, err, ratio))
    write_csv('bound_ratio_nine_points.csv', ['N', 'm', 'g', 'bound_coefficient',
                                              'measured_error_dt_1_128', 'ratio'], nine)
    ratios = {(n, m, g): r for n, m, g, _, _, r in nine}
    allr = list(ratios.values())
    c.close('Nine points: smallest ratio', min(allr), 1.51, 0.005, 'Result 1: 1.51')
    c.close('Nine points: largest ratio', max(allr), 6.84, 0.005, 'Result 1: 6.84')
    c.close('Nine points: median ratio', float(np.median(allr)), 2.01, 0.005, 'Result 1: median 2.01')
    c.close('Ratio at N=4', ratios[(4, 0.5, 1.0)], 1.81, 0.005, 'Result 1: 1.81 -> 2.31 for N = 4 -> 10')
    c.close('Ratio at N=10', ratios[(10, 0.5, 1.0)], 2.31, 0.005, 'Result 1')
    c.close('Ratio at g=0.5', ratios[(6, 0.5, 0.5)], 1.57, 0.005, 'Result 1: 1.57 -> 6.84 for g = 0.5 -> 2')
    c.close('Ratio at g=2', ratios[(6, 0.5, 2.0)], 6.84, 0.005, 'Result 1')

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
                                            'bound_coefficient': coef,
                                            'ratio_reference_point_dt_1_512_to_1_16': ratio_ref,
                                            'ratio_nine_points': {'min': min(allr), 'max': max(allr),
                                                                  'median': float(np.median(allr))},
                                            'prefactor_extraction': pref})
    return c.finish()


if __name__ == '__main__':
    run_main(main)
