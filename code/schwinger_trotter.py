"""Reference 2: classical second-order Trotter evolution (Sec. 4.1-4.2, Table 2).

Evolves |1010...> for t = 1 with U2(dt) = e^{-iH1 dt/2} e^{-iH2 dt} e^{-iH1 dt/2}
and fits operator-norm, state and (extensive) observable errors over
dt in [1/512, 1/32] for N=6, m=0.5, g=1.
"""
import numpy as np
from _common import Checks, write_json, write_csv, run_main
from schwinger_core import (split_H1_H2, exact_propagator, trotter_propagator, vacuum_state,
                            condensate_op, norm2)

N, M0, G0, T = 6, 0.5, 1.0, 1.0


def errors(N, m, g, t, dts):
    H1, H2 = split_H1_H2(N, m, g)
    H = H1 + H2
    Uex = exact_propagator(H, t)
    psi0 = vacuum_state(N)
    pex = Uex @ psi0
    O = condensate_op(N, normalised=False)
    rows = []
    for dt in dts:
        r = int(round(t / dt))
        U = trotter_propagator(H1, H2, t, r)
        p = U @ psi0
        rows.append((dt, r, norm2(U - Uex), np.linalg.norm(p - pex),
                     abs((p.conj() @ O @ p).real - (pex.conj() @ O @ pex).real)))
    return rows


def fit(dts, errs):
    s, c = np.polyfit(np.log(dts), np.log(errs), 1)
    return s, np.exp(c)


def main(quick=False):
    c = Checks('Reference 2: classical Trotter scaling (Table 2)')
    dts = [1 / 2 ** k for k in range(5, 10)]
    rows = errors(N, M0, G0, T, dts)
    write_csv('trotter_errors_N6.csv', ['dt', 'r', 'operator_norm_error', 'state_error', 'observable_error_extensive'], rows)
    res = {}
    targets = {'operator_norm': 0.389, 'state': 0.197, 'observable_extensive': 1.09}
    for i, name in enumerate(['operator_norm', 'state', 'observable_extensive']):
        e = [r[2 + i] for r in rows]
        s, pre = fit(dts, e)
        res[name] = {'exponent': s, 'prefactor': pre}
        c.close(f'{name}: fitted exponent', s, 2.000, 5e-4, 'Table 2: 2.000')
        tol = 0.005 if name != 'observable_extensive' else 0.005
        c.close(f'{name}: prefactor', pre, targets[name], tol, f'Table 2: {targets[name]}')
    write_json('trotter_fits_N6.json', res)
    return c.finish()


if __name__ == '__main__':
    run_main(main)
