"""The complete three-reference hierarchy (Sec. 6.2, Table 4, Result 5).

N=6, m=0.5, g=1, t=1.  Ref 3 = Qiskit circuit (r steps), Ref 2 = classical
three-summand Trotter, Ref 1 = exact evolution.  Also fits the Trotter exponent
(preprint: 2.003) and reports the extensive condensate error.
"""
import numpy as np
from qiskit.quantum_info import Operator
from _common import Checks, write_csv, write_json, run_main
from schwinger_core import (hamiltonian, exact_propagator, vacuum_state, condensate_op, norm2,
                            phase_aligned_distance)
from schwinger_circuit_lib import trotter_step_circuit
from schwinger_circuit import classical_three_summand

N, M0, G0, T = 6, 0.5, 1.0, 1.0
TABLE4 = {4: (3.2e-15, 2.58e-2, 7.16e-2), 8: (9.3e-15, 6.40e-3, 1.79e-2), 16: (2.0e-14, 1.60e-3, 4.46e-3),
          32: (4.7e-14, 3.99e-4, 1.12e-3), 64: (7.6e-14, 9.97e-5, 2.79e-4)}


def main(quick=False):
    c = Checks('Three-reference hierarchy (Table 4)')
    H = hamiltonian(N, M0, G0)
    Uex = exact_propagator(H, T)
    psi0 = vacuum_state(N)
    pex = Uex @ psi0
    O = condensate_op(N, normalised=False)
    rows = []
    for r, (t32, t21, tob) in TABLE4.items():
        dt = T / r
        step_c = Operator(trotter_step_circuit(N, M0, G0, dt)).data
        step_cl, _, _ = classical_three_summand(N, M0, G0, dt)
        U3 = np.linalg.matrix_power(step_c, r)
        U2 = np.linalg.matrix_power(step_cl, r)
        d32 = phase_aligned_distance(U3, U2)
        d21 = norm2(U2 - Uex)
        p = U2 @ psi0
        dob = abs((p.conj() @ O @ p).real - (pex.conj() @ O @ pex).real)
        rows.append((r, dt, d32, d21, dob))
        c.less(f'r={r}: Ref 3 vs Ref 2 (machine precision)', d32, 1e-12, f'preprint {t32:.1e}')
        c.close(f'r={r}: Ref 2 vs Ref 1', d21, t21, 0.006 * t21, f'preprint {t21:.3g}')
        c.close(f'r={r}: |d<Sigma>| (extensive)', dob, tob, 0.006 * tob, f'preprint {tob:.3g}')
    write_csv('hierarchy_table4.csv', ['r', 'dt', 'ref3_vs_ref2', 'ref2_vs_ref1', 'obs_err_extensive'], rows)
    dts = [x[1] for x in rows]
    s = np.polyfit(np.log(dts), np.log([x[3] for x in rows]), 1)[0]
    c.close('Fitted Trotter exponent (Ref 2 vs Ref 1, r = 4..64)', s, 2.003, 0.0015, 'Result 5: 2.003')
    write_json('hierarchy.json', {'fitted_exponent': s})
    return c.finish()


if __name__ == '__main__':
    run_main(main)
