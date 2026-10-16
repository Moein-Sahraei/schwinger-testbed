"""Reference 3: the gate-level Qiskit circuit (Sec. 6.1-6.3, Table 5).

Reproduces: single-step agreement between the Qiskit circuit and the classical
three-summand propagator at ~1e-15 for many (N, m, g, dt); the gate counts and
depth per Trotter step for N = 4..12 (Table 5); the R_ZZ count (N-1)(N-2)/2;
the circuit depth 513 for N=6, r=32 (Fig. 4F); the checkerboard commutators.
"""
import numpy as np
from scipy.linalg import expm
from qiskit.quantum_info import Operator
from _common import Checks, write_json, write_csv, run_main
from schwinger_core import parts, two_site_op, SP, SM, comm, norm2, phase_aligned_distance
from schwinger_circuit_lib import trotter_step_circuit, full_circuit, gate_counts

TABLE5 = {4: (11, 6, 6, 3, 3, 18), 6: (16, 10, 10, 5, 10, 35), 8: (20, 14, 14, 7, 21, 56),
          10: (24, 18, 18, 9, 36, 81), 12: (28, 22, 22, 11, 55, 110)}


def hop_part(N, bonds, a=1.0):
    dim = 2 ** N
    M = np.zeros((dim, dim), complex)
    for j, k in bonds:
        M += (two_site_op(SP, j, SM, k, N) + two_site_op(SM, j, SP, k, N)) / (2 * a)
    return M


def classical_three_summand(N, m, g, dt):
    """Classical propagator of Eq. (12) built from dense matrices (independent of Qiskit)."""
    Hh, Hm, He = parts(N, m, g)
    H1e = hop_part(N, [(n, n + 1) for n in range(0, N - 1, 2)])
    H1o = hop_part(N, [(n, n + 1) for n in range(1, N - 1, 2)])
    A, B = expm(-1j * H1e * dt / 2), expm(-1j * H1o * dt / 2)
    return A @ B @ expm(-1j * (Hm + He) * dt) @ B @ A, H1e, H1o


def main(quick=False):
    c = Checks('Reference 3: gate-level circuit')
    # checkerboard structure
    H1e = hop_part(4, [(0, 1), (2, 3)])
    H1o = hop_part(4, [(1, 2)])
    c.true('Even bonds commute exactly', np.max(abs(comm(hop_part(4, [(0, 1)]), hop_part(4, [(2, 3)])))) == 0.0)
    c.close('||[H1e, H1o]||_2 for N=4', norm2(comm(H1e, H1o)), 0.500, 1e-12, 'Sec. 6.1: 0.500')

    worst = 0.0
    params = [(N, m, g, dt) for N in (4, 6, 8) for m in (0.0, 0.5, 2.0) for g in (0.5, 1.0, 2.0)
              for dt in (0.25, 0.05)]
    if quick:
        params = [p for p in params if p[0] <= 6]
    for (N, m, g, dt) in params:
        Uc = Operator(trotter_step_circuit(N, m, g, dt)).data
        Ucl, _, _ = classical_three_summand(N, m, g, dt)
        worst = max(worst, phase_aligned_distance(Uc, Ucl))
    c.less(f'Circuit step = classical propagator up to global phase, worst over {len(params)} settings',
           worst, 1e-13, 'Result 5: ~1e-15')

    rows = []
    for N, ref in TABLE5.items():
        gc = gate_counts(N)
        got = (gc['depth'], gc['rxx'], gc['ryy'], gc['rz'], gc['rzz'], gc['total'])
        rows.append((N,) + got)
        c.true(f'Table 5 gate counts N={N} (depth, RXX, RYY, RZ, RZZ, total)', got == ref, f'got {got}, preprint {ref}')
        c.true(f'R_ZZ count = (N-1)(N-2)/2 for N={N}', gc['rzz'] == (N - 1) * (N - 2) // 2)
    write_csv('gate_counts_table5.csv', ['N', 'depth', 'rxx', 'ryy', 'rz', 'rzz', 'total'], rows)
    d = full_circuit(6, 0.5, 1.0, 1.0, 32).depth()
    c.true('Full circuit depth N=6, r=32 (with state preparation)', d == 513, f'got {d}, preprint 513')
    write_json('circuit_validation.json', {'worst_single_step_distance': worst, 'depth_N6_r32': d,
                                           'n_settings': len(params)})
    return c.finish()


if __name__ == '__main__':
    run_main(main)
