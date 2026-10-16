"""The observable gap factor G (Sec. 5, Eq. 11, Table 3, Results 2-4).

G = 2 ||O|| ||dpsi|| / |<O>_T - <O>_ex|, invariant under O -> cO.
Reproduces: the 2N normalisation artifact (Sec. 5.1), invariance G = 2.164 for
O, O/(2N), 7.3 O, constancy of G over dt in [1/512, 1/16], Table 3, and
G = 3.31 for the single-site Z_0.
"""
import numpy as np
from _common import Checks, write_json, write_csv, run_main
from schwinger_core import (split_H1_H2, exact_propagator, trotter_propagator, vacuum_state,
                            condensate_op, site_op, Z, norm2, comm)

T = 1.0
TABLE3 = [  # N, m, g, ||[H1,H2]||, state err, G   (preprint Table 3, dt = 1/128)
    (4, 0.5, 1.0, 1.609, 9.32e-6, 1.80), (6, 0.5, 1.0, 2.551, 1.20e-5, 2.16),
    (8, 0.5, 1.0, 3.813, 1.43e-5, 2.48), (10, 0.5, 1.0, 5.376, 1.62e-5, 2.76),
    (6, 0.0, 1.0, 1.644, 2.69e-6, 3.66), (6, 1.0, 1.0, 4.210, 2.68e-5, 2.56),
    (6, 2.0, 1.0, 7.672, 4.62e-5, 8.10), (6, 0.5, 0.5, 1.918, 7.79e-6, 2.28),
    (6, 0.5, 2.0, 7.277, 3.46e-5, 3.05)]


def evolve(N, m, g, dt):
    H1, H2 = split_H1_H2(N, m, g)
    psi0 = vacuum_state(N)
    pex = exact_propagator(H1 + H2, T) @ psi0
    p = trotter_propagator(H1, H2, T, int(round(T / dt))) @ psi0
    return p, pex, norm2(comm(H1, H2))


def gap(p, pex, O):
    d = np.linalg.norm(p - pex)
    do = abs((p.conj() @ O @ p).real - (pex.conj() @ O @ pex).real)
    return 2 * norm2(O) * d / do, d, do


def main(quick=False):
    c = Checks('Observable gap factor G (Table 3)')
    out = {}
    # 5.1 normalisation artifact: ratio (operator-norm error)/(observable error) changes by 2N
    art = {}
    for N in (4, 6, 8, 10):
        if quick and N == 10:
            continue
        H1, H2 = split_H1_H2(N, 0.5, 1.0)
        dt = 1 / 128
        p, pex, _ = evolve(N, 0.5, 1.0, dt)
        _, _, d_norm = gap(p, pex, condensate_op(N, True))
        _, _, d_ext = gap(p, pex, condensate_op(N, False))
        art[N] = d_ext / d_norm
        c.close(f'Normalisation artifact: ratio changes by exactly 2N, N={N}', art[N], 2 * N, 1e-9, 'Sec. 5.1')
    out['normalisation_ratio'] = art

    # 5.2 invariance
    p, pex, _ = evolve(6, 0.5, 1.0, 1 / 128)
    O = condensate_op(6, False)
    for lab, OO in [('O', O), ('O/(2N)', O / 12), ('7.3 O', 7.3 * O)]:
        c.close(f'G invariant under scaling, {lab}', gap(p, pex, OO)[0], 2.164, 5e-4, 'Sec. 5.2: 2.164')

    # Result 2: constancy in dt
    Gdt = []
    for k in range(4, 10):
        p2, pex2, _ = evolve(6, 0.5, 1.0, 1 / 2 ** k)
        Gdt.append((1 / 2 ** k, gap(p2, pex2, condensate_op(6))[0]))
    write_csv('gap_vs_dt_N6.csv', ['dt', 'G'], Gdt)
    for dt, G in Gdt:
        c.close(f'G at dt=1/{int(round(1 / dt))}', G, 2.164, 5e-4, 'Result 2: 2.164 for dt in [1/512, 1/16]')

    # Table 3
    rows = []
    for (N, m, g, cn, se, Gt) in TABLE3:
        if quick and N == 10:
            continue
        p, pex, cnorm = evolve(N, m, g, 1 / 128)
        G, d, _ = gap(p, pex, condensate_op(N))
        rows.append((N, m, g, cnorm, d, G))
        c.close(f'Table 3 ||[H1,H2]|| (N={N}, m={m}, g={g})', cnorm, cn, 1e-3)
        c.close(f'Table 3 state err (N={N}, m={m}, g={g})', d, se, 0.006 * se)
        c.close(f'Table 3 G (N={N}, m={m}, g={g})', G, Gt, 0.006)
    write_csv('gap_factor_table3.csv', ['N', 'm', 'g', 'comm_norm_H1H2', 'state_error', 'G'], rows)

    # Result 4: single-site observable
    p, pex, _ = evolve(6, 0.5, 1.0, 1 / 128)
    Gz = gap(p, pex, site_op(Z, 0, 6))[0]
    c.close('G for single-site Z_0 (N=6)', Gz, 3.31, 0.006, 'Result 4: 3.31')
    out['G_Z0'] = Gz
    write_json('gap_factor.json', out)
    return c.finish()


if __name__ == '__main__':
    run_main(main)
