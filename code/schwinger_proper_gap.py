"""The observable gap factor G (Sec. 5, Eqs. 11-12, Table 3, Results 2-4).

G = 2 ||O|| ||dpsi|| / |<O>_T - <O>_ex|, invariant under O -> cO.
Reproduces: the 2N normalisation artifact (Sec. 5.1); invariance G = 2.164 for
O, O/(2N), 7.3 O; constancy of G over dt in [1/512, 1/16] with an O(dt^2)
correction of relative coefficient ~0.025; the leading-order prediction of
Eq. (12), G = 2.163812 with E extracted at dt = 1/2048, the 2.2e-4 departure at
dt = 1/16, five-figure agreement for dt <= 1/128, and the anticommutator value
25.79; Table 3; the single-site decomposition of Result 3 (N = 4..10); and
G = 3.31 for the single-site Z_0 (Result 4).
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
LOCAL = {  # N: (||dpsi||, |d<Z_0>|, G)   (preprint Result 3 table, dt = 1/128)
    4: (9.319e-6, 7.2653e-6, 2.565), 6: (1.204e-5, 7.2637e-6, 3.315),
    8: (1.425e-5, 7.2637e-6, 3.924), 10: (1.617e-5, 7.2637e-6, 4.452)}


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


def eq12(U, Uex, psi0, O, dt, anti=False):
    """Leading-order prediction of Eq. (12) with the error operator E extracted
    from e^{iHt} U2^r = 1 - i eps E, eps = t dt^2.  anti=True replaces the
    commutator in the denominator by the anticommutator."""
    eps = T * dt ** 2
    E = 1j * (Uex.conj().T @ U - np.eye(len(U))) / eps
    Ot = Uex.conj().T @ O @ Uex                    # Heisenberg-evolved observable
    K = E @ Ot + Ot @ E if anti else E @ Ot - Ot @ E
    return 2 * norm2(O) * np.linalg.norm(E @ psi0) / abs(psi0.conj() @ K @ psi0)


def main(quick=False):
    c = Checks('Observable gap factor G (Table 3)')
    out = {}
    # 5.1 normalisation artifact: ratio (operator-norm error)/(observable error) changes by 2N
    art = {}
    for N in (4, 6, 8, 10):
        if quick and N == 10:
            continue
        p, pex, _ = evolve(N, 0.5, 1.0, 1 / 128)
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

    # Result 2: constancy in dt, and the size of the O(dt^2) correction
    H1, H2 = split_H1_H2(6, 0.5, 1.0)
    Uex = exact_propagator(H1 + H2, T)
    psi0 = vacuum_state(6)
    pex6 = Uex @ psi0
    Gdt, Upk = [], {}
    for k in range(4, 10):
        dt = 1 / 2 ** k
        U = trotter_propagator(H1, H2, T, 2 ** k)
        Upk[k] = U
        Gdt.append((dt, gap(U @ psi0, pex6, condensate_op(6))[0]))
    for dt, G in Gdt:
        c.close(f'G at dt=1/{int(round(1 / dt))}', G, 2.164, 5e-4, 'Result 2: 2.164 for dt in [1/512, 1/16]')
    dts = np.array([d for d, _ in Gdt]); Gs = np.array([g for _, g in Gdt])
    b, a = np.polyfit(dts ** 2, Gs, 1)            # G(dt) = a + b dt^2 = G0 (1 + coef dt^2)
    coef = b / a
    c.close('Relative O(dt^2) coefficient of G', coef, 0.025, 0.002, 'Result 2: ~0.025')

    # Eq. (12): E extracted once, at dt = 1/2048
    dt0 = 1 / 2048
    U0 = trotter_propagator(H1, H2, T, 2048)
    G12 = eq12(U0, Uex, psi0, O, dt0)
    G12_anti = eq12(U0, Uex, psi0, O, dt0, anti=True)
    c.close('Eq. (12) prediction, E extracted at dt = 1/2048', G12, 2.163812, 2e-6, 'Sec. 5.3: 2.163812')
    c.close('dt -> 0 limit of the measured G agrees with Eq. (12)', a, G12, 5e-6, 'fit of G = G0 + b dt^2')
    dev16 = Gdt[0][1] - G12
    c.close('Measured G departs from Eq. (12) at dt = 1/16 by', dev16, 2.2e-4, 0.1e-4, 'Sec. 5.3: 2.2e-4')
    c.true('Five-figure agreement for dt <= 1/128', all(abs(G - G12) < 1e-5 for dt, G in Gdt if dt <= 1 / 128),
           ', '.join(f'1/{int(round(1 / dt))}: {G - G12:+.1e}' for dt, G in Gdt if dt <= 1 / 128))
    c.close('Anticommutator form of the denominator', G12_anti, 25.79, 0.005, 'Sec. 5.3: 25.79 instead of 2.1638')
    # E re-extracted at each dt: E is dt-independent at leading order, so the prediction does not move
    per_dt = [(1 / 2 ** k, eq12(Upk[k], Uex, psi0, O, 1 / 2 ** k)) for k in range(4, 10)]
    c.true('Eq. (12) with E re-extracted at each dt gives the same prediction',
           all(abs(g - G12) < 3e-6 for _, g in per_dt), f'{min(g for _, g in per_dt):.7f} - {max(g for _, g in per_dt):.7f}')
    write_csv('gap_vs_dt_N6.csv', ['dt', 'G_measured', 'G_eq12_reextracted'],
              [(dt, G, g) for (dt, G), (_, g) in zip(Gdt, per_dt)])
    out['eq12'] = {'G_eq12_dt_1_2048': G12, 'G_anticommutator': G12_anti, 'fit_G0': a,
                   'relative_dt2_coefficient': coef, 'deviation_at_dt_1_16': dev16}

    # Table 3 (and the states reused for Result 3)
    rows, states = [], {}
    for (N, m, g, cn, se, Gt) in TABLE3:
        if quick and N == 10:
            continue
        p, pex, cnorm = evolve(N, m, g, 1 / 128)
        if (m, g) == (0.5, 1.0):
            states[N] = (p, pex)
        G, d, _ = gap(p, pex, condensate_op(N))
        rows.append((N, m, g, cnorm, d, G))
        c.close(f'Table 3 ||[H1,H2]|| (N={N}, m={m}, g={g})', cnorm, cn, 1e-3)
        c.close(f'Table 3 state err (N={N}, m={m}, g={g})', d, se, 0.006 * se)
        c.close(f'Table 3 G (N={N}, m={m}, g={g})', G, Gt, 0.006)
    write_csv('gap_factor_table3.csv', ['N', 'm', 'g', 'comm_norm_H1H2', 'state_error', 'G'], rows)

    # Result 3: decomposition for the single-site observable Z_0 (||O|| = 1 for every N)
    loc = []
    for N, (dpsi_t, dz_t, G_t) in LOCAL.items():
        if N not in states:
            continue
        p, pex = states[N]
        Gz, d, dz = gap(p, pex, site_op(Z, 0, N))
        Gs, _, ds = gap(p, pex, condensate_op(N, False))
        loc.append((N, d, dz, Gz, ds, Gs))
        c.close(f'Result 3: ||dpsi|| (N={N})', d, dpsi_t, 0.0006 * dpsi_t)
        c.close(f'Result 3: |d<Z_0>| (N={N})', dz, dz_t, 6e-11)
        c.close(f'Result 3: G for Z_0 (N={N})', Gz, G_t, 6e-4)
    write_csv('gap_local_observable.csv', ['N', 'state_error', 'obs_error_Z0', 'G_Z0',
                                           'obs_error_Sigma_extensive', 'G_Sigma'], loc)
    L = {r[0]: r for r in loc}
    big = [L[n][2] for n in (6, 8, 10) if n in L]
    c.true('Local-observable error identical to five figures for N = 6, 8, 10' if not quick else
           'Local-observable error identical to five figures for N = 6, 8', (max(big) - min(big)) / min(big) < 1e-5,
           f'relative spread {(max(big) - min(big)) / min(big):.1e}')
    if 4 in L and 10 in L:
        gs = L[10][1] / L[4][1]
        c.close('Result 3: state error grows N=4 -> 10 by', gs, 1.735, 6e-4, 'Result 3: 1.735')
        c.close('Result 3: G(Z_0) grows by the same factor', (L[10][3] / L[4][3]) / gs, 1.0, 5e-4, 'Result 3')
        c.close('Result 3: |d<Sigma>| (extensive) grows by', L[10][4] / L[4][4], 2.83, 0.005, 'Result 3: 2.83')
        c.close('Result 3: G(Sigma) grows by', L[10][5] / L[4][5], 1.53, 0.005, 'Result 3: 1.53')

    # Result 4: single-site observable at the reference point
    p, pex = states[6]
    Gz = gap(p, pex, site_op(Z, 0, 6))[0]
    c.close('G for single-site Z_0 (N=6)', Gz, 3.31, 0.006, 'Result 4: 3.31')
    out['G_Z0'] = Gz
    write_json('gap_factor.json', out)
    return c.finish()


if __name__ == '__main__':
    run_main(main)
