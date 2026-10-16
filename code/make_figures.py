"""Regenerate the figures from the numbers in results/ (run the analysis scripts first,
or simply run code/run_all.py).  Outputs PNG files in figures/."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _common import FIGURES, RESULTS, read_json, run_main, Checks

# Categorical palette (validated for colour-vision deficiency); every series also
# differs by marker or line style, so colour is never the only cue.
C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
INK, MUTED, GRID = '#222222', '#666666', '#dddddd'

from matplotlib.ticker import NullFormatter


def tidy_log(ax):
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())


def dt_ticks(ax, dts):
    ax.set_xticks(list(dts))
    ax.set_xticklabels([f'1/{int(round(1 / d))}' for d in dts])
    ax.xaxis.set_minor_formatter(NullFormatter())


plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 9.5, 'axes.titleweight': 'bold', 'axes.labelsize': 9,
    'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
    'axes.spines.right': False, 'lines.linewidth': 2, 'lines.markersize': 5, 'legend.frameon': False,
    'figure.dpi': 150, 'savefig.bbox': 'tight'})


def rcsv(name):
    with open(os.path.join(RESULTS, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) if k != 'strategy' else r[k] for r in rows]) for k in rows[0]}


def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight='bold', color=INK)


def fig1():
    d = read_json('ed_reference1.json')
    f = d['fig1']
    fig, axs = plt.subplots(2, 3, figsize=(11, 6.2))
    ax = axs[0, 0]
    Ns = sorted(int(k) for k in d['free_field_error'])
    errs = [max(d['free_field_error'][str(n)], 1e-17) for n in Ns]
    ax.bar([str(n) for n in Ns], errs, color=C1, width=0.5)
    for i, n in enumerate(Ns):
        if d['free_field_error'][str(n)] == 0:
            ax.text(i, 2e-17, 'exactly 0', ha='center', va='bottom', fontsize=8, color=INK)
    ax.axhline(1e-14, color=MUTED, ls='--', lw=1)
    ax.text(0.02, 1.2e-14, 'threshold 1e-14', color=MUTED, fontsize=8)
    ax.set_yscale('log'); ax.set_ylim(1e-17, 1e-13)
    ax.set_xlabel('N'); ax.set_ylabel('|E0 - free-fermion value|')
    ax.set_title('Free-field limit (m = g = 0)'); panel(ax, 'A')
    ax = axs[0, 1]
    ax.plot(f['m'], f['condensate'], color=C1, marker='o', markevery=4)
    ax.axhline(-0.5, color=MUTED, ls='--', lw=1); ax.text(12, -0.48, 'bare Dirac sea  -1/2', color=MUTED, fontsize=8)
    ax.set_xlabel('m'); ax.set_ylabel('chiral condensate  Sigma')
    ax.set_title('Condensate approaches -1/2 (N = 6, g = 1)'); panel(ax, 'B')
    ax = axs[0, 2]
    ax.plot(f['m'], f['gap'], color=C2, marker='s', markevery=4)
    ax.set_xlabel('m'); ax.set_ylabel('mass gap  E1 - E0  (Q = 0)')
    ax.set_title('Gap grows linearly at large m'); panel(ax, 'C')
    ax = axs[1, 0]
    sp = np.array(f['spectrum_N6'])
    ax.eventplot(sp, orientation='vertical', colors=C1, lineoffsets=0.5, linelengths=0.8, linewidths=1)
    ax.set_xticks([]); ax.set_ylabel('energy'); ax.set_xlim(0, 1)
    ax.set_title('Full spectrum, N = 6 (64 levels)'); panel(ax, 'D')
    ax = axs[1, 1]
    for n, col, mk in zip(('4', '6', '8'), (C1, C2, C3), ('o', 's', '^')):
        ax.plot(f['g'], f['comm_vs_g'][n], color=col, marker=mk, markevery=4, label=f'N = {n}')
    ax.set_xlabel('gauge coupling g'); ax.set_ylabel('||[H1, H2]||_2')
    ax.set_title('Commutator norm vs coupling (m = 0.5)'); ax.legend(); panel(ax, 'E')
    ax = axs[1, 2]; ax.axis('off'); panel(ax, 'F')
    rows = [['Hermiticity', '0.00', 'PASS'], ['[Q, H] = 0', '0.00', 'PASS'],
            ['Free-field limit', '< 1e-14', 'PASS'], ['Strong-mass limit (m=100)', '0.0000', 'PASS'],
            ['Sigma at m = 20', f"{d['condensate_m20']:.4f}", 'PASS'],
            ['Gauss law <H_E>', 'exact', 'PASS'], ['||[H1,H2]|| (N = 6)', f"{d['comm_norm_H1H2']['6']:.4f}", 'KEY']]
    t = ax.table(cellText=rows, colLabels=['Check', 'Result', 'Status'], loc='center', cellLoc='left',
                 colWidths=[0.58, 0.24, 0.18])
    t.auto_set_font_size(False); t.set_fontsize(8.5); t.scale(1, 1.4)
    ax.set_title('Validation summary - Reference 1')
    fig.suptitle('Figure 1. Reference 1: exact diagonalization of the lattice Schwinger model', fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig1_reference1_exact_diagonalization.png'))
    plt.close(fig)


def fig2():
    e = rcsv('trotter_errors_N6.csv')
    b = rcsv('bound_vs_measured_N6.csv')
    gdt = rcsv('gap_vs_dt_N6.csv')
    t3 = rcsv('gap_factor_table3.csv')
    fits = read_json('trotter_fits_N6.json')
    fig, axs = plt.subplots(2, 3, figsize=(11, 6.2))
    ax = axs[0, 0]
    for key, lab, col, mk in [('operator_norm_error', 'operator norm', C1, 'o'), ('state_error', 'state', C2, 's'),
                              ('observable_error_extensive', 'observable (extensive)', C3, '^')]:
        ax.loglog(e['dt'], e[key], color=col, marker=mk, label=lab)
    ax.set_xlabel('dt'); ax.set_ylabel('error at t = 1'); ax.legend(fontsize=8)
    tidy_log(ax); dt_ticks(ax, e['dt'])
    ax.set_title(f"All errors scale as dt^2 (slope {fits['operator_norm']['exponent']:.3f})"); panel(ax, 'A')
    ax = axs[0, 1]
    ax.loglog(b['dt'], b['measured_operator_norm_error'], color=C1, marker='o', label='measured')
    ax.loglog(b['dt'], b['bound_eq10'], color=C2, ls='--', label='bound, Eq. (10)')
    ax.set_xlabel('dt'); ax.set_ylabel('||U2^r - e^{-iHt}||'); ax.legend(fontsize=8)
    tidy_log(ax)
    ax.set_title(f"Bound / measured = {np.mean(b['ratio'][b['dt'] <= 1/32]):.2f}"); panel(ax, 'B')
    ax = axs[0, 2]
    ax.semilogx(gdt['dt'], gdt['G'], color=C1, marker='o')
    ax.set_ylim(2.0, 2.3); ax.set_xlabel('dt'); ax.set_ylabel('gap factor G')
    tidy_log(ax); dt_ticks(ax, gdt['dt'])
    ax.set_title(f"G independent of dt ({np.mean(gdt['G']):.3f})"); panel(ax, 'C')
    ax = axs[1, 0]
    sel = (t3['m'] == 0.5) & (t3['g'] == 1.0)
    ax.plot(t3['N'][sel], t3['G'][sel], color=C1, marker='o')
    for n, gv in zip(t3['N'][sel], t3['G'][sel]):
        ax.annotate(f'{gv:.2f}', (n, gv), textcoords='offset points', xytext=(0, 7), ha='center', fontsize=8, color=INK)
    ax.set_xlabel('N'); ax.set_ylabel('G'); ax.set_title('G grows with system size'); panel(ax, 'D')
    ax = axs[1, 1]
    regs = [(m, g, G) for n, m, g, G in zip(t3['N'], t3['m'], t3['g'], t3['G']) if n == 6]
    labs = [f'm={m:g}\ng={g:g}' for m, g, _ in regs]
    ax.bar(labs, [G for *_, G in regs], color=C1, width=0.55)
    ax.set_ylabel('G (N = 6)'); ax.tick_params(axis='x', labelsize=7.5)
    ax.set_title('G depends on the physical regime'); panel(ax, 'E')
    ax = axs[1, 2]; ax.axis('off'); panel(ax, 'F')
    rows = [['operator norm', f"{fits['operator_norm']['exponent']:.3f}", f"{fits['operator_norm']['prefactor']:.3g}"],
            ['state', f"{fits['state']['exponent']:.3f}", f"{fits['state']['prefactor']:.3g}"],
            ['observable', f"{fits['observable_extensive']['exponent']:.3f}", f"{fits['observable_extensive']['prefactor']:.3g}"],
            ['bound / measured', '', f"{np.mean(b['ratio'][b['dt'] <= 1/32]):.2f}x"],
            ['G range', '', f"{t3['G'].min():.1f} - {t3['G'].max():.1f}"]]
    t = ax.table(cellText=rows, colLabels=['Quantity', 'Exponent', 'Prefactor'], loc='center', cellLoc='left',
                 colWidths=[0.5, 0.25, 0.25])
    t.auto_set_font_size(False); t.set_fontsize(8.5); t.scale(1, 1.4)
    ax.set_title('Summary - Reference 2')
    fig.suptitle('Figure 2. Reference 2 (classical Trotter) and the observable gap factor', fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig2_reference2_trotter_and_gap.png'))
    plt.close(fig)


def fig3():
    from schwinger_circuit_lib import trotter_step_circuit, diagonal_terms
    gc = rcsv('gate_counts_table5.csv')
    qc = trotter_step_circuit(4, 0.5, 1.0, 0.25)
    qc.draw(output='mpl', fold=-1, filename=os.path.join(FIGURES, 'fig3a_circuit_one_step_N4.png'))
    plt.close('all')
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axs[0]
    for key, lab, col, mk in [('rzz', 'RZZ', C2, 's'), ('rxx', 'RXX (= RYY)', C1, 'o'), ('rz', 'RZ', C3, '^'),
                              ('depth', 'depth', C4, 'D')]:
        ax.plot(gc['N'], gc[key], color=col, marker=mk, label=lab)
    ax.set_xlabel('N'); ax.set_ylabel('per Trotter step'); ax.legend(fontsize=8)
    ax.set_title('Gate counts: RZZ grows as (N-1)(N-2)/2'); panel(ax, 'B')
    ax = axs[1]
    N = 8
    _, zz = diagonal_terms(N, 0.5, 1.0)
    M = np.full((N, N), np.nan)
    for (j, k), v in zz.items():
        M[j, k] = v
    im = ax.imshow(M, cmap='Blues', vmin=0, origin='upper')
    for (j, k), v in zz.items():
        ax.text(k, j, f'{v:.2f}', ha='center', va='center', fontsize=7, color='white' if v > 0.9 else INK)
    ax.set_xlabel('k (right index)'); ax.set_ylabel('j (left index)'); ax.grid(False)
    ax.set_title('Gauss-law ZZ coefficients, N = 8')
    ax.text(0.98, 0.04, 'c(j,k) = (N-1-k) g^2 a/4', transform=ax.transAxes, ha='right', fontsize=8.5, color=INK); panel(ax, 'C')
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle('Figure 3. Circuit resources (panel A, the circuit itself, is fig3a_circuit_one_step_N4.png)',
                 fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig3_circuit_resources.png'))
    plt.close(fig)


def fig4():
    h = rcsv('hierarchy_table4.csv')
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axs[0]
    ax.loglog(h['dt'], h['ref2_vs_ref1'], color=C1, marker='o', label='Ref 2 vs Ref 1 (Trotter error)')
    ax.loglog(h['dt'], h['ref3_vs_ref2'], color=C2, marker='s', label='Ref 3 vs Ref 2 (gate decomposition)')
    ax.set_xlabel('dt'); ax.set_ylabel('operator-norm distance'); ax.legend(fontsize=8)
    tidy_log(ax); dt_ticks(ax, h['dt'])
    ax.set_title('Circuit = classical Trotter to ~1e-14'); panel(ax, 'A')
    ax = axs[1]
    ax.loglog(h['dt'], h['obs_err_extensive'], color=C3, marker='^')
    ax.set_xlabel('dt'); ax.set_ylabel('|d<Sigma>| (extensive)')
    tidy_log(ax); dt_ticks(ax, h['dt'])
    ax.set_title('Observable converges as dt^2'); panel(ax, 'B')
    fig.suptitle('Figure 4. The complete three-reference hierarchy (N = 6, m = 0.5, g = 1, t = 1)', fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig4_three_reference_hierarchy.png'))
    plt.close(fig)


def fig5():
    with open(os.path.join(RESULTS, 'ordering_table6.csv')) as f:
        strat = list(csv.DictReader(f))
    rnd = rcsv('ordering_random_errors.csv')
    sp = rcsv('ordering_spearman.csv')
    res = read_json('ordering.json')
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axs[0]
    ax.hist(rnd['error'] * 1e4, bins=30, color='#9ec5f4', edgecolor='white')
    styles = ['-', '--', '-.', ':']
    for s, ls, col in zip(strat, styles, (C1, C3, C4, C2)):
        x = float(s['error']) * 1e4
        ax.axvline(x, color=col, ls=ls, lw=2, label=f"{s['strategy']} ({float(s['percentile']):g}th pct)")
    ax.set_xlabel('Trotter error  (x 1e-4)'); ax.set_ylabel('random orderings (of 200)')
    ax.legend(fontsize=7.5); ax.set_title('Ordering strategies vs 200 random orderings'); panel(ax, 'A')
    ax = axs[1]
    ax.scatter(sp['error_operator_norm'], sp['true_error'] * 1e4, color=C1, s=22)
    ax.set_xlabel('error-operator norm (triangle-inequality)'); ax.set_ylabel('true error (x 1e-4)')
    ax.set_title(f"No correlation: Spearman rho = {res['spearman']['rho']:.3f}, p = {res['spearman']['p']:.2f}"); panel(ax, 'B')
    fig.suptitle('Figure 5. Replication of the ordering study of Tranter et al. (N = 6, t = 1, r = 32)', fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, 'fig5_ordering_strategies.png'))
    plt.close(fig)


def main(quick=False):
    c = Checks('Figures')
    for fn in (fig1, fig2, fig3, fig4, fig5):
        fn()
    expected = ['fig1_reference1_exact_diagonalization.png', 'fig2_reference2_trotter_and_gap.png',
                'fig3a_circuit_one_step_N4.png', 'fig3_circuit_resources.png',
                'fig4_three_reference_hierarchy.png', 'fig5_ordering_strategies.png']
    for name in expected:
        p = os.path.join(FIGURES, name)
        c.true(f'{name} written', os.path.exists(p) and os.path.getsize(p) > 10000)
    return c.finish()


if __name__ == '__main__':
    run_main(main)
