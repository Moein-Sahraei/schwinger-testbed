"""Replication of the ordering strategies of Tranter et al. (Sec. 7, Table 6, Results 6-7).

N=6, m=0.5, g=1, t=1, r=32.  25 distinct Pauli terms (repeated strings merged,
identity dropped), each exponentiated separately in a symmetric second-order
formula.  Error = ||U^r - e^{-iHt}||_2.
  magnitude      : descending |coefficient| (stable sort on the term list)
  lexicographic  : sorted Pauli strings written site 0 first
  depleteGroups / equaliseGroups : on commuting sets from greedy colouring
                   (rustworkx, Degree strategy) of the anticommutation graph
  random         : 200 orderings from numpy.random.default_rng(7)
Result 7 uses the next 40 orderings from the same generator and the
triangle-inequality norm of the second-order error operator of Tranter et al.,
  (1/12) sum_{b} sum_{a<=b} sum_{c<b} (1 - delta_ab/2) ||[H_a,[H_b,H_c]]||.
"""
import numpy as np
import rustworkx as rx
from scipy.linalg import expm
from scipy.stats import spearmanr
from _common import Checks, write_json, write_csv, run_main
from schwinger_core import pauli_terms, pauli_matrix, norm2

N, M0, G0, T, R = 6, 0.5, 1.0, 1.0, 32


def anticommute(p, q):
    return sum(1 for x, y in zip(p, q) if x != 'I' and y != 'I' and x != y) % 2 == 1


def product_label(p, q):
    out = []
    for x, y in zip(p, q):
        if x == 'I':
            out.append(y)
        elif y == 'I':
            out.append(x)
        elif x == y:
            out.append('I')
        else:
            out.append(({'X', 'Y', 'Z'} - {x, y}).pop())
    return ''.join(out)


class Model:
    def __init__(self):
        terms = pauli_terms(N, M0, G0)
        self.labels = list(terms)
        self.coef = np.array([terms[l] for l in self.labels])
        self.n = len(self.labels)
        mats = [cf * pauli_matrix(l) for l, cf in zip(self.labels, self.coef)]
        H = sum(mats)
        self.Uex = expm(-1j * H * T)
        dt = T / R
        self.half = [expm(-1j * M * dt / 2) for M in mats]
        # ||[P_a,[P_b,P_c]]|| for Pauli terms is 4|c_a c_b c_c| or 0
        L = self.labels
        self.K = np.zeros((self.n,) * 3)
        for b in range(self.n):
            for cc in range(self.n):
                if anticommute(L[b], L[cc]):
                    pbc = product_label(L[b], L[cc])
                    for a in range(self.n):
                        if anticommute(L[a], pbc):
                            self.K[a, b, cc] = 4 * abs(self.coef[a] * self.coef[b] * self.coef[cc])

    def error(self, order):
        S = np.eye(2 ** N, dtype=complex)
        for i in order:
            S = self.half[i] @ S
        for i in reversed(order):
            S = self.half[i] @ S
        return norm2(np.linalg.matrix_power(S, R) - self.Uex)

    def error_operator_norm(self, order):
        s = 0.0
        for bi in range(self.n):
            for ai in range(bi + 1):
                f = 0.5 if ai == bi else 1.0
                for ci in range(bi):
                    s += f * self.K[order[ai], order[bi], order[ci]]
        return s / 12

    def groups(self):
        G = rx.PyGraph()
        G.add_nodes_from(range(self.n))
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if anticommute(self.labels[i], self.labels[j]):
                    G.add_edge(i, j, None)
        col = rx.graph_greedy_color(G, strategy=rx.ColoringStrategy.Degree)
        k = max(col.values()) + 1
        return [[i for i in range(self.n) if col[i] == cc] for cc in range(k)], G.num_edges()


def magnitude_order(mdl):
    return sorted(range(mdl.n), key=lambda i: -abs(mdl.coef[i]))


def lexicographic_order(mdl):
    return sorted(range(mdl.n), key=lambda i: mdl.labels[i][::-1])   # site 0 written first


def deplete_groups(mdl, groups):
    gs = [sorted(g, key=lambda i: -abs(mdl.coef[i])) for g in groups]
    out = []
    while any(gs):
        for g in gs:
            if g:
                out.append(g.pop(0))
    return out


def equalise_groups(mdl, groups):
    gs = [sorted(g, key=lambda i: -abs(mdl.coef[i])) for g in groups]
    out = []
    while any(gs):
        L = max(len(g) for g in gs)
        best = max((g for g in gs if len(g) == L), key=lambda g: abs(mdl.coef[g[0]]))
        out.append(best.pop(0))
    return out


def main(quick=False):
    c = Checks('Ordering strategies (Table 6, Results 6-7)')
    mdl = Model()
    groups, n_edges = mdl.groups()
    c.true('25 Pauli terms', mdl.n == 25)
    c.true('80 non-commuting pairs (density 0.267)', n_edges == 80, f'{n_edges} pairs, density {n_edges / 300:.3f}')
    rng = np.random.default_rng(7)
    rand = [list(rng.permutation(mdl.n)) for _ in range(200)]
    rand_err = np.array([mdl.error(o) for o in rand])
    write_csv('ordering_random_errors.csv', ['ordering_index', 'error'], list(enumerate(rand_err)))
    c.close('Random orderings: minimum error', rand_err.min(), 3.23e-4, 0.006e-4, 'Table 6 caption')
    c.close('Random orderings: maximum error', rand_err.max(), 1.42e-3, 0.006e-3, 'Table 6 caption')
    c.close('Random orderings: span factor', rand_err.max() / rand_err.min(), 4.40, 0.006, 'Table 6 caption')
    strategies = {'magnitude': magnitude_order(mdl), 'equaliseGroups': equalise_groups(mdl, groups),
                  'depleteGroups': deplete_groups(mdl, groups), 'lexicographic': lexicographic_order(mdl)}
    target = {'magnitude': (3.38e-4, 1.000, 0.5), 'equaliseGroups': (4.35e-4, 1.286, 2.5),
              'depleteGroups': (5.36e-4, 1.586, 11.0), 'lexicographic': (1.44e-3, 4.261, 100.0)}
    res, rows = {}, []
    e_mag = mdl.error(strategies['magnitude'])
    for name, order in strategies.items():
        e = mdl.error(order)
        pct = 100 * np.mean(rand_err <= e)
        res[name] = {'error': e, 'vs_magnitude': e / e_mag, 'percentile': pct, 'order': [mdl.labels[i] for i in order]}
        rows.append((name, e, e / e_mag, pct))
        te, tr, tp = target[name]
        c.close(f'{name}: error', e, te, 0.006 * te, 'Table 6')
        c.close(f'{name}: ratio vs magnitude', e / e_mag, tr, 0.002, 'Table 6')
        c.close(f'{name}: percentile of 200 random orderings', pct, tp, 0.01, 'Table 6')
    write_csv('ordering_table6.csv', ['strategy', 'error', 'vs_magnitude', 'percentile'], rows)
    c.close('depleteGroups / equaliseGroups', res['depleteGroups']['error'] / res['equaliseGroups']['error'],
            1.23, 0.005, 'Sec. 7: factor 1.23')

    # Result 7: next 40 orderings from the same generator
    more = [list(rng.permutation(mdl.n)) for _ in range(40)]
    true_err = [mdl.error(o) for o in more]
    op_norm = [mdl.error_operator_norm(o) for o in more]
    rho = spearmanr(op_norm, true_err)
    write_csv('ordering_spearman.csv', ['error_operator_norm', 'true_error'], list(zip(op_norm, true_err)))
    c.close('Spearman rho (error-operator norm vs true error, 40 orderings)', rho.statistic, -0.081, 0.0015, 'Result 7')
    c.close('Spearman p-value', rho.pvalue, 0.62, 0.006, 'Result 7')
    res['spearman'] = {'rho': rho.statistic, 'p': rho.pvalue}
    res['groups'] = [[mdl.labels[i] for i in g] for g in groups]
    write_json('ordering.json', res)
    return c.finish()


if __name__ == '__main__':
    run_main(main)
