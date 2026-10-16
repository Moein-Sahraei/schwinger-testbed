"""
Core definitions for the lattice Schwinger model testbed.

Kogut-Susskind staggered fermions, 1+1D QED, N sites, open boundaries,
gauge field eliminated via Gauss's law (Sahraei 2026, Eqs. 2-3):

  H = 1/(2a) sum_{n=0}^{N-2} (s+_n s-_{n+1} + h.c.)        (H_hop)
    + m/2   sum_{n=0}^{N-1} (-1)^n Z_n                      (H_mass)
    + g^2 a/2 sum_{n=0}^{N-2} L_n^2                         (H_E)
  L_n = sum_{k<=n} q_k,   q_k = (Z_k + (-1)^k)/2

Qubit convention: site n <-> qubit n, little-endian (as in Qiskit):
the dense matrix is kron(op_{N-1}, ..., op_1, op_0).
"""
import numpy as np
from scipy.linalg import expm, eigh

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma^+ = |0><1|
SM = SP.conj().T                                 # sigma^-

A_DEFAULT = 1.0  # lattice spacing used throughout the preprint


def site_op(op, n, N):
    """Single-site operator `op` on site n of an N-site chain (little-endian)."""
    out = np.array([[1.0 + 0j]])
    for k in reversed(range(N)):
        out = np.kron(out, op if k == n else I2)
    return out


def two_site_op(op1, n1, op2, n2, N):
    out = np.array([[1.0 + 0j]])
    for k in reversed(range(N)):
        if k == n1:
            out = np.kron(out, op1)
        elif k == n2:
            out = np.kron(out, op2)
        else:
            out = np.kron(out, I2)
    return out


def parts(N, m, g, a=A_DEFAULT, stagger=+1):
    """Return (H_hop, H_mass, H_E) as dense matrices built directly from
    sigma^+-, Z and the Gauss-law electric field (no Pauli expansion).
    stagger=+1 is the convention of the preprint, (-1)^k with k = 0..N-1;
    stagger=-1 flips it (convention of Nguyen et al. 2022)."""
    dim = 2 ** N
    Hhop = np.zeros((dim, dim), complex)
    for n in range(N - 1):
        Hhop += (two_site_op(SP, n, SM, n + 1, N) + two_site_op(SM, n, SP, n + 1, N)) / (2 * a)
    Hmass = sum(stagger * (-1) ** n * site_op(Z, n, N) for n in range(N)) * (m / 2)
    HE = np.zeros((dim, dim), complex)
    L = np.zeros((dim, dim), complex)
    for n in range(N - 1):
        L = L + (site_op(Z, n, N) + stagger * (-1) ** n * np.eye(dim)) / 2
        HE += L @ L
    HE *= g ** 2 * a / 2
    return Hhop, Hmass, HE


def hamiltonian(N, m, g, a=A_DEFAULT, stagger=+1):
    Hh, Hm, He = parts(N, m, g, a, stagger)
    return Hh + Hm + He


def split_H1_H2(N, m, g, a=A_DEFAULT):
    Hh, Hm, He = parts(N, m, g, a)
    return Hh, Hm + He


def charge_op(N):
    """Total charge Q = sum_k q_k."""
    dim = 2 ** N
    return sum((site_op(Z, k, N) + (-1) ** k * np.eye(dim)) / 2 for k in range(N))


def condensate_op(N, normalised=True):
    """Chiral condensate Sigma = 1/(2N) sum_n (-1)^n Z_n (normalised) or the
    extensive operator sum_n (-1)^n Z_n."""
    O = sum((-1) ** n * site_op(Z, n, N) for n in range(N))
    return O / (2 * N) if normalised else O


def vacuum_state(N):
    """Strong-coupling vacuum |1010...> (site 0 occupied by |1>): the bare Dirac sea."""
    idx = sum(1 << n for n in range(0, N, 2))
    psi = np.zeros(2 ** N, complex)
    psi[idx] = 1.0
    return psi


def comm(A, B):
    return A @ B - B @ A


def norm2(A):
    """Spectral (operator 2-) norm."""
    return np.linalg.norm(A, 2)


# ---------------------------------------------------------------- Pauli form
def pauli_terms(N, m, g, a=A_DEFAULT, drop_identity=True, tol=1e-12):
    """Pauli decomposition (Eqs. 6-8), merged by Pauli string.
    Returns dict {label: coeff}; label is a little-endian string as in
    Qiskit (rightmost character = qubit 0)."""
    terms = {}

    def add(ops, c):
        lab = ['I'] * N
        for q, p in ops:
            lab[N - 1 - q] = p
        key = ''.join(lab)
        terms[key] = terms.get(key, 0.0) + c

    for n in range(N - 1):                     # hopping, Eq. (6)
        add([(n, 'X'), (n + 1, 'X')], 1 / (4 * a))
        add([(n, 'Y'), (n + 1, 'Y')], 1 / (4 * a))
    for n in range(N):                         # mass, Eq. (7)
        add([(n, 'Z')], m / 2 * (-1) ** n)
    pref = g ** 2 * a / 2
    for n in range(N - 1):                     # electric, Eq. (8)
        c_n = sum((-1) ** k for k in range(n + 1))
        add([], pref * 0.25 * ((n + 1) + c_n ** 2))
        for j in range(n + 1):
            for k in range(j + 1, n + 1):
                add([(j, 'Z'), (k, 'Z')], pref * 0.25 * 2)
        for k in range(n + 1):
            add([(k, 'Z')], pref * 0.25 * 2 * c_n)
    out = {}
    for k, v in terms.items():
        if drop_identity and set(k) == {'I'}:
            continue
        if abs(v) > tol:
            out[k] = v
    return out


PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def pauli_matrix(label):
    out = np.array([[1.0 + 0j]])
    for ch in label:            # leftmost char = highest qubit -> first kron factor
        out = np.kron(out, PAULI[ch])
    return out


def matrix_from_terms(terms):
    lab0 = next(iter(terms))
    dim = 2 ** len(lab0)
    M = np.zeros((dim, dim), complex)
    for lab, c in terms.items():
        M += c * pauli_matrix(lab)
    return M


# ---------------------------------------------------------------- Trotter
def strang_step(H1, H2, dt):
    """Second-order symmetric product formula, Eq. (9)."""
    A = expm(-1j * H1 * dt / 2)
    return A @ expm(-1j * H2 * dt) @ A


def trotter_propagator(H1, H2, t, r):
    U = strang_step(H1, H2, t / r)
    return np.linalg.matrix_power(U, r)


def exact_propagator(H, t):
    return expm(-1j * H * t)


def phase_aligned_distance(U, V):
    """min_phi ||U - e^{i phi} V||_2 (operators compared up to a global phase)."""
    ov = np.trace(V.conj().T @ U)
    ph = ov / abs(ov) if abs(ov) > 0 else 1.0
    return norm2(U - ph * V)


def ground_state(H):
    w, v = eigh(H)
    return w[0], v[:, 0]


def report(name, ok, detail=''):
    tag = 'PASS' if ok else 'FAIL'
    print(f'  [{tag}] {name}' + (f'  ({detail})' if detail else ''))
    return ok
