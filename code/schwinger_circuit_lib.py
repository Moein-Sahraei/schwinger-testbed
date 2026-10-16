"""Gate-level Qiskit circuit for one symmetric three-summand Trotter step,
Eq. (12): e^{-i H1e dt/2} e^{-i H1o dt/2} e^{-i H2 dt} e^{-i H1o dt/2} e^{-i H1e dt/2}.

Qiskit conventions: RXX(th)=exp(-i th XX/2), RYY, RZZ likewise, RZ(th)=exp(-i th Z/2).
  exp(-i tau/(4a) (XX+YY)) = RXX(tau/(2a)) RYY(tau/(2a))
  exp(-i c tau Z)          = RZ(2 c tau)
  exp(-i c tau ZZ)         = RZZ(2 c tau)
Identity terms are dropped (global phase)."""
import numpy as np
from qiskit import QuantumCircuit
from schwinger_core import pauli_terms, A_DEFAULT


def diagonal_terms(N, m, g, a=A_DEFAULT):
    """Z and ZZ coefficients of H2 = H_mass + H_E (identity dropped, zeros removed)."""
    t = pauli_terms(N, m, g, a)
    z, zz = {}, {}
    for lab, c in t.items():
        qs = [N - 1 - i for i, ch in enumerate(lab) if ch == 'Z']
        if set(lab) <= {'I', 'Z'}:
            if len(qs) == 1:
                z[qs[0]] = c
            elif len(qs) == 2:
                zz[tuple(sorted(qs))] = c
    return dict(sorted(z.items())), dict(sorted(zz.items()))


def _hop_layer(qc, bonds, tau, a):
    for (j, k) in bonds:
        qc.rxx(tau / (2 * a), j, k)
        qc.ryy(tau / (2 * a), j, k)


def trotter_step_circuit(N, m, g, dt, a=A_DEFAULT):
    even = [(n, n + 1) for n in range(0, N - 1, 2)]
    odd = [(n, n + 1) for n in range(1, N - 1, 2)]
    z, zz = diagonal_terms(N, m, g, a)
    qc = QuantumCircuit(N, name='U(dt)')
    _hop_layer(qc, even, dt / 2, a)
    _hop_layer(qc, odd, dt / 2, a)
    for q, c in z.items():
        qc.rz(2 * c * dt, q)
    for (j, k), c in zz.items():
        qc.rzz(2 * c * dt, j, k)
    _hop_layer(qc, odd, dt / 2, a)
    _hop_layer(qc, even, dt / 2, a)
    return qc


def full_circuit(N, m, g, t, r, a=A_DEFAULT, prep=True, measure=False):
    """State preparation into |1010...> (X on even sites) + r Trotter steps."""
    qc = QuantumCircuit(N, N if measure else 0)
    if prep:
        for q in range(0, N, 2):
            qc.x(q)
    step = trotter_step_circuit(N, m, g, t / r, a)
    for _ in range(r):
        qc.compose(step, inplace=True)
    if measure:
        qc.measure(range(N), range(N))
    return qc


def gate_counts(N, m=0.5, g=1.0, dt=0.1):
    qc = trotter_step_circuit(N, m, g, dt)
    ops = qc.count_ops()
    return {'depth': qc.depth(), 'rxx': ops.get('rxx', 0), 'ryy': ops.get('ryy', 0),
            'rz': ops.get('rz', 0), 'rzz': ops.get('rzz', 0),
            'total': sum(ops.values())}
