# Schwinger testbed

**A validated three-reference quantum-circuit testbed for the lattice Schwinger model**

[![check](https://github.com/Moein-Sahraei/schwinger-testbed/actions/workflows/check.yml/badge.svg)](https://github.com/Moein-Sahraei/schwinger-testbed/actions/workflows/check.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21854391.svg)](https://doi.org/10.5281/zenodo.21854391)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Moein Sahraei, Department of Physics, Razi University, Kermanshah, Iran ·
ORCID [0009-0008-8672-3934](https://orcid.org/0009-0008-8672-3934)

This repository contains the code for the preprint
*A Validated Three-Reference Quantum-Circuit Testbed for the Lattice Schwinger Model*
(Zenodo, 2026; [10.5281/zenodo.21854391](https://doi.org/10.5281/zenodo.21854391); newer version
[10.5281/zenodo.22703884](https://doi.org/10.5281/zenodo.22703884)).

## What the project does

The lattice Schwinger model (1+1D QED, Kogut–Susskind staggered fermions, open boundaries, gauge field
eliminated through Gauss's law) is written as a spin chain and simulated with three mutually independent
references:

| Reference | Method | Role |
|---|---|---|
| 1 | Exact diagonalization / exact evolution | ground truth |
| 2 | Classical second-order (Strang) Trotter evolution | isolates the physical Trotter error |
| 3 | Gate-level Qiskit circuit (RXX, RYY, RZ, RZZ) | checks the gate decomposition |

Reference 3 agrees with Reference 2 to machine precision, so the remaining difference from Reference 1 is
genuine Trotter error. The testbed is then used on three questions from the literature:

1. **Commutator-scaling bound** (Childs et al., PRX 11, 011020, 2021): the bound overestimates the
   measured operator-norm error by a constant factor 1.74.
2. **Observable gap factor** G = 2‖O‖‖Δψ‖ / |Δ⟨O⟩|: independent of Δt (2.164) and growing with system
   size (1.80 → 2.76 for N = 4…10).
3. **Term ordering** (Tranter et al., Entropy 21, 1218, 2019): all four of their trends transfer to this
   gauge theory.

## Quick start

Python 3.11 is used for testing (GitHub Actions). From the top folder of the repository:

```bash
pip install -r requirements.txt
python code/run_all.py            # full run, about 40 s
python code/run_all.py --quick    # skips N = 10 and some parameter sweeps, about 20 s
```

Every script prints a `[PASS]` / `[FAIL]` line for each number it compares with the preprint, ends with
`SUCCESS: ...`, and exits with status 1 if any check fails. The full run makes **153 checks**.

## Scripts

All scripts run from the top folder, e.g. `python code/schwinger_ed.py`.

| Script | Preprint section | What it reproduces |
|---|---|---|
| `code/schwinger_core.py` | Eqs. (2)–(3), (6)–(9) | Hamiltonian, Pauli form, Trotter step, helper functions (library) |
| `code/schwinger_circuit_lib.py` | Eq. (12) | Qiskit circuit for one Trotter step (library) |
| `code/schwinger_ed.py` | Sec. 2.2, Table 1, Fig. 1 | Reference 1 checks; equivalence with Nguyen et al. (PRX Quantum 3, 020324) |
| `code/schwinger_pauli.py` | Sec. 2.3 | Pauli decomposition, diagonal H₂, Gauss-law ZZ coefficients |
| `code/schwinger_trotter.py` | Sec. 4.2, Table 2 | exponents and prefactors of the three error measures |
| `code/schwinger_commutator.py` | Sec. 4.3, Eq. (10), Result 1 | nested commutators, bound ratio 1.74, prefactors 1/24 and 1/12 |
| `code/schwinger_proper_gap.py` | Sec. 5, Table 3, Results 2–4 | gap factor G, its invariance and Δt independence |
| `code/schwinger_circuit.py` | Sec. 6.1, 6.3, Table 5 | circuit = classical propagator; gate counts and depth |
| `code/schwinger_hierarchy.py` | Sec. 6.2, Table 4, Result 5 | the complete three-reference hierarchy |
| `code/schwinger_ordering.py` | Sec. 7, Table 6, Results 6–7 | ordering strategies and the Spearman test |
| `code/make_figures.py` | Figs. 1–5 | regenerates the figures from `results/` |
| `code/run_all.py` | — | runs everything and prints a summary |

Numerical results are written to `results/` (JSON and CSV); figures to `figures/`.

## Key numbers: preprint vs this code

Default parameters: N = 6, m = 0.5, g = 1, a = 1, t = 1, initial state |1010…⟩.

| Quantity | Preprint | This code |
|---|---|---|
| ‖[H₁,H₂]‖₂ for N = 4, 6, 8, 10 | 1.609, 2.551, 3.813, 5.376 | 1.6093, 2.5507, 3.8132, 5.3759 |
| Free-field limit error, N = 4, 6, 8 | < 10⁻¹⁴ | 2.2×10⁻¹⁶, 0, 4.0×10⁻¹⁵ |
| Chiral condensate at m = 20 | −0.4997 | −0.49975 |
| Fitted exponents (operator, state, observable) | 2.000, 2.000, 2.000 | 2.0000, 2.0000, 2.0000 |
| Prefactors (operator, state, observable) | 0.389, 0.197, 1.09 | 0.3888, 0.1973, 1.094 |
| ‖[H₁,[H₁,H₂]]‖, ‖[H₂,[H₁,H₂]]‖ | 3.750, 6.254 | 3.7500, 6.2541 |
| Bound coefficient, Eq. (10) | 0.6774 | 0.67743 |
| **Bound / measured error** | **1.74** | **1.743** |
| Extracted prefactors c₁, c₂ | 1/24, 1/12 | 0.0416665, 0.0833334 |
| Normalisation artifact ratio, N = 4…10 | 8, 12, 16, 20 | 8, 12, 16, 20 |
| **G (chiral condensate), Δt ∈ [1/512, 1/16]** | **2.164** | **2.1638–2.1640** |
| G for N = 4, 6, 8, 10 | 1.80, 2.16, 2.48, 2.76 | 1.797, 2.164, 2.478, 2.757 |
| G range over Table 3 | 1.8–8.1 | 1.80–8.10 |
| G for single-site Z₀ | 3.31 | 3.315 |
| Ref 2 vs Ref 1, r = 4, 8, 16, 32, 64 | 2.58e-2, 6.40e-3, 1.60e-3, 3.99e-4, 9.97e-5 | 2.58e-2, 6.40e-3, 1.60e-3, 3.99e-4, 9.97e-5 |
| Ref 3 vs Ref 2 (machine precision), r = 64 | 7.6e-14 | 9.9e-14 |
| Fitted Trotter exponent (Ref 2 vs Ref 1) | 2.003 | 2.0033 |
| Gate counts per step, N = 4 … 12 (total) | 18, 35, 56, 81, 110 | 18, 35, 56, 81, 110 |
| Circuit depth, N = 6, r = 32 | 513 | 513 |
| Ordering errors: magnitude, equalise, deplete, lexicographic | 3.38e-4, 4.35e-4, 5.36e-4, 1.44e-3 | 3.380e-4, 4.347e-4, 5.360e-4, 1.440e-3 |
| Percentiles among 200 random orderings | 0.5, 2.5, 11.0, 100 | 0.5, 2.5, 11.0, 100 |
| Random-ordering span | 3.23e-4 – 1.42e-3 (×4.40) | 3.231e-4 – 1.422e-3 (×4.40) |
| Spearman ρ (error-operator norm vs true error) | −0.081 (p = 0.62) | −0.0812 (p = 0.618) |

## Figures

| | |
|---|---|
| ![Figure 1](figures/fig1_reference1_exact_diagonalization.png) | ![Figure 2](figures/fig2_reference2_trotter_and_gap.png) |
| **Fig. 1** Reference 1: exact diagonalization | **Fig. 2** Reference 2 and the gap factor G |
| ![Figure 3](figures/fig3_circuit_resources.png) | ![Figure 4](figures/fig4_three_reference_hierarchy.png) |
| **Fig. 3** Circuit resources | **Fig. 4** The three-reference hierarchy |
| ![Circuit](figures/fig3a_circuit_one_step_N4.png) | ![Figure 5](figures/fig5_ordering_strategies.png) |
| **Fig. 3A** One Trotter step, N = 4 | **Fig. 5** Ordering strategies |

## Conventions and notes on reproduction

* **Qubit order.** Site n is qubit n, little-endian as in Qiskit. The vacuum |1010…⟩ has site 0 in |1⟩.
* **Global phase.** Identity terms are dropped from the circuit, so circuit and classical propagators are
  compared up to a global phase.
* **Machine-precision values.** Distances at the 10⁻¹⁴ level depend on the order of floating-point
  operations and differ from the preprint in the last digits (e.g. 9.9×10⁻¹⁴ instead of 7.6×10⁻¹⁴ at
  r = 64). The checks require these to stay below 10⁻¹².
* **Bound ratio.** The ratio 1.743 is constant to three significant figures for Δt from 1/32 down to 1/8192,
  where the error itself spans about five decades (`results/bound_vs_measured_N6.csv`).
* **Ordering study.** Magnitude ordering uses a stable sort of the term list; the lexicographic order sorts
  Pauli strings written site 0 first; commuting groups come from greedy colouring (rustworkx, *Degree*
  strategy) of the anticommutation graph. Random orderings use `numpy.random.default_rng(7)`: the first
  200 for Table 6 and the next 40 for the Spearman test, where the error-operator norm is the
  triangle-inequality norm of the second-order error operator of Tranter et al.
* **Figures** are regenerated from the computed data; they show the same quantities as the preprint
  figures, with a lighter layout.
* All results are noiseless (no hardware noise model).

## Citation

If you use this code, please cite the preprint (see also `CITATION.cff`):

```bibtex
@misc{sahraei2026schwinger,
  author       = {Sahraei, Moein},
  title        = {A Validated Three-Reference Quantum-Circuit Testbed for the Lattice Schwinger Model},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21854391},
  note         = {Department of Physics, Razi University, Kermanshah, Iran}
}
```

## References

1. A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, S. Zhu, *Theory of Trotter Error with Commutator Scaling*,
   Phys. Rev. X **11**, 011020 (2021).
2. J. Kogut, L. Susskind, Phys. Rev. D **11**, 395 (1975).
3. A. Tranter, P. J. Love, F. Mintert, N. Wiebe, P. V. Coveney, *Ordering of Trotterization*,
   Entropy **21**, 1218 (2019).
4. N. H. Nguyen et al., *Digital Quantum Simulation of the Schwinger Model and Symmetry Protection with
   Trapped Ions*, PRX Quantum **3**, 020324 (2022).

## License

MIT © 2026 Moein Sahraei — see [LICENSE](LICENSE).
