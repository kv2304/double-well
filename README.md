# Tunnel splitting in a double-well potential

Finite-difference solver for the time-independent Schrödinger equation in a tunable double well.
It recovers the decay constant of the wavefunction under the barrier, $\kappa$, from the collapse
of the tunnel splitting with barrier width, agreeing with the analytical
$\kappa = \sqrt{2m_e(V_0-E)}/\hbar$ to 0.001 %.

The lowest two states reduce to a two-level Hamiltonian with tunnel coupling $t$ and detuning
$\varepsilon$, the same model that describes a double-quantum-dot charge qubit.

<p align="center">
<img src="figures/avoided_crossing.svg">
</p>

Tilting the wells by $\lambda$ carries the two lowest levels through an avoided crossing of
minimum gap $2t$. The dashed curves are the two-level fit
$\Delta E = \sqrt{(d\lambda)^{2}+(2t)^{2}}$; the right panel shows each state changing well.

## Results

| | measured | independent | agreement |
|---|---|---|---|
| $\kappa$ (nm⁻¹) | 2.17469 from the slope of $\ln 2t$ against $w$ | 2.17471 from $V_0-E$ | **0.001 %** |
| $2t$ (eV) | 1.7815e-03 from the two-level fit | 1.7813e-03 read at $\lambda=0$ | 0.01 % |
| $O(h^2)$ truncation | error ratio 4.0031 on halving $h$ | 4 for a three-point stencil | 0.08 % |

Derivation, further checks and figures: [method.md](method.md).

## Run

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python levels.py
```
