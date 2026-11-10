# Systems

The standard map (lifted and wrapped coordinates) and linear symplectic maps, with their tangent maps and residuals; `PairMap`, the product of a map with itself for the RLI orbit pair (block-wise tangent). The Froeschle 4D map and the chain of N coupled standard maps (open or ring), with the banded chain tangent `ChainTangent`. Hamiltonian flows for the integrators: `HarmonicOscillator`, `Pendulum`, `Kepler`, `HenonHeilesFlow` (separable, with `drift` and `kick`, their increments and Jacobians, exact flows where they exist), `TaoNonseparable`, and `CountingFlow`, which counts force evaluations.

```@autodocs
Modules = [Chaoslab.Systems]
Private = false
```
