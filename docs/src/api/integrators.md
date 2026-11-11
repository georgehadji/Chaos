# Integrators

Fixed-step integrators for Hamiltonian flows (catalog rows S1 to S9, S11 to S13, P1 to P3, G1 to G3, R1). Every explicit splitting or composition method is a `SplittingMethod`, a table of drift weights `c` and kick weights `d` run by one generic stepper: `SymplecticEuler`, `Leapfrog`, `McLachlan2`, `Ruth3`, `ForestRuth`, `PEFRL`, `CalvoSanzSerna4`, the composition builders `triple_jump` and `suzuki_jump`, `Yoshida`, `KahanLi`, `SofroniouSpaletta10`, and for near-integrable flows `SABA`, `SBAB` (from Gauss-Legendre and Gauss-Lobatto nodes), `ABA864`, `ABAH864` and the corrected `SABAC`, `SBABC`. Implicit methods for any Hamiltonian: `Gauss` (s-stage Gauss-Legendre collocation with the fixed-point solver of Hairer2008 and compensated sums) and `ImplicitMidpoint`. `TaoExtended` runs nonseparable flows explicitly in Tao's extended phase space (`extend_state`, `physical_state`). The stepping functions `advance`, `advance_n`, `advance_to`, `advance_checked`, the tangent map by the same splitting `tangent_advance`, the one-step map `IntegratorMap` as a Systems map, and `advance_ensemble` (one KernelAbstractions kernel over orbits). The surface-of-section and escape tools for a flow orbit: `section_map` with `locate_crossing` (upward crossings of `q1 = 0`, optional return-map Jacobian) and `escape_orbit` with `locate_escape` and `escape_sector`. `taylor_reference` has its methods in the TaylorIntegration extension.

The flows themselves (`HarmonicOscillator`, `Pendulum`, `Kepler`, `HenonHeilesFlow`, `TaoNonseparable`, `CountingFlow`) are in [Systems](systems.md).

```@autodocs
Modules = [Chaoslab.Integrators]
Private = false
```
