# Tangent

The propagation loop: one kernel for the CPU and the GPU, deviation matrices, reducers applied at handoff steps. A reducer may request a per-orbit `UInt32` workspace (`reducer_workspace`, `reducer_bytes`, the six-argument `reducer_step` and three-argument `reducer_final`; the `workspace` keyword of `propagate_orbit`), carried with its state across relaunches.

```@autodocs
Modules = [Chaoslab.Tangent]
Private = false
```
