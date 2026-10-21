# API

```@docs
Chaoslab
```

Every name below is exported by its layer and re-exported by `Chaoslab`, so `using Chaoslab` is enough. The layers depend only on the layers above them in this list.

## Foundation

Health errors, compensated sums, compute backend resolution, precision labels and exact conversion.

```@autodocs
Modules = [Chaoslab.Foundation]
Private = false
```

## Systems

The standard map (lifted and wrapped coordinates) and linear symplectic maps, with their tangent maps and residuals.

```@autodocs
Modules = [Chaoslab.Systems]
Private = false
```

## Tangent

The propagation loop: one kernel for the CPU and the GPU, deviation matrices, reducers applied at handoff steps.

```@autodocs
Modules = [Chaoslab.Tangent]
Private = false
```

## Diagnostics

Symplectic residuals, orthogonality loss and the other reducers that check a run.

```@autodocs
Modules = [Chaoslab.Diagnostics]
Private = false
```

## Lyapunov

Lyapunov estimators as reducers: Benettin, the QR family (Gram-Schmidt variants, Householder, Givens), singular values, two-trajectory separation.

```@autodocs
Modules = [Chaoslab.Lyapunov]
Private = false
```

## Indicators

Chaos indicators as reducers: the Smaller Alignment Index SALI, its pure functions and the frame check.

```@autodocs
Modules = [Chaoslab.Indicators]
Private = false
```

## Experiments

Preregistered campaigns: orbit panels, chunked resumable runs, the HDF5 run store and manifests.

```@autodocs
Modules = [Chaoslab.Experiments]
Private = false
```

## Analysis

The analysis of a finished study run: tables, figure data and the preregistered decisions.

```@autodocs
Modules = [Chaoslab.Analysis]
Private = false
```
