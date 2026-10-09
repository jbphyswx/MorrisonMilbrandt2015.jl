```@meta
CurrentModule = MorrisonMilbrandt2015
```

# Performance

`benchmark/runbenchmarks.jl` times [`tendencies`](@ref) for every scheme, test state, moisture
basis, floating-point type, and thermodynamics backend, and writes the times, allocations, segment
and event counts, and floors to `benchmark/results/<commit>.json`:

```
julia --project=benchmark benchmark/runbenchmarks.jl
```

Every call is free of allocations, on every combination; the test suite asserts it.

## Floor

The floor of a call is the cost of the transcendental functions it cannot avoid, at the costs the
benchmark measures on the machine it runs on:

- [`MM2015PiecewiseLinear`](@ref): two saturation evaluations (liquid and ice) for the
  [`Coefficients`](@ref).
- [`MM2015FixedT`](@ref): the same two, plus one `expm1` per segment with an active phase, one
  `log1p` per saturation event, and one `expm1` per exhaustion event.
- [`MM2015`](@ref): six saturation evaluations per integrator step, for the states at its start and
  at its two inner stages.

A host that supplies its own thermodynamics through [`coefficients`](@ref) and calls
`tendencies(scheme, k::Coefficients, Δt)` skips the two saturation evaluations.

## What the cost scales with

- [`MM2015PiecewiseLinear`](@ref): one division and a few multiplications per segment; a step has
  at most 14 segments.
- [`MM2015FixedT`](@ref): one `expm1` per segment, one `log1p` per saturation event, and one
  depletion root per exhaustion event. A depletion root is a closed-form Lambert ``W`` seed and at
  most 8 Newton corrections ([Depletion times](../depletion.md)). Steps without exhaustion cost a
  few segments; each exhaustion adds a root.
- [`MM2015`](@ref): per accepted `exprb43` step, one linearization, the ``\varphi_k`` functions of the
  ``3\times3`` Jacobian, and three parcel states: at the two inner stages and at its end, which
  becomes the start of the next step. Locating an event takes 2 to 6 further evaluations of the step
  map on the test states ([Integrator](integrator.md)). The number of steps falls as the tolerance
  loosens.

## Measured times

![Time per call against the floor](../assets/speed.png)

Times per call (filled) and floors (open), one marker shape per scheme, from the benchmark results of
the commit the figure was made from.

## First call

The package compiles [`tendencies`](@ref) and [`trajectory`](@ref) for every scheme, moisture basis,
and floating-point type with the default backend when it is precompiled, and the Thermodynamics.jl
extension does the same for its parameter sets, so the first call in a session does not compile.

## Root finders

With `RootSolversRootFinder(RootSolvers.BrentsMethod)`, [`MM2015`](@ref) is also free of
allocations. On the test states, the Brent iteration of RootSolvers.jl takes up to twice the
evaluations of the step map of the default `BrentRootFinder()` to locate an event: 8 against 4 on
`activation_cooling`.
