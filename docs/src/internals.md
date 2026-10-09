```@meta
CurrentModule = MorrisonMilbrandt2015
```

# Internals

Functions and types that the public interface calls, grouped by component. They are not exported
and can change between versions.

## Coefficients

The single conversion from an [`MM2015Problem`](@ref) in either moisture basis to the frozen
coefficients of [Appendix C](theory/appendix_c.md).

```@docs
Coefficients
coefficients
dry_fraction
specific_humidity
basis_heat_capacity
basis_heat_capacity_derivative
basis_density
basis_saturation
basis_saturation_second
total_water_derivatives
```

## Frozen-coefficient equations

```@docs
relaxation
supersaturation_tendency
equilibrium_supersaturation
frozen_evolution
condensate_increments
saturation_time
```

## Events and recorders

```@docs
active_phases
max_events
evolve
NoRecorder
SegmentRecorder
StepRecorder
```

## Frozen-coefficient schemes

```@docs
advance
next_event
linear_rates
linear_time
```

## Parcel model of `MM2015`

```@docs
parcel_state
parcel_tendencies
ParcelRHS
parcel_linearization
supersaturation_rates
below_triple
two_sum
compensated_add
scaled_error
```

## Event location in `MM2015`

```@docs
event_indicators
event_indicator_rates
indicator_resolution
event_kind
masked
EventMargin
hermite_maximum
hidden_crossing
```

## Exponential Rosenbrock step

```@docs
Vec3
Mat3
opnorm1
add_diagonal
balancing
similarity_factors
identity_scaling
PHI_TAYLOR_RADIUS
phi_taylor_terms
phi_taylor_thresholds
taylor_terms_from
taylor_coefficient
phi4_taylor
phi_double
phi_functions
exprb43_step
```
