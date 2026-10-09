```@meta
CurrentModule = MorrisonMilbrandt2015
```

# API

```@docs
MorrisonMilbrandt2015
```

## Problem

```@docs
MM2015Problem
MM2015State
MM2015Timescales
MM2015Forcing
AbstractMoistureBasis
SpecificHumidity
DryAirMixingRatio
```

## Schemes

```@docs
MM2015PiecewiseLinear
MM2015FixedT
MM2015
```

## Rates

```@docs
tendencies
validate
```

## Recorded steps

```@docs
trajectory
Trajectory
Segment
StepTrajectory
StepRecord
state_at
EventKind
ActivePhases
```

## Root finders

The keyword `root_finder` of [`MM2015`](@ref) takes one of these.

```@docs
AbstractRootFinder
BrentRootFinder
RootSolversRootFinder
bracket_root
```

## Thermodynamics backends

A backend is [`DefaultThermodynamicsBackend`](@ref) or, after `using Thermodynamics`, any
`Thermodynamics.Parameters.AbstractThermodynamicsParameters` set. A new backend implements the
constants [`R_d`](@ref), [`R_v`](@ref), [`cp_d`](@ref), [`cp_v`](@ref), [`cp_l`](@ref),
[`cp_i`](@ref), [`grav`](@ref), and [`T_triple`](@ref), and the functions
[`saturation_vapor_pressure`](@ref) and [`latent_heat`](@ref). The package derives everything else
from these.

```@docs
DefaultThermodynamicsBackend
AbstractPhase
Liquid
Ice
saturation_vapor_pressure
latent_heat
R_d
R_v
cp_d
cp_v
cp_l
cp_i
grav
T_triple
saturation
PhaseSaturation
heat_capacity
cp_m
gas_constant_air
air_density
```

## Figures

These functions have methods after `using CairoMakie`.

```@docs
plot_step_change
plot_step_change!
plot_evolution
plot_evolution!
figure_theme
series_color
phase_color
series_style
```
