# MorrisonMilbrandt2015.jl

[![CI](https://github.com/jbphyswx/MorrisonMilbrandt2015.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/jbphyswx/MorrisonMilbrandt2015.jl/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://jbphyswx.github.io/MorrisonMilbrandt2015.jl/dev/)

Condensation, evaporation, deposition, and sublimation of a homogeneous air parcel over a model time
step, after Appendix C of Morrison & Milbrandt (2015). Given the parcel's temperature, pressure,
humidities, the relaxation times of its liquid and ice, and the external forcing, the package returns
the mean phase-change rates of liquid and ice over the step, with the Wegener–Bergeron–Findeisen
exchange between them, the exhaustion of either phase, and activation from clear air.

| Scheme | Temperature | Cost per call |
|:--|:--|:--|
| `MM2015PiecewiseLinear()` | coefficients frozen at the step start; rates frozen per segment | 50–127 ns |
| `MM2015FixedT()` | coefficients frozen at the step start; the exact solution of Appendix C | 47–294 ns |
| `MM2015()` | evolves; the parcel model solved to a tolerance | 0.5–18 µs |

Costs are for one `Float64` call on the 17 test states. Every call is free of allocations.

## Example

```julia
using MorrisonMilbrandt2015

problem = MM2015Problem(
    SpecificHumidity(),
    DefaultThermodynamicsBackend(),
    MM2015State(261.0, 8.0e4, 2.07e-3, 2.0e-4, 1.0e-4),  # T [K], p [Pa], q_t, q_l, q_i [kg kg⁻¹]
    MM2015Timescales(8.0, 60.0),                        # τ_liq, τ_ice [s]
    MM2015Forcing(-5.0, 0.0, 0.0),                      # dp/dt [Pa s⁻¹], dT/dt [K s⁻¹], dq_v/dt [s⁻¹]
)
validate(problem, 300.0)
tendencies(MM2015FixedT(), problem, 300.0)  # (-6.67e-7, 1.20e-6) kg kg⁻¹ s⁻¹
tendencies(MM2015(), problem, 300.0)        # (-6.67e-7, 1.21e-6) kg kg⁻¹ s⁻¹
```

The parcel sits between ice and liquid saturation at 261 K and rises at about 0.5 m s⁻¹. Its liquid
evaporates onto the ice and is exhausted after about 76 s, so the mean liquid rate over 300 s is
`-q_l/Δt`. `trajectory` returns the same step with its segments or integrator steps recorded, and
`plot_evolution` draws it after `using CairoMakie`. After `using Thermodynamics`, Thermodynamics.jl
parameter sets serve as thermodynamics backends.

## Results

![Change over the step, between ice and liquid saturation](docs/src/assets/step_change_between_saturations.png)

Change of liquid, ice, and their sum over one step against the step, for a parcel at 261 K between
ice and liquid saturation: the liquid evaporates onto the ice and is exhausted after 26 s, after
which its change is `-q_l`. Dots are a reference solution of the parcel model; the lower panel is the
difference of each scheme from it. `plot_step_change` draws this figure after `using CairoMakie`.

![Error of each scheme against a reference solution](docs/src/assets/accuracy.png)

Left: relative error of the mean rates of each scheme against a Radau IIA solution of the parcel
model, as a function of the step, for the same state. `MM2015` stays below 10⁻⁹; `MM2015FixedT`
differs by 10⁻⁴ to 10⁻³ and `MM2015PiecewiseLinear` by up to 0.28. Right: the error of `MM2015` falls
with its relative tolerance and stays below it.

![Time per call against the floor](docs/src/assets/speed.png)

Time per call of each scheme on each test state, and the floor set by the transcendental functions
the call must evaluate.

## Reference

Morrison, H., and J. A. Milbrandt, 2015: Parameterization of cloud microphysics based on the
prediction of bulk ice particle properties. Part I: Scheme description and idealized tests.
*J. Atmos. Sci.*, **72**, 287–311, doi:[10.1175/JAS-D-14-0065.1](https://doi.org/10.1175/JAS-D-14-0065.1).
