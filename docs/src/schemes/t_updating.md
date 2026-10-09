# `MM2015`

`MM2015` solves the [parcel model](../theory/parcel_model.md) itself.
Nothing is frozen: the saturation humidities, ``\Gamma_l``, ``\Gamma_i``, ``\alpha``, ``\Delta``,
the latent heats, the heat capacity, and the density follow ``T``, ``p``, and the humidities
along the trajectory. Latent heating enters the temperature equation with ``L_v(T)`` and
``L_s(T)``.

## Method

The state ``u = (q_l, q_i, T)`` is integrated by the fourth-order exponential Rosenbrock method
described on the [Integrator](../numerics/integrator.md) page, with step-size control from its
embedded third-order solution. The relaxation times enter through the matrix exponential of the
Jacobian, so steps are limited by the change of the thermodynamic coefficients along the
trajectory, and a relaxation time much shorter than ``\Delta t`` costs no extra steps.

## Events

The active phases follow the rule in [Active phases and events](../numerics/events.md).
Liquid saturation, ice saturation, liquid exhaustion, and ice exhaustion are the events of
[`MM2015FixedT`](fixed_T.md). Because ``T`` changes, one more event exists: the crossing of the
triple point ``T = T_\mathrm{tr}``, where the ice rules change. The integrator stops at every
event: its time is the root, along the numerical solution, of the corresponding event function,
located by a bracketing root finder.

## Keywords

`MM2015{FT}(; rtol, atol_q, atol_T, max_steps, root_finder)`; `MM2015()` is `MM2015{Float64}()`.

| Keyword | Meaning | Default |
|:--|:--|:--|
| `rtol` | relative tolerance on ``q_l``, ``q_i``, ``T`` | ``\sqrt{\varepsilon}`` of `FT` |
| `atol_q` | absolute tolerance on ``q_l`` and ``q_i`` [kg kg⁻¹] | ``10^{-4}`` `rtol` |
| `atol_T` | absolute tolerance on ``T`` [K] | ``100`` `rtol` |
| `max_steps` | largest number of integrator steps; more raise an error | 10 000 |
| `root_finder` | bracketing root finder for events | `BrentRootFinder()` |

`RootSolversRootFinder(method)`, after `using RootSolvers`, locates events with a bracketing method
of RootSolvers.jl such as `RootSolvers.BrentsMethod`. An `rtol` below ``4\varepsilon`` of the
problem's floating-point type raises an error, so a `Float32` problem takes `MM2015{Float32}()`.

With the defaults, on the 17 test states (default backend, specific humidity), the mean rates
differ from a Radau IIA reference solved to a relative tolerance of ``10^{-13}`` by at most
``1.5\times10^{-8}`` of the larger rate in `Float64` (``10^{-11}`` to ``10^{-9}`` on most states)
and ``3.4\times10^{-4}`` in `Float32`. Larger tolerances take fewer steps.

## Relation to `MM2015FixedT`

Suppose a backend's saturation curves are linear in ``T`` and ``p`` with equal liquid and ice
slopes, and its latent heats, heat capacity, and density are constant. Then the frozen
coefficients are exact, and `MM2015` and `MM2015FixedT` give the same result to within the
tolerance. For real thermodynamics the two differ whenever ``T`` changes during the step.
