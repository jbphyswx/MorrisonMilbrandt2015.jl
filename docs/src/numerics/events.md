```@meta
CurrentModule = MorrisonMilbrandt2015
```

# Active phases and events

The three schemes share the rule that decides which phases are active and the list of events
that end a segment. This page states both and bounds the number of events in a step of
[`MM2015FixedT`](../schemes/fixed_T.md) and
[`MM2015PiecewiseLinear`](../schemes/piecewise_linear.md).

## Events

| Event | Condition | Ends a segment when | Schemes |
|:--|:--|:--|:--|
| liquid saturation | ``\delta = 0`` | ``q_l = 0`` | all |
| ice saturation | ``\delta_i = 0``, that is ``\delta = -\Delta`` | ``q_i = 0`` below ``T_\mathrm{tr}``; ``q_i > 0`` above ``T_\mathrm{tr}`` | all |
| liquid exhausted | ``q_l = 0`` while liquid evaporates | liquid is active | all |
| ice exhausted | ``q_i = 0`` while ice sublimates | ice is active | all |
| equilibrium | ``\delta = A_c\tau`` | a phase is active | `MM2015PiecewiseLinear` |
| triple-point crossing | ``T = T_\mathrm{tr}`` | always | `MM2015` |

A saturation crossing matters where the rate of a phase changes sign with no mass to absorb it: a
phase with zero mass forms or starts to lose what it gained, and ice with mass above
``T_\mathrm{tr}`` stops or resumes sublimating. At an event the crossed quantity is set to its
boundary value, and the active phases are re-evaluated from the new state.

## Boundary rule

A phase with zero mass that sits exactly on its saturation boundary is active when its
supersaturation increases with the phase inactive. For the frozen-coefficient dynamics

```math
f(\delta) = A_l - [l]\,\frac{\delta}{\tau_l} - [i]\,\alpha\,\frac{\delta + \Delta}{\tau_i\Gamma_i},
```

the liquid term vanishes at ``\delta = 0``, so ``f(0)`` is the same with liquid active or
inactive, and liquid forms when ``f(0) > 0``. Below ``T_\mathrm{tr}``, ice with zero mass at
``\delta = -\Delta`` forms when ``f(-\Delta) > 0``. Above ``T_\mathrm{tr}``, ice with mass at
``\delta = -\Delta`` is active when ``f(-\Delta) < 0``, that is when it is about to sublimate.
The rule also applies when a step starts on a boundary, so no segment has zero length. `MM2015`
applies the same rule with the time derivative of the full parcel model.

## Continuity and monotonicity

Hold the condensate presence (``q_l > 0`` or not, ``q_i > 0`` or not) fixed. Then:

- the liquid term switches only at ``\delta = 0``, where it is zero;
- the ice term switches only at ``\delta = -\Delta``, where it is zero, including under the
  above-freezing rule;
- each active term decreases with ``\delta``, with slope ``-1/\tau_l`` or
  ``-\alpha/(\tau_i\Gamma_i)``.

So ``f`` is continuous, piecewise linear, and non-increasing in ``\delta``. The equation
``d\delta/dt = f(\delta)`` has a unique solution, and that solution is monotone in time until the
condensate presence changes.

## Changes of direction

Presence changes in two ways.

- A phase forms from zero mass at its own boundary. Its term is zero there, so ``f`` does not
  jump.
- An active phase is exhausted. Condensate is lost only while it evaporates or sublimates, that
  is, while its supersaturation is negative. Its term then adds to ``f``, and removing it lowers
  ``f``.

``f`` therefore never jumps up. Once ``f < 0``, ``\delta`` falls toward the equilibrium of the
current dynamics, ``f`` stays non-positive along the way, and later jumps lower it further.
``\delta`` rises for a while, and then, after at most one change of direction, falls.

## Bound on the number of events

- Each of the two saturation boundaries is crossed at most once while ``\delta`` rises and once
  while it falls: at most 4 saturation events.
- A phase is exhausted at most once while ``\delta`` rises: after it re-forms at its boundary, it
  stays supersaturated for the rest of the rise. It is exhausted at most once while ``\delta``
  falls. That makes at most 4 exhaustion events for the two phases.

A step of `MM2015FixedT` has at most 8 events, so at most 9 segments.

`MM2015PiecewiseLinear` moves ``\delta`` along the same monotone path with linear pieces and stops
it at the equilibrium of the current dynamics. After an equilibrium event ``\delta`` stays
constant, so only an exhaustion can end the next segment: at most one equilibrium event more than
the exhaustion events, so at most 5. A step of `MM2015PiecewiseLinear` has at most 13 events.
Exceeding either bound is an error.

For `MM2015` the coefficients change along the trajectory, and the argument above gives no bound;
the integrator limits the number of steps.

## Thresholds

The keyword `thresholds` of each scheme takes [`Thresholds`](@ref). At the start of each segment of
`MM2015FixedT` and `MM2015PiecewiseLinear`, liquid or ice below `x_min` returns to the vapor, which
raises ``\delta`` and ``\delta_i`` by its mass, and a supersaturation below `δ_min` over liquid or
`δ_i_min` over ice in magnitude is set to zero, where the boundary rule decides the phase. Neither
is an event. `MM2015` applies the thresholds at the start of each integrator step
([Integrator](integrator.md#Thresholds)). The defaults round nothing.

## Clear air

A segment of `MM2015FixedT` or `MM2015PiecewiseLinear` that starts without condensate ends the step
when both supersaturations stay negative to the end of the step at the rate ``A_l``: no phase can
form, and the rates over the remaining time are zero. Above ``T_\mathrm{tr}`` only ``\delta``
counts, because ice does not form there.

## Triple point

``\Delta = 0`` at ``T_\mathrm{tr}``, where both saturation boundaries coincide. The parcel counts
as below the triple point when ``T < T_\mathrm{tr}``, or when ``T = T_\mathrm{tr}`` and the
temperature tendency without latent heating, ``\dot T_\mathrm{ext} + \dot p/(\rho\,c_{pm})``, is
negative. `MM2015FixedT` and `MM2015PiecewiseLinear` hold ``T`` fixed and decide once at the start
of the step. `MM2015` decides again at every triple-point crossing.
