```@meta
CurrentModule = MorrisonMilbrandt2015
```

# `MM2015FixedT`

`MM2015FixedT()` is Appendix C of Morrison & Milbrandt (2015): it solves the
[parcel model](../theory/parcel_model.md) with every coefficient frozen at its value at the start
of the step ([derivation](../theory/appendix_c.md)). Within the step it is exact for that
frozen-coefficient problem.

## Dynamics

The scheme carries the supersaturations over liquid, ``\delta``, and over ice, ``\delta_i``, as two
states. With ``\Delta`` fixed, ``\delta_i = \delta + \Delta`` and both change at the rate

```math
\frac{d\delta}{dt} = \frac{d\delta_i}{dt} = f(\delta) = A_l - [l]\,\frac{\delta}{\tau_l}
  - [i]\,\alpha\,\frac{\delta_i}{\tau_i\,\Gamma_i},
```

which is C1 with the inactive phase's term removed from C2 and C4. The same equation serves liquid
only, ice only, and both phases.

## Segments

Between events the set of active phases is fixed, and both supersaturations relax at the rate

```math
\tau^{-1} = \frac{[l]}{\tau_l} + [i]\,\frac{\alpha}{\tau_i\Gamma_i}
```

toward their own equilibria:

```math
\frac{d\delta}{dt} = A_c - \frac{\delta}{\tau}, \quad A_c = A_l - [i]\,\frac{\Delta\,\alpha}{\tau_i\Gamma_i},
\qquad
\frac{d\delta_i}{dt} = A_{c,i} - \frac{\delta_i}{\tau}, \quad A_{c,i} = A_l + [l]\,\frac{\Delta}{\tau_l}.
```

The segment solution is C5 for each. The liquid increment integrates
``S_l = [l]\,\delta/(\tau_l\Gamma_l)`` and the ice increment integrates
``S_i = [i]\,\delta_i/(\tau_i\Gamma_i)``, which gives C6 and C7
([Frozen coefficients](../theory/appendix_c.md#Frozen-coefficients)). Near ice saturation
``\delta_i`` is small against ``\delta`` and ``\Delta``; carried as its own state, it keeps its
relative precision. With neither phase active, ``\tau = \infty`` and both supersaturations change
linearly at the rate ``A_l``.

## Events

A segment ends at the earliest of the following times, or at the end of the step.

| Event | Condition | Time from the segment start |
|:--|:--|:--|
| liquid saturation | ``\delta = 0`` | ``\tau\ln\dfrac{\delta_0 - A_c\tau}{-A_c\tau}`` |
| ice saturation | ``\delta_i = 0`` | ``\tau\ln\dfrac{\delta_{i0} - A_{c,i}\tau}{-A_{c,i}\tau}`` |
| liquid exhausted | ``q_l(t) = 0`` | smallest positive root, closed form ([Depletion](../depletion.md)) |
| ice exhausted | ``q_i(t) = 0`` | smallest positive root, closed form ([Depletion](../depletion.md)) |

A saturation time exists when the argument of the logarithm exceeds 1; with ``\tau = \infty`` it
is ``-\delta_0/A_l`` (``-\delta_{i0}/A_l`` for ice) when positive. With ``x`` the supersaturation
of the phase (``\delta`` for liquid, ``\delta_i`` for ice) and ``A`` its forcing (``A_c`` or
``A_{c,i}``), the exhaustion time is the smallest positive ``t`` with

```math
A\,t + (x_0 - A\tau)\left(1 - e^{-t/\tau}\right) = -\frac{q_{c0}\,\tau_c\Gamma_c}{\tau}.
```

At an event, the crossed quantity is set to its boundary value, the active phases are
re-evaluated by the rule in [Active phases and events](../numerics/events.md), and the next
segment starts from that state. The same page shows that a step has at most eight events, and
states the thresholds and the exit for clear air that apply at each segment start.

## Coefficients from the host model

The scheme reads the problem only through its [`Coefficients`](@ref). A host model that already
holds the saturation humidities, latent heats, heat capacity, and density passes them to
[`coefficients`](@ref) as [`ThermodynamicInputs`](@ref), together with ``\delta`` and
``\delta_i`` when it holds them more precisely than the difference of the vapor and the saturation
humidity, and calls `tendencies(MM2015FixedT(), k, Δt)`.

The coefficients split the change of the saturation humidity into its temperature, pressure, and
total-water parts, so the temperature derivatives are taken at constant pressure and total water,
and the pressure part follows as ``\partial_p x_{sl} = -x_{sl}/(p - e_{sl})``. A derivative of
``q_s = e_s/(\rho R_v T)`` at constant density is smaller by the factor
``(1 - R_v T/L)(p - e_s)/p``, 4 to 8 % on the test states, and with it the rates of the test states
change by up to 11 %.

## Result

The mean rates are the summed condensate increments of all segments divided by ``\Delta t``.

## Accuracy

The result is exact for the frozen-coefficient problem. Its difference from the parcel model
comes from the frozen coefficients alone, listed in
[What the frozen coefficients leave out](../theory/appendix_c.md#What-the-frozen-coefficients-leave-out).
