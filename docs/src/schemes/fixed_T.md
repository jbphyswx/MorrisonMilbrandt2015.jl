# `MM2015FixedT`

`MM2015FixedT()` is Appendix C of Morrison & Milbrandt (2015): it solves the
[parcel model](../theory/parcel_model.md) with every coefficient frozen at its value at the start
of the step ([derivation](../theory/appendix_c.md)). Within the step it is exact for that
frozen-coefficient problem.

## Dynamics

The scheme tracks one supersaturation, ``\delta`` over liquid. Ice sees ``\delta + \Delta`` with
``\Delta`` fixed. For every set of active phases,

```math
\frac{d\delta}{dt} = f(\delta) = A_l - [l]\,\frac{\delta}{\tau_l}
  - [i]\,\alpha\,\frac{\delta + \Delta}{\tau_i\,\Gamma_i},
```

which is C1 with the inactive phase's term removed from C2 and C4. The same liquid-frame
equation serves liquid only, ice only, and both phases.

## Segments

Between events the set of active phases is fixed and ``f`` is linear in ``\delta``. With

```math
\tau^{-1} = \frac{[l]}{\tau_l} + [i]\,\frac{\alpha}{\tau_i\Gamma_i}, \qquad
A_c = A_l - [i]\,\frac{\Delta\,\alpha}{\tau_i\Gamma_i},
```

the segment solution is C5 for ``\delta`` and the integrated C6 and C7 for ``q_l`` and ``q_i``
(see [Frozen coefficients](../theory/appendix_c.md#Frozen-coefficients)). With neither phase
active, ``\tau = \infty`` and ``\delta`` changes linearly at the rate ``A_l``.

## Events

A segment ends at the earliest of the following times, or at the end of the step.

| Event | Condition | Time from the segment start |
|:--|:--|:--|
| liquid saturation | ``\delta = 0`` | ``\tau\ln\dfrac{\delta_0 - A_c\tau}{-A_c\tau}`` |
| ice saturation | ``\delta = -\Delta`` | ``\tau\ln\dfrac{\delta_0 - A_c\tau}{-\Delta - A_c\tau}`` |
| liquid exhausted | ``q_l(t) = 0`` | smallest positive root, closed form ([Depletion](../depletion.md)) |
| ice exhausted | ``q_i(t) = 0`` | smallest positive root, closed form ([Depletion](../depletion.md)) |

A saturation time exists when the argument of the logarithm exceeds 1; with ``\tau = \infty`` it
is ``(\delta_\mathrm{target} - \delta_0)/A_l`` when positive. The exhaustion times solve

```math
\bar A\,t + (\delta_0 - A_c\tau)\left(1 - e^{-t/\tau}\right) = -\frac{q_{c0}\,\tau_c\Gamma_c}{\tau},
\qquad
\bar A = \begin{cases} A_c & \text{liquid} \\ A_c + \Delta/\tau & \text{ice} \end{cases}
```

for the smallest positive ``t``.

At an event, the crossed quantity is set to its boundary value, the active phases are
re-evaluated by the rule in [Active phases and events](../numerics/events.md), and the next
segment starts from that state. The same page shows that a step has at most eight events.

## Result

The mean rates are the summed condensate increments of all segments divided by ``\Delta t``.

## Accuracy

The result is exact for the frozen-coefficient problem. Its difference from the parcel model
comes from the frozen coefficients alone, listed in
[What the frozen coefficients leave out](../theory/appendix_c.md#What-the-frozen-coefficients-leave-out).
