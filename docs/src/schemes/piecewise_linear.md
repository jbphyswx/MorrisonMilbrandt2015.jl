# `MM2015PiecewiseLinear`

`MM2015PiecewiseLinear()` uses the dynamics, active phases, and events of
[`MM2015FixedT`](fixed_T.md) and advances each segment with rates frozen at the segment start.
It is a forward-Euler step of C1 from each event to the next.

## Segments

At the start ``t_s`` of a segment, with supersaturations ``\delta_s`` over liquid and
``\delta_{i,s}`` over ice, the rates are

```math
S_l = [l]\,\frac{\delta_s}{\tau_l\Gamma_l}, \qquad
S_i = [i]\,\frac{\delta_{i,s}}{\tau_i\Gamma_i},
```

and they stay constant over the segment. The vapor budget with constant rates gives both
supersaturations the same constant slope,

```math
\delta(t) = \delta_s + f(\delta_s)\,(t - t_s), \qquad
\delta_i(t) = \delta_{i,s} + f(\delta_s)\,(t - t_s), \qquad
q_l(t) = q_l(t_s) + S_l\,(t - t_s), \qquad
q_i(t) = q_i(t_s) + S_i\,(t - t_s),
```

with ``f`` from [`MM2015FixedT`](fixed_T.md). Because ``f`` is linear with slope ``-1/\tau``,
``\delta`` reaches the equilibrium ``\delta_\mathrm{eq} = A_c\tau`` at ``t_s + \tau``.

## Events

The events of [`MM2015FixedT`](fixed_T.md), with linear times, and one more:

| Event | Condition | Time from the segment start |
|:--|:--|:--|
| liquid saturation | ``\delta = 0`` | ``-\delta_s / f(\delta_s)`` when positive |
| ice saturation | ``\delta_i = 0`` | ``-\delta_{i,s} / f(\delta_s)`` when positive |
| liquid exhausted | ``q_l = 0`` | ``-q_l(t_s)/S_l`` when positive |
| ice exhausted | ``q_i = 0`` | ``-q_i(t_s)/S_i`` when positive |
| equilibrium | ``\delta = \delta_\mathrm{eq}`` | ``\tau`` |

After the equilibrium event, ``\delta`` stays at ``\delta_\mathrm{eq}`` and the rates take their
equilibrium values until another event or the end of the step. ``\delta_\mathrm{eq} = 0`` is an
ordinary state.

## Result

The mean rates are the summed condensate increments of all segments divided by ``\Delta t``.

## Properties

- The rates converge to those of [`MM2015FixedT`](fixed_T.md) at first order as ``\Delta t \to 0``.
- In a step with one active phase and no other event, the linear segment to the equilibrium and the
  exponential approach of C5 cover the same area of ``\delta - \delta_\mathrm{eq}``,
  ``(\delta_s - \delta_\mathrm{eq})\,\tau``. The increments of the two schemes then differ by
  ``(\delta_s - \delta_\mathrm{eq})\,\tau\,e^{-\Delta t/\tau}`` divided by ``\tau_c\Gamma_c``, which
  vanishes for ``\Delta t \gg \tau``.
- A phase that forms at its saturation boundary starts its segment at ``\delta_s = 0`` or
  ``\delta_{i,s} = 0``, so its rate stays zero until the equilibrium event, a time ``\tau`` later.
