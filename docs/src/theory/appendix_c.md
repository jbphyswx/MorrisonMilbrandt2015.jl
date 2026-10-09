# Appendix C of Morrison & Milbrandt (2015)

Morrison, H., and J. A. Milbrandt, 2015: Parameterization of cloud microphysics based on the
prediction of bulk ice particle properties. Part I: Scheme description and idealized tests.
*J. Atmos. Sci.*, **72**, 287–311, Appendix C, section b.

## Published equations

The paper writes mixing ratios: ``q`` is the vapor mixing ratio, ``q_{sl}`` and ``q_{si}`` the
saturation mixing ratios over liquid and ice, and ``\delta = q - q_{sl}``.

```math
\frac{d\delta}{dt} = \frac{dq}{dt} - \frac{dq_{sl}}{dT}\frac{dT}{dt} = A_c - \frac{\delta}{\tau}
\tag{C1}
```

```math
\tau^{-1} = \tau_c^{-1} + \tau_r^{-1}
  + \left(1 + \frac{L_s}{c_p}\frac{dq_{sl}}{dT}\right)\frac{\tau_i^{-1}}{\Gamma_i}
\tag{C2}
```

```math
\Gamma_i = 1 + \frac{L_s}{c_p}\frac{dq_{si}}{dT}
\tag{C3}
```

```math
A_c = \left(\frac{\partial q}{\partial t}\right)_\mathrm{mix}
  - \frac{q_{sl}\,\rho g w}{p - e_s}
  - \frac{dq_{sl}}{dT}\left[\left(\frac{\partial T}{\partial t}\right)_\mathrm{rad}
  + \left(\frac{\partial T}{\partial t}\right)_\mathrm{mix} - \frac{w g}{c_p}\right]
  - \frac{q_{sl} - q_{si}}{\tau_i\,\Gamma_l}\left(1 + \frac{L_s}{c_p}\frac{dq_{sl}}{dT}\right)
\tag{C4}
```

```math
\delta(t) = A_c\tau + (\delta_{t=0} - A_c\tau)\,e^{-t/\tau}
\tag{C5}
```

```math
\mathrm{QCCON} = A_c\frac{\tau}{\tau_c\Gamma_l}
  + (\delta_{t=0} - A_c\tau)\frac{\tau}{\Delta t\,\tau_c\Gamma_l}\left(1 - e^{-\Delta t/\tau}\right)
\tag{C6}
```

```math
\mathrm{QICON} = A_c\frac{\tau}{\tau_i\Gamma_i}
  + (\delta_{t=0} - A_c\tau)\frac{\tau}{\Delta t\,\tau_i\Gamma_i}\left(1 - e^{-\Delta t/\tau}\right)
  + \frac{q_{sl} - q_{si}}{\tau_i\Gamma_i}
\tag{C7}
```

``\Gamma_l`` is C3 with ``L_s \to L_v`` and ``q_{si} \to q_{sl}``. The paper caps ``\tau`` at
``10^8`` s when no hydrometeors are present. The rain rate is C6 with ``\tau_c \to \tau_r``;
this package combines cloud and rain into one liquid category with
``\tau_l^{-1} = \tau_c^{-1} + \tau_r^{-1}``.

## Derivation from the parcel model

Differentiate ``\delta = q_v - q_{sl}(T, p, q_t)`` along the [parcel model](parcel_model.md) and
substitute its temperature equation. The terms free of the phase-change rates collect into
the external forcing ``A_l``, and the latent heating of the rates gives one factor
per phase:

```math
\frac{d\delta}{dt} = A_l - \Gamma_l\,S_l - \alpha\,S_i,
```

```math
A_l = F_v - \partial_T q_{sl}\left(\dot T_\mathrm{ext} + \frac{\dot p}{\rho\,c_{pm}}\right)
  - \partial_p q_{sl}\,\dot p - \partial_{q_t} q_{sl}\,F_v .
```

This identity holds at every instant; it involves no approximation. With ``\dot p = -\rho g w``,
the second and third terms of ``A_l`` are the second and third terms of C4. The last term exists
only in specific humidity, because ``q_{sl} = (1 - q_t)\,r_{sl}`` depends on ``q_t`` and the
saturation mixing ratio ``r_{sl}`` does not.

Inserting the rates ``S_l = [l]\,\delta/(\tau_l\Gamma_l)`` and
``S_i = [i]\,(\delta + \Delta)/(\tau_i\Gamma_i)`` gives C1 with

```math
\tau^{-1} = \frac{[l]}{\tau_l} + [i]\,\frac{\alpha}{\tau_i\,\Gamma_i}, \qquad
A_c = A_l - [i]\,\frac{\Delta\,\alpha}{\tau_i\,\Gamma_i}.
```

With both phases active this is C2 and C4. For liquid, the ``\Gamma_l`` of the rate and the
``\Gamma_l`` from latent heating cancel, leaving ``1/\tau_l``. For ice, the rate carries
``1/\Gamma_i`` and the latent heating carries ``\alpha``. C2 contains both factors, so in the
paper the ``\Gamma`` in each rate and the latent heating in ``dT/dt`` are separate terms, and
the parcel model keeps both.

## The C4 denominator

The derivation gives ``\tau_i\Gamma_i`` in the last term of C4; the paper prints
``\tau_i\Gamma_l``. The two agree only when ``\Gamma_l = \Gamma_i``. C2 and C7 both use
``\Gamma_i``, and only ``\Gamma_i`` reproduces C1. This package uses ``\Gamma_i``.

The check below compares ``d\delta/dt`` of the full parcel model, computed by a centered
difference along its trajectory, against ``A_c - \delta/\tau`` from the derived coefficients
and from the printed C4. Thermodynamics.jl 1.3 supplies the thermodynamic functions.

| state | T [K] | p [hPa] | phases | derived C4, relative error | printed C4, relative error |
|:--|--:|--:|:--|--:|--:|
| warm updraft, w ≈ 1 m s⁻¹ | 285 | 900 | liquid | −1.9e−9 | −1.9e−9 |
| mixed phase | 261 | 800 | liquid + ice | −4.1e−10 | −1.2e−2 |
| mixed phase, moistening and cooling | 261 | 800 | liquid + ice | −2.0e−10 | −1.1e−2 |
| above freezing | 275 | 850 | liquid + ice | −2.0e−9 | +1.6e−1 |
| ice only, supersaturated over ice | 225 | 250 | ice | +1.2e−9 | +1.9e−2 |

The derived form agrees to the truncation error of the difference. The printed form is off by
1 to 16 % whenever ice is active.

## Frozen coefficients

Appendix C holds ``A_c``, ``\tau``, ``\Gamma_l``, ``\Gamma_i``, ``\alpha``, and ``\Delta`` at their
values at the start of the interval. C1 is then linear with constant coefficients, and its
solution is C5. The condensate follows by integrating the rates:

```math
q_l(t) = q_{l0} + \frac{[l]}{\tau_l\Gamma_l}\left[A_c\tau\,t
  + (\delta_0 - A_c\tau)\,\tau\left(1 - e^{-t/\tau}\right)\right],
```

```math
q_i(t) = q_{i0} + \frac{[i]}{\tau_i\Gamma_i}\left[(A_c\tau + \Delta)\,t
  + (\delta_0 - A_c\tau)\,\tau\left(1 - e^{-t/\tau}\right)\right].
```

Dividing the increments by ``t = \Delta t`` gives C6 and C7. The term ``\Delta/(\tau_i\Gamma_i)``
in C7 converts the liquid-frame supersaturation to the ice supersaturation that drives
deposition. With neither phase active, ``\tau = \infty`` and ``\delta(t) = \delta_0 + A_l\,t``.

## What the frozen coefficients leave out

Freezing the coefficients linearizes ``q_{sl}(T, p)`` and ``q_{si}(T, p)`` about the start of the
interval and holds ``\Gamma``, ``L``, ``c_{pm}``, and ``\rho`` fixed. It also holds ``\Delta``
fixed, while the parcel model gives

```math
\frac{d\Delta}{dt} = (\partial_T q_{sl} - \partial_T q_{si})\,\frac{dT}{dt}
  + (\partial_p q_{sl} - \partial_p q_{si})\,\dot p
  + (\partial_{q_t} q_{sl} - \partial_{q_t} q_{si})\,F_v .
```

Keeping ``\Delta`` free turns the linearization into a coupled system for ``\delta`` and
``\delta_i``:

```math
\frac{d}{dt}\begin{pmatrix}\delta\\ \delta_i\end{pmatrix}
= \begin{pmatrix}A_l\\ A_i\end{pmatrix}
- \begin{pmatrix} [l]/\tau_l & [i]\,\alpha/(\tau_i\Gamma_i) \\
  [l]\,\alpha'/(\tau_l\Gamma_l) & [i]/\tau_i \end{pmatrix}
  \begin{pmatrix}\delta\\ \delta_i\end{pmatrix},
\qquad \alpha' = 1 + \frac{L_v}{c_{pm}}\,\partial_T q_{si},
```

where ``A_i`` is ``A_l`` written with the ice saturation quantities. With both phases active,
the determinant of the matrix is

```math
\frac{\Gamma_l\Gamma_i - \alpha\alpha'}{\tau_l\tau_i\Gamma_l\Gamma_i}, \qquad
\Gamma_l\Gamma_i - \alpha\alpha' = -\frac{(L_s - L_v)(\partial_T q_{sl} - \partial_T q_{si})}{c_{pm}} .
```

Below ``T_\mathrm{tr}``, ``\partial_T q_{sl} > \partial_T q_{si}``, so the determinant is negative
and one mode of the coupled linearization grows slowly: the fusion heat of the
liquid-to-ice conversion raises ``T``, which widens ``\Delta``. For the mixed-phase state in the
table the two decay rates are 0.142 s⁻¹ and −1.5 × 10⁻⁶ s⁻¹; for the ice-only state,
1.7 × 10⁻³ s⁻¹ and −1.5 × 10⁻¹¹ s⁻¹. Above ``T_\mathrm{tr}`` both are positive (0.150 s⁻¹ and
1.5 × 10⁻⁴ s⁻¹ for the 275 K state).

[`MM2015FixedT`](../schemes/fixed_T.md) solves the frozen-coefficient problem exactly.
[`MM2015`](../schemes/t_updating.md) solves the parcel model itself and keeps every term above.
