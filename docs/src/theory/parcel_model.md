# Parcel model

Every scheme in this package approximates one problem: the phase changes of a single
homogeneous air parcel that holds water vapor, liquid, and ice, over a time step ``\Delta t``.
This page states that problem. [Appendix C](appendix_c.md) derives Morrison & Milbrandt
(2015) Appendix C from it, and each scheme page states which approximation that scheme makes.

The equations on this page use specific humidity. [Moisture bases](moisture_bases.md) gives the
dry-air mixing-ratio form.

## Variables and forcing

| Symbol | Meaning | Unit |
|:--|:--|:--|
| ``T`` | temperature | K |
| ``p`` | pressure | Pa |
| ``q_t,\ q_l,\ q_i`` | total water, liquid, and ice specific humidity | kg kg⁻¹ |
| ``q_v = q_t - q_l - q_i`` | vapor specific humidity | kg kg⁻¹ |
| ``\tau_l,\ \tau_i`` | supersaturation relaxation times of liquid and of ice | s |
| ``\dot p`` | pressure tendency | Pa s⁻¹ |
| ``\dot T_\mathrm{ext}`` | temperature tendency from radiation and mixing | K s⁻¹ |
| ``F_v`` | vapor tendency from mixing | s⁻¹ |

Liquid combines cloud water and rain, ``\tau_l^{-1} = \tau_c^{-1} + \tau_r^{-1}`` (paper C2).
The three forcings are the external terms of the paper's C4 and are held constant over the step.
A hydrostatic parcel rising at speed ``w`` has ``\dot p = -\rho g w``. Liquid and ice change
only through condensation, evaporation, deposition, and sublimation; the vapor forcing changes
``q_t`` at the rate ``F_v``.

## Saturation

With ``\varepsilon = R_d/R_v`` and saturation vapor pressures ``e_l(T)`` over liquid and
``e_i(T)`` over ice,

```math
q_{sl} = (1 - q_t)\,\frac{\varepsilon\, e_l}{p - e_l}, \qquad
q_{si} = (1 - q_t)\,\frac{\varepsilon\, e_i}{p - e_i}.
```

The supersaturations over liquid and over ice, and their difference, are

```math
\delta = q_v - q_{sl}, \qquad \delta_i = q_v - q_{si}, \qquad
\Delta = q_{sl} - q_{si}, \qquad \delta_i = \delta + \Delta .
```

``\Delta > 0`` below the triple-point temperature ``T_\mathrm{tr}``, ``\Delta < 0`` above it, and
``\Delta = 0`` at ``T_\mathrm{tr}``, where ``e_l = e_i``.

A thermodynamics backend supplies ``e_l``, ``e_i``, the latent heats ``L_v(T)`` and
``L_s(T)``, and the heat capacities. The backend must satisfy the Clausius–Clapeyron relation
with its own latent heats and Kirchhoff's relation with its own heat capacities:

```math
\frac{d\ln e_l}{dT} = \frac{L_v(T)}{R_v T^2}, \qquad
\frac{d\ln e_i}{dT} = \frac{L_s(T)}{R_v T^2}, \qquad
\frac{dL_v}{dT} = c_{pv} - c_l, \qquad
\frac{dL_s}{dT} = c_{pv} - c_i .
```

The partial derivatives of the liquid saturation humidity follow; the ice ones are the same with
``e_i`` and ``L_s``:

```math
\partial_T q_{sl} = q_{sl}\,\frac{p}{p - e_l}\,\frac{L_v}{R_v T^2}, \qquad
\partial_p q_{sl} = -\frac{q_{sl}}{p - e_l}, \qquad
\partial_{q_t} q_{sl} = -\frac{q_{sl}}{1 - q_t}.
```

## Psychrometric factors

The paper's C3 defines the psychrometric corrections of liquid and ice, and C2 and C4 contain a
cross factor:

```math
\Gamma_l = 1 + \frac{L_v}{c_{pm}}\,\partial_T q_{sl}, \qquad
\Gamma_i = 1 + \frac{L_s}{c_{pm}}\,\partial_T q_{si}, \qquad
\alpha = 1 + \frac{L_s}{c_{pm}}\,\partial_T q_{sl},
```

with the moist isobaric heat capacity
``c_{pm} = c_{pd}(1 - q_t) + c_{pv} q_v + c_l q_l + c_i q_i``.

## Active phases

A phase exchanges mass with the vapor when it is active. The indicators ``[l]`` and ``[i]`` are 1
for an active phase and 0 for an inactive one.

- Liquid is active when ``q_l > 0`` or ``\delta > 0``.
- Below ``T_\mathrm{tr}``, ice is active when ``q_i > 0`` or ``\delta_i > 0``.
- Above ``T_\mathrm{tr}``, ice is active only when ``q_i > 0`` and ``\delta_i < 0``: existing ice
  sublimates, and no deposition occurs. Melting is a separate process.
- A phase with zero mass that sits exactly on its saturation boundary (``\delta = 0`` with
  ``q_l = 0``, or ``\delta_i = 0`` with ``q_i = 0``) is active when its supersaturation increases
  with the phase inactive. [Active phases and events](../numerics/events.md) shows that this
  rule makes the dynamics continuous across the boundary.
- At ``T = T_\mathrm{tr}`` exactly, the parcel counts as below the triple point when the
  temperature tendency without latent heating, ``\dot T_\mathrm{ext} + \dot p/(\rho\,c_{pm})``, is
  negative.

## Rates

The instantaneous condensation and deposition rates are the paper's C6 and C7 before the time
average:

```math
S_l = [l]\,\frac{\delta}{\tau_l\,\Gamma_l}, \qquad
S_i = [i]\,\frac{\delta_i}{\tau_i\,\Gamma_i}.
```

## Equations of motion

```math
\frac{dq_l}{dt} = S_l, \qquad
\frac{dq_i}{dt} = S_i, \qquad
c_{pm}\,\frac{dT}{dt} = \frac{\dot p}{\rho} + c_{pm}\,\dot T_\mathrm{ext}
  + L_v(T)\,S_l + L_s(T)\,S_i,
```

with ``p(t) = p_0 + \dot p\,t``, ``q_t(t) = q_{t0} + F_v\,t``, and the ideal-gas density
``\rho = p / (R_m T)``, ``R_m = R_d(1 - q_t) + R_v q_v``.

### Energy equation

For a closed parcel the first law is ``dh = dp/\rho``, with the moist enthalpy

```math
h = c_{pm}\,(T - T_0) + q_v\,L_{v0} - q_i\,L_{f0},
```

where ``L_{v0}`` and ``L_{f0} = L_{s0} - L_{v0}`` are the latent heats at the reference
temperature ``T_0``. Differentiating ``h`` at constant ``q_t`` and using
``dq_v = -dq_l - dq_i`` gives

```math
c_{pm}\,dT = \frac{dp}{\rho}
  + \left[L_{v0} + (c_{pv} - c_l)(T - T_0)\right] dq_l
  + \left[L_{s0} + (c_{pv} - c_i)(T - T_0)\right] dq_i
  = \frac{dp}{\rho} + L_v(T)\,dq_l + L_s(T)\,dq_i ,
```

so the latent heats in the temperature equation are the temperature-dependent ones that also
appear in ``\Gamma_l``, ``\Gamma_i``, and the saturation curves. External heating adds
``c_{pm}\dot T_\mathrm{ext}``. The vapor forcing ``F_v`` enters at the parcel temperature; any
temperature change that mixing causes is part of ``\dot T_\mathrm{ext}``, as in the paper's C4.

## Result of a step

A scheme returns the mean phase-change rates over the step,

```math
\bar S_l = \frac{q_l(\Delta t) - q_l(0)}{\Delta t}, \qquad
\bar S_i = \frac{q_i(\Delta t) - q_i(0)}{\Delta t}.
```

A step of length zero returns ``(0, 0)``.
