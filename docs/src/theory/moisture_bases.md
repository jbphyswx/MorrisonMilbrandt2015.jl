# Moisture bases

A problem states its moisture in one of two bases, and the returned rates use the same basis.

| Basis | Type | Definition |
|:--|:--|:--|
| specific humidity | `SpecificHumidity()` | ``q_x = m_x / m`` |
| dry-air mixing ratio | `DryAirMixingRatio()` | ``r_x = m_x / m_d`` |

Here ``m_x`` is the mass of water species ``x``, ``m_d`` the dry-air mass, and ``m`` the total
mass of the parcel (dry air, vapor, and condensate). With ``q_d = m_d/m = 1 - q_t``,

```math
q_x = q_d\, r_x, \qquad q_d = \frac{1}{1 + r_t}.
```

Morrison & Milbrandt (2015) write mixing ratios; the [parcel model](parcel_model.md) page writes
specific humidity.

## The parcel model in each basis

The model has the same form in both bases once every quantity is normalized by the same mass:

| Quantity | Specific humidity | Dry-air mixing ratio |
|:--|:--|:--|
| saturation over liquid | ``q_{sl} = q_d\,\varepsilon e_l/(p - e_l)`` | ``r_{sl} = \varepsilon e_l/(p - e_l)`` |
| heat capacity | ``c_{pm}`` (per unit moist-air mass) | ``c_p^{(d)} = c_{pm}/q_d`` (per unit dry-air mass) |
| density in ``\dot p/\rho`` | ``\rho`` | ``\rho_d = q_d\,\rho`` |
| total-water dependence of saturation | ``\partial_{q_t} q_{sl} = -q_{sl}/q_d`` | ``\partial_{r_t} r_{sl} = 0`` |

Because ``\partial_T q_{sl} = q_d\,\partial_T r_{sl}``, the products that form the psychrometric
factors are the same in both bases,

```math
\frac{L}{c_{pm}}\,\partial_T q_{s} = \frac{L}{c_p^{(d)}}\,\partial_T r_{s},
```

so ``\Gamma_l``, ``\Gamma_i``, ``\alpha``, and ``\tau`` are identical in both bases. The
supersaturation and the forcing scale with the basis: ``\delta^{(q)} = q_d\,\delta^{(r)}``.

## Closed parcel

Without vapor forcing, ``q_t`` and ``q_d`` stay constant, so every rate converts by the constant
factor ``q_d``,

```math
S^{(q)} = q_d\, S^{(r)},
```

and a problem stated in either basis gives the same physical result.

## Vapor forcing

A vapor forcing changes ``q_t``, and with it ``q_d``. The statement "liquid and ice change only
through phase change" then refers to the basis of the problem:

- in the mixing-ratio basis, ``r_l`` and ``r_i`` are unchanged by the forcing;
- in the specific-humidity basis, ``q_l`` and ``q_i`` are unchanged by the forcing.

The two describe parcels that differ by the change of the normalizing mass, a relative difference
of order ``F_v\,\Delta t`` in the condensate. Each problem is solved in its own basis, and the
converted rates of the two agree to that order.

The vapor forcings of the two bases relate by differentiating ``r_v = q_v/q_d``:

```math
\frac{dr_v}{dt} = \frac{1}{q_d}\frac{dq_v}{dt} + \frac{q_v}{q_d^2}\frac{dq_t}{dt}.
```
