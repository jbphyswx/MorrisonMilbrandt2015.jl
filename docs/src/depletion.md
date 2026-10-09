```@meta
CurrentModule = MorrisonMilbrandt2015
```

# Condensate exhaustion

A phase holding condensate ``q_c`` under exponential relaxation is exhausted at the smallest positive
root of

```math
f(t) = \bar A\,t + K\left(1 - e^{-t/\tau}\right) - c,
\qquad K = \delta_0 - A_c\tau, \qquad c = -\frac{q_c\,\tau_c\,\Gamma}{\tau},
```

with ``\bar A = A_c`` without WBF and ``\bar A = A_c + (q_{sl} - q_{si})/\tau`` with it. In
``\hat u = t/\tau``, with ``u = \bar A\tau`` and ``\beta = \tau f'(0)``,

```math
g(\hat u) = u\,p(\hat u) + \beta\,e(\hat u) - c,
\qquad e = 1 - e^{-\hat u}, \qquad p = \hat u - e .
```

## Root structure

``g'' = -K e^{-\hat u}`` has one sign, so ``g`` has at most one turning point and the roots follow from
signs.

- If ``u`` and ``K`` have equal signs, or ``\beta`` is zero or has the sign of ``u``, ``g`` is monotone on
  ``\hat u > 0`` and a root exists exactly when ``c`` and ``u`` have equal signs.
- If ``u`` and ``K`` have opposite signs and ``\beta`` has the sign of ``K``, ``g`` turns at
  ``\hat u_\star = -\ln(1 - \beta/K) > 0``. With the depth ``D = g(\hat u_\star)``,

```math
g(\hat u_\star + s) = D + u\left(e^{-s} - 1 + s\right),
```

so the roots are ``\hat u_\star + s_\pm``, where ``e^{-s} - 1 + s = \rho = -D/u``. No root exists when
``D`` has the sign of ``u``. The first root lies before ``\hat u_\star`` when ``g(0) = -c`` has the sign
of ``u``, and after it when ``-c`` has the opposite sign. When ``|D|`` is below the rounding floor, the
two roots coincide at ``\hat u_\star`` to working precision.

## Seeds and certification

Each root starts from a closed form:

- ``s_\pm = 1 + \rho + W_{0,-1}\!\left(-e^{-1-\rho}\right)``, and its series
  ``s_\pm = \pm h + h^2/6 \pm h^3/36 + h^4/270 \pm h^5/4320`` in ``h = \sqrt{2\rho}`` near the turning
  point;
- ``\hat u = s_0 + W_0\!\left(\pm e^{E}\right)`` with ``s_0 = (c - K)/u`` and
  ``E = \ln|K/u| - s_0`` in the monotone case, with the asymptotic series of ``W_0`` where ``e^E``
  overflows;
- the smallest positive root of the quadratic Taylor model ``\beta\hat u - K\hat u^2/2 - c`` when its
  cubic truncation error is below ``\sqrt\varepsilon``.

Newton's method on ``g`` polishes the seed until
``|g| \le 32\varepsilon\,(|u\,p| + |\beta\,e| + |c|)``, with ``p`` and ``e`` evaluated to a few ulps. An
iterate that leaves the bracket ``(0, \hat u_\star)``, ``(\hat u_\star, s_0)``, ``(0, s_0)`` or
``(0, \infty)`` restarts at the bracket end where ``g\,g'' > 0``, from which Newton converges
monotonically; ``g(s_0) = -K e^{-s_0}`` makes ``s_0`` such an end. One more Newton correction is applied
when the Kantorovich condition ``|g''|\,|g| \le g'^2/2`` holds.

[`MM2015FixedT`](@ref) uses these depletion times. [`MM2015`](@ref) locates exhaustion as an event
along its numerical solution ([Integrator](numerics/integrator.md)).

## Functions

```@docs
get_t_out_of_q_no_WBF
get_t_out_of_q_WBF
_mm_smallest_depletion_time
_mm_depl_eval
_mm_depl_newton
_depletion_quadratic_seed
_depletion_phi_inverse
_depletion_pe
fast_lambertw0
fast_lambertwm1
fast_lambertwm1_from_ln
```
