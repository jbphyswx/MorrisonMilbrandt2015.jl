# Integrator

[`MM2015`](../schemes/t_updating.md) integrates the [parcel model](../theory/parcel_model.md)

```math
u' = F(t, u), \qquad u = (q_l,\ q_i,\ T),
```

where ``F`` is smooth between events because the set of active phases is fixed there. ``F``
depends on time explicitly through ``p(t)`` and ``q_t(t)``.

## Exponential Rosenbrock method

The method is `exprb43` of Hochbruck, Ostermann & Schweitzer (2009), *SIAM J. Numer. Anal.*
**47**, 786–803, doi:10.1137/080717717: fourth order, three stages, with an embedded third-order
solution for step-size control. Each step linearizes ``F`` at ``(t_n, u_n)``:

```math
J_n = \frac{\partial F}{\partial u}(t_n, u_n), \qquad
v_n = \frac{\partial F}{\partial t}(t_n, u_n), \qquad
g_n(t, u) = F(t, u) - J_n u - v_n t,
```

and uses the entire functions

```math
\varphi_0(z) = e^z, \qquad
\varphi_k(z) = \int_0^1 e^{(1-\sigma)z}\,\frac{\sigma^{k-1}}{(k-1)!}\,d\sigma, \qquad
\varphi_k(z) = \frac{\varphi_{k-1}(z) - \varphi_{k-1}(0)}{z}.
```

In the non-autonomous form (the paper's eqs. 6.5a and 6.5b) with step ``h``,

```math
U_{ni} = u_n + c_i h\,\varphi_1(c_i h J_n)\,F(t_n, u_n) + c_i^2 h^2\,\varphi_2(c_i h J_n)\,v_n
  + h\sum_{j=2}^{i-1} a_{ij}(h J_n)\,D_{nj},
```

```math
u_{n+1} = u_n + h\,\varphi_1(h J_n)\,F(t_n, u_n) + h^2\,\varphi_2(h J_n)\,v_n
  + h\sum_{i=2}^{3} b_i(h J_n)\,D_{ni},
\qquad D_{nj} = g_n(t_n + c_j h, U_{nj}) - g_n(t_n, u_n).
```

The `exprb43` coefficients are ``c = (0, \tfrac12, 1)``, ``a_{32} = \varphi_1``, and

| | ``b_2`` | ``b_3`` |
|:--|:--|:--|
| fourth order | ``16\varphi_3 - 48\varphi_4`` | ``-2\varphi_3 + 12\varphi_4`` |
| embedded third order | ``16\varphi_3`` | ``-2\varphi_3`` |

The linear part ``J_n u`` is integrated exactly by the matrix functions, and it carries the
relaxation rates ``1/\tau_l`` and ``\alpha/(\tau_i\Gamma_i)``. A relaxation time far shorter than
the step therefore places no limit on the step.

## Jacobian

The entries of ``J_n`` follow by differentiating the parcel model. With the liquid rate
``S_l = \delta/(\tau_l\Gamma_l)``,

```math
\frac{\partial S_l}{\partial q_l} = -\frac{1}{\tau_l\Gamma_l} - \frac{S_l}{\Gamma_l}\frac{\partial\Gamma_l}{\partial q_l},
\qquad
\frac{\partial S_l}{\partial T} = -\frac{\partial_T q_{sl}}{\tau_l\Gamma_l} - \frac{S_l}{\Gamma_l}\,\partial_T\Gamma_l,
\qquad
\partial_T\Gamma_l = \frac{1}{c_{pm}}\left(\frac{dL_v}{dT}\,\partial_T q_{sl} + L_v\,\partial^2_T q_{sl}\right),
```

and ``\partial\Gamma_l/\partial q_l = -(\Gamma_l - 1)\,(c_l - c_{pv})/c_{pm}``, because
``c_{pm}`` depends on the condensate. The ice entries are analogous. The temperature row
differentiates ``\dot p/(\rho c_{pm}) + \dot T_\mathrm{ext} + (L_v S_l + L_s S_i)/c_{pm}``,
including ``\partial\rho/\partial T = -\rho/T`` and ``\partial\rho/\partial q_l = \rho R_v/R_m``.
The second derivative of the saturation humidity is

```math
\partial^2_T q_{sl} = \partial_T q_{sl}\left[\beta + \frac{e_l'}{p - e_l} + \frac{L_v'}{L_v} - \frac{2}{T}\right],
\qquad \beta = \frac{p}{p - e_l}\frac{L_v}{R_v T^2}, \qquad e_l' = e_l\,\frac{L_v}{R_v T^2}.
```

## Spectrum

The eigenvalues of the ``3\times3`` Jacobian on four test states (Thermodynamics.jl 1.3) show one
fast mode near ``-1/\tau`` and two slow modes. With both phases active, one slow mode is slightly
positive.

| state | phases | eigenvalues [s⁻¹] |
|:--|:--|:--|
| warm updraft, 285 K, ``\tau_l = 5`` s | liquid | −0.200, −2.2e−5, 0 |
| mixed phase, 261 K, ``\tau_l = 8``, ``\tau_i = 60`` s | liquid + ice | −0.141, −3.4e−5, +3.3e−5 |
| above freezing, 275 K, ``\tau_l = 10``, ``\tau_i = 20`` s | liquid + ice | −0.150, −1.7e−4, +2.8e−5 |
| ice only, 225 K, ``\tau_i = 600`` s | ice | −1.7e−3, −1.7e−8, 0 |

The zero eigenvalue belongs to the inactive phase. The matrix functions must therefore be accurate
for arguments near zero, of both signs, and far into the left half-plane.

## Matrix functions

The Jacobian mixes kg kg⁻¹ and K. Its temperature row, with entries of order
``L/(c_{pm}\tau)``, dominates ``\lVert hJ_n\rVert_1``, which on the test states exceeds the
magnitude of the eigenvalues by a factor of about ``10^3``. The package first forms
``B = D^{-1}(hJ_n/2)D``, with ``D`` a diagonal of powers of 2 that brings the off-diagonal norms of
each row and column within a factor 4 of each other (Parlett & Reinsch 1969, *Numer. Math.* **13**,
293–304), and keeps ``D = I`` when ``B`` does not have the smaller 1-norm. The similarity is exact,
and ``\varphi_k(D B D^{-1}) = D\,\varphi_k(B)\,D^{-1}``. ``D`` depends only on the ratios of the
entries of ``J_n``, so it is computed once per linearization and serves every step length tried
from it; the iteration for each linearization starts from the ``D`` of the previous one.

It then evaluates the Taylor polynomial of ``\varphi_4`` at ``Z = B/2^s``, with ``s`` the smallest
integer that gives ``\lVert Z\rVert_1 \le 4``, to the least degree ``m`` with
``\lVert Z\rVert_1^m\,4!/(m+4)! \le \varepsilon``, read from a table of thresholds in `Float32` and
`Float64`. The polynomial is evaluated in blocks of four powers (Paterson & Stockmeyer 1973,
*SIAM J. Comput.* **2**, 60–66). The package then obtains ``\varphi_3`` to ``\varphi_0`` from
``\varphi_{k-1}(Z) = Z\varphi_k(Z) + I/(k-1)!``, and doubles the argument ``s`` times with
(Skaflestad & Wright 2009, *Appl. Numer. Math.* **59**, 783–799)

```math
\varphi_k(2Z) = \frac{1}{2^k}\left[\varphi_0(Z)\,\varphi_k(Z)
  + \sum_{j=1}^{k} \frac{\varphi_j(Z)}{(k-j)!}\right].
```

The functions at ``hJ_n/2`` serve the second stage, and one more doubling gives those at
``hJ_n``. For the mixed-phase test state with ``h = 60`` s, balancing lowers the 1-norm of
``hJ_n`` from ``1.6\times10^4`` to 11 and the doublings from 12 to 2.

Against 512-bit references, the error of each of ``\varphi_1`` to ``\varphi_4``, relative to its
largest entry, stays below ``2.7\,(\kappa_k + \varepsilon)``, where the componentwise condition
number ``\kappa_k`` is the summed change of ``\varphi_k``, relative to its largest entry, under a
relative change ``\varepsilon`` of each entry of the argument. The matrices tested include
Jacobians of test states, matrices with a zero row, a complex pair, or a zero eigenvalue, and
norms up to ``7\times10^5``.

## Step-size control

With ``\hat u_{n+1}`` the embedded solution, ``u_{n+1} - \hat u_{n+1} = h\,\varphi_4(hJ_n)\,(12D_{n3} - 48D_{n2})``,
and the error of a step is

```math
\mathrm{err} = \max_c \frac{\lvert u_{n+1,c} - \hat u_{n+1,c}\rvert}
  {\mathrm{atol}_c + \mathrm{rtol}\,\max(\lvert u_{n,c}\rvert, \lvert u_{n+1,c}\rvert)},
```

with ``\mathrm{atol}_c`` equal to `atol_q` for ``q_l`` and ``q_i`` and `atol_T` for ``T``. A
step is accepted when ``\mathrm{err} \le 1``. The next step size is
``h \min\!\left(5,\ \max\!\left(0.2,\ 0.9\,\mathrm{err}^{-1/4}\right)\right)``. The first
trial step covers the remaining interval. An `rtol` below ``4\varepsilon`` of the problem's
floating-point type is an error.

The state is accumulated with compensated summation: each addition ``u_n + \Delta u`` is split
exactly into a rounded sum and its rounding error (Knuth's TwoSum), and the error is added to the
next increment. Temperature increments smaller than half a unit in the last place of ``T``
(``1.5\times10^{-5}`` K at 280 K in `Float32`) therefore accumulate.

## Events

The event functions follow the active phases: ``-q_l`` for active liquid and ``\delta`` for
inactive liquid; ``-q_i`` for active ice; ``\delta_i`` for inactive ice below ``T_\mathrm{tr}`` and
for active ice above it, ``-\delta_i`` for inactive ice with mass above it; and the signed
distance past ``T_\mathrm{tr}``. An event happens where one of them becomes positive.

After an accepted step, the event functions positive at its end have fired. The event time
``s^\star \in (0, h]`` is a zero of ``s \mapsto \max_k g_k(\Phi_s(u_n))/r_k`` over the fired event
functions, where ``\Phi_s`` is the `exprb43` step of length ``s`` from ``u_n`` with the same
``J_n``, ``v_n``, and ``D``, and ``r_k`` is the rounding of ``g_k`` at the end of the step. A
bracketing root finder narrows the sign change until its width is below ``16\varepsilon\,(t_n + h)``
or its positive end is within ``r_k`` of zero, and returns the positive end; across a run of equal
zeros it doubles its step. When an event function that did not fire is positive at ``s^\star``, or
its cubic Hermite interpolant on ``[0, s^\star]`` rises above zero, the step is halved.

At ``s^\star``, each event function that fired must not exceed the sum of its absolute tolerance,
its rounding resolution, and twice its rate times the bracket width, or the step raises an error.
The rounding resolution covers the rounding of the state at both ends of the step and the response
of the rates, over the step, to the rounding of ``T`` and of the humidities; ``r_k`` is its part at
the end of the step. The step is then taken to ``s^\star``, an exhausted phase is set to zero mass,
and the active phases are re-evaluated.

When no event function is positive at the end of a step, the step is halved if the cubic Hermite
interpolant of an event function through its values and time derivatives at the two ends rises
above zero, so a turning point cannot hide two crossings inside one step.

The root finder is the keyword `root_finder`: `BrentRootFinder()` (Brent–Dekker iteration with a
bisection safeguard) by default, or `RootSolversRootFinder(method)` with a bracketing method of
RootSolvers.jl after `using RootSolvers`. Both stop on the same width and value criteria. From an
exact zero, where RootSolvers.jl stops, the bracket is narrowed by the same doubling steps and then
by bisection.

The number of steps is bounded by the keyword `max_steps`. Exceeding it raises an error.
