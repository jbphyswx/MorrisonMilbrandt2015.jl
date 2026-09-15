# Mathematical contract

This package implements the homogeneous condensation and deposition subsystem in
Morrison and Milbrandt (2015), Appendix C. It is not the complete P3 microphysics
scheme. The piecewise-linear limiter is a separate approximation.

## Moisture coordinates

Let \(m_d\) be dry-air mass and \(m\) total moist-air mass. For any water species \(x\),

\[
r_x = \frac{m_x}{m_d},\qquad
q_x = \frac{m_x}{m},\qquad
q_d = \frac{m_d}{m}=1-q_t,\qquad
q_x=q_d r_x.
\]

For phase changes in a closed parcel, \(q_t\) and \(q_d\) are constant. Changing
between \(r\) and \(q\) is therefore only a change of coordinates:

\[
\dot q_x=q_d\dot r_x,\qquad
c_p^{(d)}=\frac{c_p^{(m)}}{q_d}.
\]

The psychrometric products are invariant:

\[
\frac{L}{c_p^{(d)}}\frac{\partial r_s}{\partial T}
=
\frac{L}{c_p^{(m)}}\frac{\partial q_s}{\partial T}.
\]

When external forcing changes total water, the constant-factor rate conversion is
not valid. Differentiating \(r_x=q_x/q_d\) gives

\[
\dot r_x =
\frac{\dot q_x}{q_d}
+\frac{q_x\dot q_t}{q_d^2},\qquad
\dot q_x =
q_d\dot r_x-\frac{q_x}{q_d}\dot q_t.
\]

For saturation vapor pressure \(e_s(T)\),

\[
r_s=\frac{\epsilon e_s}{p-e_s},\qquad
q_s=q_d r_s.
\]

Thus a specific-humidity saturation value depends on the parcel normalization.
The identity \(q_s=e_s/(\rho R_vT)\) holds when the parcel density is evaluated
with vapor partial pressure \(e_s\); away from saturation, inserting the actual
off-equilibrium parcel density defines a different saturation coordinate. This
package uses the paper's pressure-based \(r_s(T,p)\) and transforms it explicitly.
The pressure-only expression
\(\epsilon e_s/[p-(1-\epsilon)e_s]\) instead assumes a condensate-free saturated
parcel and is not the general conversion of the paper's \(r_s\).

## Appendix C in one mass basis

The following derivation uses \(x\) for either \(r\) or \(q\), provided every
specific heat, saturation value, derivative, and source uses that same mass basis.
Define

\[
\delta=x_v-x_{sl},\qquad \Delta_s=x_{sl}-x_{si},
\]

\[
\Gamma_l=1+\frac{L_v}{c_p}\partial_Tx_{sl},\qquad
\Gamma_i=1+\frac{L_s}{c_p}\partial_Tx_{si},
\]

and

\[
\alpha_{i|l}=1+\frac{L_s}{c_p}\partial_Tx_{sl}.
\]

The instantaneous liquid and ice rates are

\[
C_l=\frac{\delta}{\tau_l\Gamma_l},\qquad
C_i=\frac{\delta+\Delta_s}{\tau_i\Gamma_i}.
\]

Here \(\tau_l^{-1}=\tau_c^{-1}+\tau_r^{-1}\) when cloud and rain are represented by
one aggregated liquid category.

Let \(A_{\mathrm{ext},l}\) be the liquid-frame supersaturation tendency from
external vapor, pressure, radiation, and mixing, excluding condensation,
deposition, and their latent heating. Substitution of \(C_l\) and \(C_i\) into the
vapor and temperature budgets gives

\[
\frac{1}{\tau}
=\frac{1}{\tau_l}
+\frac{\alpha_{i|l}}{\tau_i\Gamma_i},
\]

\[
A_c=A_{\mathrm{ext},l}
-\frac{\Delta_s\alpha_{i|l}}{\tau_i\Gamma_i},
\qquad
\dot\delta=A_c-\frac{\delta}{\tau}.
\]

The solution is

\[
\delta(t)=A_c\tau+(\delta_0-A_c\tau)e^{-t/\tau}.
\]

Averaging the instantaneous rates over \([0,\Delta t]\) gives

\[
\overline C_l =
\frac{1}{\Gamma_l}\left[
\frac{A_c\tau}{\tau_l}
+\frac{(\delta_0-A_c\tau)\tau}{\Delta t\,\tau_l}
\left(1-e^{-\Delta t/\tau}\right)\right],
\]

\[
\overline C_i =
\frac{1}{\Gamma_i}\left[
\frac{A_c\tau}{\tau_i}
+\frac{(\delta_0-A_c\tau)\tau}{\Delta t\,\tau_i}
\left(1-e^{-\Delta t/\tau}\right)
+\frac{\Delta_s}{\tau_i}
\right].
\]

The printed denominator \(\tau_i\Gamma_l\) in MM2015 equation C4 is inconsistent
with C1, C2, and C7 unless \(\Gamma_l=\Gamma_i\). The reconstruction above requires
\(\tau_i\Gamma_i\). No published erratum was found; the package treats the
reconstructed identity as its mathematical contract and tests it directly. It
does not alter latent heats or saturation derivatives to force the two
psychrometric factors to agree.

## External forcing

The public temperature forcing is
\(\dot T_{\mathrm{external}}=(\partial T/\partial t)_{\mathrm{rad}}+
(\partial T/\partial t)_{\mathrm{mix}}\), as in C4. It excludes adiabatic and
microphysical latent-heating tendencies.

The pressure tendency is supplied independently. A hydrostatic parcel adapter may
form \(\dot p=-\rho g w\). For the T-updating extension, the external heat source is
mapped on each frozen-coefficient segment to

\[
\dot\theta_{li,\mathrm{external}}=
\left(\frac{\partial\theta_{li}}{\partial T}\right)_{p,x}
\dot T_{\mathrm{external}}.
\]

Pressure and latent effects then enter through the updated pressure, condensate,
and the nonequilibrium \(p\)-\(\theta_{li}\)-\(x\) inversion; they are not added
again to \(\dot T_{\mathrm{external}}\).
