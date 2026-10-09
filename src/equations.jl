# Appendix C of Morrison & Milbrandt (2015) on frozen coefficients, in the moisture basis of the problem.

"""
    relaxation(k, a) -> (inv_τ, τ, A_δ, A_δi)

Inverse relaxation time `1/τ` (C2), `τ` (`Inf` with neither phase active), and the forcings of
`dδ/dt = A_δ − δ/τ` (C1, C4 with `Γ_i`) and `dδ_i/dt = A_δi − δ_i/τ` for the active phases `a`.
"""
@inline function relaxation(k::Coefficients{FT}, a::ActivePhases) where {FT}
    inv_τ = (a.liquid ? k.k_l : zero(FT)) + (a.ice ? k.k_i : zero(FT))
    τ = a.liquid ? (a.ice ? k.τ_both : k.τ_l) : (a.ice ? k.τ_ice_only : FT(Inf))
    A_δ = k.A_l - (a.ice ? k.k_iΔ : zero(FT))
    A_δi = k.A_l + (a.liquid ? k.k_lΔ : zero(FT))
    return inv_τ, τ, A_δ, A_δi
end

"""Time derivative of `δ` and of `δ_i` (C1) at the supersaturations `δ`, `δ_i` with the active phases `a`."""
@inline function supersaturation_tendency(k::Coefficients{FT}, a::ActivePhases, δ, δ_i) where {FT}
    return k.A_l - (a.liquid ? k.k_l * δ : zero(FT)) - (a.ice ? k.k_i * δ_i : zero(FT))
end

"""Equilibrium supersaturations `(A_δ τ, A_δi τ)` of the active phases `a`, which must include one phase."""
@inline function equilibrium_supersaturations(k::Coefficients, a::ActivePhases)
    _, τ, A_δ, A_δi = relaxation(k, a)
    return A_δ * τ, A_δi * τ
end

"""
    frozen_evolution(δ_0, δ_i0, inv_τ, τ, A_δ, A_δi, t) -> (δ, δ_i, I, I_i)

Supersaturations after time `t` (C5) and their integrals `I` and `I_i` over `[0, t]`.
"""
@inline function frozen_evolution(δ_0::FT, δ_i0::FT, inv_τ::FT, τ::FT, A_δ::FT, A_δi::FT, t::FT) where {FT}
    iszero(inv_τ) && return δ_0 + A_δ * t, δ_i0 + A_δi * t, δ_0 * t + A_δ * t^2 / 2, δ_i0 * t + A_δi * t^2 / 2
    δ_eq, δ_ieq = A_δ * τ, A_δi * τ
    decay = -expm1(-inv_τ * t)
    return δ_0 + (δ_eq - δ_0) * decay, δ_i0 + (δ_ieq - δ_i0) * decay,
           δ_eq * t + (δ_0 - δ_eq) * decay * τ, δ_ieq * t + (δ_i0 - δ_ieq) * decay * τ
end

"""Liquid and ice increments (C6, C7) from the integrals `I` of `δ` and `I_i` of `δ_i`."""
@inline function condensate_increments(k::Coefficients{FT}, a::ActivePhases, I::FT, I_i::FT) where {FT}
    return a.liquid ? k.r_l * I : zero(FT), a.ice ? k.r_i * I_i : zero(FT)
end

"""
    saturation_time(x_0, inv_τ, τ, A) -> t

First time at which a supersaturation following `dx/dt = A − x/τ` from `x_0` reaches zero (the
inverse of C5), or `Inf`.
"""
@inline function saturation_time(x_0::FT, inv_τ::FT, τ::FT, A::FT) where {FT}
    if iszero(inv_τ)
        t = -x_0 / A
        return t > zero(FT) ? t : FT(Inf)
    end
    r = -x_0 / (A * τ)
    return r > zero(FT) ? log1p(r) * τ : FT(Inf)
end
