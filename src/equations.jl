"""Form the C3 liquid and ice psychrometric factors without altering inputs."""
function prepare_coefficients(
    g::FT,
    L_i::FT,
    L_l::FT,
    c_p::FT,
    T_freeze::FT,
    dqsl_dT::FT,
    dqsi_dT::FT,
    e_sl::FT,
    e_si::FT,
    q_sl::FT,
    q_si::FT,
    q_liq::FT,
    q_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    dqvdt::FT,
    dTdt::FT,
) where {FT}
    Γ_l = one(FT) + (L_l / c_p) * dqsl_dT
    Γ_i = one(FT) + (L_i / c_p) * dqsi_dT
    return (;
        g,
        L_i,
        L_l,
        c_p,
        e_sl,
        e_si,
        dqsl_dT,
        dqsi_dT,
        q_sl,
        q_si,
        q_liq,
        q_ice,
        T_freeze,
        δ_0,
        δ_0i,
        Γ_l,
        Γ_i,
        dqvdt,
        dTdt,
    )
end

"""Combined liquid+ice relaxation time, Morrison & Milbrandt (2015) Eq. C2."""
τ_func(τ_liq::FT, τ_ice::FT, L_i::FT, c_p::FT, dqsl_dT::FT, Γ_i::FT) where {FT} =
    one(FT) / (one(FT) / τ_liq + (one(FT) + (L_i / c_p) * dqsl_dT) * ((one(FT) / τ_ice) / Γ_i))

"""Combined liquid+ice relaxation time; identical to the reconstructed C2."""
τ_func_combined(τ_liq::FT, τ_ice::FT, L_i::FT, c_p::FT, dqsl_dT::FT, Γ_i::FT) where {FT} =
    τ_func(τ_liq, τ_ice, L_i, c_p, dqsl_dT, Γ_i)

"""External supersaturation forcing with **no** WBF term (used by MM2015PiecewiseLinear and as A_c when WBF is off)."""
A_c_func_no_WBF(q_sl::FT, g::FT, w::FT, c_p::FT, e_sl::FT, dqsl_dT::FT, dqvdt::FT, dTdt::FT, p::FT, ρ::FT) where {FT} =
    dqvdt - (q_sl * ρ * g * w) / (p - e_sl) - dqsl_dT * (dTdt - (w * g) / c_p)

"""
Reconstructed C4 forcing in a reference saturation frame.

`τ_other`, `Γ_other`, and `L_other` belong to the other condensate phase. The
printed C4 denominator `Γ_l` is inconsistent with C1, C2, and C7; reconstruction
of the supersaturation budget requires `Γ_other` (`Γ_i` in the liquid frame).
"""
A_c_func(
    τ_other::FT,
    Γ_other::FT,
    q_sl::FT,
    q_si::FT,
    g::FT,
    w::FT,
    c_p::FT,
    e_sl::FT,
    L_i::FT,
    dqsl_dT::FT,
    dqvdt::FT,
    dTdt::FT,
    p::FT,
    ρ::FT,
) where {FT} =
    A_c_func_no_WBF(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ) -
    (q_sl - q_si) / (τ_other * Γ_other) * (1 + (L_i / c_p) * dqsl_dT)

function A_c_func_with_and_without_WBF(
    τ_other::FT,
    Γ_other::FT,
    q_sl::FT,
    q_si::FT,
    g::FT,
    w::FT,
    c_p::FT,
    e_sl::FT,
    L_i::FT,
    dqsl_dT::FT,
    dqvdt::FT,
    dTdt::FT,
    p::FT,
    ρ::FT,
) where {FT}
    A_c_no_WBF = A_c_func_no_WBF(q_sl, g, w, c_p, e_sl, dqsl_dT, dqvdt, dTdt, p, ρ)
    A_c =
        A_c_no_WBF -
        (q_sl - q_si) / (τ_other * Γ_other) * (1 + (L_i / c_p) * dqsl_dT)
    return (; A_c, A_c_no_WBF)
end

"""Time at which `δ(t) = value` under C5, or `Inf` if that value is unreachable in the future."""
function t_δ_hit_value(value::FT, δ_0::FT, A_c::FT, τ::FT) where {FT}
    logand = (δ_0 - A_c * τ) / (value - A_c * τ)
    return ((logand ≤ 1) || !isfinite(logand)) ? FT(Inf) : τ * log(logand)
end

"""Linear time for `x` to reach 0 at rate `dxdt`, or `Inf` if it never does."""
@inline function t_out_of_x(x::FT, dxdt::FT) where {FT}
    (iszero(dxdt) || iszero(x)) && return FT(Inf)
    t = x / -dxdt
    return (t < zero(FT)) ? FT(Inf) : t
end

"""C5 supersaturation (expm1 form)."""
function δ_func(δ_0::FT, A_c::FT, τ::FT, t::FT) where {FT}
    term = -expm1(-t / τ)
    return δ_0 + (A_c * τ - δ_0) * term
end

"""EP C5 written as `A_c τ + (δ_0 − A_c τ) e^{−t/τ}`."""
function δ_func_exponential(A_c::FT, τ::FT, δ_0::FT, Δt::FT) where {FT}
    term_1 = (δ_0 - A_c * τ)
    term_2 = exp(-Δt / τ)
    prod = iszero(term_2) ? zero(FT) : (term_1 * term_2)
    return A_c * τ + prod
end

const δ_func_EP = δ_func_exponential

function dδ_func_exponential(A_c::FT, τ::FT, δ_0::FT, Δt::FT) where {FT}
    term = -expm1(-Δt / τ)
    prod_1 = iszero(term) ? zero(FT) : ((A_c * τ) * term)
    prod_2 = iszero(term) ? zero(FT) : (-δ_0 * term)
    return iszero(term) ? zero(FT) : (prod_1 + prod_2)
end

"""One-species mean source over Δt (C6), τ_species = τ, no WBF.  Γ placement."""
function S_func_indiv(A_c::FT, τ::FT, δ_0::FT, Δt::FT, Γ::FT) where {FT}
    if isfinite(τ)
        term_1 = (δ_0 - A_c * τ) / Δt
        term_2 = -expm1(-Δt / τ)
        prod_ = iszero(term_2) ? zero(FT) : (term_1 * term_2)
        return (A_c + prod_) / Γ
    else
        return zero(FT)
    end
end

"""Two-species mean source for one condensate (C6), no WBF. EP Γ placement."""
function S_func_no_WBF(A_c::FT, τ::FT, τ_c::FT, δ_0::FT, Δt::FT, Γ::FT) where {FT}
    if !isfinite(τ_c)
        return FT(0)
    else
        term_1 = (δ_0 - A_c * τ) * (τ / τ_c) / Δt
        if isinf(term_1)
            term_1Δt = (δ_0 - A_c * τ) * (τ / τ_c)
            isinf(term_1Δt) && return floatmax(FT) * sign(term_1Δt)
            term_2Δt = -expm1(-Δt / τ)
            prodΔt = iszero(term_2Δt) ? zero(FT) : term_1Δt * term_2Δt
            alt = (((A_c * τ * Δt) / τ_c + prodΔt) / Γ) / Δt
            return isinf(alt) ? floatmax(FT) * sign(alt) : alt
        end
        term_2 = -expm1(-Δt / τ)
        prod = iszero(term_2) ? zero(FT) : (term_1 * term_2)
        S = A_c * τ / τ_c + prod
        return S / Γ
    end
end

"""Ice source including the C7 liquid-to-ice saturation-frame addend."""
S_func_WBF(
    A_c::FT,
    τ::FT,
    τ_ice::FT,
    δ_0::FT,
    Δt::FT,
    Γ::FT,
    q_sl::FT,
    q_si::FT,
) where {FT} = S_func_no_WBF(A_c, τ, τ_ice, δ_0, Δt, Γ) + ((q_sl - q_si) / (τ_ice * Γ))

const S_ql_func = S_func_no_WBF
const S_qi_func = S_func_WBF
const S_ql_func_indiv = S_func_indiv
const S_qi_func_indiv = S_func_indiv

"""
Instantaneous C6 rates (not the mean over `[0, t]`).

`q_c(t) = q_c0 + S̄(t) t` with `S̄` the C6 mean, so `q_c'(t) = S_inst(t)` off the
`max(q,0)` floor. Liquid frame when both phases are active:

`S_ql = δ(t)/(τ_liq Γ_l)`, `S_qi = δ(t)/(τ_ice Γ_i) + (q_sl−q_si) wbf /(τ_ice Γ_i)`,
`δ(t) = A_c τ + (δ_0 − A_c τ) e^{−t/τ}`.
"""
function c6_instantaneous_rates(
    δ0,
    q_sl,
    q_si,
    τ_liq,
    τ_ice,
    Γ_l,
    Γ_i,
    A_c,
    τ,
    la::Bool,
    ia::Bool,
    t,
)
    z = zero(δ0)
    if !isfinite(τ)
        return z, z
    elseif la && ia
        δt = δ_func_EP(A_c, τ, δ0, t)
        S_ql = δt / (τ_liq * Γ_l)
        S_qi = δt / (τ_ice * Γ_i) + ((q_sl - q_si) / (τ_ice * Γ_i))
        return S_ql, S_qi
    elseif la
        δt = δ_func_EP(A_c, τ, δ0, t)
        return δt / (τ_liq * Γ_l), z
    elseif ia
        δi0 = δ0 + (q_sl - q_si)
        δit = δ_func_EP(A_c, τ, δi0, t)
        return z, δit / (τ_ice * Γ_i)
    else
        return z, z
    end
end

"""Weighted average of two sub-step mean rates over the full `Δt`."""
@inline function resolve_S_S_addit(
    S_ql::FT,
    S_qi::FT,
    Δt_S::FT,
    S_ql_addit::FT,
    S_qi_addit::FT,
    Δt_S_addit::FT,
    Δt::FT,
) where {FT}
    return (S_ql * Δt_S) / Δt + (S_ql_addit * Δt_S_addit) / Δt,
    (S_qi * Δt_S) / Δt + (S_qi_addit * Δt_S_addit) / Δt
end

"""Argmin of a few times; a numerical 0 is bumped to `eps` so a true event is not skipped."""
function find_min_t(x1::FT, rest::Vararg{FT}) where {FT}
    xs = (x1, rest...)
    min_t = x1
    i_min_t = 1
    for (i, val) in enumerate(xs)
        if val < min_t
            min_t = val
            i_min_t = i
        end
    end
    if min_t == zero(FT)
        best = FT(typemax(FT))
        for val in xs
            if val != zero(FT) && val < best
                best = val
            end
        end
        min_t = best == FT(typemax(FT)) ? eps(FT) : min(best, eps(FT))
    end
    return min_t, i_min_t
end
