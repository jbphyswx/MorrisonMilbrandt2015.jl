# Closed-form condensate depletion time. See docs/src/depletion.md.
# f(t) = Ā·t + K·(1 − e^{−t/τ}) − c ; smallest positive root, or Inf if none.

@inline function _mm_depl_eval(û::FT, u::FT, K::FT, c::FT, β::FT) where {FT}
    if û < FT(0.05)
        p = û^2 / 2 - û^3 / 6 + û^4 / 24 - û^5 / 120
        g = β * û - K * p - c
        gp = β - K * (û - û^2 / 2 + û^3 / 6 - û^4 / 24)
        terms = abs(u * û) + abs(K) * û + abs(c)
        gp_floor = abs(β) + abs(K) * û
    else
        em = -expm1(-û)
        Kem = K * em
        g = u * û + Kem - c
        gp = β - Kem
        terms = abs(u * û) + abs(Kem) + abs(c)
        gp_floor = abs(β) + abs(Kem)
    end
    return g, gp, FT(64) * eps(FT) * terms, FT(64) * eps(FT) * gp_floor
end

@inline function _mm_depl_take(best::FT, û::FT, u::FT, K::FT, c::FT, β::FT, τ::FT) where {FT}
    (isfinite(û) && û > zero(FT)) || return best
    g, gp, _, gpf = _mm_depl_eval(û, u, K, c, β)
    (isfinite(gp) && abs(gp) > gpf) && (û -= g / gp)
    (isfinite(û) && û > zero(FT)) || return best
    g, gp, _, gpf = _mm_depl_eval(û, u, K, c, β)
    (isfinite(gp) && abs(gp) > gpf) && (û -= g / gp)
    (isfinite(û) && û > zero(FT)) || return best
    g, _, gfloor, _ = _mm_depl_eval(û, u, K, c, β)
    return abs(g) ≤ gfloor ? min(best, τ * û) : best
end

"""
    _mm_smallest_depletion_time(Ā, K, τ, c, β) -> FT

Smallest positive root of `Ā·t + K·(1 − e^{−t/τ}) − c`, or `Inf` if none.
`β = τ·f′(0)` must be supplied from the inputs, not as `u + K`.
"""
function _mm_smallest_depletion_time(Ā::FT, K::FT, τ::FT, c::FT, β::FT) where {FT}
    iszero(c) && return zero(FT)
    best = FT(Inf)
    u = Ā * τ

    if iszero(Ā)
        if !iszero(K)
            r = c / K
            (zero(FT) < r < one(FT)) && (best = -τ * log1p(-r))
        end
        return best
    end
    if iszero(K)
        t = c / Ā
        t > zero(FT) && (best = t)
        return best
    end

    disc = β^2 - FT(2) * K * c
    if disc ≥ zero(FT)
        den = β + (β ≥ zero(FT) ? sqrt(disc) : -sqrt(disc))
        if !iszero(den)
            û_qs = FT(2) * c / den
            (zero(FT) < û_qs < FT(0.05)) && (best = _mm_depl_take(best, û_qs, u, K, c, β, τ))
            best < FT(1e-3) * τ && return best
            û_ql = den / K
            (zero(FT) < û_ql < FT(0.05)) && (best = _mm_depl_take(best, û_ql, u, K, c, β, τ))
        end
    end

    s0 = (c - K) / u
    E = log(abs(K / u)) - s0
    sgn = sign(K) * sign(Ā)

    if sgn < zero(FT)
        d = -one(FT) - E
        if d > FT(1e-7)
            if E > FT(-690)
                best = _mm_depl_take(best, fast_lambertwm1_from_ln(E) + s0, u, K, c, β, τ)
                isfinite(best) || (best = _mm_depl_take(best, fast_lambertw0(-exp(E)) + s0, u, K, c, β, τ))
            else
                L1 = E
                L2 = log(-E)
                Wm1 =
                    L2 / L1 + L2 * (L2 - 2) / (2 * L1^2) + L2 * (2 * L2^2 - 9 * L2 + 6) / (6 * L1^3) +
                    L2 * (-12 + 36 * L2 - 22 * L2^2 + 3 * L2^3) / (12 * L1^4)
                best = _mm_depl_take(best, log(abs(K / u) / (-E)) + Wm1, u, K, c, β, τ)
                isfinite(best) || (best = _mm_depl_take(best, s0, u, K, c, β, τ))
            end
        elseif d > -FT(64) * eps(FT)
            û_st = log(-K / u)
            r2d = sqrt(FT(2) * max(d, zero(FT)))
            best = _mm_depl_take(best, û_st - r2d * (one(FT) - r2d / FT(6)), u, K, c, β, τ)
            isfinite(best) || (best = _mm_depl_take(best, û_st + r2d * (one(FT) + r2d / FT(6)), u, K, c, β, τ))
        end
    else
        if E > FT(500)
            L1 = E
            L2 = log(E)
            W0mE =
                -L2 + L2 / L1 + L2 * (L2 - 2) / (2 * L1^2) + L2 * (2 * L2^2 - 9 * L2 + 6) / (6 * L1^3) +
                L2 * (3 * L2^3 - 22 * L2^2 + 36 * L2 - 12) / (12 * L1^4) +
                L2 * (12 * L2^4 - 125 * L2^3 + 350 * L2^2 - 300 * L2 + 60) / (60 * L1^5)
            best = _mm_depl_take(best, log(abs(K / u)) + W0mE, u, K, c, β, τ)
        else
            best = _mm_depl_take(best, fast_lambertw0(exp(E)) + s0, u, K, c, β, τ)
        end
    end

    return best
end

"""
    get_t_out_of_q_no_WBF(δ_0, A_c, τ, τ_c, q_c, Γ) -> FT

Smallest positive time at which condensate `q_c` is exhausted under exponential (no-WBF)
relaxation, or `Inf` if the phase never exhausts.
"""
function get_t_out_of_q_no_WBF(δ_0::FT, A_c::FT, τ::FT, τ_c::FT, q_c::FT, Γ::FT) where {FT}
    c = -q_c * (τ_c * Γ / τ)
    return _mm_smallest_depletion_time(A_c, δ_0 - A_c * τ, τ, c, δ_0)
end
const get_t_out_of_q_liq = get_t_out_of_q_no_WBF
const get_t_out_of_q_ice_no_WBF = get_t_out_of_q_no_WBF

"""
    get_t_out_of_q_WBF(δ_0, A_c, τ, τ_c, q_ice, Γ, q_sl, q_si) -> FT

Ice exhaustion time including the C7 liquid-to-ice saturation-frame term.
"""
function get_t_out_of_q_WBF(
    δ_0::FT,
    A_c::FT,
    τ::FT,
    τ_c::FT,
    q_ice::FT,
    Γ::FT,
    q_sl::FT,
    q_si::FT,
) where {FT}
    B = (q_sl - q_si) / τ
    c = -q_ice * (τ_c * Γ / τ)
    return _mm_smallest_depletion_time(FT(A_c + B), FT(δ_0 - A_c * τ), τ, c, FT(δ_0 + B * τ))
end
const get_t_out_of_q_ice = get_t_out_of_q_WBF
