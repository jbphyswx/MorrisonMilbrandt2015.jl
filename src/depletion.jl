# Condensate depletion time: the smallest positive root of f(t) = Ā·t + K·(1 − e^{−t/τ}) − c, or Inf.
# With û = t/τ, u = Āτ and β = τ·f′(0): g(û) = u·p + β·e − c, where e = 1 − e^{−û} and p = û − e; g″ = −K·e^{−û}.

"""
    _depletion_pe(û) -> (p, e)

`p = û − (1 − e^{−û})` and `e = 1 − e^{−û}` for `û ≥ 0`, each with relative error below `4ε`.
"""
@inline function _depletion_pe(û::FT) where {FT <: Union{Float32, Float64}}
    if û < FT(0.5)
        # Taylor series of p through û¹⁵; the first omitted term is below 7e-18 relative at û = 0.5
        p = û^2 * evalpoly(û, (inv(FT(2)), -inv(FT(6)), inv(FT(24)), -inv(FT(120)), inv(FT(720)), -inv(FT(5040)),
            inv(FT(40320)), -inv(FT(362880)), inv(FT(3628800)), -inv(FT(39916800)), inv(FT(479001600)),
            -inv(FT(6227020800)), inv(FT(87178291200)), -inv(FT(1307674368000))))
        return p, û - p
    end
    e = -expm1(-û)
    return û - e, e
end

function _depletion_pe(û::FT) where {FT <: AbstractFloat}
    if û < FT(0.5)
        # Taylor series of p with compensated summation
        term = û^2 / 2
        p = term
        compensation = zero(FT)
        n = 2
        while abs(term) > eps(FT) * p
            n += 1
            term *= -û / n
            y = term - compensation
            s = p + y
            compensation = (s - p) - y
            p = s
        end
        return p, û - p
    end
    e = -expm1(-û)
    return û - e, e
end

"""
    _mm_depl_eval(û, u, K, c, β) -> (g, g′, g″, g_floor)

`g(û) = u·p + β·e − c`, `g′(û) = β − K·e`, `g″(û) = K·(e − 1)`, and
`g_floor = 32ε·(|u·p| + |β·e| + |c|)`, below which `g` is zero to working precision.
"""
@inline function _mm_depl_eval(û::FT, u::FT, K::FT, c::FT, β::FT) where {FT}
    p, e = _depletion_pe(û)
    up = u * p
    βe = β * e
    return up + βe - c, β - K * e, K * (e - 1), FT(32) * eps(FT) * (abs(up) + abs(βe) + abs(c))
end

"""Largest number of Newton corrections applied to a closed-form depletion seed."""
const MAX_DEPLETION_CORRECTIONS = 8

@noinline throw_depletion_not_converged() =
    error("depletion root: $MAX_DEPLETION_CORRECTIONS Newton corrections did not reach the rounding floor")

"""
    _mm_depl_newton(û, lo, hi, restart, u, K, c, β) -> FT

Newton corrections of the closed-form seed `û` of the root of `g` in `(lo, hi)`, stopped at the
rounding floor, at most `MAX_DEPLETION_CORRECTIONS`. An iterate outside `(lo, hi)` is replaced by
`restart`, a point where `g·g″ > 0`, from which Newton converges monotonically to the root. The first
iterate within the floor receives one more Newton correction when the Kantorovich condition
`|g″|·|g| ≤ g′²/2` holds.
"""
@inline function _mm_depl_newton(û::FT, lo::FT, hi::FT, restart::FT, u::FT, K::FT, c::FT, β::FT) where {FT}
    lo < û < hi || (û = restart)
    for _ in 1:MAX_DEPLETION_CORRECTIONS
        g, g′, g″, g_floor = _mm_depl_eval(û, u, K, c, β)
        Δ = g / g′
        û_next = û - Δ
        if abs(g) ≤ g_floor
            return abs(g″ * Δ) ≤ abs(g′) / 2 && lo < û_next < hi ? û_next : û
        end
        û = lo < û_next < hi ? û_next : restart
    end
    throw_depletion_not_converged()
end

"""
    _depletion_quadratic_seed(K, c, β) -> FT

Smallest positive root of the quadratic model `β·û − (K/2)·û² − c` of `g` when its cubic truncation
error is below `√ε` relative, else `NaN`.
"""
@inline function _depletion_quadratic_seed(K::FT, c::FT, β::FT) where {FT}
    disc = β^2 - 2 * K * c
    disc ≥ zero(FT) || return FT(NaN)
    q = β + copysign(sqrt(disc), β)
    û_a = 2 * c / q
    û_b = q / K
    û = û_a > zero(FT) && !(zero(FT) < û_b < û_a) ? û_a : û_b
    return û > zero(FT) && abs(K) * û^2 ≤ 6 * sqrt(eps(FT)) * abs(β - K * û) ? û : FT(NaN)
end

"""
    _depletion_phi_inverse(ρ, σ) -> FT

Solution `s` of `e^{−s} − 1 + s = ρ` with `sign(s) = σ = ±1`, for `ρ > 0`.
"""
@inline function _depletion_phi_inverse(ρ::FT, σ::Int) where {FT}
    if ρ ≤ FT(0.005)
        # inversion series in h = √(2ρ); relative error below 6e-10 at ρ = 0.005
        h = sqrt(2 * ρ)
        return σ * h + h^2 / 6 + σ * h^3 / 36 + h^4 / 270 + σ * h^5 / 4320
    end
    W = σ < 0 ? fast_lambertwm1_from_ln(-1 - ρ) : fast_lambertw0(-exp(-1 - ρ))
    return 1 + ρ + W
end

# E = ln|z| beyond which e^E overflows in FT; W₀(e^E) then comes from its asymptotic series
@inline _depletion_E_max(::Type{FT}) where {FT} = min(FT(500), log(floatmax(FT)) - 8)

"""
    _mm_smallest_depletion_time(Ā, K, τ, c, β) -> FT

Smallest positive root of `f(t) = Ā·t + K·(1 − e^{−t/τ}) − c`, or `Inf` if none, with
`β = τ·f′(0)` taken from the caller's inputs (`δ_0`, or `δ_0 + Bτ` with WBF). The root is certified
by its residual: `|g| ≤ g_floor` from [`_mm_depl_eval`](@ref).
"""
function _mm_smallest_depletion_time(Ā::FT, K::FT, τ::FT, c::FT, β::FT) where {FT}
    iszero(c) && return zero(FT)
    if iszero(Ā)
        r = c / K
        return zero(FT) < r < one(FT) ? -τ * log1p(-r) : FT(Inf)
    end
    if iszero(K)
        t = c / Ā
        return t > zero(FT) ? t : FT(Inf)
    end
    u = Ā * τ
    s0 = (c - K) / u  # g(s0) = −K·e^{−s0}
    if (u > zero(FT)) != (K > zero(FT)) && !iszero(β) && (β > zero(FT)) == (K > zero(FT))
        # g turns at û_st > 0, where e^{−û_st} = −u/K = 1 − β/K and g(û_st + s) = D + u·(e^{−s} − 1 + s)
        r = β / K
        û_st = r < FT(0.5) ? -log1p(-r) : log(-K / u)
        D, _, _, D_floor = _mm_depl_eval(û_st, u, K, c, β)
        abs(D) ≤ D_floor && return τ * û_st
        (D > zero(FT)) == (u > zero(FT)) && return FT(Inf)
        ρ = -D / u
        if (c > zero(FT)) == (u > zero(FT))
            return τ * _mm_depl_newton(û_st + _depletion_phi_inverse(ρ, 1), û_st, s0, s0, u, K, c, β)
        end
        û = _depletion_quadratic_seed(K, c, β)
        isnan(û) && (û = û_st + _depletion_phi_inverse(ρ, -1))
        return τ * _mm_depl_newton(û, zero(FT), û_st, zero(FT), u, K, c, β)
    end
    (c > zero(FT)) == (u > zero(FT)) || return FT(Inf)
    converges_from_above = (K > zero(FT)) != (u > zero(FT))
    û = _depletion_quadratic_seed(K, c, β)
    if isnan(û)
        lnKu = log(abs(K / u))
        E = lnKu - s0
        û = if converges_from_above
            s0 + fast_lambertw0(-exp(E))
        elseif E > _depletion_E_max(FT)
            L1 = E
            L2 = log(E)
            W0mE =
                -L2 + L2 / L1 + L2 * (L2 - 2) / (2 * L1^2) + L2 * (2 * L2^2 - 9 * L2 + 6) / (6 * L1^3) +
                L2 * (3 * L2^3 - 22 * L2^2 + 36 * L2 - 12) / (12 * L1^4) +
                L2 * (12 * L2^4 - 125 * L2^3 + 350 * L2^2 - 300 * L2 + 60) / (60 * L1^5)
            lnKu + W0mE
        else
            s0 + fast_lambertw0(exp(E))
        end
    end
    hi = converges_from_above ? s0 : FT(Inf)
    return τ * _mm_depl_newton(û, zero(FT), hi, converges_from_above ? s0 : zero(FT), u, K, c, β)
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

"""
    get_t_out_of_q_WBF(δ_0, A_c, τ, τ_c, q_ice, Γ, q_sl, q_si) -> FT

Smallest positive time at which ice `q_ice` is exhausted when its rate includes the C7 term in
`q_sl − q_si`, or `Inf` if the ice never exhausts.
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
    return _mm_smallest_depletion_time(A_c + B, δ_0 - A_c * τ, τ, c, δ_0 + B * τ)
end
