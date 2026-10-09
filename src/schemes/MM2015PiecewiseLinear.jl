# MM2015PiecewiseLinear: forward-Euler segments of C1 with rates frozen at each segment start.

"""
    linear_segment(k, a, δ, δ_i) -> (S_l, S_i, slope, τ)

Liquid and ice rates at `(δ, δ_i)`, the common slope of `δ` and `δ_i`, and the time `τ` at which the
linear path reaches the equilibrium of the active phases `a`. At that equilibrium the slope is zero and
`τ = Inf`.
"""
@inline function linear_segment(k::Coefficients{FT}, a::ActivePhases, δ::FT, δ_i::FT) where {FT}
    S_l = a.liquid ? k.r_l * δ : zero(FT)
    S_i = a.ice ? k.r_i * δ_i : zero(FT)
    slope = supersaturation_tendency(k, a, δ, δ_i)
    inv_τ, τ, A_δ, _ = relaxation(k, a)
    iszero(inv_τ) && return S_l, S_i, slope, FT(Inf)
    return δ == A_δ * τ ? (S_l, S_i, zero(FT), FT(Inf)) : (S_l, S_i, slope, τ)
end

"""Time for `x` to reach zero at the rate `slope`, or `Inf`."""
@inline function linear_time(x::FT, slope::FT) where {FT}
    t = -x / slope
    return t > zero(FT) ? t : FT(Inf)
end

"""
    advance(::MM2015PiecewiseLinear, k, a, δ, δ_i, t) -> (δ, δ_i, Δx_l, Δx_i)

Supersaturations after time `t` and the liquid and ice increments, with the rates at `(δ, δ_i)` held fixed.
"""
@inline function advance(::MM2015PiecewiseLinear, k::Coefficients{FT}, a::ActivePhases, δ::FT, δ_i::FT, t::FT) where {FT}
    S_l, S_i, slope, _ = linear_segment(k, a, δ, δ_i)
    return δ + slope * t, δ_i + slope * t, S_l * t, S_i * t
end

"""
    segment(::MM2015PiecewiseLinear, k, a, δ, δ_i, x_l, x_i, remaining) -> (t, kind, δ, δ_i, Δx_l, Δx_i)

Duration of the segment that starts at `(δ, δ_i, x_l, x_i)` with the active phases `a` and ends at its
first event or after `remaining`, the event (`EndOfStep` for the end of the step), and the state change
over it.
"""
@inline function segment(::MM2015PiecewiseLinear, k::Coefficients{FT}, a::ActivePhases, δ::FT, δ_i::FT, x_l::FT, x_i::FT, remaining::FT) where {FT}
    S_l, S_i, slope, τ = linear_segment(k, a, δ, δ_i)
    t, kind = remaining, EndOfStep
    τ < t && ((t, kind) = (τ, Equilibrium))
    if x_l > 0
        if S_l < 0
            t_l = -x_l / S_l
            t_l < t && ((t, kind) = (t_l, LiquidExhausted))
        end
    else
        t_l = linear_time(δ, slope)
        t_l < t && ((t, kind) = (t_l, LiquidSaturation))
    end
    if a.ice && x_i > 0 && S_i < 0
        t_i = -x_i / S_i
        t_i < t && ((t, kind) = (t_i, IceExhausted))
    end
    if k.below_triple ? iszero(x_i) : x_i > 0
        t_i = linear_time(δ_i, slope)
        t_i < t && ((t, kind) = (t_i, IceSaturation))
    end
    return t, kind, δ + slope * t, δ_i + slope * t, S_l * t, S_i * t
end
