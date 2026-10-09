"""Event that ends a segment."""
@enum EventKind::UInt8 begin
    EndOfStep
    LiquidSaturation
    IceSaturation
    LiquidExhausted
    IceExhausted
    Equilibrium
    TripleCrossing
end

"""Phases that exchange mass with the vapor during a segment."""
struct ActivePhases
    liquid::Bool
    ice::Bool
end

"""
    active_phases(k, δ, δ_i, x_l, x_i) -> ActivePhases

Liquid is active when `x_l > 0` or `δ > 0`. Below the triple point, ice is active when `x_i > 0` or
`δ_i > 0`; above it, when `x_i > 0` and `δ_i < 0`. On its saturation boundary, a phase with zero mass
is active when its supersaturation rises with that phase inactive, and ice with mass above the triple
point is active when `δ_i` falls.
"""
@inline function active_phases(k::Coefficients, δ, δ_i, x_l, x_i)
    ice_interior = k.below_triple ? (x_i > 0 || δ_i > 0) : (x_i > 0 && δ_i < 0)
    liquid = x_l > 0 || δ > 0 ||
             (iszero(x_l) && iszero(δ) && supersaturation_tendency(k, ActivePhases(false, ice_interior), δ, δ_i) > 0)
    ice = if k.below_triple
        ice_interior || (iszero(x_i) && iszero(δ_i) && supersaturation_tendency(k, ActivePhases(liquid, false), δ, δ_i) > 0)
    else
        ice_interior || (x_i > 0 && iszero(δ_i) && supersaturation_tendency(k, ActivePhases(liquid, false), δ, δ_i) < 0)
    end
    return ActivePhases(liquid, ice)
end

"""Largest number of events in one step of a frozen-coefficient scheme."""
max_events(::MM2015FixedT) = 8
max_events(::MM2015PiecewiseLinear) = 13

"""
Whether a parcel without condensate stays below saturation over liquid, and below the triple point
over ice, for the time `t`; with no phase active, both supersaturations change at the rate `A_l`.
"""
@inline function stays_subsaturated(k::Coefficients, δ, δ_i, t)
    liquid = δ < 0 && δ + k.A_l * t < 0
    return k.below_triple ? liquid && δ_i < 0 && δ_i + k.A_l * t < 0 : liquid
end

"""
    evolve(scheme, k, Δt, recorder) -> (Δx_l, Δx_i)

Liquid and ice increments over `[0, Δt]` of the frozen-coefficient dynamics of `k`, advanced by
`scheme` from each event to the next, with the state rounded by the thresholds of `scheme` at the start
of each segment. Each segment is passed to `recorder`.
"""
@inline function evolve(scheme::FrozenCoefficientScheme, k::Coefficients{FT}, Δt::FT, recorder) where {FT}
    (; x_min, δ_min, δ_i_min) = scheme.thresholds
    x_min, δ_min, δ_i_min = FT(x_min), FT(δ_min), FT(δ_i_min)
    δ, δ_i, x_l, x_i = k.δ, k.δ_i, k.x_l, k.x_i
    Δx_l = Δx_i = zero(FT)
    t = zero(FT)
    for _ in 0:max_events(scheme)
        if zero(FT) < x_l < x_min
            δ, δ_i, Δx_l, x_l = δ + x_l, δ_i + x_l, Δx_l - x_l, zero(FT)
        end
        if zero(FT) < x_i < x_min
            δ, δ_i, Δx_i, x_i = δ + x_i, δ_i + x_i, Δx_i - x_i, zero(FT)
        end
        abs(δ) < δ_min && (δ = zero(FT))
        abs(δ_i) < δ_i_min && (δ_i = zero(FT))
        remaining = Δt - t
        if iszero(x_l) && iszero(x_i) && stays_subsaturated(k, δ, δ_i, remaining)
            record!(recorder, t, remaining, ActivePhases(false, false), x_l, x_i, k.T, δ, δ_i, EndOfStep)
            return Δx_l, Δx_i
        end
        a = active_phases(k, δ, δ_i, x_l, x_i)
        s, kind, δ_next, δ_i_next, dx_l, dx_i = segment(scheme, k, a, δ, δ_i, x_l, x_i, remaining)
        kind == LiquidSaturation && (δ_next = zero(FT))
        kind == IceSaturation && (δ_i_next = zero(FT))
        kind == Equilibrium && ((δ_next, δ_i_next) = equilibrium_supersaturations(k, a))
        kind == LiquidExhausted && (dx_l = -x_l)
        kind == IceExhausted && (dx_i = -x_i)
        record!(recorder, t, s, a, x_l, x_i, k.T, δ, δ_i, kind)
        δ, δ_i = δ_next, δ_i_next
        x_l = kind == LiquidExhausted ? zero(FT) : x_l + dx_l
        x_i = kind == IceExhausted ? zero(FT) : x_i + dx_i
        Δx_l += dx_l
        Δx_i += dx_i
        t += s
        kind == EndOfStep && return Δx_l, Δx_i
    end
    throw_too_many_events(scheme)
end

@noinline throw_too_many_events(scheme) =
    throw(ErrorException("MorrisonMilbrandt2015: more than $(max_events(scheme)) events in one step of $(nameof(typeof(scheme)))"))
