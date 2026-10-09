# MM2015FixedT: segments of the exact frozen-coefficient solution (C5–C7), exhaustion times from Lambert W.

"""
    advance(::MM2015FixedT, k, a, δ, δ_i, t) -> (δ, δ_i, Δx_l, Δx_i)

Supersaturations after time `t` and the liquid and ice increments, for the active phases `a`.
"""
@inline function advance(::MM2015FixedT, k::Coefficients{FT}, a::ActivePhases, δ::FT, δ_i::FT, t::FT) where {FT}
    inv_τ, τ, A_δ, A_δi = relaxation(k, a)
    δ_t, δ_i_t, I, I_i = frozen_evolution(δ, δ_i, inv_τ, τ, A_δ, A_δi, t)
    return δ_t, δ_i_t, condensate_increments(k, a, I, I_i)...
end

"""
    segment(::MM2015FixedT, k, a, δ, δ_i, x_l, x_i, remaining) -> (t, kind, δ, δ_i, Δx_l, Δx_i)

Duration of the segment that starts at `(δ, δ_i, x_l, x_i)` with the active phases `a` and ends at its
first event or after `remaining`, the event (`EndOfStep` for the end of the step), and the state change
over it.
"""
@inline function segment(::MM2015FixedT, k::Coefficients{FT}, a::ActivePhases, δ::FT, δ_i::FT, x_l::FT, x_i::FT, remaining::FT) where {FT}
    inv_τ, τ, A_δ, A_δi = relaxation(k, a)
    t, kind = remaining, EndOfStep
    if x_l > 0
        if δ < 0 || A_δ < 0
            t_l = get_t_out_of_q_no_WBF(δ, A_δ, τ, k.τ_l, x_l, k.Γ_l)
            t_l < t && ((t, kind) = (t_l, LiquidExhausted))
        end
    else
        t_l = saturation_time(δ, inv_τ, τ, A_δ)
        t_l < t && ((t, kind) = (t_l, LiquidSaturation))
    end
    if a.ice && x_i > 0 && (δ_i < 0 || A_δi < 0)
        t_i = get_t_out_of_q_no_WBF(δ_i, A_δi, τ, k.τ_i, x_i, k.Γ_i)
        t_i < t && ((t, kind) = (t_i, IceExhausted))
    end
    if k.below_triple ? iszero(x_i) : x_i > 0
        t_i = saturation_time(δ_i, inv_τ, τ, A_δi)
        t_i < t && ((t, kind) = (t_i, IceSaturation))
    end
    δ_t, δ_i_t, I, I_i = frozen_evolution(δ, δ_i, inv_τ, τ, A_δ, A_δi, t)
    return t, kind, δ_t, δ_i_t, condensate_increments(k, a, I, I_i)...
end
