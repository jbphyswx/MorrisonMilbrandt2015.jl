"""Vapor specific humidity at which linear liquid and ice sources balance `dδdt_no_S`."""
function get_qv_eq_point(q_sl::FT, q_si::FT, τ_liq::FT, τ_ice::FT; dδdt_no_S::FT = FT(0)) where {FT}
    if isinf(τ_liq)
        if isinf(τ_ice)
            return FT(Inf)
        else
            return τ_ice * dδdt_no_S + (q_si - q_sl) + q_sl
        end
    elseif isinf(τ_ice)
        return τ_liq * dδdt_no_S + q_sl
    end
    qv_eq =
        q_sl * (τ_ice / (τ_liq + τ_ice)) +
        q_si * (τ_liq / (τ_liq + τ_ice)) +
        dδdt_no_S * (τ_liq * τ_ice) / (τ_liq + τ_ice)
    if qv_eq == q_sl
        error("Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast (qv_eq indistinguishable from q_sl).")
    elseif qv_eq == q_si
        error("Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast (qv_eq indistinguishable from q_si).")
    end
    return qv_eq
end

"""`(δ_eq, δi_eq)` at which linear sources balance `dδdt_no_S`."""
function get_δ_eq_point(q_sl::FT, q_si::FT, τ_liq::FT, τ_ice::FT; dδdt_no_S::FT = FT(0)) where {FT}
    if isinf(τ_liq)
        if isinf(τ_ice)
            return FT(Inf), FT(Inf)
        else
            δ_eq = τ_ice * dδdt_no_S + (q_si - q_sl)
            δi_eq = τ_ice * dδdt_no_S
            return δ_eq, δi_eq
        end
    elseif isinf(τ_ice)
        δ_eq = τ_liq * dδdt_no_S
        δi_eq = τ_liq * dδdt_no_S + (q_sl - q_si)
        return δ_eq, δi_eq
    end
    δ_eq = (q_si - q_sl) * (τ_liq / (τ_liq + τ_ice)) + dδdt_no_S * (τ_liq * τ_ice) / (τ_liq + τ_ice)
    δi_eq = (q_sl - q_si) * (τ_ice / (τ_liq + τ_ice)) + dδdt_no_S * (τ_liq * τ_ice) / (τ_liq + τ_ice)
    if iszero(δi_eq)
        error("Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast (δi_eq = 0).")
    end
    if iszero(δ_eq)
        error("Timescales τ_liq = $τ_liq and/or τ_ice = $τ_ice are too fast (δ_eq = 0).")
    end
    return δ_eq, δi_eq
end

function calculate_next_standard_milestone_time_given_milestones(
    S_ql::FT,
    S_qi::FT,
    dδdt::FT,
    δ_0::FT,
    δ_0i::FT,
    δ_eq::FT,
    δi_eq::FT,
    q_liq::FT,
    q_ice::FT;
    dδdt_is_full_tendency::Bool = false,
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
)::Tuple{FT, MilestoneType} where {FT}
    t_out_of_liq = t_out_of_x(q_liq, S_ql)
    t_out_of_ice = t_out_of_x(q_ice, S_qi)

    if !at_δ_eq_point
        dδ = dδdt_is_full_tendency ? dδdt : dδdt - (S_ql * Γ_l + S_qi * Γ_i)
        t_hit_liq_sat = t_out_of_x(δ_0, dδ)
        t_hit_ice_sat = t_out_of_x(δ_0i, dδ)
        if allow_δ_eq_point && is_same_supersaturation_regime(δ_eq, δ_0, δ_0i)
            t_hit_eq_point = min(t_out_of_x(δ_0 - δ_eq, dδ), t_out_of_x(δ_0i - δi_eq, dδ))
        else
            t_hit_eq_point = FT(Inf)
        end
    else
        t_hit_eq_point = FT(Inf)
        t_hit_liq_sat = FT(Inf)
        t_hit_ice_sat = FT(Inf)
    end

    min_t, i_min_t = findmin((t_out_of_liq, t_out_of_ice, t_hit_eq_point, t_hit_liq_sat, t_hit_ice_sat))
    milestone = if isinf(min_t)
        NotAtSupersaturationMilestone
    elseif i_min_t == 1
        OutOfLiquidMilestone
    elseif i_min_t == 2
        OutOfIceMilestone
    elseif i_min_t == 3
        AtSupersaturationStationaryPointMilestone
    elseif i_min_t == 4
        AtSaturationOverLiquidMilestone
    else
        AtSaturationOverIceMilestone
    end
    return min_t, milestone
end

function calculate_next_standard_milestone_time(
    regime::AbstractSaturationRegime,
    q_sl::FT,
    q_si::FT,
    q_liq::FT,
    q_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    below_freezing::Bool,
    τ_liq::FT,
    τ_ice::FT;
    dδdt_no_S::FT = FT(0),
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
    at_δ_eq_point::Bool = false,
    allow_δ_eq_point::Bool = true,
)::Tuple{FT, MilestoneType, FT, FT, FT, FT} where {FT}
    if below_freezing
        if regime isa Supersaturated
            S_ql = δ_0 / (Γ_l * τ_liq)
            S_qi = δ_0i / (Γ_i * τ_ice)
        elseif regime isa WBF
            S_ql = (q_liq > FT(0)) ? δ_0 / (Γ_l * τ_liq) : FT(0)
            τ_liq = (q_liq > FT(0)) ? τ_liq : FT(Inf)
            S_qi = δ_0i / (Γ_i * τ_ice)
        elseif regime isa Subsaturated
            S_ql = (q_liq > FT(0)) ? δ_0 / (Γ_l * τ_liq) : FT(0)
            τ_liq = (q_liq > FT(0)) ? τ_liq : FT(Inf)
            S_qi = (q_ice > FT(0)) ? δ_0i / (Γ_i * τ_ice) : FT(0)
            τ_ice = (q_ice > FT(0)) ? τ_ice : FT(Inf)
        else
            error("Unknown regime: $regime")
        end
    else
        if regime isa Supersaturated
            S_ql = δ_0 / (Γ_l * τ_liq)
            S_qi = FT(0)
            τ_ice = FT(Inf)
        elseif regime isa WBF
            S_ql = δ_0 / (Γ_l * τ_liq)
            S_qi = (q_ice > FT(0)) ? δ_0i / (Γ_i * τ_ice) : FT(0)
            τ_ice = (q_ice > FT(0)) ? τ_ice : FT(Inf)
        elseif regime isa Subsaturated
            S_ql = (q_liq > FT(0)) ? δ_0 / (Γ_l * τ_liq) : FT(0)
            τ_liq = (q_liq > FT(0)) ? τ_liq : FT(Inf)
            S_qi = (q_ice > FT(0)) ? δ_0i / (Γ_i * τ_ice) : FT(0)
            τ_ice = (q_ice > FT(0)) ? τ_ice : FT(Inf)
        else
            error("Unknown regime: $regime")
        end
    end

    δ_eq, δi_eq = get_δ_eq_point(q_sl, q_si, τ_liq, τ_ice; dδdt_no_S = dδdt_no_S)
    dδdt = dδdt_no_S - (S_ql * Γ_l + S_qi * Γ_i)
    min_t, milestone = calculate_next_standard_milestone_time_given_milestones(
        S_ql,
        S_qi,
        dδdt,
        δ_0,
        δ_0i,
        δ_eq,
        δi_eq,
        q_liq,
        q_ice;
        dδdt_is_full_tendency = true,
        Γ_l = Γ_l,
        Γ_i = Γ_i,
        at_δ_eq_point = at_δ_eq_point,
        allow_δ_eq_point = allow_δ_eq_point,
    )
    return min_t, milestone, S_ql, S_qi, δ_eq, δi_eq
end

function get_new_saturation_regime_type_from_milestone(
    milestone::MilestoneType,
    regime::AbstractSaturationRegime,
    old_δ_0::FT,
    old_δ_0i::FT,
) where {FT}
    if milestone ∈ (
        NotAtSupersaturationMilestone,
        OutOfLiquidMilestone,
        OutOfIceMilestone,
        AtSupersaturationStationaryPointMilestone,
    )
        if regime isa Supersaturated
            return Supersaturated
        elseif regime isa WBF
            return WBF
        elseif regime isa Subsaturated
            return Subsaturated
        else
            error("invalid regime type $regime")
        end
    elseif milestone == AtSaturationOverLiquidMilestone
        if is_below_freezing(regime)
            if regime isa Supersaturated
                return WBF
            elseif regime isa WBF
                return Supersaturated
            elseif regime isa Subsaturated
                iszero(old_δ_0i) || error("Below freezing, cannot hit liquid sat from Subsaturated unless δ_0i was zero")
                return Supersaturated
            else
                error("invalid regime type $regime")
            end
        else
            if regime isa Supersaturated
                iszero(old_δ_0i) || error("Above freezing, cannot hit liquid sat from Supersaturated unless δ_0i was zero")
                return Subsaturated
            elseif regime isa WBF
                return Subsaturated
            elseif regime isa Subsaturated
                return WBF
            else
                error("invalid regime type $regime")
            end
        end
    elseif milestone == AtSaturationOverIceMilestone
        if is_below_freezing(regime)
            if regime isa Supersaturated
                iszero(old_δ_0) || error("Below freezing, cannot hit ice sat from Supersaturated unless δ_0 was zero")
                return Subsaturated
            elseif regime isa WBF
                return Subsaturated
            elseif regime isa Subsaturated
                return WBF
            else
                error("invalid regime type $regime")
            end
        else
            if regime isa Supersaturated
                return WBF
            elseif regime isa WBF
                return Supersaturated
            elseif regime isa Subsaturated
                iszero(old_δ_0) || error("Above freezing, cannot hit ice sat from Subsaturated unless δ_0 was zero")
                return Supersaturated
            else
                error("invalid regime type $regime")
            end
        end
    else
        error("unknown milestone $milestone")
    end
end

get_new_saturation_regime_enum_type_from_milestone(
    milestone::MilestoneType,
    regime::AbstractSaturationRegime,
    old_δ_0::FT,
    old_δ_0i::FT,
) where {FT} = get_saturation_regime_enum_type(
    get_new_saturation_regime_type_from_milestone(milestone, regime, old_δ_0, old_δ_0i),
)

function step(
    regime::AbstractSaturationRegime,
    Δt::FT,
    q_liq::FT,
    q_ice::FT,
    δ_0::FT,
    δ_0i::FT,
    δ_eq::FT,
    δi_eq::FT,
    q_sl::FT,
    q_si::FT,
    S_ql::FT,
    S_qi::FT,
    milestone::MilestoneType = NotAtSupersaturationMilestone;
    at_δ_eq_point::Bool = false,
    dδdt_no_S::FT = FT(0),
    Γ_l::FT = FT(1),
    Γ_i::FT = FT(1),
)::Tuple{FT, FT, FT, FT, SaturationRegimeEnumTypes} where {FT}
    old_δ_0 = δ_0
    old_δ_0i = δ_0i

    if milestone == NotAtSupersaturationMilestone
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        if !at_δ_eq_point
            dδ = -(S_ql * Γ_l * Δt) - (S_qi * Γ_i * Δt)
            δ_0 += dδ + (dδdt_no_S * Δt)
            δ_0i += dδ + (dδdt_no_S * Δt)
        end
    elseif milestone == OutOfLiquidMilestone
        if !at_δ_eq_point
            dδ = +q_liq * Γ_l
            δ_0 += dδ
            δ_0i += dδ
        end
        q_liq = FT(0)
        q_ice += S_qi * Δt
        if !at_δ_eq_point
            dδ = -(S_qi * Γ_i * Δt)
            δ_0 += dδ + (dδdt_no_S * Δt)
            δ_0i += dδ + (dδdt_no_S * Δt)
        end
    elseif milestone == OutOfIceMilestone
        q_liq += S_ql * Δt
        if !at_δ_eq_point
            dδ = -(S_ql * Γ_l * Δt)
            δ_0 += dδ + (dδdt_no_S * Δt)
            δ_0i += dδ + (dδdt_no_S * Δt)
            dδ = +q_ice * Γ_i
            δ_0 += dδ
            δ_0i += dδ
        end
        q_ice = FT(0)
    elseif milestone == AtSupersaturationStationaryPointMilestone
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        δ_0i = δi_eq
        δ_0 = δ_eq
    elseif milestone == AtSaturationOverLiquidMilestone
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        δ_0i = (δ_0i - δ_0)
        δ_0 = FT(0)
    elseif milestone == AtSaturationOverIceMilestone
        q_liq += S_ql * Δt
        q_ice += S_qi * Δt
        δ_0 = (δ_0 - δ_0i)
        δ_0i = FT(0)
    else
        error("Unknown milestone type: $milestone")
    end

    new_regime_enum_type = get_new_saturation_regime_enum_type_from_milestone(milestone, regime, old_δ_0, old_δ_0i)
    return q_liq, q_ice, δ_0, δ_0i, new_regime_enum_type
end

"""
    morrison_milbrandt_2015_piecewise_linear(g, L_i, L_l, c_p, T_freeze, dqsl_dT, dqsi_dT, e_sl, e_si,
                                     ρ, p, T, w, τ_liq, τ_ice, q_tot, q_liq, q_ice, q_vap_eq_liq, q_vap_eq_ice, Δt; opts) -> (S_ql, S_qi)

Piecewise-linear supersaturation limiter: `S = δ / (τ Γ)` until the next milestone
(out of condensate, hit saturation, or hit a δ-stationary point). External forcing enters
as `dδdt_no_S = A_c` without the WBF term.
"""
function morrison_milbrandt_2015_piecewise_linear(
    g::FT,
    L_i::FT,
    L_l::FT,
    c_p::FT,
    T_freeze::FT,
    dqsl_dT::FT,
    dqsi_dT::FT,
    e_sl::FT,
    e_si::FT,
    ρ::FT,
    p::FT,
    T::FT,
    w::FT,
    τ_liq::FT,
    τ_ice::FT,
    q_tot::FT, q_liq::FT, q_ice::FT,
    q_vap_eq_liq::FT, q_vap_eq_ice::FT,
    Δt::FT;
    opts::MM2015PiecewiseLinearOpts{FT} = MM2015PiecewiseLinearOpts{FT}(),
    dqvdt::FT = zero(FT),
    dTdt::FT = zero(FT),
) where {FT}
    iszero(Δt) && return (zero(FT), zero(FT))

    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    δ_0 = limit_δ(q_vap - q_vap_eq_liq, q_vap, q_liq, q_ice)
    δ_0i = limit_δ(q_vap - q_vap_eq_ice, q_vap, q_liq, q_ice)
    P = prepare_coefficients(
        g, L_i, L_l, c_p, T_freeze, dqsl_dT, dqsi_dT, e_sl, e_si,
        q_vap_eq_liq, q_vap_eq_ice, q_liq, q_ice, δ_0, δ_0i, dqvdt, dTdt,
    )
    q_sl, q_si = P.q_sl, P.q_si
    Γ_l, Γ_i = P.Γ_l, P.Γ_i
    q_liq, q_ice = P.q_liq, P.q_ice
    δ_0, δ_0i = P.δ_0, P.δ_0i
    dδdt_no_S = A_c_func_no_WBF(q_sl, P.g, w, P.c_p, P.e_sl, P.dqsl_dT, dqvdt, dTdt, p, ρ)

    S_ql = FT(0)
    S_qi = FT(0)
    Δt_left = Δt
    dt = FT(0)
    last_milestone = NotAtSupersaturationMilestone
    depth = 0
    below_freezing = T < P.T_freeze

    if iszero(δ_0) || iszero(δ_0i)
        if iszero(δ_0)
            last_milestone = AtSaturationOverLiquidMilestone
        elseif iszero(δ_0i)
            last_milestone = AtSaturationOverIceMilestone
        end
        dδdt_0 = get_dδdt_0(δ_0, δ_0i, q_liq, q_ice, τ_liq, τ_ice, dδdt_no_S, below_freezing)
        regime = get_saturation_regime(δ_0, δ_0i, q_liq, q_ice, below_freezing; dδdt = dδdt_0)
    else
        regime = get_saturation_regime(δ_0, δ_0i, q_liq, q_ice, below_freezing)
    end

    while (Δt_left > FT(0)) && (depth < opts.max_events)
        milestone_t, milestone, S_ql_addit, S_qi_addit, δ_eq, δi_eq = calculate_next_standard_milestone_time(
            regime,
            q_sl,
            q_si,
            q_liq,
            q_ice,
            δ_0,
            δ_0i,
            below_freezing,
            τ_liq,
            τ_ice;
            dδdt_no_S = dδdt_no_S,
            Γ_l = Γ_l,
            Γ_i = Γ_i,
            at_δ_eq_point = (last_milestone == AtSupersaturationStationaryPointMilestone),
        )
        dt_here = min(Δt_left, milestone_t)
        if !(milestone_t < Δt_left)
            return resolve_S_S_addit(
                S_ql,
                S_qi,
                dt,
                S_ql_addit,
                S_qi_addit,
                Δt_left,
                Δt,
            )
        end
        at_δ_eq_point = last_milestone == AtSupersaturationStationaryPointMilestone && ((δ_eq == δ_0) || (δ_eq == δ_0i))

        q_liq, q_ice, δ_0, δ_0i, new_regime_enum_type = step(
            regime,
            dt_here,
            q_liq,
            q_ice,
            δ_0,
            δ_0i,
            δ_eq,
            δi_eq,
            q_sl,
            q_si,
            S_ql_addit,
            S_qi_addit,
            (milestone_t < Δt_left) ? milestone : NotAtSupersaturationMilestone;
            at_δ_eq_point = at_δ_eq_point,
            dδdt_no_S = dδdt_no_S,
            Γ_l = Γ_l,
            Γ_i = Γ_i,
        )
        regime = add_regime_parameters(new_regime_enum_type, q_liq, q_ice, below_freezing)
        S_ql, S_qi = resolve_S_S_addit(S_ql, S_qi, dt, S_ql_addit, S_qi_addit, dt_here, dt + dt_here)
        dt += dt_here
        Δt_left -= dt_here
        if !(milestone == NotAtSupersaturationMilestone) && (milestone == last_milestone)
            error("Hit the same MM2015PiecewiseLinear milestone twice: $milestone")
        end
        last_milestone = milestone
        depth += 1
    end
    Δt_left > 0 && error("MM2015PiecewiseLinear solver exceeded max_events = $(opts.max_events) with Δt_left = $Δt_left")
    return S_ql, S_qi
end



function morrison_milbrandt_2015(
    ::MM2015PiecewiseLinear,
    args...;
    kwargs...,
)
    return morrison_milbrandt_2015_piecewise_linear(args...; kwargs...)
end

function morrison_milbrandt_2015(::MM2015PiecewiseLinear, inputs::MM2015Inputs{FT}; kwargs...) where {FT}
    return morrison_milbrandt_2015_piecewise_linear(
        inputs.g, inputs.L_i, inputs.L_l, inputs.c_p, inputs.T_freeze,
        inputs.dqsl_dT, inputs.dqsi_dT, inputs.e_sl, inputs.e_si,
        inputs.ρ, inputs.p, inputs.T, inputs.w, inputs.τ_liq, inputs.τ_ice,
        inputs.q_tot, inputs.q_liq, inputs.q_ice, inputs.q_vap_eq_liq, inputs.q_vap_eq_ice, inputs.Δt;
        kwargs...,
    )
end