"""
Shared analytic EP segment: derive active phases, step to the earliest event, return
updated `(δ, q_liq, q_ice)` and mean rates over `dt`.
"""
function ep_next_segment(
    δ::FT,
    q_liq::FT,
    q_ice::FT,
    q_sl::FT,
    q_si::FT,
    τ_liq::FT,
    τ_ice::FT,
    Γ_l::FT,
    Γ_i::FT,
    A_nwL::FT,
    A_nwI::FT,
    A_wbf::FT,
    τ_both::FT,
    BF::Bool,
    Δt_left::FT,
) where {FT}
    cofs = q_sl - q_si
    δ_i = δ + cofs
    la0 = (q_liq > 0) || (δ > 0) || (iszero(δ) && A_nwL > 0)
    ia0 = if BF
        (q_ice > 0) || (δ_i > 0) || (iszero(δ_i) && A_nwI > 0)
    else
        q_ice > 0 && (δ_i < 0 || (iszero(δ_i) && A_nwI < 0))
    end
    A0 = (la0 && ia0) ? A_wbf : (la0 ? A_nwL : (ia0 ? A_nwI : A_nwL))
    τ0 = (la0 && ia0) ? τ_both : (la0 ? τ_liq : (ia0 ? τ_ice : FT(Inf)))
    deq0 = A0 * τ0
    δ_eq0 = (ia0 && !la0) ? (deq0 - cofs) : deq0
    δi_eq0 = (ia0 && !la0) ? deq0 : (deq0 + cofs)
    la = la0 || (δ ≥ 0 && isfinite(τ0) && δ_eq0 > 0)
    ia = if BF
        ia0 || (δ_i ≥ 0 && isfinite(τ0) && δi_eq0 > 0)
    else
        ia0 || (q_ice > 0 && iszero(δ_i) && isfinite(τ0) && δi_eq0 < 0)
    end

    A_c = (la && ia) ? A_wbf : (la ? A_nwL : (ia ? A_nwI : A_nwL))
    τ = (la && ia) ? τ_both : (la ? τ_liq : (ia ? τ_ice : FT(Inf)))
    relax_ice = ia && !la

    if relax_ice
        x = δ_i
        t_liqsat = t_δ_hit_value(cofs, x, A_c, τ)
        t_icesat = t_δ_hit_value(zero(FT), x, A_c, τ)
    else
        x = δ
        t_liqsat = isfinite(τ) ? t_δ_hit_value(zero(FT), x, A_c, τ) : t_out_of_x(x, A_c)
        t_icesat = isfinite(τ) ? t_δ_hit_value(-cofs, x, A_c, τ) : t_out_of_x(x + cofs, A_c)
    end
    δ_eq = A_c * τ
    δi_eq = relax_ice ? δ_eq : δ_eq + cofs
    t_liqout = FT(Inf)
    t_iceout = FT(Inf)
    if la && q_liq > 0 && (δ < 0 || δ_eq < 0)
        t_liqout = get_t_out_of_q_liq(δ, A_c, τ, τ_liq, q_liq, Γ_l)
    end
    if ia && q_ice > 0 && (δ_i < 0 || δi_eq < 0)
        t_iceout = if la && ia
            get_t_out_of_q_ice(δ, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si)
        else
            get_t_out_of_q_ice_no_WBF(δ_i, A_c, τ, τ_ice, q_ice, Γ_i)
        end
    end

    min_t, iev = find_min_t(t_liqsat, t_icesat, t_liqout, t_iceout)
    dt = min(min_t, Δt_left)
    terminal = !(min_t < Δt_left)

    if iszero(dt)
        S_ql_leg = zero(FT)
        S_qi_leg = zero(FT)
    elseif la && ia
        S_ql_leg = S_ql_func(A_c, τ, τ_liq, δ, dt, Γ_l)
        S_qi_leg = S_qi_func(A_c, τ, τ_ice, δ, dt, Γ_i, q_sl, q_si)
    elseif la
        S_ql_leg = S_ql_func_indiv(A_c, τ_liq, δ, dt, Γ_l)
        S_qi_leg = zero(FT)
    elseif ia
        S_ql_leg = zero(FT)
        S_qi_leg = S_qi_func_indiv(A_c, τ_ice, δ_i, dt, Γ_i)
    else
        S_ql_leg = zero(FT)
        S_qi_leg = zero(FT)
    end
    if !terminal && iev == 3
        S_ql_leg = -q_liq / dt
    elseif !terminal && iev == 4
        S_qi_leg = -q_ice / dt
    end

    if !isfinite(τ)
        δ += A_c * dt
    elseif relax_ice
        δ_i_new = δ_func_EP(A_c, τ, δ_i, dt)
        δ = δ_i_new - cofs
    else
        δ = δ_func_EP(A_c, τ, δ, dt)
    end
    q_liq = max(q_liq + S_ql_leg * dt, zero(FT))
    q_ice = max(q_ice + S_qi_leg * dt, zero(FT))
    if !terminal
        if iev == 1
            δ = zero(FT)
        elseif iev == 2
            δ = -cofs
        elseif iev == 3
            q_liq = zero(FT)
        else
            q_ice = zero(FT)
        end
    end
    return (; dt, iev, terminal, S_ql_leg, S_qi_leg, δ, q_liq, q_ice, A_c, τ, la, ia)
end

"""
    morrison_milbrandt_2015(::MM2015FixedT, g, L_i, …, Δt; opts) -> (S_ql, S_qi)

Frozen-T event-driven Appendix C solver. Rates are specific-humidity means over `Δt`.
"""
function morrison_milbrandt_2015(
    ::MM2015FixedT,
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
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    q_vap_eq_liq::FT,
    q_vap_eq_ice::FT,
    Δt::FT;
    opts::MM2015FixedTOpts{FT} = MM2015FixedTOpts{FT}(),
) where {FT}
    iszero(Δt) && return (zero(FT), zero(FT))

    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    δ = limit_δ(q_vap - q_vap_eq_liq, q_vap, q_liq, q_ice)
    δ_0i = limit_δ(q_vap - q_vap_eq_ice, q_vap, q_liq, q_ice)
    P = prepare_coefficients(
        g, L_i, L_l, c_p, T_freeze, dqsl_dT, dqsi_dT, e_sl, e_si,
        q_vap_eq_liq, q_vap_eq_ice, q_liq, q_ice, δ, δ_0i, opts.dqvdt, opts.dTdt,
    )
    q_sl, q_si = P.q_sl, P.q_si
    Γ_l, Γ_i = P.Γ_l, P.Γ_i
    BF = T < P.T_freeze
    A_nwL = A_c_func_no_WBF(q_sl, P.g, w, P.c_p, P.e_sl, P.dqsl_dT, P.dqvdt, P.dTdt, p, ρ)
    A_nwI = A_c_func_no_WBF(q_si, P.g, w, P.c_p, P.e_si, P.dqsi_dT, P.dqvdt, P.dTdt, p, ρ)
    A_wbf = A_c_func_with_and_without_WBF(
        τ_ice, Γ_i, q_sl, q_si, P.g, w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT, P.dqvdt, P.dTdt, p, ρ,
    ).A_c
    τ_both = τ_func_combined(τ_liq, τ_ice, P.L_i, P.c_p, P.dqsl_dT, Γ_i)

    S_ql_mass = zero(FT)
    S_qi_mass = zero(FT)
    Δt_left = Δt

    for _ in 1:opts.max_events
        Δt_left ≤ zero(FT) && break
        seg = ep_next_segment(
            δ, q_liq, q_ice, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i,
            A_nwL, A_nwI, A_wbf, τ_both, BF, Δt_left,
        )
        δ = seg.δ
        q_liq = seg.q_liq
        q_ice = seg.q_ice
        S_ql_mass += seg.S_ql_leg * seg.dt
        S_qi_mass += seg.S_qi_leg * seg.dt
        Δt_left -= seg.dt
        seg.terminal && break
    end
    Δt_left > 0 && error("EP solver exceeded max_events = $(opts.max_events) with Δt_left = $Δt_left")
    return (S_ql_mass / Δt, S_qi_mass / Δt)
end

function morrison_milbrandt_2015(::MM2015FixedT, inputs::MM2015Inputs{FT}; kwargs...) where {FT}
    return morrison_milbrandt_2015(
        MM2015FixedT(),
        inputs.g, inputs.L_i, inputs.L_l, inputs.c_p, inputs.T_freeze,
        inputs.dqsl_dT, inputs.dqsi_dT, inputs.e_sl, inputs.e_si,
        inputs.ρ, inputs.p, inputs.T, inputs.w, inputs.τ_liq, inputs.τ_ice,
        inputs.q_tot, inputs.q_liq, inputs.q_ice, inputs.q_vap_eq_liq, inputs.q_vap_eq_ice, inputs.Δt;
        kwargs...,
    )
end
