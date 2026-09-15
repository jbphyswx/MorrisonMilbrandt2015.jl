"""
C6 mean condensation / deposition rates over `dt` at frozen coefficients.

When both phases are active this is Appendix C liquid-frame C6 (`S_ql_func` /
`S_qi_func`). Single-species C6 when only one phase is active.
"""
@inline function c6_rates(
    δ,
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
    dt,
)
    δi = δ + (q_sl - q_si)
    if iszero(dt)
        z = zero(δ)
        return z, z
    elseif la && ia
        S_ql = S_ql_func(A_c, τ, τ_liq, δ, dt, Γ_l)
        S_qi = S_qi_func(A_c, τ, τ_ice, δ, dt, Γ_i, q_sl, q_si)
        return S_ql, S_qi
    elseif la
        return S_ql_func_indiv(A_c, τ_liq, δ, dt, Γ_l), zero(δ)
    elseif ia
        return zero(δ), S_qi_func_indiv(A_c, τ_ice, δi, dt, Γ_i)
    else
        z = zero(δ)
        return z, z
    end
end

"""
Apply frozen C6 over `t`, then NonEquilibrium `pθq` at the new hydrostatic `p`.

`p(t) = p_0 − w ρ_0 g t` (adiabatic cooling is this `p` in the invert; do not add
`−w g / c_p` again). Latent heat is the condensate in `θ_li`; do not add `L S / c_p`
again. Host `dTdt` evolves `θ_li` as `θ_li(t) = θ_li(0) + (∂θ_li/∂T) dTdt t`.
"""
@inline function residual_state_at(
    thermo,
    t,
    θ_li0,
    p0,
    ρ0,
    g,
    w,
    q_tot0,
    q_liq0,
    q_ice0,
    dqvdt,
    S_ql,
    S_qi,
    dθdT0,
    dTdt,
)
    q_liq = max(q_liq0 + S_ql * t, zero(q_liq0))
    q_ice = max(q_ice0 + S_qi * t, zero(q_ice0))
    q_tot = q_tot0 + dqvdt * t
    p = p0 - w * ρ0 * g * t
    θ_li = θ_li_evolved(θ_li0, dθdT0, dTdt, t)
    T = air_temperature_noneq_pθq(thermo, p, θ_li, q_tot, q_liq, q_ice)
    sc = thermo_scalars(thermo, T, p, q_tot, q_liq, q_ice)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    return (;
        T,
        p,
        ρ = sc.ρ,
        q_tot,
        q_liq,
        q_ice,
        q_vap,
        q_sl = sc.q_sl,
        q_si = sc.q_si,
        F_δ = q_vap - sc.q_sl,
        F_δi = q_vap - sc.q_si,
        F_T = T - sc.T_freeze,
        sc,
        θ_li,
    )
end

@inline function residual_state_at(
    thermo,
    t,
    θ_li0,
    p0,
    ρ0,
    g,
    w,
    q_tot0,
    q_liq0,
    q_ice0,
    dqvdt,
    δ,
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
    dθdT0,
    dTdt,
)
    S_ql, S_qi = c6_rates(
        δ, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i, A_c, τ, la, ia, t,
    )
    st = residual_state_at(
        thermo, t, θ_li0, p0, ρ0, g, w, q_tot0, q_liq0, q_ice0, dqvdt, S_ql, S_qi, dθdT0, dTdt,
    )
    return merge(st, (; S_ql, S_qi))
end

"""`q = max(q0 + S̄ t, 0)`; `q' = S_inst` above the floor, else `0`."""
@inline function floored_q_rate(q0, S_mean, S_inst, t)
    u = q0 + S_mean * t
    return u > zero(u) ? S_inst : zero(u)
end

"""
Analytic `F′` on the residual trajectory.

`p' = −w ρ₀ g`, `θ' = (∂θ_li/∂T)_0 dTdt`, `q_tot' = dqvdt`,
`q_liq' = S_{ql,inst}` (or `0` on the floor), `T'` from [`dT_noneq_dt`](@ref),
`F_δ' = q_vap' − q_sl'`, `F_T' = T'`.
"""
@inline function residual_deriv_at(
    thermo,
    t,
    θ_li0,
    p0,
    ρ0,
    g,
    w,
    q_tot0,
    q_liq0,
    q_ice0,
    dqvdt,
    δ,
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
    dθdT0,
    dTdt,
)
    st = residual_state_at(
        thermo, t, θ_li0, p0, ρ0, g, w, q_tot0, q_liq0, q_ice0, dqvdt,
        δ, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i, A_c, τ, la, ia,
        dθdT0, dTdt,
    )
    S_ql_inst, S_qi_inst = c6_instantaneous_rates(
        δ, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i, A_c, τ, la, ia, t,
    )
    dql = floored_q_rate(q_liq0, st.S_ql, S_ql_inst, t)
    dqi = floored_q_rate(q_ice0, st.S_qi, S_qi_inst, t)
    dqt = dqvdt
    dqv = dqt - dql - dqi
    dpdt = -w * ρ0 * g
    dθdt = dθdT0 * dTdt
    dT = dT_noneq_dt(
        thermo, st.T, st.p, st.θ_li, st.q_tot, st.q_liq, st.q_ice,
        dpdt, dθdt, dqt, dql, dqi,
    )
    dqsl = dq_vap_sat_liq_dt(thermo, st.T, st.p, st.q_tot, st.q_sl, dT, dpdt, dqt)
    dqsi = dq_vap_sat_ice_dt(thermo, st.T, st.p, st.q_tot, st.q_si, dT, dpdt, dqt)
    return (;
        dF_δ = dqv - dqsl,
        dF_δi = dqv - dqsi,
        dF_T = dT,
    )
end

"""Active phases for one frozen-coefficient segment (same rules as [`ep_next_segment`](@ref))."""
function mm2015_active_phases(
    δ::FT,
    q_liq::FT,
    q_ice::FT,
    q_sl::FT,
    q_si::FT,
    τ_liq::FT,
    τ_ice::FT,
    A_nwL::FT,
    A_nwI::FT,
    A_wbf::FT,
    τ_both::FT,
    BF::Bool,
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
    return la, ia
end

"""
    mm2015_next_segment(...) -> named tuple

One residual T-updating segment: Lambert for condensate out; [`leftmost_root`](@ref)
for saturation and `T_triple`. Ice-frame C5 (`A_wbfI`, `τI`) is only a `t_C5` knot
for `F_δi`. Mass is liquid-frame C6 when both phases are active.

Returned `iev`: 1 liquid sat, 2 ice sat, 3 liquid out, 4 ice out, 5 freeze, 0 end of `Δt_left`.
"""
function mm2015_next_segment(
    thermo,
    solver,
    δ::FT,
    q_tot::FT,
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
    A_wbfL::FT,
    A_wbfI::FT,
    τL::FT,
    τI::FT,
    BF::Bool,
    Δt_left::FT,
    θ_li0::FT,
    dθdT0::FT,
    dTdt::FT,
    p0::FT,
    ρ0::FT,
    T0::FT,
    g::FT,
    w::FT,
    dqvdt::FT,
    T_triple::FT,
) where {FT}
    la, ia = mm2015_active_phases(δ, q_liq, q_ice, q_sl, q_si, τ_liq, τ_ice, A_nwL, A_nwI, A_wbfL, τL, BF)
    if la && ia
        A_c = A_wbfL
        τ = τL
    elseif la
        A_c = A_nwL
        τ = τ_liq
    elseif ia
        A_c = A_nwI
        τ = τ_ice
    else
        A_c = A_nwL
        τ = FT(Inf)
    end

    cofs = q_sl - q_si
    δi = δ + cofs
    δ_eq = A_c * τ
    δi_eq = (ia && !la) ? δ_eq : (isfinite(τ) ? δ_eq + cofs : FT(Inf))

    t_liqout = FT(Inf)
    t_iceout = FT(Inf)
    if la && q_liq > 0 && (δ < 0 || δ_eq < 0)
        t_liqout = get_t_out_of_q_liq(δ, A_c, τ, τ_liq, q_liq, Γ_l)
    end
    if ia && q_ice > 0 && (δi < 0 || δi_eq < 0)
        t_iceout = if la && ia
            get_t_out_of_q_ice(δ, A_c, τ, τ_ice, q_ice, Γ_i, q_sl, q_si)
        else
            get_t_out_of_q_ice_no_WBF(δi, A_c, τ, τ_ice, q_ice, Γ_i)
        end
    end

    function state_at(t)
        return residual_state_at(
            thermo, t, θ_li0, p0, ρ0, g, w, q_tot, q_liq, q_ice, dqvdt,
            δ, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i, A_c, τ, la, ia,
            dθdT0, dTdt,
        )
    end
    function deriv_at(t)
        return residual_deriv_at(
            thermo, t, θ_li0, p0, ρ0, g, w, q_tot, q_liq, q_ice, dqvdt,
            δ, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i, A_c, τ, la, ia,
            dθdT0, dTdt,
        )
    end

    t_c5_liq = isfinite(τ) ? t_δ_hit_value(zero(FT), δ, A_c, τ) : t_out_of_x(δ, A_c)
    t_c5_ice = if la && ia
        t_δ_hit_value(zero(FT), δi, A_wbfI, τI)
    elseif ia
        isfinite(τ) ? t_δ_hit_value(zero(FT), δi, A_c, τ) : t_out_of_x(δi, A_c)
    else
        FT(Inf)
    end

    fastest_scale = min(
        Δt_left,
        isfinite(τ) ? τ : Δt_left,
        isfinite(τ_liq) ? τ_liq : Δt_left,
        isfinite(τ_ice) ? τ_ice : Δt_left,
    )
    lo = max(nextfloat(zero(FT)), eps(FT) * fastest_scale)
    hi = Δt_left
    time_atol = FT(64) * eps(FT) * max(abs(hi), one(FT))
    q_scale = max(abs(δ), abs(δi), abs(q_sl), abs(q_si), floatmin(FT))
    q_residual_atol = FT(64) * eps(FT) * q_scale
    q_derivative_atol = q_residual_atol / max(abs(hi), one(FT))
    T_residual_atol = FT(64) * eps(FT) * max(abs(T0), one(FT))
    T_derivative_atol = T_residual_atol / max(abs(hi), one(FT))
    t_liqsat = if !iszero(δ)
        leftmost_root(
            t -> state_at(t).F_δ,
            t -> deriv_at(t).dF_δ,
            lo,
            hi,
            solver;
            t_c5 = t_c5_liq,
            time_atol,
            residual_atol = q_residual_atol,
            derivative_atol = q_derivative_atol,
        )
    else
        FT(Inf)
    end
    t_icesat = if (BF || q_ice > zero(FT)) && !iszero(δi)
        leftmost_root(
            t -> state_at(t).F_δi,
            t -> deriv_at(t).dF_δi,
            lo,
            hi,
            solver;
            t_c5 = t_c5_ice,
            time_atol,
            residual_atol = q_residual_atol,
            derivative_atol = q_derivative_atol,
        )
    else
        FT(Inf)
    end
    t_freeze = if !iszero(T0 - T_triple)
        leftmost_root(
            t -> state_at(t).F_T,
            t -> deriv_at(t).dF_T,
            lo,
            hi,
            solver;
            time_atol,
            residual_atol = T_residual_atol,
            derivative_atol = T_derivative_atol,
        )
    else
        FT(Inf)
    end

    min_t, iev = find_min_t(t_liqsat, t_icesat, t_liqout, t_iceout, t_freeze)
    dt = min(min_t, Δt_left)
    terminal = !(min_t < Δt_left)
    if terminal
        iev = 0
    end

    st = state_at(dt)
    S_ql_leg, S_qi_leg = st.S_ql, st.S_qi
    q_liq_n, q_ice_n = st.q_liq, st.q_ice
    δ_n, δi_n, T_n, sc_n = st.F_δ, st.F_δi, st.T, st.sc
    if !terminal && iev == 3
        S_ql_leg = iszero(dt) ? zero(FT) : -q_liq / dt
        q_liq_n = zero(FT)
    elseif !terminal && iev == 4
        S_qi_leg = iszero(dt) ? zero(FT) : -q_ice / dt
        q_ice_n = zero(FT)
    elseif !terminal && iev == 1
        δ_n = zero(FT)
        δi_n = st.q_sl - st.q_si
    elseif !terminal && iev == 2
        δi_n = zero(FT)
        δ_n = st.q_si - st.q_sl
    elseif !terminal && iev == 5
        T_n = T_triple
        sc_n = thermo_scalars(thermo, T_n, st.p, st.q_tot, q_liq_n, q_ice_n)
        q_vap_n = vapor_specific_humidity(st.q_tot, q_liq_n, q_ice_n)
        δ_n = q_vap_n - sc_n.q_sl
        δi_n = q_vap_n - sc_n.q_si
    end
    return (;
        dt,
        iev,
        terminal,
        S_ql_leg,
        S_qi_leg,
        δ = δ_n,
        q_tot = st.q_tot,
        q_liq = q_liq_n,
        q_ice = q_ice_n,
        T = T_n,
        p = st.p,
        ρ = st.ρ,
        q_sl = sc_n.q_sl,
        q_si = sc_n.q_si,
        F_δ = δ_n,
        F_δi = δi_n,
        F_T = T_n - T_triple,
        A_c,
        τ,
        la,
        ia,
        sc = sc_n,
    )
end

"""
    morrison_milbrandt_2015(::MM2015, g, L_i, …, Δt; thermo, solver, opts) -> (S_ql, S_qi)

T-updating Appendix C. Mass on a segment uses frozen liquid-frame C6. Saturation and
freeze are the soonest residual roots after a NonEquilibrium `pθq` update.

`thermo` implements [`air_temperature_noneq_pθq`](@ref), [`liquid_ice_pottemp`](@ref),
[`dθ_li_dT`](@ref), and [`thermo_scalars`](@ref). Defaults to
[`DefaultThermodynamicsBackend`](@ref).
"""
function morrison_milbrandt_2015(
    ::MM2015,
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
    thermo = DefaultThermodynamicsBackend(),
    solver = default_mm2015_solver,
    opts::MM2015Opts{FT} = MM2015Opts{FT}(),
) where {FT}
    iszero(Δt) && return (zero(FT), zero(FT))

    θ_li = liquid_ice_pottemp(thermo, T, p, q_tot, q_liq, q_ice)
    dθdT0 = dθ_li_dT(thermo, T, p, q_tot, q_liq, q_ice)
    q_vap = vapor_specific_humidity(q_tot, q_liq, q_ice)
    δ = limit_δ(q_vap - q_vap_eq_liq, q_vap, q_liq, q_ice)
    δ_0i = limit_δ(q_vap - q_vap_eq_ice, q_vap, q_liq, q_ice)
    dqvdt = opts.dqvdt
    dTdt = opts.dTdt

    P = prepare_coefficients(
        g, L_i, L_l, c_p, T_freeze, dqsl_dT, dqsi_dT, e_sl, e_si,
        q_vap_eq_liq, q_vap_eq_ice, q_liq, q_ice, δ, δ_0i, dqvdt, dTdt,
    )
    q_sl, q_si = P.q_sl, P.q_si
    Γ_l, Γ_i = P.Γ_l, P.Γ_i
    T_triple = P.T_freeze
    BF = T < T_triple || (T == T_triple && (dTdt - w * P.g / P.c_p) < zero(FT))

    A_nwL = A_c_func_no_WBF(q_sl, P.g, w, P.c_p, P.e_sl, P.dqsl_dT, dqvdt, dTdt, p, ρ)
    A_nwI = A_c_func_no_WBF(q_si, P.g, w, P.c_p, P.e_si, P.dqsi_dT, dqvdt, dTdt, p, ρ)
    A_wbfL = A_c_func(
        τ_ice, Γ_i, q_sl, q_si, P.g, w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT, dqvdt, dTdt, p, ρ,
    )
    A_wbfI = A_c_func(
        τ_liq, Γ_l, q_si, q_sl, P.g, w, P.c_p, P.e_si, P.L_l, P.dqsi_dT, dqvdt, dTdt, p, ρ,
    )
    τL = τ_func_combined(τ_liq, τ_ice, P.L_i, P.c_p, P.dqsl_dT, Γ_i)
    τI = τ_func_combined(τ_ice, τ_liq, P.L_l, P.c_p, P.dqsi_dT, Γ_l)

    S_ql_mass = zero(FT)
    S_qi_mass = zero(FT)
    Δt_left = Δt

    for _ in 1:opts.max_events
        Δt_left ≤ zero(FT) && break
        seg = mm2015_next_segment(
            thermo, solver, δ, q_tot, q_liq, q_ice, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i,
            A_nwL, A_nwI, A_wbfL, A_wbfI, τL, τI, BF, Δt_left,
            θ_li, dθdT0, dTdt, p, ρ, T, P.g, w, dqvdt, T_triple,
        )
        S_ql_mass += seg.S_ql_leg * seg.dt
        S_qi_mass += seg.S_qi_leg * seg.dt
        Δt_left -= seg.dt
        q_tot = seg.q_tot
        q_liq = seg.q_liq
        q_ice = seg.q_ice
        T = seg.T
        p = seg.p
        ρ = seg.ρ
        δ = seg.δ
        θ_li = liquid_ice_pottemp(thermo, T, p, q_tot, q_liq, q_ice)
        dθdT0 = dθ_li_dT(thermo, T, p, q_tot, q_liq, q_ice)
        sc = seg.sc
        P = prepare_coefficients(
            sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze, sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
            sc.q_sl, sc.q_si, q_liq, q_ice, δ, δ + (sc.q_sl - sc.q_si), dqvdt, dTdt,
        )
        q_sl, q_si = P.q_sl, P.q_si
        Γ_l, Γ_i = P.Γ_l, P.Γ_i
        T_triple = P.T_freeze
        BF = (!seg.terminal && seg.iev == 5) ? !BF :
             (T < T_triple || (T == T_triple && (dTdt - w * P.g / P.c_p) < zero(FT)))
        A_nwL = A_c_func_no_WBF(q_sl, P.g, w, P.c_p, P.e_sl, P.dqsl_dT, dqvdt, dTdt, p, ρ)
        A_nwI = A_c_func_no_WBF(q_si, P.g, w, P.c_p, P.e_si, P.dqsi_dT, dqvdt, dTdt, p, ρ)
        A_wbfL = A_c_func(
            τ_ice, Γ_i, q_sl, q_si, P.g, w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT, dqvdt, dTdt, p, ρ,
        )
        A_wbfI = A_c_func(
            τ_liq, Γ_l, q_si, q_sl, P.g, w, P.c_p, P.e_si, P.L_l, P.dqsi_dT, dqvdt, dTdt, p, ρ,
        )
        τL = τ_func_combined(τ_liq, τ_ice, P.L_i, P.c_p, P.dqsl_dT, Γ_i)
        τI = τ_func_combined(τ_ice, τ_liq, P.L_l, P.c_p, P.dqsi_dT, Γ_l)
        seg.terminal && break
    end
    Δt_left > 0 && error("MM2015 solver exceeded max_events = $(opts.max_events) with Δt_left = $Δt_left")
    return (S_ql_mass / Δt, S_qi_mass / Δt)
end

function morrison_milbrandt_2015(::MM2015, inputs::MM2015Inputs{FT}; kwargs...) where {FT}
    return morrison_milbrandt_2015(
        MM2015(),
        inputs.g, inputs.L_i, inputs.L_l, inputs.c_p, inputs.T_freeze,
        inputs.dqsl_dT, inputs.dqsi_dT, inputs.e_sl, inputs.e_si,
        inputs.ρ, inputs.p, inputs.T, inputs.w, inputs.τ_liq, inputs.τ_ice,
        inputs.q_tot, inputs.q_liq, inputs.q_ice, inputs.q_vap_eq_liq, inputs.q_vap_eq_ice, inputs.Δt;
        kwargs...,
    )
end
