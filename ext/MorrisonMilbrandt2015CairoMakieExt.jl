module MorrisonMilbrandt2015CairoMakieExt

using CairoMakie: CairoMakie as CM
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

const _IEV_LABEL = Dict(1 => "liq sat", 2 => "ice sat", 3 => "liq out", 4 => "ice out", 5 => "T_triple")
const _QL = CM.RGB(0.13, 0.55, 0.13)
const _QI = CM.RGB(0.12, 0.25, 0.70)
const _SUM = CM.RGB(0.80, 0.15, 0.15)
const _QL0 = CM.RGB(0.05, 0.28, 0.05)
const _QI0 = CM.RGB(0.05, 0.12, 0.40)
const _EVT = CM.RGB(0.55, 0.05, 0.55)

_save(path, fig) = (CM.save(path, fig; px_per_unit = 2); fig)

_with_dt(inputs, Δt) = MM2015.MM2015Inputs(inputs; Δt = convert(typeof(inputs.Δt), Δt))

function _ep_opts(opts::MM2015.MM2015Opts{FT}) where {FT}
    return MM2015.MM2015FixedTOpts{FT}(;
        dqvdt = opts.dqvdt,
        dTdt = opts.dTdt,
        max_events = opts.max_events,
    )
end

function _pl_opts(opts::MM2015.MM2015Opts{FT}) where {FT}
    return MM2015.MM2015PiecewiseLinearOpts{FT}(; max_events = opts.max_events)
end

"""Initial-phase flags and frozen C6 coefficients (same start as the solvers)."""
function _c6_start(inputs; opts)
    q_vap = MM2015.vapor_specific_humidity(inputs.q_tot, inputs.q_liq, inputs.q_ice)
    δ = MM2015.limit_δ(q_vap - inputs.q_vap_eq_liq, q_vap, inputs.q_liq, inputs.q_ice)
    δi = MM2015.limit_δ(q_vap - inputs.q_vap_eq_ice, q_vap, inputs.q_liq, inputs.q_ice)
    P = MM2015.prepare_coefficients(
        inputs.g, inputs.L_i, inputs.L_l, inputs.c_p, inputs.T_freeze,
        inputs.dqsl_dT, inputs.dqsi_dT, inputs.e_sl, inputs.e_si,
        inputs.q_vap_eq_liq, inputs.q_vap_eq_ice, inputs.q_liq, inputs.q_ice, δ, δi, opts.dqvdt, opts.dTdt,
    )
    BF = inputs.T < P.T_freeze
    la = (inputs.q_liq > 0) || (δ > 0)
    ia = BF ? ((inputs.q_ice > 0) || (δi > 0)) : (inputs.q_ice > 0 && δi < 0)
    A_nw = MM2015.A_c_func_no_WBF(P.q_sl, P.g, inputs.w, P.c_p, P.e_sl, P.dqsl_dT, opts.dqvdt, opts.dTdt, inputs.p, inputs.ρ)
    A_wbf = MM2015.A_c_func(
        inputs.τ_ice, P.Γ_i, P.q_sl, P.q_si, P.g, inputs.w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT,
        opts.dqvdt, opts.dTdt, inputs.p, inputs.ρ,
    )
    τ_both = MM2015.τ_func_combined(inputs.τ_liq, inputs.τ_ice, P.L_i, P.c_p, P.dqsl_dT, P.Γ_i)
    A_c = (la && ia) ? A_wbf : A_nw
    τ = (la && ia) ? τ_both : (la ? inputs.τ_liq : (ia ? inputs.τ_ice : Inf))
    S0_ql, S0_qi = MM2015.c6_instantaneous_rates(
        δ, P.q_sl, P.q_si, inputs.τ_liq, inputs.τ_ice, P.Γ_l, P.Γ_i, A_c, τ, la, ia, zero(δ),
    )
    return (; δ, δi, S0_ql, S0_qi, P, la, ia)
end

function _dt_grid(inputs; n::Int = 96, tmin = 1e-4, tmax = 1e2)
    FT = typeof(inputs.Δt)
    return FT.(exp10.(range(log10(tmin), log10(tmax); length = n)))
end

"""Mean rates vs `Δt` for T-updating MM2015, frozen-T EP, and piecewise-linear."""
function _rate_sweep(inputs; thermo, opts, Δts)
    FT = typeof(inputs.Δt)
    n = length(Δts)
    S_ql = Vector{FT}(undef, n)
    S_qi = Vector{FT}(undef, n)
    S_ql_ep = Vector{FT}(undef, n)
    S_qi_ep = Vector{FT}(undef, n)
    S_ql_pl = Vector{FT}(undef, n)
    S_qi_pl = Vector{FT}(undef, n)
    eopts = _ep_opts(opts)
    popts = _pl_opts(opts)
    for (i, Δt) in enumerate(Δts)
        inp = _with_dt(inputs, Δt)
        S_ql[i], S_qi[i] = MM2015.morrison_milbrandt_2015(MM2015.MM2015(), inp; thermo, opts)
        S_ql_ep[i], S_qi_ep[i] = MM2015.morrison_milbrandt_2015(MM2015.MM2015FixedT(), inp; opts = eopts)
        S_ql_pl[i], S_qi_pl[i] = MM2015.morrison_milbrandt_2015(MM2015.MM2015PiecewiseLinear(), inp; opts = popts)
    end
    return (; Δts, S_ql, S_qi, S_ql_ep, S_qi_ep, S_ql_pl, S_qi_pl)
end

function _liq_out_dt(Δts, S_ql, q_liq)
    q_liq ≤ 0 && return nothing
    for (Δt, S) in zip(Δts, S_ql)
        if S * Δt ≤ -0.99 * q_liq
            return Δt
        end
    end
    return nothing
end

function _annotate_events!(ax, events, y)
    seen = Set{Int}()
    for ev in events
        CM.vlines!(ax, [ev.t]; color = (_EVT, 0.65), linestyle = :dash)
        if !(ev.iev in seen)
            push!(seen, ev.iev)
            CM.text!(
                ax, ev.t, y;
                text = get(_IEV_LABEL, ev.iev, "iev=$(ev.iev)"),
                align = (:center, :bottom),
                fontsize = 12,
                color = _EVT,
            )
        end
    end
    return ax
end

# ---------------------------------------------------------------------------
# Frozen-coeff samples on MM2015 segments (closed form between milestones)
# ---------------------------------------------------------------------------

function _sample_segment!(
    t, q_liq, q_ice, q_vap, q_sl, q_si, δ, δi, T, p, ρ,
    thermo, t0, dt, n,
    θ_li0, p0, ρ0, g, w, q_tot0, q_liq0, q_ice0, dqvdt,
    δ0, q_sl0, q_si0, τ_liq, τ_ice, Γ_l, Γ_i, A_c, τ, la, ia, dθdT0, dTdt,
    skip_t0::Bool,
)
    n = max(2, n)
    i0 = skip_t0 ? 2 : 1
    for i in i0:n
        tloc = dt * (i - 1) / (n - 1)
        st = MM2015.residual_state_at(
            thermo, tloc, θ_li0, p0, ρ0, g, w, q_tot0, q_liq0, q_ice0, dqvdt,
            δ0, q_sl0, q_si0, τ_liq, τ_ice, Γ_l, Γ_i, A_c, τ, la, ia, dθdT0, dTdt,
        )
        push!(t, t0 + tloc)
        push!(q_liq, st.q_liq)
        push!(q_ice, st.q_ice)
        push!(q_vap, st.q_vap)
        push!(q_sl, st.q_sl)
        push!(q_si, st.q_si)
        push!(δ, st.F_δ)
        push!(δi, st.F_δi)
        push!(T, st.T)
        push!(p, st.p)
        push!(ρ, st.ρ)
    end
    return nothing
end

function _march(inputs::MM2015.MM2015Inputs{FT}; thermo, solver, opts, n_per_segment) where {FT}
    q_tot, q_liq, q_ice = inputs.q_tot, inputs.q_liq, inputs.q_ice
    p, ρ, T, w = inputs.p, inputs.ρ, inputs.T, inputs.w
    τ_liq, τ_ice, Δt = inputs.τ_liq, inputs.τ_ice, inputs.Δt
    q_tot ≥ 1 && error("q_tot ≥ 1")
    θ_li = MM2015.liquid_ice_pottemp(thermo, T, p, q_tot, q_liq, q_ice)
    dθdT0 = MM2015.dθ_li_dT(thermo, T, p, q_tot, q_liq, q_ice)
    q_vap = MM2015.vapor_specific_humidity(q_tot, q_liq, q_ice)
    δ = MM2015.limit_δ(q_vap - inputs.q_vap_eq_liq, q_vap, q_liq, q_ice)
    δ_0i = MM2015.limit_δ(q_vap - inputs.q_vap_eq_ice, q_vap, q_liq, q_ice)
    dqvdt, dTdt = opts.dqvdt, opts.dTdt
    P = MM2015.prepare_coefficients(
        inputs.g, inputs.L_i, inputs.L_l, inputs.c_p, inputs.T_freeze,
        inputs.dqsl_dT, inputs.dqsi_dT, inputs.e_sl, inputs.e_si,
        inputs.q_vap_eq_liq, inputs.q_vap_eq_ice, q_liq, q_ice, δ, δ_0i, dqvdt, dTdt,
    )
    q_sl, q_si, Γ_l, Γ_i, T_triple = P.q_sl, P.q_si, P.Γ_l, P.Γ_i, P.T_freeze
    BF = T < T_triple
    A_nwL = MM2015.A_c_func_no_WBF(q_sl, P.g, w, P.c_p, P.e_sl, P.dqsl_dT, dqvdt, dTdt, p, ρ)
    A_nwI = MM2015.A_c_func_no_WBF(q_si, P.g, w, P.c_p, P.e_si, P.dqsi_dT, dqvdt, dTdt, p, ρ)
    A_wbfL = MM2015.A_c_func(τ_ice, Γ_i, q_sl, q_si, P.g, w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT, dqvdt, dTdt, p, ρ)
    A_wbfI = MM2015.A_c_func(τ_liq, Γ_l, q_si, q_sl, P.g, w, P.c_p, P.e_si, P.L_l, P.dqsi_dT, dqvdt, dTdt, p, ρ)
    τL = MM2015.τ_func_combined(τ_liq, τ_ice, P.L_i, P.c_p, P.dqsl_dT, Γ_i)
    τI = MM2015.τ_func_combined(τ_ice, τ_liq, P.L_l, P.c_p, P.dqsi_dT, Γ_l)

    t = FT[]; q_liq_t = FT[]; q_ice_t = FT[]; q_vap_t = FT[]
    q_sl_t = FT[]; q_si_t = FT[]; δ_t = FT[]; δi_t = FT[]
    T_t = FT[]; p_t = FT[]; ρ_t = FT[]
    events = @NamedTuple{t::FT, iev::Int}[]
    t_abs = zero(Δt)
    Δt_left = Δt
    skip_t0 = false
    T_triple_out = T_triple

    for _ in 1:opts.max_events
        Δt_left ≤ zero(FT) && break
        seg = MM2015.mm2015_next_segment(
            thermo, solver, δ, q_tot, q_liq, q_ice, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i,
            A_nwL, A_nwI, A_wbfL, A_wbfI, τL, τI, BF, Δt_left,
            θ_li, dθdT0, dTdt, p, ρ, T, P.g, w, dqvdt, T_triple,
        )
        _sample_segment!(
            t, q_liq_t, q_ice_t, q_vap_t, q_sl_t, q_si_t, δ_t, δi_t, T_t, p_t, ρ_t,
            thermo, t_abs, seg.dt, n_per_segment,
            θ_li, p, ρ, P.g, w, q_tot, q_liq, q_ice, dqvdt,
            δ, q_sl, q_si, τ_liq, τ_ice, Γ_l, Γ_i, seg.A_c, seg.τ, seg.la, seg.ia, dθdT0, dTdt,
            skip_t0,
        )
        q_liq_t[end] = seg.q_liq
        q_ice_t[end] = seg.q_ice
        q_vap_t[end] = MM2015.vapor_specific_humidity(seg.q_tot, seg.q_liq, seg.q_ice)
        q_sl_t[end] = seg.q_sl
        q_si_t[end] = seg.q_si
        δ_t[end] = seg.F_δ
        δi_t[end] = seg.F_δi
        T_t[end] = seg.T
        p_t[end] = seg.p
        ρ_t[end] = seg.ρ
        if !seg.terminal && !iszero(seg.iev)
            push!(events, (; t = t_abs + seg.dt, iev = Int(seg.iev)))
        end
        t_abs += seg.dt
        Δt_left -= seg.dt
        skip_t0 = true
        q_tot, q_liq, q_ice = seg.q_tot, seg.q_liq, seg.q_ice
        T, p, ρ, δ = seg.T, seg.p, seg.ρ, seg.δ
        θ_li = MM2015.liquid_ice_pottemp(thermo, T, p, q_tot, q_liq, q_ice)
        dθdT0 = MM2015.dθ_li_dT(thermo, T, p, q_tot, q_liq, q_ice)
        sc = seg.sc
        P = MM2015.prepare_coefficients(
            sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze, sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
            sc.q_sl, sc.q_si, q_liq, q_ice, δ, δ + (sc.q_sl - sc.q_si), dqvdt, dTdt,
        )
        q_sl, q_si, Γ_l, Γ_i, T_triple = P.q_sl, P.q_si, P.Γ_l, P.Γ_i, P.T_freeze
        T_triple_out = T_triple
        BF = T < T_triple
        A_nwL = MM2015.A_c_func_no_WBF(q_sl, P.g, w, P.c_p, P.e_sl, P.dqsl_dT, dqvdt, dTdt, p, ρ)
        A_nwI = MM2015.A_c_func_no_WBF(q_si, P.g, w, P.c_p, P.e_si, P.dqsi_dT, dqvdt, dTdt, p, ρ)
        A_wbfL = MM2015.A_c_func(τ_ice, Γ_i, q_sl, q_si, P.g, w, P.c_p, P.e_sl, P.L_i, P.dqsl_dT, dqvdt, dTdt, p, ρ)
        A_wbfI = MM2015.A_c_func(τ_liq, Γ_l, q_si, q_sl, P.g, w, P.c_p, P.e_si, P.L_l, P.dqsi_dT, dqvdt, dTdt, p, ρ)
        τL = MM2015.τ_func_combined(τ_liq, τ_ice, P.L_i, P.c_p, P.dqsl_dT, Γ_i)
        τI = MM2015.τ_func_combined(τ_ice, τ_liq, P.L_l, P.c_p, P.dqsi_dT, Γ_l)
        seg.terminal && break
    end
    Δt_left > 0 && error("plot march exceeded max_events = $(opts.max_events)")
    return (;
        t, q_liq = q_liq_t, q_ice = q_ice_t, q_vap = q_vap_t, q_sl = q_sl_t, q_si = q_si_t,
        δ = δ_t, δi = δi_t, T = T_t, p = p_t, ρ = ρ_t, events, T_triple = T_triple_out,
        q_liq0 = inputs.q_liq, q_ice0 = inputs.q_ice,
    )
end

function _ep_q_vs_t(inputs, t; opts)
    eopts = _ep_opts(opts)
    q_liq = similar(t)
    q_ice = similar(t)
    for (i, ti) in enumerate(t)
        if iszero(ti)
            q_liq[i] = inputs.q_liq
            q_ice[i] = inputs.q_ice
            continue
        end
        S_ql, S_qi = MM2015.morrison_milbrandt_2015(MM2015.MM2015FixedT(), _with_dt(inputs, ti); opts = eopts)
        q_liq[i] = max(inputs.q_liq + S_ql * ti, zero(ti))
        q_ice[i] = max(inputs.q_ice + S_qi * ti, zero(ti))
    end
    return q_liq, q_ice
end

function _traj(inputs; thermo, solver, opts, n_per_segment)
    return _march(inputs; thermo, solver, opts, n_per_segment)
end

function MM2015.plot_rates(
    inputs::MM2015.MM2015Inputs;
    thermo = MM2015.DefaultThermodynamicsBackend(),
    opts = MM2015.MM2015Opts{typeof(inputs.Δt)}(),
    Δts = _dt_grid(inputs),
    size = (1100, 460),
    title = nothing,
)
    sw = _rate_sweep(inputs; thermo, opts, Δts)
    st = _c6_start(inputs; opts)
    floor_l = -inputs.q_liq ./ sw.Δts
    floor_i = -inputs.q_ice ./ sw.Δts
    t_out = _liq_out_dt(sw.Δts, sw.S_ql, inputs.q_liq)
    ttl = if title !== nothing
        title
    else
        "mean sources vs Δt   T=$(inputs.T) K,  q_liq=$(inputs.q_liq),  q_ice=$(inputs.q_ice),  τ_liq=$(inputs.τ_liq) s,  τ_ice=$(inputs.τ_ice) s"
    end
    ys = vcat(sw.S_ql, sw.S_qi, sw.S_ql_ep, sw.S_qi_ep, sw.S_ql_pl, sw.S_qi_pl, [st.S0_ql, st.S0_qi, zero(st.S0_ql)])
    ymin, ymax = extrema(ys)
    pad = 0.18 * (ymax - ymin)

    fig = CM.Figure(size = size, fontsize = 13)
    ax = CM.Axis(
        fig[1, 1];
        xlabel = "Δt [s]",
        ylabel = "mean S [s⁻¹]",
        title = ttl,
        xscale = log10,
        xminorticksvisible = true,
        xminorgridvisible = true,
        xgridvisible = true,
        ygridvisible = true,
    )
    CM.hlines!(ax, [0.0]; color = (:black, 0.45), linestyle = :dash, linewidth = 1)
    CM.hlines!(ax, [st.S0_ql]; color = _QL0, linestyle = :dot, linewidth = 1.6, label = "S_ql (t = 0 C6)")
    CM.hlines!(ax, [st.S0_qi]; color = _QI0, linestyle = :dot, linewidth = 1.6, label = "S_qi (t = 0 C6)")
    CM.lines!(ax, sw.Δts, floor_l; color = _QL0, linewidth = 1.5, label = "−q_liq / Δt")
    CM.lines!(ax, sw.Δts, floor_i; color = _QI0, linewidth = 1.5, label = "−q_ice / Δt")
    CM.lines!(ax, sw.Δts, sw.S_ql_pl; color = _QL, linestyle = :dot, linewidth = 1.6, label = "S_ql piecewise-linear")
    CM.lines!(ax, sw.Δts, sw.S_qi_pl; color = _QI, linestyle = :dot, linewidth = 1.6, label = "S_qi piecewise-linear")
    CM.lines!(ax, sw.Δts, sw.S_ql_pl .+ sw.S_qi_pl; color = _SUM, linestyle = :dot, linewidth = 1.1, label = "sum piecewise-linear")
    CM.lines!(ax, sw.Δts, sw.S_ql_ep; color = _QL, linestyle = :dash, linewidth = 2, label = "S_ql EP (frozen T)")
    CM.lines!(ax, sw.Δts, sw.S_qi_ep; color = _QI, linestyle = :dash, linewidth = 2, label = "S_qi EP (frozen T)")
    CM.lines!(ax, sw.Δts, sw.S_ql_ep .+ sw.S_qi_ep; color = _SUM, linestyle = :dash, linewidth = 1.2, label = "sum EP")
    CM.lines!(ax, sw.Δts, sw.S_ql; color = _QL, linewidth = 2.6, label = "S_ql MM2015 (T-updating)")
    CM.lines!(ax, sw.Δts, sw.S_qi; color = _QI, linewidth = 2.6, label = "S_qi MM2015 (T-updating)")
    CM.lines!(ax, sw.Δts, sw.S_ql .+ sw.S_qi; color = _SUM, linewidth = 1.6, label = "sum MM2015")
    CM.vlines!(ax, [sw.Δts[1], sw.Δts[end]]; color = (_EVT, 0.45), linestyle = :dash)
    if t_out !== nothing
        CM.vlines!(ax, [t_out]; color = (_EVT, 0.9), linestyle = :dash, linewidth = 1.5)
        CM.text!(ax, t_out, ymax; text = "liq out", align = (:center, :bottom), fontsize = 12, color = _EVT)
    end
    CM.ylims!(ax, ymin - pad, ymax + pad)
    CM.Legend(fig[1, 2], ax; framevisible = false, labelsize = 11, rowgap = 2)
    CM.colgap!(fig.layout, 8)
    return fig
end

function _evolution_figure(traj, inputs; opts, size, title)
    q_ep_l, q_ep_i = _ep_q_vs_t(inputs, traj.t; opts)
    st = _c6_start(inputs; opts)
    fig = CM.Figure(size = size, fontsize = 13)
    axq = CM.Axis(fig[1, 1]; ylabel = "q [kg kg⁻¹]", title, xticklabelsvisible = false)
    axδ = CM.Axis(fig[2, 1]; ylabel = "δ [kg kg⁻¹]", xlabel = "t [s]")
    CM.linkxaxes!(axq, axδ)

    CM.hlines!(axq, [0.0]; color = (:black, 0.4), linestyle = :dash, linewidth = 1, label = "q = 0")
    CM.hlines!(axq, [traj.q_liq0]; color = (_QL, 0.35), linestyle = :dot, label = "q_liq(0)")
    CM.hlines!(axq, [traj.q_ice0]; color = (_QI, 0.35), linestyle = :dot, label = "q_ice(0)")
    CM.lines!(axq, traj.t, q_ep_l; color = _QL, linestyle = :dash, linewidth = 1.8, label = "q_liq EP (frozen T)")
    CM.lines!(axq, traj.t, q_ep_i; color = _QI, linestyle = :dash, linewidth = 1.8, label = "q_ice EP (frozen T)")
    CM.lines!(axq, traj.t, traj.q_liq; color = _QL, linewidth = 2.4, label = "q_liq MM2015")
    CM.lines!(axq, traj.t, traj.q_ice; color = _QI, linewidth = 2.4, label = "q_ice MM2015")
    ymaxq = maximum((maximum(traj.q_liq), maximum(traj.q_ice), maximum(q_ep_l), maximum(q_ep_i)))
    _annotate_events!(axq, traj.events, ymaxq)
    CM.ylims!(axq, -0.05 * ymaxq, 1.12 * ymaxq)

    CM.hlines!(axδ, [0.0]; color = (:black, 0.4), linestyle = :dash, linewidth = 1, label = "saturation")
    CM.hlines!(axδ, [st.δ]; color = (_QL, 0.45), linestyle = :dot, label = "δ(0)")
    CM.hlines!(axδ, [st.δi]; color = (_QI, 0.45), linestyle = :dot, label = "δi(0)")
    CM.lines!(axδ, traj.t, traj.δ; color = _QL, linewidth = 2.2, label = "δ = q_vap − q_sl")
    CM.lines!(axδ, traj.t, traj.δi; color = _QI, linewidth = 2.2, label = "δi = q_vap − q_si")
    δabs = maximum((maximum(abs, traj.δ), maximum(abs, traj.δi), abs(st.δ), abs(st.δi)))
    _annotate_events!(axδ, traj.events, 0.92 * δabs)
    CM.ylims!(axδ, -1.15 * δabs, 1.15 * δabs)

    CM.Legend(fig[1:2, 2], axq; framevisible = false, labelsize = 11, rowgap = 2)
    CM.Legend(fig[1:2, 3], axδ; framevisible = false, labelsize = 11, rowgap = 2)
    return fig
end

function MM2015.plot_evolution(
    traj;
    inputs = nothing,
    opts = MM2015.MM2015Opts{Float64}(),
    size = (1000, 720),
    title = "T-updating C6 (solid) vs frozen-T EP (dashed)",
)
    inputs === nothing && error("plot_evolution(traj) needs `inputs` so EP expected q(t) can be overlaid")
    return _evolution_figure(traj, inputs; opts, size, title)
end

function MM2015.plot_condensate(traj; inputs = nothing, opts = MM2015.MM2015Opts{Float64}(), size = (900, 400), title = "condensate")
    fig = CM.Figure(size = size, fontsize = 13)
    ax = CM.Axis(fig[1, 1]; xlabel = "t [s]", ylabel = "q [kg kg⁻¹]", title)
    CM.hlines!(ax, [0.0]; color = (:black, 0.4), linestyle = :dash, label = "q = 0")
    if inputs !== nothing
        q_ep_l, q_ep_i = _ep_q_vs_t(inputs, traj.t; opts)
        CM.lines!(ax, traj.t, q_ep_l; color = _QL, linestyle = :dash, linewidth = 1.8, label = "q_liq EP")
        CM.lines!(ax, traj.t, q_ep_i; color = _QI, linestyle = :dash, linewidth = 1.8, label = "q_ice EP")
    end
    CM.lines!(ax, traj.t, traj.q_liq; color = _QL, linewidth = 2.4, label = "q_liq MM2015")
    CM.lines!(ax, traj.t, traj.q_ice; color = _QI, linewidth = 2.4, label = "q_ice MM2015")
    ymaxq = maximum((maximum(traj.q_liq), maximum(traj.q_ice)))
    _annotate_events!(ax, traj.events, ymaxq)
    CM.Legend(fig[1, 2], ax; framevisible = false, labelsize = 11)
    return fig
end

function MM2015.plot_supersaturation(traj; inputs = nothing, opts = MM2015.MM2015Opts{Float64}(), size = (900, 400), title = "supersaturation")
    fig = CM.Figure(size = size, fontsize = 13)
    ax = CM.Axis(fig[1, 1]; xlabel = "t [s]", ylabel = "δ [kg kg⁻¹]", title)
    CM.hlines!(ax, [0.0]; color = (:black, 0.4), linestyle = :dash, label = "saturation")
    if inputs !== nothing
        st = _c6_start(inputs; opts)
        CM.hlines!(ax, [st.δ]; color = (_QL, 0.45), linestyle = :dot, label = "δ(0)")
        CM.hlines!(ax, [st.δi]; color = (_QI, 0.45), linestyle = :dot, label = "δi(0)")
    end
    CM.lines!(ax, traj.t, traj.δ; color = _QL, linewidth = 2.2, label = "δ (liquid)")
    CM.lines!(ax, traj.t, traj.δi; color = _QI, linewidth = 2.2, label = "δi (ice)")
    ymaxδ = maximum((maximum(abs, traj.δ), maximum(abs, traj.δi)))
    _annotate_events!(ax, traj.events, 0.9 * ymaxδ)
    CM.Legend(fig[1, 2], ax; framevisible = false, labelsize = 11)
    return fig
end

function MM2015.plot_specific_humidities(traj; inputs = nothing, opts = nothing, size = (900, 400), title = "vapor and saturation")
    fig = CM.Figure(size = size, fontsize = 13)
    ax = CM.Axis(fig[1, 1]; xlabel = "t [s]", ylabel = "q [kg kg⁻¹]", title)
    CM.lines!(ax, traj.t, traj.q_vap; color = _SUM, linewidth = 2.2, label = "q_vap")
    CM.lines!(ax, traj.t, traj.q_sl; color = _QL, linestyle = :dash, linewidth = 1.8, label = "q_sl")
    CM.lines!(ax, traj.t, traj.q_si; color = _QI, linestyle = :dash, linewidth = 1.8, label = "q_si")
    ymax = maximum((maximum(traj.q_vap), maximum(traj.q_sl), maximum(traj.q_si)))
    _annotate_events!(ax, traj.events, ymax)
    CM.Legend(fig[1, 2], ax; framevisible = false, labelsize = 11)
    return fig
end

function MM2015.plot_temperature(traj; inputs = nothing, opts = nothing, size = (900, 360), title = "temperature")
    fig = CM.Figure(size = size, fontsize = 13)
    ax = CM.Axis(fig[1, 1]; xlabel = "t [s]", ylabel = "T [K]", title)
    CM.lines!(ax, traj.t, traj.T; color = :black, linewidth = 2.2, label = "T")
    Tmin, Tmax = extrema(traj.T)
    span = Tmax - Tmin
    pad = max(0.08 * (span == 0 ? one(span) : span), 0.02)
    CM.ylims!(ax, Tmin - pad, Tmax + pad)
    if Tmin - pad ≤ traj.T_triple ≤ Tmax + pad
        CM.hlines!(ax, [traj.T_triple]; color = (:royalblue, 0.8), linestyle = :dash, label = "T_triple")
    else
        ax.ylabel = "T [K]   (T_triple = $(round(traj.T_triple; digits = 2)))"
    end
    _annotate_events!(ax, traj.events, Tmax)
    CM.Legend(fig[1, 2], ax; framevisible = false, labelsize = 11)
    return fig
end

function MM2015.plot_rates(path::AbstractString, inputs::MM2015.MM2015Inputs; kwargs...)
    return _save(path, MM2015.plot_rates(inputs; kwargs...))
end

for name in (:plot_evolution, :plot_condensate, :plot_supersaturation, :plot_specific_humidities, :plot_temperature)
    @eval begin
        function MM2015.$name(path::AbstractString, traj; kwargs...)
            return _save(path, MM2015.$name(traj; kwargs...))
        end
        function MM2015.$name(
            path::AbstractString,
            inputs::MM2015.MM2015Inputs{FT};
            thermo = MM2015.DefaultThermodynamicsBackend(),
            solver = MM2015.default_mm2015_solver,
            opts::MM2015.MM2015Opts{FT} = MM2015.MM2015Opts{FT}(),
            n_per_segment::Int = 48,
            kwargs...,
        ) where {FT}
            traj = _traj(inputs; thermo, solver, opts, n_per_segment)
            return MM2015.$name(path, traj; inputs, opts, kwargs...)
        end
    end
end

end # module
