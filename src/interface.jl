"""
    validate(problem, Δt)

Opt-in validation for setup and debugging. The hot [`tendencies`](@ref) path does
not call this function.
"""
function validate(problem::MM2015Problem, Δt)
    st = problem.state
    ts = problem.timescales
    forcing = problem.forcing
    all(isfinite, (st.T, st.p, st.x_tot, st.x_liq, st.x_ice, ts.τ_liq, ts.τ_ice, Δt)) ||
        throw(ArgumentError("state, timescales, and Δt must be finite"))
    Δt ≥ zero(Δt) || throw(ArgumentError("Δt must be nonnegative"))
    st.T > zero(st.T) || throw(ArgumentError("temperature must be positive"))
    st.p > zero(st.p) || throw(ArgumentError("pressure must be positive"))
    ts.τ_liq > zero(ts.τ_liq) || throw(ArgumentError("τ_liq must be positive"))
    ts.τ_ice > zero(ts.τ_ice) || throw(ArgumentError("τ_ice must be positive"))
    st.x_liq ≥ zero(st.x_liq) || throw(ArgumentError("liquid must be nonnegative"))
    st.x_ice ≥ zero(st.x_ice) || throw(ArgumentError("ice must be nonnegative"))
    st.x_tot ≥ st.x_liq + st.x_ice ||
        throw(ArgumentError("total water must contain liquid and ice"))
    all(isfinite, (forcing.dpdt, forcing.dTdt_external, forcing.dx_vap_dt, forcing.dx_tot_dt)) ||
        throw(ArgumentError("forcing must be finite"))
    return nothing
end

@inline function _specific_state(::SpecificHumidity, st::MM2015State)
    return st, one(st.x_tot) - st.x_tot
end

@inline function _specific_state(::DryAirMixingRatio, st::MM2015State)
    q_d = inv(one(st.x_tot) + st.x_tot)
    qstate = MM2015State(st.T, st.p, q_d * st.x_tot, q_d * st.x_liq, q_d * st.x_ice)
    return qstate, q_d
end

"""
Prepare the old scalar kernel arguments from one thermodynamic source of truth.

This helper is internal during the API transition. Saturation and material
properties come from `thermo`; callers supply only state, forcing, and relaxation
times.
"""
function _scalar_inputs(problem::MM2015Problem{FT, Basis}, Δt::FT) where {FT, Basis}
    qstate, q_d = _specific_state(problem.basis, problem.state)
    sc = thermo_scalars(
        problem.thermo,
        qstate.T,
        qstate.p,
        qstate.x_tot,
        qstate.x_liq,
        qstate.x_ice,
    )
    if problem.basis isa SpecificHumidity
        xstate = problem.state
        x_sl, x_si = sc.q_sl, sc.q_si
        dxsl_dT, dxsi_dT = sc.dqsl_dT, sc.dqsi_dT
        c_p = sc.c_p
    else
        xstate = problem.state
        x_sl, x_si = sc.q_sl / q_d, sc.q_si / q_d
        dxsl_dT, dxsi_dT = sc.dqsl_dT / q_d, sc.dqsi_dT / q_d
        c_p = sc.c_p / q_d
    end
    w = -problem.forcing.dpdt / (sc.ρ * sc.g)
    inputs = MM2015Inputs(
        sc.g,
        sc.L_i,
        sc.L_l,
        c_p,
        sc.T_freeze,
        dxsl_dT,
        dxsi_dT,
        sc.e_sl,
        sc.e_si,
        sc.ρ,
        xstate.p,
        xstate.T,
        w,
        problem.timescales.τ_liq,
        problem.timescales.τ_ice,
        xstate.x_tot,
        xstate.x_liq,
        xstate.x_ice,
        x_sl,
        x_si,
        Δt,
    )
    return inputs
end

@inline function _scheme_options(
    ::MM2015PiecewiseLinear,
    forcing::MM2015Forcing{FT},
) where {FT}
    return MM2015PiecewiseLinearOpts{FT}()
end

@inline function _scheme_options(::MM2015FixedT, forcing::MM2015Forcing{FT}) where {FT}
    return MM2015FixedTOpts{FT}(;
        dqvdt = forcing.dx_vap_dt,
        dTdt = forcing.dTdt_external,
    )
end

@inline function _scheme_options(::MM2015, forcing::MM2015Forcing{FT}) where {FT}
    return MM2015Opts{FT}(;
        dqvdt = forcing.dx_vap_dt,
        dTdt = forcing.dTdt_external,
    )
end

"""
    tendencies(scheme, problem, Δt) -> (S_liq, S_ice)

Mean homogeneous phase-change tendencies over `Δt`, in the same moisture basis
as `problem.state`. External forcing affects the trajectory but is not included
in the returned phase-change tendencies.
"""
function tendencies(
    scheme::Union{MM2015PiecewiseLinear, MM2015FixedT},
    problem::MM2015Problem{FT},
    Δt::FT,
) where {FT}
    inputs = _scalar_inputs(problem, Δt)
    opts = _scheme_options(scheme, problem.forcing)
    if scheme isa MM2015PiecewiseLinear
        return morrison_milbrandt_2015(
            scheme,
            inputs;
            opts,
            dqvdt = problem.forcing.dx_vap_dt,
            dTdt = problem.forcing.dTdt_external,
        )
    end
    return morrison_milbrandt_2015(scheme, inputs; opts)
end

function tendencies(
    scheme::MM2015,
    problem::MM2015Problem{FT, SpecificHumidity},
    Δt::FT,
) where {FT}
    inputs = _scalar_inputs(problem, Δt)
    opts = _scheme_options(scheme, problem.forcing)
    return morrison_milbrandt_2015(scheme, inputs; thermo = problem.thermo, opts)
end

function tendencies(
    scheme::MM2015,
    problem::MM2015Problem{FT, DryAirMixingRatio},
    Δt::FT,
) where {FT}
    qstate, q_d = _specific_state(problem.basis, problem.state)
    rforcing = problem.forcing
    q_tot_dt = q_d^2 * rforcing.dx_tot_dt
    q_vap = qstate.x_tot - qstate.x_liq - qstate.x_ice
    q_vap_dt = q_d * rforcing.dx_vap_dt - q_d * q_vap * rforcing.dx_tot_dt
    qforcing = MM2015Forcing(
        rforcing.dpdt,
        rforcing.dTdt_external,
        q_vap_dt,
        q_tot_dt,
    )
    qproblem = MM2015Problem(
        SpecificHumidity(),
        problem.thermo,
        qstate,
        problem.timescales,
        qforcing,
    )
    q_rates = tendencies(scheme, qproblem, Δt)
    return q_rates[1] / q_d, q_rates[2] / q_d
end
