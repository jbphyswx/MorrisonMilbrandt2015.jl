# =============================================== #
# Schemes
# =============================================== #

abstract type AbstractMM2015Scheme end

"""Piecewise-linear supersaturation limiter `S = δ / (τ Γ)` until a milestone."""
struct MM2015PiecewiseLinear <: AbstractMM2015Scheme end

"""Appendix C exponential with residual T-updating (NonEquilibrium `pθq`). Host `dTdt` evolves `θ_li`; adiabatic cooling is hydrostatic `p(t)` in the invert."""
struct MM2015 <: AbstractMM2015Scheme end

"""Appendix C exponential at frozen T (event-driven C5)."""
struct MM2015FixedT <: AbstractMM2015Scheme end

# =============================================== #
# Options
# =============================================== #

abstract type AbstractMM2015Options end

"""Options for [`MM2015PiecewiseLinear`](@ref)."""
Base.@kwdef struct MM2015PiecewiseLinearOpts{FT} <: AbstractMM2015Options
    max_events::Int = 15
end

"""Options for [`MM2015FixedT`](@ref)."""
Base.@kwdef struct MM2015FixedTOpts{FT} <: AbstractMM2015Options
    dqvdt::FT = zero(FT)
    dTdt::FT = zero(FT)
    max_events::Int = 12
end

"""Options for [`MM2015`](@ref) (same knobs as EP; saturation/freeze events are residual roots)."""
Base.@kwdef struct MM2015Opts{FT} <: AbstractMM2015Options
    dqvdt::FT = zero(FT)
    dTdt::FT = zero(FT)
    max_events::Int = 12
end

# =============================================== #
# Public problem description
# =============================================== #

"""Mass normalization used by the moisture coordinates in an [`MM2015Problem`](@ref)."""
abstract type AbstractMoistureBasis end

"""Water mass per unit dry-air mass."""
struct DryAirMixingRatio <: AbstractMoistureBasis end

"""Water mass per unit total moist-air mass."""
struct SpecificHumidity <: AbstractMoistureBasis end

"""
    MM2015State(T, p, x_tot, x_liq, x_ice)

Parcel thermodynamic and moisture state. Every `x_*` uses the problem's declared
moisture basis.
"""
struct MM2015State{FT}
    T::FT
    p::FT
    x_tot::FT
    x_liq::FT
    x_ice::FT
end

function MM2015State(T, p, x_tot, x_liq, x_ice)
    values = promote(T, p, x_tot, x_liq, x_ice)
    return MM2015State(values...)
end

"""Liquid and ice supersaturation-relaxation times [s]."""
struct MM2015Timescales{FT}
    τ_liq::FT
    τ_ice::FT
end

function MM2015Timescales(τ_liq, τ_ice)
    values = promote(τ_liq, τ_ice)
    return MM2015Timescales(values...)
end

"""
    MM2015Forcing(dpdt, dTdt_external, dx_vap_dt, dx_tot_dt)

Nonmicrophysical parcel forcing in the declared moisture basis. `dTdt_external`
contains radiation/mixing heating only; adiabatic and latent tendencies are
computed by the solver. `dx_vap_dt` and `dx_tot_dt` are separate so coordinate
transformations can retain the normalization tendency.
"""
struct MM2015Forcing{FT}
    dpdt::FT
    dTdt_external::FT
    dx_vap_dt::FT
    dx_tot_dt::FT
end

function MM2015Forcing(dpdt, dTdt_external, dx_vap_dt, dx_tot_dt = dx_vap_dt)
    values = promote(dpdt, dTdt_external, dx_vap_dt, dx_tot_dt)
    return MM2015Forcing(values...)
end

"""Complete homogeneous parcel problem consumed by [`tendencies`](@ref)."""
struct MM2015Problem{FT, Basis <: AbstractMoistureBasis, Thermo}
    basis::Basis
    thermo::Thermo
    state::MM2015State{FT}
    timescales::MM2015Timescales{FT}
    forcing::MM2015Forcing{FT}
end



# =============================================== #
# Saturation Regimes & Milestones
# =============================================== #



abstract type AbstractSaturationRegime end

has_ice(r::AbstractSaturationRegime) = r.has_ice
has_liquid(r::AbstractSaturationRegime) = r.has_liq
is_below_freezing(r::AbstractSaturationRegime) = r.below_freezing

struct Supersaturated <: AbstractSaturationRegime
    has_liq::Bool
    has_ice::Bool
    below_freezing::Bool
end

struct WBF <: AbstractSaturationRegime
    has_liq::Bool
    has_ice::Bool
    below_freezing::Bool
end

WBF(has_liq::Bool, has_ice::Bool) = WBF(has_liq, has_ice, true)

struct Subsaturated <: AbstractSaturationRegime
    has_liq::Bool
    has_ice::Bool
    below_freezing::Bool
end

Supersaturated(q_liq::FT, q_ice::FT, below_freezing::Bool) where {FT} =
    Supersaturated(q_liq > zero(FT), q_ice > zero(FT), below_freezing)
WBF(q_liq::FT, q_ice::FT, below_freezing::Bool) where {FT} =
    WBF(q_liq > zero(FT), q_ice > zero(FT), below_freezing)
Subsaturated(q_liq::FT, q_ice::FT, below_freezing::Bool) where {FT} =
    Subsaturated(q_liq > zero(FT), q_ice > zero(FT), below_freezing)

@enum SaturationRegimeEnumTypes begin
    SupersaturatedEnumType
    WBFEnumType
    SubsaturatedEnumType
end
get_saturation_regime_type(::Val{SupersaturatedEnumType}) = Supersaturated
get_saturation_regime_type(::Val{WBFEnumType}) = WBF
get_saturation_regime_type(::Val{SubsaturatedEnumType}) = Subsaturated
get_saturation_regime_enum_type(::Type{<:Supersaturated}) = SupersaturatedEnumType
get_saturation_regime_enum_type(::Type{<:WBF}) = WBFEnumType
get_saturation_regime_enum_type(::Type{<:Subsaturated}) = SubsaturatedEnumType

function add_regime_parameters(
    regime_enum_type::SaturationRegimeEnumTypes,
    q_liq::FT,
    q_ice::FT,
    below_freezing::Bool,
) where {FT}
    if regime_enum_type == SupersaturatedEnumType
        return Supersaturated(q_liq, q_ice, below_freezing)
    elseif regime_enum_type == WBFEnumType
        return WBF(q_liq, q_ice, below_freezing)
    else
        return Subsaturated(q_liq, q_ice, below_freezing)
    end
end

@enum MilestoneType begin
    OutOfLiquidMilestone
    OutOfIceMilestone
    AtSaturationOverLiquidMilestone
    AtSaturationOverIceMilestone
    AtSupersaturationStationaryPointMilestone
    NotAtSupersaturationMilestone
end

δi_from_δ(δ::FT, q_sl::FT, q_si::FT) where {FT} = δ + (q_sl - q_si)
δ_from_δi(δi::FT, q_sl::FT, q_si::FT) where {FT} = δi - (q_sl - q_si)

"""Initial regime type; at a sat boundary, `dδdt` breaks the tie (growth unless `dδdt < 0`)."""
function get_saturation_regime_type(δ::FT, δi::FT, BF::Bool; dδdt::FT = FT(0)) where {FT}
    if BF
        if δ < FT(0)
            if δi > FT(0)
                return WBF
            elseif iszero(δi)
                return (dδdt < FT(0)) ? Subsaturated : WBF
            else
                return Subsaturated
            end
        else
            if iszero(δ)
                return (dδdt < FT(0)) ? WBF : Supersaturated
            end
            return Supersaturated
        end
    else
        if δi < FT(0)
            if δ > FT(0)
                return WBF
            elseif iszero(δ)
                return (dδdt < FT(0)) ? Subsaturated : WBF
            else
                return Subsaturated
            end
        else
            if iszero(δi)
                return (dδdt < FT(0)) ? WBF : Supersaturated
            end
            return Supersaturated
        end
    end
end
get_saturation_regime_type(δ::FT, δi::FT, T::FT, T_freeze::FT) where {FT} =
    get_saturation_regime_type(δ, δi, T < T_freeze)
get_saturation_regime_type(q_vap::FT, q_sl::FT, q_si::FT, T::FT, T_freeze::FT; dδdt::FT = FT(0)) where {FT} =
    get_saturation_regime_type(q_vap - q_sl, q_vap - q_si, T < T_freeze; dδdt = dδdt)

function get_saturation_regime_enum_type(δ::FT, δi::FT, BF::Bool; dδdt::FT = FT(0)) where {FT}
    if BF
        if δ < FT(0)
            if δi > FT(0)
                return WBFEnumType
            elseif iszero(δi)
                return (dδdt < FT(0)) ? SubsaturatedEnumType : WBFEnumType
            else
                return SubsaturatedEnumType
            end
        else
            if iszero(δ)
                return (dδdt < FT(0)) ? WBFEnumType : SupersaturatedEnumType
            end
            return SupersaturatedEnumType
        end
    else
        if δi < FT(0)
            if δ > FT(0)
                return WBFEnumType
            elseif iszero(δ)
                return (dδdt < FT(0)) ? SubsaturatedEnumType : WBFEnumType
            else
                return SubsaturatedEnumType
            end
        else
            if iszero(δi)
                return (dδdt < FT(0)) ? WBFEnumType : SupersaturatedEnumType
            end
            return SupersaturatedEnumType
        end
    end
end
get_saturation_regime_enum_type(δ::FT, δi::FT, T::FT, T_freeze::FT) where {FT} =
    get_saturation_regime_enum_type(δ, δi, T < T_freeze)
get_saturation_regime_enum_type(q_vap::FT, q_sl::FT, q_si::FT, T::FT, T_freeze::FT; dδdt::FT = FT(0)) where {FT} =
    get_saturation_regime_enum_type(q_vap - q_sl, q_vap - q_si, T < T_freeze; dδdt = dδdt)

function get_saturation_regime(
    δ::FT,
    δi::FT,
    q_liq::FT,
    q_ice::FT,
    BF::Bool;
    dδdt::FT = FT(0),
) where {FT}
    regime_enum = get_saturation_regime_enum_type(δ, δi, BF; dδdt = dδdt)
    return add_regime_parameters(regime_enum, q_liq, q_ice, BF)
end

get_saturation_regime(
    q_vap::FT,
    q_liq::FT,
    q_ice::FT,
    q_sl::FT,
    q_si::FT,
    below_freezing::Bool;
    dδdt::FT = FT(0),
) where {FT} = get_saturation_regime(q_vap - q_sl, q_vap - q_si, q_liq, q_ice, below_freezing; dδdt = dδdt)

get_saturation_regime(
    q_vap::FT,
    q_liq::FT,
    q_ice::FT,
    q_sl::FT,
    q_si::FT,
    T::FT,
    T_freeze::FT;
    dδdt::FT = FT(0),
) where {FT} = get_saturation_regime(q_vap, q_liq, q_ice, q_sl, q_si, T < T_freeze; dδdt = dδdt)

get_saturation_regime(δ::FT, δi::FT, q_liq::FT, q_ice::FT, T::FT, T_freeze::FT; dδdt::FT = FT(0)) where {FT} =
    get_saturation_regime(δ, δi, q_liq, q_ice, T < T_freeze; dδdt = dδdt)

"""Initial `dδ/dt` from linear rates (MM2015PiecewiseLinear), used to break an exact-saturation tie."""
function get_dδdt_0(
    δ_0::FT,
    δ_0i::FT,
    q_liq::FT,
    q_ice::FT,
    τ_liq::FT,
    τ_ice::FT,
    dδdt_no_S::FT,
    below_freezing::Bool,
) where {FT}
    if below_freezing
        if δ_0 ≥ FT(0)
            S_ql_δ = (δ_0 / τ_liq)
            S_qi_δ = (δ_0i / τ_ice)
        elseif δ_0i ≥ FT(0)
            S_ql_δ = (q_liq > FT(0)) ? (δ_0 / τ_liq) : FT(0)
            S_qi_δ = (δ_0i / τ_ice)
        else
            S_ql_δ = (q_liq > FT(0)) ? (δ_0 / τ_liq) : FT(0)
            S_qi_δ = (q_ice > FT(0)) ? (δ_0i / τ_ice) : FT(0)
        end
    else
        if δ_0i ≥ FT(0)
            S_ql_δ = δ_0 / τ_liq
            S_qi_δ = FT(0)
        elseif δ_0 ≥ FT(0)
            S_ql_δ = δ_0 / τ_liq
            S_qi_δ = (q_ice > FT(0)) ? (δ_0i / τ_ice) : FT(0)
        else
            S_ql_δ = (q_liq > FT(0)) ? (δ_0 / τ_liq) : FT(0)
            S_qi_δ = (q_ice > FT(0)) ? (δ_0i / τ_ice) : FT(0)
        end
    end
    return dδdt_no_S - (S_ql_δ + S_qi_δ)
end

comparable_sign(x::FT, y::FT) where {FT} = (sign(x) == sign(y)) || iszero(x) || iszero(y)

function is_same_supersaturation_regime(δ_candidate::FT, δ_0::FT, δ_0i::FT) where {FT}
    return comparable_sign(δ_candidate, δ_0) && comparable_sign(δ_candidate + (δ_0i - δ_0), δ_0i)
end

"""Internal, prepared scalar inputs shared by the three solver kernels."""
struct MM2015Inputs{FT}
    g::FT
    L_i::FT
    L_l::FT
    c_p::FT
    T_freeze::FT
    dqsl_dT::FT
    dqsi_dT::FT
    e_sl::FT
    e_si::FT
    ρ::FT
    p::FT
    T::FT
    w::FT
    τ_liq::FT
    τ_ice::FT
    q_tot::FT
    q_liq::FT
    q_ice::FT
    q_vap_eq_liq::FT
    q_vap_eq_ice::FT
    Δt::FT
end

"""Internal copy constructor used by diagnostics to replace the timestep."""
function MM2015Inputs(inputs::MM2015Inputs{FT}; Δt::FT = inputs.Δt) where {FT}
    return MM2015Inputs(
        inputs.g,
        inputs.L_i,
        inputs.L_l,
        inputs.c_p,
        inputs.T_freeze,
        inputs.dqsl_dT,
        inputs.dqsi_dT,
        inputs.e_sl,
        inputs.e_si,
        inputs.ρ,
        inputs.p,
        inputs.T,
        inputs.w,
        inputs.τ_liq,
        inputs.τ_ice,
        inputs.q_tot,
        inputs.q_liq,
        inputs.q_ice,
        inputs.q_vap_eq_liq,
        inputs.q_vap_eq_ice,
        Δt,
    )
end

