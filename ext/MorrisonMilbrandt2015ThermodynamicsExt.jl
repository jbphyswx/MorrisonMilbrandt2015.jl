module MorrisonMilbrandt2015ThermodynamicsExt

using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015
using Thermodynamics: Thermodynamics as TD, Parameters as TP


const TD_VERSION = pkgversion(TD)

# AbstractThermodynamicsParameters was introduced in TD 0.12.6.
const APS = if TD_VERSION ≥ v"0.12.6"
    TP.AbstractThermodynamicsParameters
else
    TP.ThermodynamicsParameters
end

MM2015.R_d(ps::APS, ::Type{FT}) where {FT} = FT(TP.R_d(ps))
MM2015.R_v(ps::APS, ::Type{FT}) where {FT} = FT(TP.R_v(ps))
MM2015.cp_d(ps::APS, ::Type{FT}) where {FT} = FT(TP.cp_d(ps))
MM2015.cp_v(ps::APS, ::Type{FT}) where {FT} = FT(TP.cp_v(ps))
MM2015.cp_l(ps::APS, ::Type{FT}) where {FT} = FT(TP.cp_l(ps))
MM2015.cp_i(ps::APS, ::Type{FT}) where {FT} = FT(TP.cp_i(ps))
MM2015.grav(ps::APS, ::Type{FT}) where {FT} = FT(TP.grav(ps))
MM2015.T_triple(ps::APS, ::Type{FT}) where {FT} = FT(TP.T_triple(ps))

_td_phase(::MM2015.Liquid) = TD.Liquid()
_td_phase(::MM2015.Ice) = TD.Ice()

MM2015.saturation_vapor_pressure(ps::APS, T::FT, phase::MM2015.AbstractPhase) where {FT} =
    FT(TD.saturation_vapor_pressure(ps, T, _td_phase(phase)))

MM2015.latent_heat(ps::APS, T::FT, ::MM2015.Liquid) where {FT} = FT(TD.latent_heat_vapor(ps, T))
MM2015.latent_heat(ps::APS, T::FT, ::MM2015.Ice) where {FT} = FT(TD.latent_heat_sublim(ps, T))

# tendencies and trajectory of every scheme, basis, and float type with parameter sets of either float type
for FT in (Float32, Float64), FT_parameters in (Float32, Float64), Basis in (MM2015.SpecificHumidity, MM2015.DryAirMixingRatio)
    P = MM2015.MM2015Problem{FT, Basis, TP.ThermodynamicsParameters{FT_parameters}}
    for S in (MM2015.MM2015PiecewiseLinear{Float64}, MM2015.MM2015FixedT{Float64}, MM2015.MM2015{FT, MM2015.BrentRootFinder})
        precompile(MM2015.tendencies, (S, P, FT))
        precompile(MM2015.trajectory, (S, P, FT))
    end
end

end # module
