"""
    generate_assets.jl

Write the figures `docs/src/assets` holds.

    julia --project=docs/generate_assets docs/generate_assets/generate_assets.jl
"""

using CairoMakie: CairoMakie
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

const ASSETS = joinpath(dirname(@__DIR__), "src", "assets")

function _wbf_parcel(::Type{FT} = Float64) where {FT}
    thermo = MM2015.DefaultThermodynamicsBackend()
    T, p = FT(261), FT(8.0e4)
    q_liq, q_ice = FT(2.0e-4), FT(1.0e-4)
    ε = MM2015.molmass_ratio(thermo, FT)
    e_sl = MM2015.saturation_vapor_pressure_liq(thermo, T)
    e_si = MM2015.saturation_vapor_pressure_ice(thermo, T)
    r_sl = ε * e_sl / (p - e_sl)
    r_si = ε * e_si / (p - e_si)
    r_vap = (r_sl + r_si) / 2
    q_d = (one(FT) - q_liq - q_ice) / (one(FT) + r_vap)
    q_vap = q_d * r_vap
    q_tot = q_vap + q_liq + q_ice
    ρ = MM2015.air_density(thermo, T, p, q_tot, q_liq, q_ice)
    sc = MM2015.thermo_scalars(thermo, T, p, q_tot, q_liq, q_ice)
    inputs = MM2015.MM2015Inputs(
        sc.g, sc.L_i, sc.L_l, sc.c_p, sc.T_freeze,
        sc.dqsl_dT, sc.dqsi_dT, sc.e_sl, sc.e_si,
        ρ, p, T, FT(0.3), FT(8), FT(12), q_tot, q_liq, q_ice, sc.q_sl, sc.q_si, FT(40),
    )
    return inputs, thermo, MM2015.MM2015Opts{FT}()
end

function generate_assets()
    mkpath(ASSETS)
    inputs, thermo, opts = _wbf_parcel()
    MM2015.plot_rates(joinpath(ASSETS, "mm2015_rates.png"), inputs; thermo, opts)
    MM2015.plot_evolution(joinpath(ASSETS, "mm2015_evolution.png"), inputs; thermo, opts)
    MM2015.plot_condensate(joinpath(ASSETS, "mm2015_condensate.png"), inputs; thermo, opts)
    MM2015.plot_supersaturation(joinpath(ASSETS, "mm2015_supersaturation.png"), inputs; thermo, opts)
    MM2015.plot_specific_humidities(joinpath(ASSETS, "mm2015_specific_humidities.png"), inputs; thermo, opts)
    MM2015.plot_temperature(joinpath(ASSETS, "mm2015_temperature.png"), inputs; thermo, opts)
    println("wrote ", join(sort(readdir(ASSETS)), ", "), " to $ASSETS")
    return nothing
end

generate_assets()
