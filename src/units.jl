"""
Morrison & Milbrandt (2015) write mixing ratios. This package's public rates are
**specific humidities** .

With total specific humidity `q_tot` conserved under phase change:

- specific humidity → mixing ratio: `q / (1 - q_tot)`
- mixing ratio → specific humidity: `q_mix * (1 - q_tot)`

The same factor converts a phase-change *rate* when `q_tot` is constant.
"""

specific_humidity_to_mixing_ratio(q::FT, q_tot::FT) where {FT} = q / (1 - q_tot)
mixing_ratio_to_specific_humidity(q_mix::FT, q_tot::FT) where {FT} = q_mix * (1 - q_tot)

"""
    limit_δ(δ, q_vap, q_liq, q_ice)

Zero a supersaturation that cannot be realized by changing condensate or vapor at this
floating-point precision.
"""
function limit_δ(δ::FT, q_vap::FT, q_liq::FT, q_ice::FT) where {FT}
    if δ > FT(0)
        if (δ < (nextfloat(q_liq) - q_liq)) ||
           (δ < (nextfloat(q_ice) - q_ice)) ||
           (-δ > (prevfloat(q_vap) - q_vap))
            return FT(0)
        end
    elseif δ < FT(0)
        if (δ > (prevfloat(q_liq) - q_liq)) ||
           (δ > (prevfloat(q_ice) - q_ice)) ||
           (-δ < (nextfloat(q_vap) - q_vap))
            return FT(0)
        end
    end
    return δ
end
