# tendencies and trajectory of every scheme, basis, and float type with the default backend
for FT in (Float32, Float64)
    for S in (MM2015PiecewiseLinear{Float64}, MM2015FixedT{Float64})
        precompile(tendencies, (S, Coefficients{FT}, FT))
        precompile(trajectory, (S, Coefficients{FT}, FT))
    end
    for Basis in (SpecificHumidity, DryAirMixingRatio)
        P = MM2015Problem{FT, Basis, DefaultThermodynamicsBackend}
        for S in (MM2015PiecewiseLinear{Float64}, MM2015FixedT{Float64}, MM2015{FT, BrentRootFinder})
            precompile(tendencies, (S, P, FT))
            precompile(trajectory, (S, P, FT))
        end
    end
end
