# Depletion times and Lambert W against BigFloat references at a precision above the tested type.

using Test: Test
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :TestHelpers) || include(joinpath(@__DIR__, "test_helpers.jl"))
using .TestHelpers: allocated

const PRECISION = 256

"""BigFloat precision of the references for inputs of type `FT`."""
reference_precision(FT) = FT === BigFloat ? 2 * PRECISION : PRECISION

"""Exact f of the Float inputs in BigFloat."""
exact_f(t, Ā, K, τ, c) = Ā * t + K * (-expm1(-t / τ)) - c

"""Root of `f` in `(lo, hi)` by bisection, geometric while `hi/lo > 4`."""
function bisect(f, lo, hi)
    flo = f(lo)
    for _ in 1:(4 * precision(BigFloat))
        mid = (lo > 0 && hi / lo > 4) ? sqrt(lo * hi) : (lo + hi) / 2
        (mid == lo || mid == hi) && break
        fm = f(mid)
        if (fm > 0) == (flo > 0)
            lo, flo = mid, fm
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

"""Turning point of the exact `f`, or `Inf`."""
function turning_point(Ā, K, τ)
    iszero(K) && return big(Inf)
    r = -Ā * τ / K
    return 0 < r < 1 ? -τ * log(r) : big(Inf)
end

"""Smallest positive root of the exact equation, or `Inf`; `f` is monotone on each side of its turning point."""
function reference_root(Ā::BigFloat, K::BigFloat, τ::BigFloat, c::BigFloat)
    f(t) = exact_f(t, Ā, K, τ, c)
    start = τ * big(1e-300)
    function first_root(a, b)
        fa = f(a)
        if isfinite(b)
            return ((fa > 0) != (f(b) > 0)) ? bisect(f, a, b) : nothing
        end
        lo, t = a, max(2 * a, τ)
        while t < big(1e300)
            ((fa > 0) != (f(t) > 0)) && return bisect(f, lo, t)
            lo, t = t, 2 * t
        end
        return nothing
    end
    tstar = turning_point(Ā, K, τ)
    pieces = isfinite(tstar) ? ((start, tstar), (tstar, big(Inf))) : ((start, big(Inf)),)
    for (a, b) in pieces
        root = first_root(a, b)
        root === nothing || return root
    end
    return big(Inf)
end

"""Size of f at `t` under relative perturbations of the inputs (δ_0, A, B, c)."""
function input_scale(t, δ0, A, B, τ, c)
    û = t / τ
    em = -expm1(-û)
    return abs(δ0) * em + abs(A) * τ * (û - em) + abs(B) * τ * û + abs(c)
end

"""
`:ok`, `:missing`, `:spurious`, or `:inaccurate`. A root is correct when it matches the reference to
`rtol` or when `f` stays within the backward-error floor between it and the reference. `Inf` is correct
when no root exists or when the roots are a pair that a floor-sized perturbation removes.
"""
function classify(sol, ref, Ā, K, τ, c, δ0, A, B; ε, rtol)
    within_floor(t) = abs(exact_f(t, Ā, K, τ, c)) ≤ 64 * ε * input_scale(t, δ0, A, B, τ, c)
    tstar = turning_point(Ā, K, τ)
    if !isfinite(ref)
        return (!isfinite(sol) || within_floor(big(sol))) ? :ok : :spurious
    end
    if !isfinite(sol)
        pair = isfinite(tstar) && ref < tstar && (c > 0) != (Ā > 0)
        return (pair && within_floor(tstar)) ? :ok : :missing
    end
    abs(big(sol) - ref) ≤ rtol * ref && return :ok
    lo, hi = minmax(big(sol), ref)
    return (within_floor(big(sol)) && (!(lo < tstar < hi) || within_floor(tstar))) ? :ok : :inaccurate
end

function check_no_wbf(δ0::FT, A::FT, τ::FT, q_c::FT, Γ::FT; rtol) where {FT}
    sol = MM2015.get_t_out_of_q_no_WBF(δ0, A, τ, τ, q_c, Γ)
    ε = eps(FT)
    return setprecision(BigFloat, reference_precision(FT)) do
        δ0b, Ab, τb = big(δ0), big(A), big(τ)
        K = δ0b - Ab * τb
        c = -big(q_c) * (τb * big(Γ) / τb)
        ref = reference_root(Ab, K, τb, c)
        return classify(sol, ref, Ab, K, τb, c, δ0b, Ab, big(0); ε, rtol), sol, ref
    end
end

function check_wbf(δ0::FT, A::FT, τ::FT, q_c::FT, Γ::FT, q_sl::FT, q_si::FT; rtol) where {FT}
    sol = MM2015.get_t_out_of_q_WBF(δ0, A, τ, τ, q_c, Γ, q_sl, q_si)
    ε = eps(FT)
    return setprecision(BigFloat, reference_precision(FT)) do
        δ0b, Ab, τb = big(δ0), big(A), big(τ)
        B = (big(q_sl) - big(q_si)) / τb
        K = δ0b - Ab * τb
        c = -big(q_c) * (τb * big(Γ) / τb)
        ref = reference_root(Ab + B, K, τb, c)
        return classify(sol, ref, Ab + B, K, τb, c, δ0b, Ab, B; ε, rtol), sol, ref
    end
end

"""Deterministic SplitMix64 uniform in [0, 1)."""
function splitmix(seed::UInt64)
    state = Ref(seed)
    return function ()
        state[] = state[] * 0x2545F4914F6CDD1D + 0x9E3779B97F4A7C15
        z = state[]
        z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
        z = z ⊻ (z >> 31)
        return Float64(z >> 11) / Float64(UInt64(1) << 53)
    end
end

"""`(q_c*, t*)` at which the two depletion roots coalesce, or `nothing`."""
function tangent_qc(A, δ0, τ, Γ)
    K = δ0 - A * τ
    iszero(K) && return nothing
    r = -A * τ / K
    (0 < r < 1) || return nothing
    tstar = -τ * log(r)
    cstar = A * tstar + K * (-expm1(-tstar / τ))
    qcstar = -cstar / Γ
    return qcstar > 0 ? (qcstar, tstar) : nothing
end

counts(results) = (; (s => count(==(s), results) for s in (:ok, :missing, :spurious, :inaccurate))...)

const RTOL = Dict(Float64 => 1e-13, Float32 => 1e-5, BigFloat => 1e-74)

Test.@testset "Depletion times and Lambert W" begin
    setprecision(BigFloat, PRECISION) do
        Test.@testset "$FT p = û − (1 − e^{−û}) and e = 1 − e^{−û}" for FT in (Float32, Float64, BigFloat)
            ûs = FT.(exp10.(range(-8, log10(50); length = 400)))
            append!(ûs, (prevfloat(FT(0.5)), FT(0.5), nextfloat(FT(0.5)), prevfloat(one(FT)), one(FT), nextfloat(one(FT))))
            for û in ûs
                p, e = MM2015._depletion_pe(û)
                p_ref, e_ref = setprecision(BigFloat, reference_precision(FT)) do
                    e_big = -expm1(-big(û))
                    return big(û) - e_big, e_big
                end
                Test.@test abs(p - p_ref) ≤ 4 * eps(FT) * p_ref
                Test.@test abs(e - e_ref) ≤ 4 * eps(FT) * e_ref
            end
        end

        Test.@testset "$FT no-WBF grid" for FT in (Float64, Float32)
            magnitudes = exp10.(range(-9, -1; length = 9))
            signed = vcat(-reverse(magnitudes), 0.0, magnitudes)
            results = Symbol[]
            for τ in (0.5, 8.67, 100.0), A in signed, δ0 in signed, q_c in exp10.(range(FT === Float32 ? -20 : -30, -2; length = 8))
                push!(results, first(check_no_wbf(FT(δ0), FT(A), FT(τ), FT(q_c), FT(1.25); rtol = RTOL[FT])))
            end
            c = counts(results)
            Test.@test c.missing == 0
            Test.@test c.spurious == 0
            Test.@test c.inaccurate == 0
        end

        Test.@testset "$FT WBF grid" for FT in (Float64, Float32)
            magnitudes = exp10.(range(-9, -1; length = 7))
            signed = vcat(-reverse(magnitudes), 0.0, magnitudes)
            results = Symbol[]
            q_sl = FT(5e-3)
            for τ in (0.5, 8.67, 100.0), A in signed, δ0 in signed, q_c in exp10.(range(FT === Float32 ? -20 : -30, -2; length = 5)), B in (2e-4, 3e-3)
                q_si = q_sl - FT(B * τ)
                push!(results, first(check_wbf(FT(δ0), FT(A), FT(τ), FT(q_c), FT(1.25), q_sl, q_si; rtol = RTOL[FT])))
            end
            c = counts(results)
            Test.@test c.missing == 0
            Test.@test c.spurious == 0
            Test.@test c.inaccurate == 0
        end

        Test.@testset "$FT fuzz over wide log-uniform ranges" for (FT, n) in ((Float64, 2000), (Float32, 2000), (BigFloat, 200))
            u01 = splitmix(0x9E3779B97F4A7C15)
            sgn() = u01() < 0.5 ? -1.0 : 1.0
            results = Symbol[]
            for _ in 1:n
                τ = 10.0^(-1 + 6 * u01())
                A = sgn() * 10.0^(-9 + 6 * u01())
                δ0 = sgn() * 10.0^(-9 + 6 * u01())
                q_c = 10.0^((FT === Float32 ? -20 : -30) + (FT === Float32 ? 18 : 28) * u01())
                Γ = 1 + 2 * u01()
                push!(results, first(check_no_wbf(FT(δ0), FT(A), FT(τ), FT(q_c), FT(Γ); rtol = RTOL[FT])))
            end
            c = counts(results)
            Test.@test c.missing == 0
            Test.@test c.spurious == 0
            Test.@test c.inaccurate == 0
        end

        Test.@testset "$FT near-tangent coalescing roots at slow relaxation" for FT in (Float64, Float32)
            failures = []
            n = 0
            for τ in (1e3, 1e5, 1e7), A in (1e-6, 3.5e-6, 1e-5), δ0 in (-1e-5, -2e-5, -5e-5)
                tg = tangent_qc(A, δ0, τ, 1.58)
                tg === nothing && continue
                for mult in (0.5, 0.9, 0.99, 0.999, 0.9999, 1 - 1e-6, 1 - 1e-9, 1 - 1e-12, 1 + 1e-12, 1 + 1e-9, 1 + 1e-6, 1.0001, 1.01, 1.1)
                    status, sol, ref = check_no_wbf(FT(δ0), FT(A), FT(τ), FT(tg[1] * mult), FT(1.58); rtol = RTOL[FT])
                    n += 1
                    status == :ok || push!(failures, (; status, τ, A, δ0, mult, sol, ref = Float64(ref)))
                end
            end
            isempty(failures) || @info "near-tangent failures" FT failures
            Test.@test n > 0
            Test.@test isempty(failures)
        end

        Test.@testset "named cases" begin
            status, sol, _ = check_no_wbf(-0.0007554169251637422, 9.659082385976545e-7, 8.668880477547646, 3.944304526105059e-31, 1.2458987067041472; rtol = 1e-13)
            Test.@test status == :ok
            Test.@test 0 < sol < 1e-20
            for (δ0, A, τ, q_c, Γ) in ((-0.00387, -1.4116656937688102e-6, 1.0, 1.3e-4, 2.79), (-0.003, -5.689866029018293e-6, 0.5, 5e-4, 2.0))
                Test.@test first(check_no_wbf(δ0, A, τ, q_c, Γ; rtol = 1e-13)) == :ok
            end
            q_sl, q_si, τ = 0.00390625, 0.001953125, 0.5
            δ0 = -(q_sl - q_si)
            for A in (-0.001, -0.02, 0.001), q_c in (1e-18, 1e-10, 1e-4)
                Test.@test first(check_wbf(δ0, A, τ, q_c, 1.25, q_sl, q_si; rtol = 1e-13)) == :ok
            end
            for A in (1e-6, 1e-3), δ0 in (1e-4, 6e-3), τ in (0.5, 8.67), q_c in (1e-6, 1e-3)
                status, sol, _ = check_no_wbf(δ0, A, τ, q_c, 1.25; rtol = 1e-13)
                Test.@test status == :ok
                Test.@test isinf(sol)
            end
        end

        Test.@testset "$FT Lambert W at piece boundaries, branch point, and extremes" for FT in (Float64, Float32, BigFloat)
            ref_w0(z) = setprecision(BigFloat, reference_precision(FT)) do
                z == 0 ? big(0) : bisect(w -> w * exp(w) - z, z < 0 ? big(-1) : big(0), z < 0 ? big(0) : max(big(1), log1p(big(z))))
            end
            ref_wm1(z) = setprecision(BigFloat, reference_precision(FT)) do
                bisect(w -> -(w * exp(w) - z), min(big(-3), 2 * log(-big(z))), big(-1))
            end
            # Float64 pieces lose relative accuracy toward W₀ = 0 (up to 51ε at |z| = 0.01); the iterative path is bounded by the condition number 1/|1 + W|
            tolerance = Dict(Float64 => 64 * eps(Float64), Float32 => Float64(eps(Float32)), BigFloat => 8 * eps(BigFloat))[FT]
            bound(w) = tolerance * abs(w) * (FT === BigFloat ? max(1, inv(abs(1 + w))) : 1)
            inv_e = exp(-big(1))
            w0_bounds = (2.5498939065034735716, 43.613924462669367895, 598.45353371878276946, 8049.4919850757619109,
                111124.95412121781420, 1.5870429812082297112e6, 2.3414708401875459509e7, 3.5576474308009965225e8,
                5.5501716296163627854e9, 8.8674704839657775331e10, 1.4477791865272902816e12, 2.4111458632511851931e13,
                4.0897036442600845564e14, 7.0555901476789972402e15, 1.2366607557976727287e17, 2.1999373487930999775e18,
                3.9685392198344016155e19)
            zs0 = FT[]
            for b in w0_bounds, s in (1 - 1e-12, 1.0, 1 + 1e-12)
                push!(zs0, FT(big(b * s) - inv_e))
            end
            append!(zs0, FT.((-inv_e + big(10.0)^-k for k in (15, 12, 9, 6, 3))))
            append!(zs0, FT.((0.01 - 1e-12, 0.01, 0.01 + 1e-12, -0.01, -0.01 - 1e-12, 1e-8, -1e-8, 0.5, 3.0, 1e3, 1e30)))
            append!(zs0, FT.((-0.3, -0.3054, -0.2, -0.1, -0.05, -0.020, -0.010356505380278628, 0.015, 0.05, 0.2, 1.0, 2.0)))
            for z in filter(z -> z > -inv_e, zs0)
                w = ref_w0(z)
                Test.@test abs(MM2015.fast_lambertw0(z) - w) ≤ max(bound(w), tolerance * big(floatmin(Float32)))
            end
            wm1_bounds = (-0.3542913309442164, -0.18872688282289434049, -0.060497597226958343647, -0.017105334740676008194,
                -0.0045954962127943706433, -0.0012001610672197724173, -0.00030728805932191499844, -0.000077447159838062184354,
                -4.5808119698158173174e-17)
            zsm = FT[]
            for b in wm1_bounds, s in (1 - 1e-12, 1.0, 1 + 1e-12)
                push!(zsm, FT(b * s))
            end
            append!(zsm, FT.((-inv_e + big(10.0)^-k for k in (15, 12, 9, 6, 3))))
            append!(zsm, FT.((-0.3, -1e-5, -1e-30)))
            for z in filter(z -> -inv_e < z < 0, zsm)
                w = ref_wm1(z)
                Test.@test abs(MM2015.fast_lambertwm1(z) - w) ≤ bound(w)
            end
            for l in (-1.0, -1 - 1e-12, -1 - 1e-7, -1.001, -1.037635735532508, -1.03763573553251, -1.04, -1.5, -9.46, -9.48,
                -37.6, -37.65, -180.09, -180.1, -700.0, -745.0, -746.0, -1e4, -1e8)
                lnmz = FT(l)
                w = setprecision(() -> ref_wm1(-exp(big(lnmz))), BigFloat, reference_precision(FT))
                Test.@test abs(MM2015.fast_lambertwm1_from_ln(lnmz) - w) ≤ bound(w)
            end
        end

        Test.@testset "$FT inferred and allocation-free" for FT in (Float32, Float64)
            args_nw = (FT(-1e-5), FT(2e-6), FT(8), FT(8), FT(1e-6), FT(1.25))
            args_wbf = (args_nw..., FT(5e-3), FT(4.8e-3))
            Test.@test Test.@inferred(MM2015.get_t_out_of_q_no_WBF(args_nw...)) isa FT
            Test.@test Test.@inferred(MM2015.get_t_out_of_q_WBF(args_wbf...)) isa FT
            Test.@test allocated(MM2015.get_t_out_of_q_no_WBF, args_nw...) == 0
            Test.@test allocated(MM2015.get_t_out_of_q_WBF, args_wbf...) == 0
            for (f, z) in ((MM2015.fast_lambertw0, FT(0.3)), (MM2015.fast_lambertwm1, FT(-0.1)), (MM2015.fast_lambertwm1_from_ln, FT(-50)))
                Test.@test Test.@inferred(f(z)) isa FT
                Test.@test allocated(f, z) == 0
            end
        end
    end
end
