#=
Principled tests for the MM2015 exponential-relaxation depletion-time root solvers
`get_t_out_of_q_no_WBF` / `get_t_out_of_q_WBF` (in ../morrison_milbrandt_2015_style.jl).

Both solvers return the SMALLEST POSITIVE root t of
    no-WBF:  f(t) = A_c·t          + (δ_0 − A_c·τ)(1 − e^{−t/τ}) − c = 0
    WBF:     f(t) = (A_c + B)·t     + (δ_0 − A_c·τ)(1 − e^{−t/τ}) − c = 0
(or +∞ if no positive root exists), where c = −q_c·(τ_c·Γ/τ) and, for WBF, B = (q_sl−q_si)/τ·CF.
Physically t is the first time the evaporating/subliming condensate `q_c` is exhausted.

The reference is INDEPENDENT of the solver and correct by construction. Writing the equation as
    f(t) = Ā·t + K·(1 − e^{−t/τ}) − c ,   Ā = linear coeff (A_c, or A_c+B),  K = δ_0 − A_c·τ ,
we have f'(t) = Ā + (K/τ)·e^{−t/τ}, which is MONOTONIC in t, so f has at most one turning point and
hence at most two positive roots. `_smallest_positive_root` locates the (at most one) turning point
analytically, then brackets + bisects each monotone piece — this cannot miss the smallest root and
cannot return a spurious one. It is a test reference, so its speed is irrelevant.

Evaluation note: for t ≪ τ the terms Ā·t and K(1−e^{−t/τ}) nearly cancel (their sum is β·t/τ + O(t²)
with β = τ·f'(0) possibly ≪ |Ā·τ|), so the reference evaluates f in the cancellation-free SERIES form
with β supplied EXACTLY from the inputs (β = δ_0, resp. δ_0 + B·τ) — same accurate-evaluation identity
the production code documents in docs/src/depletion.md. The root-FINDING (bracket + bisect) shares nothing
with the closed-form solver under test.

Run via `test/runtests.jl`.
=#

using Test
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

const FT = Float64

# ---------------------------------------------------------------------------------------------------
# The equation (cancellation-free evaluation) and the independent reference root-finder
# ---------------------------------------------------------------------------------------------------
function _f(t, Ā, K, τ, c, β)
    û = t / τ
    if û < 0.05   # series form in exact β (Ā·t and K(1−e^{−t/τ}) cancel to β·û + O(û²) here)
        p = û^2 / 2 - û^3 / 6 + û^4 / 24 - û^5 / 120
        return β * û - K * p - c
    end
    return Ā * t + K * (-expm1(-t / τ)) - c
end

# bisection on [lo,hi] with a known sign change; returns the root to machine precision
function _bisect(Ā, K, τ, c, β, lo, hi)
    flo = _f(lo, Ā, K, τ, c, β)
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        fm = _f(mid, Ā, K, τ, c, β)
        if (fm > 0) == (flo > 0)
            lo = mid
            flo = fm
        else
            hi = mid
        end
    end
    (lo + hi) / 2
end

# smallest root of f on a single monotone piece [a,b] (b may be Inf); returns the root or `nothing`
function _root_on_piece(Ā, K, τ, c, β, a, b)
    fa = _f(a, Ā, K, τ, c, β)
    fa == 0 && a > 0 && return a
    if isfinite(b)
        fb = _f(b, Ā, K, τ, c, β)
        ((fa > 0) != (fb > 0)) && return _bisect(Ā, K, τ, c, β, a, b)
        return nothing
    else
        t = max(a, τ)
        t = t == 0 ? τ : t
        for _ in 1:400
            t *= 2
            ft = _f(t, Ā, K, τ, c, β)
            ((fa > 0) != (ft > 0)) && return _bisect(Ā, K, τ, c, β, a, t)
            t > 1e300 && break
        end
        return nothing
    end
end

"""Smallest positive root of Ā·t + K·(1−e^{−t/τ}) − c (β = τ·f'(0) exact), or `Inf` if none."""
function _smallest_positive_root(Ā, K, τ, c, β)
    # turning point t*: f'(t*) = Ā + (K/τ)e^{−t*/τ} = 0 ⇒ e^{−t*/τ} = −Āτ/K ∈ (0,1)
    tstar = Inf
    if K != 0
        r = -Ā * τ / K
        (0 < r < 1) && (tstar = -τ * log(r))
    end
    if isfinite(tstar) && tstar > 0
        r1 = _root_on_piece(Ā, K, τ, c, β, 0.0, tstar)   # piece [0, t*]
        !isnothing(r1) && return r1
        r2 = _root_on_piece(Ā, K, τ, c, β, tstar, Inf)   # piece [t*, ∞)
        return isnothing(r2) ? Inf : r2
    else
        r = _root_on_piece(Ā, K, τ, c, β, 0.0, Inf)       # single monotone piece
        return isnothing(r) ? Inf : r
    end
end

# Tangency constructor: given (A, δ0, τ, τc, Γ) with an interior turning point t* (a minimum of
# bracket(t)=A·t+K(1−e^{−t/τ})), the two depletion roots COALESCE when c = bracket(t*) — i.e. at a
# specific q_c*. For q_c < q_c* there are two (increasingly close) roots; q_c = q_c* is the double
# root; q_c > q_c* has none. Sweeping q_c across q_c* walks the near-tangent manifold that a
# fixed-decade grid never lands on (it is measure-zero in raw parameter space). Returns (q_c*, t*)
# or `nothing` if there is no interior turning point / q_c* ≤ 0.
function _tangent_qc(A, δ0, τ, τc, Γ)
    K = δ0 - A * τ
    K == 0 && return nothing
    r = -A * τ / K
    (0 < r < 1) || return nothing
    tstar = -τ * log(r)
    cstar = A * tstar + K * (-expm1(-tstar / τ))     # = bracket(t*); tangency at c = cstar
    qcstar = -cstar * τ / (τc * Γ)                    # c = −q_c·(τc·Γ/τ)  ⇒  q_c* = −cstar·τ/(τc·Γ)
    return qcstar > 0 ? (qcstar, tstar) : nothing
end

# relative agreement between solver `s` and reference `o` (both may be Inf = "no root")
function _agree(s, o; rtol = 1e-6, atol = 1e-12)
    (isinf(o) || o > 1e290) && return !isfinite(s) || s > 1e12        # reference: no (finite) root
    isfinite(s) || return false                                      # reference has a root, solver said none
    abs(s - o) ≤ rtol * o + atol
end

# ---------------------------------------------------------------------------------------------------
# Coverage grids (both signs & many magnitudes of every parameter, incl. the degenerate limits)
# ---------------------------------------------------------------------------------------------------
const ΔBOTH = vcat(-exp10.(range(-9, -1, length = 17)), FT(0), exp10.(range(-9, -1, length = 17)))
const A_BOTH = vcat(-exp10.(range(-9, -1, length = 17)), FT(0), exp10.(range(-9, -1, length = 17)))
const TAUS = FT[0.5, 1.0, 8.67, 100.0]
const QCS = exp10.(range(FT(-31), FT(-2), length = 20))         # 30 decades incl. the tiny-q_c crash regime

@testset "MM2015 depletion-time root solvers" begin

    @testset "no-WBF: equals the reference's smallest root (or Inf) across the grid" begin
        nfail = 0
        worst = 0.0
        example = nothing
        for τ in TAUS, A in A_BOTH, δ0 in ΔBOTH, qc in QCS
            Γ = FT(1.25)
            τc = τ
            c = -qc * (τc * Γ / τ)
            ref = _smallest_positive_root(A, δ0 - A * τ, τ, c, δ0)          # β = δ0 exactly
            sol = MM2015.get_t_out_of_q_no_WBF(FT(δ0), FT(A), FT(τ), FT(τc), FT(qc), FT(Γ))
            if !_agree(sol, ref)
                nfail += 1
                isnothing(example) && (example = (; δ0, A, τ, qc, ref, sol))
            elseif isfinite(ref) && ref > 0
                worst = max(worst, abs(sol - ref) / ref)
            end
        end
        !isnothing(example) && @info "no-WBF first mismatch" example
        @test nfail == 0
        @test worst < 1e-6
    end

    @testset "WBF: equals the reference's smallest root (or Inf) across the grid" begin
        nfail = 0
        worst = 0.0
        example = nothing
        q_sl = FT(5e-3)
        for τ in TAUS, A in A_BOTH, δ0 in ΔBOTH, qc in QCS[1:2:end], B in FT[2e-4, 1e-3, 3e-3]
            Γ = FT(1.25)
            τc = τ
            c = -qc * (τc * Γ / τ)
            q_si = q_sl - B * τ
            # the equation the solver is ASKED to solve is defined by its actual inputs (q_sl, q_si):
            # it reconstructs B = (q_sl−q_si)/τ, which differs from the synthetic B above at the last
            # ulp (q_si round-trip). The reference must use the same effective B or it scores the
            # solver against a slightly different equation (visible only at β ≈ 0 with tiny q_c).
            B_eff = (q_sl - q_si) / τ
            ref = _smallest_positive_root(A + B_eff, δ0 - A * τ, τ, c, δ0 + B_eff * τ)
            sol = MM2015.get_t_out_of_q_WBF(FT(δ0), FT(A), FT(τ), FT(τc), FT(qc), FT(Γ), q_sl, q_si)
            if !_agree(sol, ref)
                nfail += 1
                isnothing(example) && (example = (; δ0, A, B, τ, qc, ref, sol))
            elseif isfinite(ref) && ref > 0
                worst = max(worst, abs(sol - ref) / ref)
            end
        end
        !isnothing(example) && @info "WBF first mismatch" example
        @test nfail == 0
        @test worst < 1e-6
    end

    # -------- named edge cases (each historically produced a wrong root; keep them explicit) --------

    @testset "RF09 crash regression: tiny q_liq, subsaturated (was the spurious ~790 s root)" begin
        δ0, A, τ, qc, Γ = -0.0007554169251637422, 9.659082385976545e-7, 8.668880477547646,
            8.668880477547646, 3.944304526105059e-31, 1.2458987067041472
        c = -qc * (τ * Γ / τ)
        sol = MM2015.get_t_out_of_q_no_WBF(FT(δ0), FT(A), FT(τ), FT(τ), FT(qc), FT(Γ))
        ref = _smallest_positive_root(A, δ0 - A * τ, τ, c, δ0)
        @test _agree(sol, ref)
        @test 0 < sol < 1e-20              # the near-0 depletion root (~5.6e-27), NOT ~790
    end

    @testset "weak-forcing overflow (E≳2500) and near-overflow (E≈709)" begin
        for (δ0, A, τ, qc, Γ) in (
            (-0.00387, -1.4116656937688102e-6, 1.0, 1.3e-4, 2.79),
            (-0.003, -5.689866029018293e-6, 0.5, 5e-4, 2.0),
        )
            c = -qc * (τ * Γ / τ)
            sol = MM2015.get_t_out_of_q_no_WBF(FT(δ0), FT(A), FT(τ), FT(τ), FT(qc), FT(Γ))
            ref = _smallest_positive_root(A, δ0 - A * τ, τ, c, δ0)
            @test _agree(sol, ref)
        end
    end

    @testset "WBF β = δ0 + Bτ = 0 exactly (linear term vanishes; root is quadratic-dominant)" begin
        # binary-exact values so β = 0 survives the wrapper's B = (q_sl−q_si)/τ reconstruction exactly:
        # q_sl = 2⁻⁸, q_si = 2⁻⁹ ⇒ q_sl−q_si = 2⁻⁹ exact; τ = 0.5 ⇒ B = 2⁻⁸, Bτ = 2⁻⁹; δ0 = −2⁻⁹ ⇒ β ≡ 0
        q_sl, q_si, τ = 0.00390625, 0.001953125, 0.5
        B = (q_sl - q_si) / τ
        δ0 = -B * τ
        for A in (-0.001, -0.02, 0.001), qc in (1e-18, 1e-10, 1e-4)
            Γ = 1.25
            c = -qc * (τ * Γ / τ)
            sol = MM2015.get_t_out_of_q_WBF(FT(δ0), FT(A), FT(τ), FT(τ), FT(qc), FT(Γ), FT(q_sl), FT(q_si))
            ref = _smallest_positive_root(A + B, δ0 - A * τ, τ, c, δ0 + B * τ)
            @test _agree(sol, ref)
        end
    end

    @testset "no depletion root ⇒ Inf (supersaturated growth, condensate never exhausts)" begin
        # δ0 > 0 and equilibrium A_c·τ ≥ 0: the phase grows, so f(t) > 0 ∀ t>0 ⇒ reference = Inf
        for A in (1e-6, 1e-3), δ0 in (1e-4, 6e-3), τ in (0.5, 8.67), qc in (1e-6, 1e-3)
            Γ = 1.25
            c = -qc * (τ * Γ / τ)
            ref = _smallest_positive_root(A, δ0 - A * τ, τ, c, δ0)
            sol = MM2015.get_t_out_of_q_no_WBF(FT(δ0), FT(A), FT(τ), FT(τ), FT(qc), FT(Γ))
            @test _agree(sol, ref)
        end
    end

    # -------- structural coverage of the fragile regimes (constructed, not decade-sampled) --------

    @testset "near-tangent coalescing roots at stiff τ (the class the fixed-decade grid misses)" begin
        # Walk q_c across the tangency q_c* at stiff τ: the two depletion roots crowd together
        # (discriminant β²−2Kc → 0). The reference tracks the smallest root; the solver must agree.
        # A spurious Inf on the two-root side (q_c < q_c*) is the missed-near-tangent-root bug.
        nfail = 0
        example = nothing
        for τ in (FT(1e3), FT(1e5), FT(1e6), FT(1e7)),
            A in (FT(1e-6), FT(3.5e-6), FT(1e-5)),
            δ0 in (FT(-1e-5), FT(-2e-5), FT(-5e-5))

            Γ = FT(1.58)
            τc = τ
            tg = _tangent_qc(A, δ0, τ, τc, Γ)
            isnothing(tg) && continue
            qcstar = tg[1]
            for mult in (FT(0.5), FT(0.9), FT(0.99), FT(0.999), FT(0.9999), FT(1.0001), FT(1.01), FT(1.1))
                qc = qcstar * mult
                (qc > 0) || continue
                c = -qc * (τc * Γ / τ)
                ref = _smallest_positive_root(A, δ0 - A * τ, τ, c, δ0)
                sol = MM2015.get_t_out_of_q_no_WBF(FT(δ0), FT(A), FT(τ), FT(τc), FT(qc), FT(Γ))
                if !_agree(sol, ref)
                    nfail += 1
                    isnothing(example) && (example = (; A, δ0, τ, qc, qcstar, mult, ref, sol))
                end
            end
        end
        !isnothing(example) && @info "near-tangent first mismatch" example
        @test nfail == 0
    end

    @testset "deterministic fuzz vs reference over wide log-uniform ranges (incl. stiff τ)" begin
        # dependency-free reproducible PRNG (SplitMix64) — no Random stdlib dep needed; UInt64
        # arithmetic wraps in Julia. Differential test: solver must match the independent reference
        # everywhere, which occasionally lands near tangencies the structured grid also targets.
        st = Ref(UInt64(0x9E3779B97F4A7C15))
        u01 = function ()
            st[] = st[] * 0x2545F4914F6CDD1D + 0x9E3779B97F4A7C15
            z = st[]
            z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
            z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
            z = z ⊻ (z >> 31)
            return Float64(z >> 11) / Float64(UInt64(1) << 53)   # 53-bit mantissa in [0,1)
        end
        sgn = () -> (u01() < 0.5 ? -1.0 : 1.0)
        nfail = 0
        example = nothing
        for _ in 1:5000
            τ = FT(10.0)^FT(-1 + 9 * u01())                      # 1e-1 .. 1e8
            A = FT(sgn() * 10.0^(-9 + 6 * u01()))                # ±(1e-9 .. 1e-3)
            δ0 = FT(sgn() * 10.0^(-9 + 6 * u01()))
            qc = FT(10.0)^FT(-30 + 28 * u01())                   # 1e-30 .. 1e-2
            Γ = FT(1 + 2 * u01())
            τc = τ
            c = -qc * (τc * Γ / τ)
            ref = _smallest_positive_root(A, δ0 - A * τ, τ, c, δ0)
            sol = MM2015.get_t_out_of_q_no_WBF(δ0, A, τ, τc, qc, Γ)
            if !_agree(sol, ref)
                nfail += 1
                isnothing(example) && (example = (; A, δ0, τ, qc, Γ, ref, sol))
            end
        end
        !isnothing(example) && @info "fuzz first mismatch" example
        @test nfail == 0
    end
end

@testset "Lambert branches preserve input precision" begin
    for T in (Float32, BigFloat)
        for z0 in (-0.3, -1e-4, 1e-8, 0.2, 10.0)
            z = T(z0)
            w = MM2015.fast_lambertw0(z)
            tolerance = T === Float32 ? T(2e-6) : T(128) * eps(T)
            @test isapprox(w * exp(w), z; rtol = tolerance, atol = tolerance * max(abs(z), one(T)))
        end
        for z0 in (-0.3, -0.1, -1e-4, -1e-20)
            z = T(z0)
            w = MM2015.fast_lambertwm1(z)
            tolerance = T === Float32 ? T(3e-6) : T(256) * eps(T)
            @test w ≤ -one(T)
            @test isapprox(w * exp(w), z; rtol = tolerance, atol = tolerance * abs(z))
        end
    end
end
