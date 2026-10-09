# Exponential Rosenbrock method exprb43 of Hochbruck, Ostermann & Schweitzer (2009), SIAM J. Numer. Anal. 47, 786–803,
# in its non-autonomous form (their eqs. 6.5a, 6.5b), for systems of three unknowns.

"""3-vector."""
const Vec3{FT} = NTuple{3, FT}

"""3×3 matrix, stored by columns."""
const Mat3{FT} = NTuple{9, FT}

@inline function mat_mul(A::Mat3, B::Mat3)
    return ntuple(Val(9)) do k
        i, j = mod1(k, 3), (k - 1) ÷ 3
        muladd(A[i + 6], B[3j + 3], muladd(A[i + 3], B[3j + 2], A[i] * B[3j + 1]))
    end
end

@inline mat_vec(A::Mat3, x::Vec3) = (
    muladd(A[7], x[3], muladd(A[4], x[2], A[1] * x[1])),
    muladd(A[8], x[3], muladd(A[5], x[2], A[2] * x[1])),
    muladd(A[9], x[3], muladd(A[6], x[2], A[3] * x[1])),
)

"""Largest absolute column sum of `A`."""
@inline opnorm1(A::Mat3) =
    max(abs(A[1]) + abs(A[2]) + abs(A[3]), abs(A[4]) + abs(A[5]) + abs(A[6]), abs(A[7]) + abs(A[8]) + abs(A[9]))

"""All-ones factors: the identity similarity."""
@inline identity_scaling(::Type{FT}) where {FT} = (ntuple(_ -> one(FT), Val(9)), ntuple(_ -> one(FT), Val(9)))

"""
    balancing(A, start = identity_scaling(FT)) -> (to, from)

Entrywise factors `to[i, j] = d_j/d_i` and `from[i, j] = d_i/d_j` of the diagonal similarity
`B = D⁻¹ A D = A .* to`, `D = diag(d)` with powers of 2 that bring the off-diagonal 1-norms of each
row and column of `B` within a factor 4 of each other (Parlett & Reinsch 1969, Numer. Math. 13,
293–304); `f(A) = f(B) .* from` for a function `f` of the matrix. `D = I` when `B` does not have the
smaller 1-norm. The balance of `A` also balances every multiple of `A`. The iteration starts from
the factors `start`, such as those of a nearby matrix.
"""
@inline function balancing(A::Mat3{FT}, (to, from)::NTuple{2, Mat3{FT}} = identity_scaling(FT)) where {FT}
    d, d_inv = (one(FT), to[4], to[7]), (one(FT), from[4], from[7])
    changed = true
    while changed
        changed = false
        for i in 1:3
            j, k = mod1(i + 1, 3), mod1(i + 2, 3)
            # off-diagonal norms of column and row i of B, with B[p, q] = A[p, q] d_q/d_p
            c = (abs(A[3i - 3 + j]) * d_inv[j] + abs(A[3i - 3 + k]) * d_inv[k]) * d[i]
            r = (abs(A[3j - 3 + i]) * d[j] + abs(A[3k - 3 + i]) * d[k]) * d_inv[i]
            (iszero(c) || iszero(r) || !isfinite(c + r)) && continue
            n = (exponent(r) - exponent(c)) >> 1
            iszero(n) && continue
            f, f_inv = ldexp(one(FT), n), ldexp(one(FT), -n)
            c * f + r * f_inv < FT(0.95) * (c + r) || continue
            d = Base.setindex(d, d[i] * f, i)
            d_inv = Base.setindex(d_inv, d_inv[i] * f_inv, i)
            changed = true
        end
    end
    return similarity_factors(A, d, d_inv)
end

"""`(to, from)` of [`balancing`](@ref) from the diagonal `d` and its inverse, or ones unless they lower `‖A‖₁`."""
@inline function similarity_factors(A::Mat3{FT}, d::Vec3{FT}, d_inv::Vec3{FT}) where {FT}
    to = ntuple(m -> d[(m - 1) ÷ 3 + 1] * d_inv[mod1(m, 3)], Val(9))
    from = ntuple(m -> d[mod1(m, 3)] * d_inv[(m - 1) ÷ 3 + 1], Val(9))
    all(isfinite, to) && all(isfinite, from) && opnorm1(A .* to) < opnorm1(A) || return identity_scaling(FT)
    return to, from
end

"""`M` with `c` added to its diagonal."""
@inline add_diagonal(M::Mat3, c) = (M[1] + c, M[2], M[3], M[4], M[5] + c, M[6], M[7], M[8], M[9] + c)

"""Largest `‖Z‖₁` at which the Taylor series of `φ₄` is summed."""
const PHI_TAYLOR_RADIUS = 4

"""Number of Taylor terms of `φ₄` that reach rounding for an argument of norm `ν`: the least `m` with `ν^m 4!/(m + 4)! ≤ eps`."""
@inline function phi_taylor_terms(ν::FT) where {FT}
    m, power, denominator = 0, one(FT), one(FT)
    while power > eps(FT) * denominator
        m += 1
        power *= ν
        denominator *= m + 4
    end
    return m
end

"""Thresholds `θ_m = (eps (m + 4)!/4!)^{1/m}`, `m = 1, 2, …`, up to the first above `PHI_TAYLOR_RADIUS`: `m` terms reach rounding for norms up to `θ_m`."""
function phi_taylor_thresholds(::Type{FT}) where {FT}
    thresholds = FT[]
    while isempty(thresholds) || last(thresholds) ≤ PHI_TAYLOR_RADIUS
        m = length(thresholds) + 1
        push!(thresholds, FT((big(eps(FT)) * factorial(big(m + 4)) / 24)^(big(1) / m)))
    end
    return Tuple(thresholds)
end

const PHI_TAYLOR_THRESHOLDS_FLOAT32 = phi_taylor_thresholds(Float32)
const PHI_TAYLOR_THRESHOLDS_FLOAT64 = phi_taylor_thresholds(Float64)

@inline phi_taylor_terms(ν::Float32) = taylor_terms_from(PHI_TAYLOR_THRESHOLDS_FLOAT32, ν)
@inline phi_taylor_terms(ν::Float64) = taylor_terms_from(PHI_TAYLOR_THRESHOLDS_FLOAT64, ν)

"""Least `m` with `ν ≤ thresholds[m]`, or the number of thresholds."""
@inline function taylor_terms_from(thresholds::NTuple{N, FT}, ν::FT) where {N, FT}
    for m in 1:N
        ν ≤ thresholds[m] && return m
    end
    return N
end

"""`(c_j, c_{j-1})` of `c_j = 4!/(j + 4)!` from `c = c_j` for `j ≤ m`, and `(0, c)` for `j > m`."""
@inline taylor_coefficient(j::Int, m::Int, c) = j > m ? (zero(c), c) : (c, c * (j + 4))

"""
    phi4_taylor(Z, m) -> P

`P = Σ_{j=0}^{m} Z^j 4!/(j + 4)!`, so that `φ₄(Z) ≈ P/24`, by Paterson–Stockmeyer evaluation in
blocks of four powers (Paterson & Stockmeyer 1973, SIAM J. Comput. 2, 60–66).
"""
@inline function phi4_taylor(Z::Mat3{FT}, m::Int) where {FT}
    Z² = mat_mul(Z, Z)
    Z³ = mat_mul(Z², Z)
    Z⁴ = mat_mul(Z², Z²)
    denominator = one(FT)
    for i in 5:(m + 4)
        denominator *= i
    end
    c = inv(denominator)
    r = m ÷ 4
    P = ntuple(_ -> zero(FT), Val(9))
    for k in r:-1:0
        c₃, c = taylor_coefficient(4k + 3, m, c)
        c₂, c = taylor_coefficient(4k + 2, m, c)
        c₁, c = taylor_coefficient(4k + 1, m, c)
        c₀, c = taylor_coefficient(4k, m, c)
        block = add_diagonal(muladd.(c₃, Z³, muladd.(c₂, Z², c₁ .* Z)), c₀)
        P = k == r ? block : mat_mul(P, Z⁴) .+ block
    end
    return P
end

"""
    phi_functions(A, scaling = balancing(A)) -> (φ₀, φ₁, φ₂, φ₃, φ₄)

`φ_k(A) = Σ_j A^j / (j + k)!` for a 3×3 matrix `A`: with `B = A .* to` for `(to, from) = scaling`
from [`balancing`](@ref), a Taylor series of `φ₄` at `B / 2^s`, with `‖B‖₁ / 2^s ≤ 4`, followed by `s`
doublings [`phi_double`](@ref), and `φ_k(A) = φ_k(B) .* from`.
"""
function phi_functions(A::Mat3{FT}, (to, from)::NTuple{2, Mat3{FT}} = balancing(A)) where {FT}
    B = A .* to
    norm = opnorm1(B)
    s = norm > PHI_TAYLOR_RADIUS ? ceil(Int, log2(norm / PHI_TAYLOR_RADIUS)) : 0
    Z = B .* exp2(-FT(s))
    φ₄ = phi4_taylor(Z, phi_taylor_terms(norm * exp2(-FT(s)))) .* inv(FT(24))
    φ₃ = add_diagonal(mat_mul(Z, φ₄), inv(FT(6)))
    φ₂ = add_diagonal(mat_mul(Z, φ₃), inv(FT(2)))
    φ₁ = add_diagonal(mat_mul(Z, φ₂), one(FT))
    φ₀ = add_diagonal(mat_mul(Z, φ₁), one(FT))
    φ = (φ₀, φ₁, φ₂, φ₃, φ₄)
    for _ in 1:s
        φ = phi_double(φ)
    end
    return map(φ_k -> φ_k .* from, φ)
end

"""
    phi_double(φ(Z)) -> φ(2Z)

`φ_k(2Z) = 2^{−k} [φ₀(Z) φ_k(Z) + Σ_{j=1}^{k} φ_j(Z) / (k − j)!]` (Skaflestad & Wright 2009,
Appl. Numer. Math. 59, 783–799).
"""
@inline function phi_double((φ₀, φ₁, φ₂, φ₃, φ₄)::NTuple{5, Mat3{FT}}) where {FT}
    return (
        mat_mul(φ₀, φ₀),
        (mat_mul(φ₀, φ₁) .+ φ₁) ./ 2,
        (mat_mul(φ₀, φ₂) .+ φ₁ .+ φ₂) ./ 4,
        (mat_mul(φ₀, φ₃) .+ φ₁ ./ 2 .+ φ₂ .+ φ₃) ./ 8,
        (mat_mul(φ₀, φ₄) .+ φ₁ .* inv(FT(6)) .+ φ₂ ./ 2 .+ φ₃ .+ φ₄) ./ 16,
    )
end

"""
    exprb43_step(rhs, t, u, h, F, J, v, scaling) -> (Δu, error)

Increment `Δu` of one exprb43 step of length `h` for `u′ = rhs(t, u)` from `(t, u)`, with
`F = rhs(t, u)`, the Jacobian `J = ∂rhs/∂u`, `v = ∂rhs/∂t` at `(t, u)`, and `scaling = balancing(J)`.
`error` is the difference from the embedded third-order solution, `h φ₄(hJ) (12 D₃ − 48 D₂)`.
"""
@inline function exprb43_step(rhs::R, t::FT, u::Vec3{FT}, h::FT, F::Vec3{FT}, J::Mat3{FT}, v::Vec3{FT}, scaling::NTuple{2, Mat3{FT}}) where {R, FT}
    half = phi_functions(J .* (h / 2), scaling)
    _, φ₁, φ₂, φ₃, φ₄ = phi_double(half)
    ΔU₂ = (h / 2) .* mat_vec(half[2], F) .+ (h^2 / 4) .* mat_vec(half[3], v)
    D₂ = rhs(t + h / 2, u .+ ΔU₂) .- F .- mat_vec(J, ΔU₂) .- (h / 2) .* v
    base = h .* mat_vec(φ₁, F) .+ h^2 .* mat_vec(φ₂, v)
    ΔU₃ = base .+ h .* mat_vec(φ₁, D₂)
    D₃ = rhs(t + h, u .+ ΔU₃) .- F .- mat_vec(J, ΔU₃) .- h .* v
    φ₃D = mat_vec(φ₃, 16 .* D₂ .- 2 .* D₃)
    φ₄D = mat_vec(φ₄, 12 .* D₃ .- 48 .* D₂)
    return base .+ h .* (φ₃D .+ φ₄D), h .* φ₄D
end
