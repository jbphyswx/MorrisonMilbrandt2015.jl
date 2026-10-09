"""Root finder that brackets the events of [`MM2015`](@ref)."""
abstract type AbstractRootFinder end

"""Brent–Dekker iteration with a bisection safeguard."""
struct BrentRootFinder <: AbstractRootFinder end

"""
    RootSolversRootFinder(method)

Root finder of RootSolvers.jl with the bracketing method `method`, a type such as
`RootSolvers.BrentsMethod`. Requires `using RootSolvers`.
"""
struct RootSolversRootFinder{M} <: AbstractRootFinder end

RootSolversRootFinder(method::Type) = RootSolversRootFinder{method}()

"""
    bracket_root(finder, f, lo, hi, f_lo, f_hi, tolerance) -> (lo, hi)

Bracket around a sign change of `f`, from `f(lo) = f_lo ≤ 0` and `f(hi) = f_hi > 0`, of width at most
`tolerance + 4 eps max(|lo|, |hi|)`. The returned points keep those signs.
"""
function bracket_root end

function bracket_root(::BrentRootFinder, f::F, lo::FT, hi::FT, f_lo::FT, f_hi::FT, tolerance::FT) where {F, FT}
    a, b, c = lo, hi, hi
    fa, fb, fc = f_lo, f_hi, f_hi
    d = e = b - a
    plateau = false
    for _ in 1:(4 * precision(FT))
        if (fb > 0) == (fc > 0)
            c, fc = a, fa
            d = e = b - a
        end
        if abs(fc) < abs(fb)
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        end
        tol = 2 * eps(FT) * abs(b) + tolerance / 2
        m = (c - b) / 2
        abs(m) ≤ tol && return fb > 0 ? (c, b) : (b, c)
        if iszero(fb)
            # a root at b: step toward c by tol, doubling the step across a run of zeros, at most to the midpoint
            plateau |= iszero(fa)
            d = e = plateau ? copysign(min(2 * abs(d), abs(m)), m) : copysign(tol, m)
        elseif abs(e) ≥ tol && abs(fa) > abs(fb)
            s = fb / fa
            if a == c
                p, q = 2 * m * s, 1 - s
            else
                q, r = fa / fc, fb / fc
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)
            end
            p > 0 ? (q = -q) : (p = -p)
            if 2 * p < min(3 * m * q - abs(tol * q), abs(e * q))
                e, d = d, p / q
            else
                d = e = m
            end
        else
            d = e = m
        end
        a, fa = b, fb
        b += abs(d) > tol ? d : copysign(tol, m)
        fb = f(b)
    end
    throw_unbracketed(tolerance, 4 * precision(FT))
end

@noinline throw_unbracketed(tolerance, iterations) =
    throw(ErrorException("MorrisonMilbrandt2015: event root not bracketed to $tolerance in $iterations iterations"))

function bracket_root(::RootSolversRootFinder, _, _, _, _, _, _)
    throw(ErrorException("RootSolversRootFinder requires `using RootSolvers`"))
end
