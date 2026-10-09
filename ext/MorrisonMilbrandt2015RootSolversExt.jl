module MorrisonMilbrandt2015RootSolversExt

using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015
using RootSolvers: RootSolvers as RS

narrow_enough(a, b, tolerance) = abs(b - a) ≤ tolerance + 4 * eps(typeof(a)) * max(abs(a), abs(b))

"""Stops at the bracket width of `bracket_root`."""
struct EventTolerance{FT} <: RS.AbstractTolerance{FT}
    width::FT
end
(tol::EventTolerance)(a, b, _) = narrow_enough(a, b, tol.width)

"""From `f(a) = 0`, step toward `b` by the tolerance, doubling each step up to the midpoint, until `f(a) ≠ 0`."""
function gallop_off_zero(f::F, a::FT, b::FT, fa::FT, fb::FT, tolerance::FT) where {F, FT}
    step = zero(FT)
    while iszero(fa) && !narrow_enough(a, b, tolerance)
        step = max(2 * step, tolerance + 4 * eps(FT) * max(abs(a), abs(b)))
        m = a + copysign(min(step, abs(b - a) / 2), b - a)
        fm = f(m)
        fm > 0 ? ((b, fb) = (m, fm)) : ((a, fa) = (m, fm))
    end
    return a, b, fa, fb
end

function MM2015.bracket_root(::MM2015.RootSolversRootFinder{M}, f::F, lo::FT, hi::FT, f_lo::FT, f_hi::FT, tolerance::FT) where {M, F, FT}
    # RootSolvers needs f(lo) < 0
    lo, hi, f_lo, f_hi = gallop_off_zero(f, lo, hi, f_lo, f_hi, tolerance)
    iszero(f_lo) && return lo, hi
    scale = max(abs(f_lo), abs(f_hi))
    g(s) = f(s) / scale
    solution = RS.find_zero(g, M, lo, hi, f_lo / scale, f_hi / scale, RS.TwoPointSolution(), EventTolerance(tolerance), 4 * precision(FT))
    solution.converged ||
        throw(ErrorException("MorrisonMilbrandt2015: $M did not converge on the event root (bracket ($(solution.x0), $(solution.x1)))"))
    a, fa, b, fb = solution.y0 ≤ 0 ? (solution.x0, solution.y0, solution.x1, solution.y1) : (solution.x1, solution.y1, solution.x0, solution.y0)
    if fb ≤ 0
        # an exact zero with the positive end dropped
        a, fa = iszero(fa) ? (a, fa) : (b, fb)
        b, fb = hi, f_hi / scale
    end
    fa ≤ 0 < fb || throw(ErrorException("MorrisonMilbrandt2015: $M returned no sign change of the event function"))
    a, b, fa, fb = gallop_off_zero(g, a, b, fa, fb, tolerance)
    while !narrow_enough(a, b, tolerance)
        m = (a + b) / 2
        g_m = g(m)
        g_m > 0 ? ((b, fb) = (m, g_m)) : ((a, fa) = (m, g_m))
    end
    return a, b
end

end
