module MorrisonMilbrandt2015RootSolversExt

using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015
using RootSolvers: RootSolvers as RS

function _rs_bracketed(f, method, a::FT, b::FT) where {FT}
    sol = RS.find_zero(f, method)
    t = FT(sol.root)
    lo, hi = min(a, b), max(a, b)
    return (sol.converged && lo ≤ t ≤ hi) ? t : FT(Inf)
end

function (::MM2015.MM2015Solver{MM2015.MM2015RSBrentSolverMethod})(
    f,
    bracket::NTuple{2, FT};
    max_iter::Int = 80,
    xatol::FT = sqrt(eps(FT)),
    fatol::FT = sqrt(eps(FT)),
    rtol::FT = sqrt(eps(FT)),
) where {FT}
    a, b = bracket
    fa, fb = f(a), f(b)
    (fa * fb > zero(FT)) && return FT(Inf)
    return _rs_bracketed(f, RS.BrentsMethod{FT}(a, b), a, b)
end

function (::MM2015.MM2015Solver{MM2015.MM2015RSRegulaFalsiSolverMethod})(
    f,
    bracket::NTuple{2, FT};
    max_iter::Int = 80,
    xatol::FT = sqrt(eps(FT)),
    fatol::FT = sqrt(eps(FT)),
    rtol::FT = sqrt(eps(FT)),
) where {FT}
    a, b = bracket
    fa, fb = f(a), f(b)
    (fa * fb > zero(FT)) && return FT(Inf)
    return _rs_bracketed(f, RS.RegulaFalsiMethod{FT}(a, b), a, b)
end

function (::MM2015.MM2015Solver{MM2015.MM2015RSBisectionSolverMethod})(
    f,
    bracket::NTuple{2, FT};
    max_iter::Int = 80,
    xatol::FT = sqrt(eps(FT)),
    fatol::FT = sqrt(eps(FT)),
    rtol::FT = sqrt(eps(FT)),
) where {FT}
    a, b = bracket
    fa, fb = f(a), f(b)
    (fa * fb > zero(FT)) && return FT(Inf)
    return _rs_bracketed(f, RS.BisectionMethod{FT}(a, b), a, b)
end

"""Newton AD from the bracket midpoint. The root is rejected if it leaves `[a, b]`."""
function (::MM2015.MM2015Solver{MM2015.MM2015RSNewtonADSolverMethod})(
    f,
    bracket::NTuple{2, FT};
    max_iter::Int = 80,
    xatol::FT = sqrt(eps(FT)),
    fatol::FT = sqrt(eps(FT)),
    rtol::FT = sqrt(eps(FT)),
) where {FT}
    a, b = bracket
    x0 = (a + b) / 2
    return _rs_bracketed(f, RS.NewtonsMethodAD{FT}(x0), a, b)
end

end # module
