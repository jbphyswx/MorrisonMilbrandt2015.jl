#=

Scalar root finders for residual T-updating ([`MM2015`](@ref)).

Phase A ([`leftmost_root`](@ref)) finds the soonest root by splitting at turning
points of `F`, using the caller’s analytic `f'`. Phase B is a bracketed iterator
inside one pair. Unbracketed methods (Secant, Newton) are never the event search.

`using RootSolvers` adds Phase B methods via MorrisonMilbrandt2015RootSolversExt.

=#

abstract type AbstractMM2015SolverMethod end

abstract type AbstractMM2015BracketedSolverMethod <: AbstractMM2015SolverMethod end
abstract type AbstractMM2015UnbracketedSolverMethod <: AbstractMM2015SolverMethod end

struct MM2015BisectionSolverMethod <: AbstractMM2015BracketedSolverMethod end
struct MM2015FalsePositionSolverMethod <: AbstractMM2015BracketedSolverMethod end
struct MM2015BrentSolverMethod <: AbstractMM2015BracketedSolverMethod end
struct MM2015NewtonSolverMethod <: AbstractMM2015UnbracketedSolverMethod end
struct MM2015SecantSolverMethod <: AbstractMM2015UnbracketedSolverMethod end

"""Phase B via RootSolvers `BrentsMethod`. Requires `using RootSolvers`."""
struct MM2015RSBrentSolverMethod <: AbstractMM2015BracketedSolverMethod end
"""Phase B via RootSolvers `RegulaFalsiMethod`. Requires `using RootSolvers`."""
struct MM2015RSRegulaFalsiSolverMethod <: AbstractMM2015BracketedSolverMethod end
"""Phase B via RootSolvers `BisectionMethod`. Requires `using RootSolvers`."""
struct MM2015RSBisectionSolverMethod <: AbstractMM2015BracketedSolverMethod end
"""Phase B via RootSolvers `NewtonsMethodAD` from the bracket midpoint. Requires `using RootSolvers`."""
struct MM2015RSNewtonADSolverMethod <: AbstractMM2015UnbracketedSolverMethod end

abstract type AbstractMM2015Solver{M <: AbstractMM2015SolverMethod} end

"""Callable scalar Phase B solver. Not a moisture scheme."""
struct MM2015Solver{M <: AbstractMM2015SolverMethod} <: AbstractMM2015Solver{M} end

const default_mm2015_solver = MM2015Solver{MM2015BisectionSolverMethod}()

function _rs_solver_needs_rootsolvers(name::String)
    error("$name requires `using RootSolvers`")
end

function (::MM2015Solver{MM2015RSBrentSolverMethod})(f, bracket; kwargs...)
    _rs_solver_needs_rootsolvers("MM2015RSBrentSolverMethod")
end
function (::MM2015Solver{MM2015RSRegulaFalsiSolverMethod})(f, bracket; kwargs...)
    _rs_solver_needs_rootsolvers("MM2015RSRegulaFalsiSolverMethod")
end
function (::MM2015Solver{MM2015RSBisectionSolverMethod})(f, bracket; kwargs...)
    _rs_solver_needs_rootsolvers("MM2015RSBisectionSolverMethod")
end
function (::MM2015Solver{MM2015RSNewtonADSolverMethod})(f, bracket; kwargs...)
    _rs_solver_needs_rootsolvers("MM2015RSNewtonADSolverMethod")
end

function (::MM2015Solver{MM2015NewtonSolverMethod})(f, x0)
    error("MM2015NewtonSolverMethod is not the residual event search (Phase B is bracketed)")
end

function (::MM2015Solver{MM2015BisectionSolverMethod})(
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
    c = a
    for _ in 1:max_iter
        c = (a + b) / 2
        fc = f(c)
        if abs(fc) ≤ fatol
            return c
        end
        if fa * fc ≤ zero(FT)
            b, fb = c, fc
        else
            a, fa = c, fc
        end
    end
    return FT(Inf)
end

function (::MM2015Solver{MM2015FalsePositionSolverMethod})(
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
    iszero(fb - fa) && return FT(Inf)
    c = a
    for _ in 1:max_iter
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        if abs(fc) ≤ fatol
            return c
        end
        if fa * fc ≤ zero(FT)
            b, fb = c, fc
        else
            a, fa = c, fc
        end
        iszero(fb - fa) && return c
    end
    return FT(Inf)
end

"""Brent (bisection + inverse quadratic). The interpolant is discarded if it leaves `[a, b]`."""
function (::MM2015Solver{MM2015BrentSolverMethod})(
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
    if abs(fa) < abs(fb)
        a, b, fa, fb = b, a, fb, fa
    end
    c, fc = a, fa
    mflag = true
    d = b - a
    for _ in 1:max_iter
        if abs(fb) ≤ fatol
            return b
        end
        if fa != fc && fb != fc
            s =
                (a * fb * fc) / ((fa - fb) * (fa - fc)) +
                (b * fa * fc) / ((fb - fa) * (fb - fc)) +
                (c * fa * fb) / ((fc - fa) * (fc - fb))
        else
            s = b - fb * (b - a) / (fb - fa)
        end
        δ = xatol + rtol * abs(b)
        lo_s = (FT(3) * a + b) / FT(4)
        hi_s = b
        if lo_s > hi_s
            lo_s, hi_s = hi_s, lo_s
        end
        cond_out = (s < lo_s) || (s > hi_s)
        cond_slow = mflag ? abs(s - b) ≥ abs(b - c) / 2 : abs(s - b) ≥ abs(c - d) / 2
        cond_tiny = mflag ? abs(b - c) < δ : abs(c - d) < δ
        if cond_out || cond_slow || cond_tiny
            s = (a + b) / 2
            mflag = true
        else
            mflag = false
        end
        fs = f(s)
        d = c
        c, fc = b, fb
        if fa * fs ≤ zero(FT)
            b, fb = s, fs
        else
            a, fa = s, fs
        end
        if abs(fa) < abs(fb)
            a, b, fa, fb = b, a, fb, fa
        end
    end
    return FT(Inf)
end

function (::MM2015Solver{MM2015SecantSolverMethod})(
    f,
    guess::NTuple{2, FT};
    max_iter::Int = 80,
    atol::FT = sqrt(eps(FT)),
    rtol::FT = sqrt(eps(FT)),
) where {FT}
    x0, x1 = guess
    f0, f1 = f(x0), f(x1)
    for _ in 1:max_iter
        iszero(f1 - f0) && return x1
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        f2 = f(x2)
        if abs(f2) ≤ atol || abs(x2 - x1) ≤ atol + rtol * abs(x2)
            return x2
        end
        x0, f0 = x1, f1
        x1, f1 = x2, f2
    end
    return x1
end

function (::MM2015Solver{MM2015SecantSolverMethod})(f, x0::FT; kwargs...) where {FT}
    shift = max(abs(x0) * eps(FT)^(FT(1) / 3), eps(FT))
    return MM2015Solver{MM2015SecantSolverMethod}()(f, (x0, x0 + shift); kwargs...)
end

@inline function _certified_root(
    f,
    a::FT,
    b::FT,
    fa::FT,
    fb::FT,
    solver,
    time_atol::FT,
    residual_atol::FT,
) where {FT}
    abs(fa) ≤ residual_atol && return a
    abs(fb) ≤ residual_atol && return b
    fa * fb > zero(FT) && return FT(Inf)
    t = solver(
        f,
        (a, b);
        xatol = time_atol,
        fatol = residual_atol,
        rtol = FT(8) * eps(FT),
    )
    if isfinite(t) && a ≤ t ≤ b && abs(f(t)) ≤ residual_atol
        return t
    end
    return FT(Inf)
end

function _scan_root_interval(
    f,
    df,
    lo::FT,
    hi::FT,
    solver,
    time_atol::FT,
    residual_atol::FT,
    derivative_atol::FT,
    subdivisions::Int,
) where {FT}
    a = lo
    fa = f(a)
    dfa = df(a)
    fhi = f(hi)
    direct = _certified_root(
        f, lo, hi, fa, fhi, solver, time_atol, residual_atol,
    )
    isfinite(direct) && return direct
    width = (hi - lo) / subdivisions
    for i in 1:subdivisions
        b = i == subdivisions ? hi : lo + i * width
        fb = i == subdivisions ? fhi : f(b)
        dfb = df(b)
        if dfa * dfb ≤ zero(FT)
            turn = solver(
                df,
                (a, b);
                xatol = time_atol,
                fatol = derivative_atol,
                rtol = FT(8) * eps(FT),
            )
            if isfinite(turn) && a < turn < b && abs(df(turn)) ≤ derivative_atol
                fturn = f(turn)
                root = _certified_root(
                    f, a, turn, fa, fturn, solver, time_atol, residual_atol,
                )
                isfinite(root) && return root
                root = _certified_root(
                    f, turn, b, fturn, fb, solver, time_atol, residual_atol,
                )
                isfinite(root) && return root
            end
        end
        root = _certified_root(
            f, a, b, fa, fb, solver, time_atol, residual_atol,
        )
        isfinite(root) && return root
        a, fa, dfa = b, fb, dfb
    end
    return FT(Inf)
end

"""
    leftmost_root(f, df, lo, hi, solver; t_c5=Inf, ...)

Leftmost certified root detected by a deterministic, allocation-free subdivision
of `(lo, hi]`. Each subdivision is split at a certified derivative root when one
is bracketed. `t_c5` is an additional exact partition boundary for fast C5
transients. The finite subdivision count is an explicit resolution contract; it
is not a proof that arbitrarily close even-multiplicity roots cannot be missed.
"""
function leftmost_root(
    f,
    df,
    lo::FT,
    hi::FT,
    solver::AbstractMM2015Solver = default_mm2015_solver;
    t_c5::FT = FT(Inf),
    time_atol::FT = FT(8) * eps(FT) * max(abs(hi), one(FT)),
    residual_atol::FT = FT(64) * eps(FT),
    derivative_atol::FT = FT(64) * eps(FT),
    subdivisions::Int = 8,
) where {FT}
    lo ≥ hi && return FT(Inf)
    flo = f(lo)
    if abs(flo) ≤ residual_atol
        lo = nextfloat(lo)
        lo ≥ hi && return FT(Inf)
    end

    if isfinite(t_c5) && lo < t_c5 < hi
        root = _scan_root_interval(
            f, df, lo, t_c5, solver, time_atol, residual_atol,
            derivative_atol, subdivisions,
        )
        isfinite(root) && return root
        lo = t_c5
    end
    return _scan_root_interval(
        f, df, lo, hi, solver, time_atol, residual_atol,
        derivative_atol, subdivisions,
    )
end
