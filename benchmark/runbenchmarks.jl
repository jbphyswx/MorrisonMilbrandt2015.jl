# Time, allocations, and event counts of every scheme on the parcel corpus, against a floor of the
# transcendental work each call must do. Writes benchmark/results/<commit>.json.

using BenchmarkTools: BenchmarkTools
using ClimaParams: ClimaParams
using Dates: Dates
using RootSolvers: RootSolvers
using Thermodynamics: Thermodynamics as TD
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(@__DIR__, "..", "test", "corpus.jl"))

const FLOAT_TYPES = (Float64, Float32)
const BACKENDS = (("default", FT -> MM2015.DefaultThermodynamicsBackend()), ("Thermodynamics.jl", FT -> TD.Parameters.ThermodynamicsParameters(FT)))
const BASES = (MM2015.SpecificHumidity(), MM2015.DryAirMixingRatio())
const SCHEMES = (
    ("MM2015PiecewiseLinear", FT -> MM2015.MM2015PiecewiseLinear()),
    ("MM2015FixedT", FT -> MM2015.MM2015FixedT()),
    ("MM2015", FT -> MM2015.MM2015{FT}()),
)

"""Minimum time [ns] of `f(args...)` over BenchmarkTools samples, with its memory and allocation count."""
function measure(f, args...)
    f(args...)
    benchmark = BenchmarkTools.@benchmarkable $f($(args)...) seconds = 0.05
    BenchmarkTools.tune!(benchmark)
    estimate = minimum(BenchmarkTools.run(benchmark))
    return (; time_ns = estimate.time, memory_bytes = estimate.memory, allocations = estimate.allocs)
end

"""Sum of `f` over `xs`, so that no call is elided."""
function sum_over(f::F, xs) where {F}
    s = zero(eltype(xs))
    for x in xs
        s += f(x)
    end
    return s
end

"""Throughput cost [ns per call] of `f` over the inputs `xs`."""
primitive_cost(f, xs) = measure(sum_over, f, xs).time_ns / length(xs)

"""Evenly spread inputs on `[lo, hi]` in type `FT`."""
spread(::Type{FT}, lo, hi; n = 4096) where {FT} = FT.(range(lo, hi; length = n))

function primitive_costs(::Type{FT}) where {FT}
    return (;
        exp = primitive_cost(exp, spread(FT, -5, 5)),
        log = primitive_cost(log, spread(FT, 0.8, 1.2)),
        expm1 = primitive_cost(expm1, spread(FT, -2, 0)),
        log1p = primitive_cost(log1p, spread(FT, 0, 2)),
        sqrt = primitive_cost(sqrt, spread(FT, 0, 2)),
        inv = primitive_cost(inv, spread(FT, 0.5, 2)),
        fast_lambertw0 = primitive_cost(MM2015.fast_lambertw0, spread(FT, -0.3, 10)),
        fast_lambertwm1 = primitive_cost(MM2015.fast_lambertwm1, spread(FT, -0.36, -1e-3)),
    )
end

"""Cost [ns] of one `saturation` evaluation of `thermo`, over temperatures of the corpus."""
function saturation_cost(thermo, ::Type{FT}) where {FT}
    Ts = spread(FT, 225, 295)
    p = FT(8e4)
    return primitive_cost(T -> MM2015.saturation(thermo, T, p, MM2015.Liquid()).dr_dT, Ts)
end

"""Number of segments or steps, events, and the transcendental floor [ns] of one call."""
function work(scheme::Union{MM2015.MM2015PiecewiseLinear, MM2015.MM2015FixedT}, problem, Δt, costs, c_sat)
    segments = MM2015.trajectory(scheme, problem, Δt).segments
    events = count(s -> s.event != MM2015.EndOfStep, segments)
    floor = 2 * c_sat  # coefficients: saturation over liquid and ice
    if scheme isa MM2015.MM2015FixedT
        for s in segments
            (s.active.liquid || s.active.ice) && (floor += costs.expm1)
            s.event in (MM2015.LiquidSaturation, MM2015.IceSaturation) && (floor += costs.log1p)
            s.event in (MM2015.LiquidExhausted, MM2015.IceExhausted) && (floor += costs.expm1)
        end
    end
    return (; segments = length(segments), events, floor_ns = floor)
end

function work(scheme::MM2015.MM2015, problem, Δt, costs, c_sat)
    steps = MM2015.trajectory(scheme, problem, Δt).segments
    events = count(s -> s.event != MM2015.EndOfStep, steps)
    # each step: the state at its start and at its two inner stages, each with saturation over liquid and ice
    return (; segments = length(steps), events, floor_ns = length(steps) * 6 * c_sat)
end

"""Costs [ns] of the pieces of one `MM2015` step at the start of `problem`."""
function step_components(problem::MM2015.MM2015Problem{FT}, Δt) where {FT}
    u = (problem.state.x_liq, problem.state.x_ice, problem.state.T)
    t = zero(FT)
    s = MM2015.parcel_state(problem, t, u)
    below = MM2015.below_triple(problem, s)
    a = MM2015.active_phases(problem, s, below)
    (; F, J, v) = MM2015.parcel_linearization(problem, a, s)
    scaling = MM2015.balancing(J)
    rhs = MM2015.ParcelRHS(problem, a)
    h = Δt / 8
    return (;
        parcel_state = measure(MM2015.parcel_state, problem, t, u).time_ns,
        parcel_linearization = measure(MM2015.parcel_linearization, problem, a, s).time_ns,
        balancing = measure(MM2015.balancing, J).time_ns,
        rhs = measure(rhs, t, u).time_ns,
        phi_functions = measure(MM2015.phi_functions, J .* (h / 2), scaling).time_ns,
        exprb43_step = measure(MM2015.exprb43_step, rhs, t, u, h, Tuple(F), J, v, scaling).time_ns,
        event_indicators = measure(MM2015.event_indicators, problem, a, below, s).time_ns,
        event_indicator_rates = measure(MM2015.event_indicator_rates, problem, a, below, s).time_ns,
    )
end

function collect_results(io)
    results = Dict{String, Any}()
    for FT in FLOAT_TYPES
        costs = primitive_costs(FT)
        results["primitives_$FT"] = costs
        println(io, "primitives $FT [ns/call]: ", costs)
        for (bname, make_thermo) in BACKENDS
            thermo = make_thermo(FT)
            c_sat = saturation_cost(thermo, FT)
            results["saturation_$(FT)_$bname"] = c_sat
            println(io, "saturation $FT $bname: $(round(c_sat; digits = 2)) ns")
            problem, Δt = ParcelCorpus.corpus_problem(ParcelCorpus.corpus_case(:wbf), thermo; FT)
            components = step_components(problem, Δt)
            results["mm2015_step_components_$(FT)_$bname"] = components
            println(io, "MM2015 step components $FT $bname (wbf) [ns]: ", components)
            for basis in BASES, (sname, make_scheme) in SCHEMES
                scheme = make_scheme(FT)
                key = "$(sname)_$(FT)_$(bname)_$(nameof(typeof(basis)))"
                cases = Dict{String, Any}()
                for case in ParcelCorpus.CORPUS
                    problem, Δt = ParcelCorpus.corpus_problem(case, thermo; basis, FT)
                    timing = measure(MM2015.tendencies, scheme, problem, Δt)
                    counts = work(scheme, problem, Δt, costs, c_sat)
                    cases[string(case.name)] = merge(timing, counts, (; ratio_to_floor = timing.time_ns / counts.floor_ns))
                    println(io, rpad(key, 72), " ", rpad(case.name, 30),
                        "time ", lpad(round(timing.time_ns; digits = 1), 10), " ns  alloc ", timing.memory_bytes,
                        "  segments ", counts.segments, "  events ", counts.events,
                        "  floor ", round(counts.floor_ns; digits = 1), "  ratio ", round(timing.time_ns / counts.floor_ns; digits = 2))
                end
                results[key] = cases
            end
        end
    end
    finder = MM2015.MM2015(; root_finder = MM2015.RootSolversRootFinder(RootSolvers.BrentsMethod))
    thermo = MM2015.DefaultThermodynamicsBackend()
    root_finders = Dict{String, Any}()
    for case in ParcelCorpus.CORPUS
        problem, Δt = ParcelCorpus.corpus_problem(case, thermo)
        root_finders[string(case.name)] = (;
            brent_ns = measure(MM2015.tendencies, MM2015.MM2015(), problem, Δt).time_ns,
            rootsolvers_ns = measure(MM2015.tendencies, finder, problem, Δt).time_ns,
        )
    end
    results["mm2015_root_finders_Float64_default_SpecificHumidity"] = root_finders
    println(io, "MM2015 root finders (Float64, default, SpecificHumidity) [ns]: ", root_finders)
    return results
end

json_escape(s::AbstractString) = replace(s, '\\' => "\\\\", '"' => "\\\"", '\n' => "\\n", '\r' => "\\r", '\t' => "\\t")

function write_json(io, value, indent::Int = 0)
    pad, nextpad = " "^indent, " "^(indent + 2)
    if value === nothing
        print(io, "null")
    elseif value isa Bool
        print(io, value ? "true" : "false")
    elseif value isa Number
        print(io, isfinite(value) ? repr(Float64(value)) : "\"$(value)\"")
    elseif value isa AbstractString || value isa Symbol
        print(io, '"', json_escape(string(value)), '"')
    elseif value isa AbstractVector || value isa Tuple
        print(io, '[')
        for (i, item) in enumerate(value)
            i > 1 && print(io, ", ")
            write_json(io, item, indent)
        end
        print(io, ']')
    elseif value isa NamedTuple || value isa AbstractDict
        entries = sort!(collect(pairs(value)); by = entry -> string(first(entry)))
        print(io, "{\n")
        for (i, (key, item)) in enumerate(entries)
            print(io, nextpad, '"', json_escape(string(key)), "\": ")
            write_json(io, item, indent + 2)
            i < length(entries) && print(io, ',')
            print(io, '\n')
        end
        print(io, pad, '}')
    else
        write_json(io, string(value), indent)
    end
end

function main(io = stdout)
    repo = normpath(joinpath(@__DIR__, ".."))
    commit = readchomp(`git -C $repo rev-parse --short=7 HEAD`)
    results = collect_results(io)
    report = (;
        commit,
        worktree_clean = isempty(readchomp(`git -C $repo status --short`)),
        generated_at_utc = string(Dates.now(Dates.UTC)),
        julia_version = string(VERSION),
        cpu = Sys.cpu_info()[1].model,
        threads = Threads.nthreads(),
        versions = (;
            MorrisonMilbrandt2015 = string(pkgversion(MM2015)),
            Thermodynamics = string(pkgversion(TD)),
            ClimaParams = string(pkgversion(ClimaParams)),
            RootSolvers = string(pkgversion(RootSolvers)),
            BenchmarkTools = string(pkgversion(BenchmarkTools)),
        ),
        floor = "saturation evaluations the scheme must make (2 per coefficient set; 6 per MM2015 step) at the measured cost of `saturation`, plus for MM2015FixedT one expm1 per segment with an active phase, one log1p per saturation event, and one expm1 per exhaustion event",
        results,
    )
    path = joinpath(@__DIR__, "results", "$commit.json")
    mkpath(dirname(path))
    open(path, "w") do file
        write_json(file, report)
        println(file)
    end
    println(io, path)
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
