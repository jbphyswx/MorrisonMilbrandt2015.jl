using BenchmarkTools: BenchmarkTools
using Dates: Dates
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

include(joinpath(@__DIR__, "..", "test", "fixture.jl"))

const BASELINE_CASES = (
    idle = mm2015_build_state(; freezing = :AF, regime = :sub, q_liq = 0.0, q_ice = 0.0),
    liquid_growth = mm2015_build_state(;
        freezing = :AF, regime = :super, q_liq = 1e-4, q_ice = 0.0, τ_liq = 8.0,
    ),
    ice_growth = mm2015_build_state(;
        freezing = :BF, regime = :super, q_liq = 0.0, q_ice = 1e-4, τ_ice = 8.0,
    ),
    wbf = mm2015_build_state(;
        freezing = :BF, regime = :wbf, q_liq = 1e-4, q_ice = 1e-4,
    ),
    liquid_depletion = mm2015_build_state(;
        freezing = :AF, regime = :sub, q_liq = 1e-4, q_ice = 0.0, Δt = 30.0,
    ),
    ice_depletion = mm2015_build_state(;
        freezing = :BF, regime = :sub, q_liq = 0.0, q_ice = 1e-4, Δt = 30.0,
    ),
    activation = mm2015_build_state(;
        freezing = :AF, regime = :sub, q_liq = 0.0, q_ice = 0.0,
        dqvdt = 2e-6, Δt = 30.0,
    ),
    freeze_crossing = mm2015_build_state(;
        freezing = :AF, regime = :sub, q_liq = 0.0, q_ice = 0.0,
        T = MM2015_T_FREEZE + 0.3, w = 1.0, Δt = 80.0,
    ),
    multi_event = mm2015_build_state(;
        freezing = :BF, regime = :wbf, q_liq = 2e-5, q_ice = 3e-5,
        τ_liq = 2.0, τ_ice = 5.0, dqvdt = -2e-7, Δt = 60.0,
    ),
)

const BASELINE_SCHEMES = (
    piecewise_linear = MM2015.MM2015PiecewiseLinear(),
    ep = MM2015.MM2015FixedT(),
    t_updating = MM2015.MM2015(),
)

json_escape(s::AbstractString) =
    replace(s, '\\' => "\\\\", '"' => "\\\"", '\n' => "\\n", '\r' => "\\r", '\t' => "\\t")

function write_json(io, value, indent::Int = 0)
    pad = " "^indent
    nextpad = " "^(indent + 2)
    if value === nothing
        print(io, "null")
    elseif value isa Bool
        print(io, value ? "true" : "false")
    elseif value isa Number
        print(io, isfinite(value) ? repr(value) : "\"$(value)\"")
    elseif value isa AbstractString
        print(io, '"', json_escape(value), '"')
    elseif value isa AbstractVector || value isa Tuple
        print(io, '[')
        for (i, item) in enumerate(value)
            i > 1 && print(io, ", ")
            write_json(io, item, indent)
        end
        print(io, ']')
    elseif value isa NamedTuple || value isa AbstractDict
        pairs_value = collect(pairs(value))
        print(io, "{\n")
        for (i, (key, item)) in enumerate(pairs_value)
            print(io, nextpad, '"', json_escape(string(key)), "\": ")
            write_json(io, item, indent + 2)
            i < length(pairs_value) && print(io, ',')
            print(io, '\n')
        end
        print(io, pad, '}')
    else
        write_json(io, string(value), indent)
    end
end

function benchmark_case(scheme, state)
    output = try
        mm2015_call_sources(scheme, state)
    catch err
        return (;
            status = "error",
            error = sprint(showerror, err),
            output = nothing,
            event_count = nothing,
            minimum_time_ns = nothing,
            memory_bytes = nothing,
            allocations = nothing,
        )
    end
    trial = run(
        BenchmarkTools.@benchmarkable(mm2015_call_sources($scheme, $state));
        samples = 200,
        evals = 1,
        seconds = 0.25,
    )
    estimate = minimum(trial)
    return (;
        status = "ok",
        error = nothing,
        output = collect(output),
        # The current public API does not expose event counts. Preserve that fact
        # rather than modifying production code before this baseline is frozen.
        event_count = nothing,
        minimum_time_ns = estimate.time,
        memory_bytes = estimate.memory,
        allocations = estimate.allocs,
    )
end

function collect_baseline()
    results = Dict{String, Any}()
    for (scheme_name, scheme) in pairs(BASELINE_SCHEMES)
        cases = Dict{String, Any}()
        for (case_name, state) in pairs(BASELINE_CASES)
            cases[string(case_name)] = benchmark_case(scheme, state)
        end
        results[string(scheme_name)] = cases
    end
    repo = normpath(joinpath(@__DIR__, ".."))
    commit = readchomp(`git -C $repo log -1 --format=%H`)
    worktree = readchomp(`git -C $repo status --short`)
    return (;
        schema_version = 1,
        generated_at_utc = string(Dates.now(Dates.UTC)),
        julia_version = string(VERSION),
        commit,
        worktree_status = worktree,
        event_count_note = "null because the pre-repair public API did not expose event counts",
        results,
    )
end

output_path = joinpath(@__DIR__, "baselines", "pre_repair.json")
mkpath(dirname(output_path))
open(output_path, "w") do io
    write_json(io, collect_baseline())
    println(io)
end
println(output_path)
