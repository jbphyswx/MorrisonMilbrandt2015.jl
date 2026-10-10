"""
    generate_assets.jl

Write the figures in `docs/src/assets`, each checked against the package values it shows.

    julia --project=docs/generate_assets docs/generate_assets/generate_assets.jl
"""

using CairoMakie: CairoMakie as CM
using JSON: JSON
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

repository() = dirname(dirname(@__DIR__))
isdefined(@__MODULE__, :ParcelCorpus) || include(joinpath(repository(), "test", "corpus.jl"))
isdefined(@__MODULE__, :ParcelReference) || include(joinpath(repository(), "test", "reference", "parcel_ode.jl"))

# text, chrome, and diverging tokens of the reference palette
TOKENS = (;
    ink = "#0b0b0b", secondary = "#52514e", muted = "#898781", axis = "#c3c2b7", surface = "#fcfcfb",
    positive = "#2a78d6", negative = "#e34948", neutral = "#f0efec",
    ordinal = ("#86b6ef", "#2a78d6", "#104281"),
)

"""Throw unless `condition` holds."""
check(condition, message) = condition || error("figure check failed: $message")

"""Ticks at `10^n` for `n` in `exponents`, labeled as powers of ten."""
decade_ticks(exponents) = (exp10.(exponents), [CM.rich("10", CM.superscript(string(n))) for n in exponents])

"""Save `figure` as `docs/src/assets/<name>` and return the name."""
function save_asset(name, figure)
    CM.save(joinpath(repository(), "docs", "src", "assets", name), figure; px_per_unit = 2)
    return name
end

backend() = MM2015.DefaultThermodynamicsBackend()

"""
A problem at `T` [K], 800 hPa, and `RH_l = r_v / r_sl`, with 0.1 g kg⁻¹ of liquid and of ice, relaxation
times of 10 s (liquid) and 50 s (ice), and the forcing `w` [m s⁻¹], `dTdt` [K s⁻¹], `dq_vap_dt` [s⁻¹],
in the corpus conventions.
"""
function regime_problem(T, RH_l; q_l = 1e-4, q_i = 1e-4, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0)
    case = (; name = :regime, T, p = 8.0e4, humidity = (:RH_l, RH_l), q_l, q_i, τ_l = 10.0, τ_i = 50.0, w, dTdt, dq_vap_dt, Δt = 1.0)
    return first(ParcelCorpus.corpus_problem(case, backend()))
end

"""Regimes at the start of the step: file name, title, temperature [K], `r_v/r_sl`, ice [kg kg⁻¹], and the step [s] of the evolution figure."""
REGIMES = (
    (; name = "supersaturated", title = "supersaturated over liquid and ice", T = 261.0, RH_l = 1.02, q_i = 1e-4, Δt = 300.0),
    (; name = "between_saturations", title = "between ice and liquid saturation", T = 261.0, RH_l = 0.95, q_i = 1e-4, Δt = 300.0),
    (; name = "subsaturated", title = "subsaturated over liquid and ice", T = 261.0, RH_l = 0.85, q_i = 1e-4, Δt = 300.0),
    (; name = "above_triple_point", title = "above the triple point, between the saturations", T = 276.0, RH_l = 1.01, q_i = 1e-4, Δt = 300.0),
    (; name = "warm_supersaturated", title = "warm, supersaturated over liquid", T = 295.0, RH_l = 1.05, q_i = 0.0, Δt = 60.0),
)

problem_of(regime; kwargs...) = regime_problem(regime.T, regime.RH_l; q_i = regime.q_i, kwargs...)
regime_named(name) = only(filter(r -> r.name == name, REGIMES))

SCHEMES = (MM2015.MM2015(), MM2015.MM2015FixedT(), MM2015.MM2015PiecewiseLinear())
SCHEME_LABELS = ("MM2015", "MM2015FixedT", "MM2015PiecewiseLinear")

"""Steps of the step-change figures [s]: from 1 s to an hour, and the steps with a reference solution."""
STEPS = exp10.(range(0.0, log10(3600.0); length = 200))
REFERENCE_STEPS = exp10.(range(0.0, log10(3600.0); length = 13))

reference_solution(problem, Δt) =
    ParcelReference.solve(ParcelReference.Nonlinear(problem), Δt; rtol = 1e-11, atol_q = 1e-17, atol_T = 1e-9)

"""Check that the default `MM2015` matches the reference rates `S_ref` of `problem` over `Δt` to ten times its tolerances."""
function check_against_reference(problem, Δt, S_ref, label)
    scheme = MM2015.MM2015()
    S = MM2015.tendencies(scheme, problem, Δt)
    bound = 10 * (scheme.rtol * maximum(abs, S_ref) + scheme.atol_q / Δt)
    check(maximum(abs.(Tuple(S) .- S_ref)) ≤ bound, "$label: MM2015 against the reference at Δt = $Δt s")
end

"""Liquid and ice changes of the reference solution of `problem` over each step in `Δts`, each checked against `MM2015`."""
function reference_changes(problem, Δts, label)
    model = ParcelReference.Nonlinear(problem)
    return map(Δts) do Δt
        S_ref = ParcelReference.rates(model, reference_solution(problem, Δt), Δt)
        check_against_reference(problem, Δt, S_ref, label)
        (S_ref[1] * Δt, S_ref[2] * Δt)
    end
end

"""Instantaneous liquid and ice rates at the start of `problem` [kg kg⁻¹ s⁻¹]."""
function initial_rates(problem)
    k = MM2015.coefficients(problem)
    a = MM2015.active_phases(k, k.δ, k.δ_i, k.x_l, k.x_i)
    return (a.liquid ? k.r_l * k.δ : 0.0, a.ice ? k.r_i * k.δ_i : 0.0)
end

# where each phase grows or shrinks at the start of a step
function figure_regime_map()
    Ts = range(225.0, 300.0; length = 301)
    RHs = range(0.75, 1.15; length = 201)
    ice = [sign(initial_rates(regime_problem(T, RH))[2]) for T in Ts, RH in RHs]
    p = 8.0e4
    ice_saturation = [MM2015.saturation(backend(), T, p, MM2015.Ice()).r / MM2015.saturation(backend(), T, p, MM2015.Liquid()).r for T in Ts]
    T_tr = MM2015.T_triple(backend(), Float64)
    for regime in REGIMES
        problem = problem_of(regime)
        S = MM2015.tendencies(MM2015.MM2015FixedT(), problem, 1e-3)
        check(sign.(Tuple(S)) == sign.(initial_rates(problem)), "regime map: rates at the start of $(regime.name)")
    end
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (760, 520))
        axis = CM.Axis(figure[2, 1]; xlabel = "T [K]", ylabel = "r_v / r_sl")
        classes = [(TOKENS.negative, 0.30), (TOKENS.neutral, 1.0), (TOKENS.positive, 0.30)]
        CM.heatmap!(axis, Ts, RHs, ice; colormap = CM.cgrad([CM.Makie.to_color(class) for class in classes]; categorical = true), colorrange = (-1, 1))
        CM.vlines!(axis, [T_tr]; color = TOKENS.axis, linewidth = 1)
        CM.hlines!(axis, [1.0]; color = TOKENS.ink, linewidth = 1.5)
        CM.lines!(axis, Ts, ice_saturation; color = TOKENS.secondary, linewidth = 1.5)
        CM.text!(axis, 226.0, 1.0; text = "liquid saturation", align = (:left, :bottom), color = TOKENS.ink, fontsize = 12)
        label_index = searchsortedfirst(Ts, 264.0)
        CM.text!(axis, Ts[label_index], ice_saturation[label_index]; text = "ice saturation", align = (:left, :top), color = TOKENS.secondary, fontsize = 12, offset = (6, 0))
        CM.text!(axis, T_tr, 0.76; text = " T_tr", align = (:left, :bottom), color = TOKENS.secondary, fontsize = 12)
        labels = (
            (240.0, 1.10, "liquid and ice grow"),
            (243.0, 0.86, "liquid evaporates, ice grows"),
            (255.0, 0.775, "liquid and ice shrink"),
            (282.0, 1.13, "liquid grows"),
            (286.0, 1.04, "liquid grows,\nice sublimates"),
            (286.0, 0.85, "liquid and ice shrink"),
        )
        for (T, RH, text) in labels
            CM.text!(axis, T, RH; text, align = (:center, :center), color = TOKENS.ink, fontsize = 12)
        end
        CM.scatter!(axis, [r.T for r in REGIMES], [r.RH_l for r in REGIMES]; color = TOKENS.ink, markersize = 9, strokecolor = TOKENS.surface, strokewidth = 2)
        CM.limits!(axis, first(Ts), last(Ts), first(RHs), last(RHs))
        CM.Legend(
            figure[1, 1],
            [[CM.PolyElement(; color = (c, a)) for (c, a) in reverse(classes)], [CM.MarkerElement(; marker = :circle, color = TOKENS.ink, markersize = 9)]],
            [["ice grows", "ice unchanged", "ice sublimates"], ["regime states"]],
            [nothing, nothing];
            halign = :left,
        )
        figure
    end
end

# what the psychrometric factor is
function figure_psychrometric()
    Ts = range(220.0, 310.0; length = 181)
    pressures = (5.0e4, 7.0e4, 9.0e4)
    state(T, p) = first(ParcelCorpus.corpus_problem(
        (; name = :state, T, p, humidity = (:RH_l, 1.0), q_l = 0.0, q_i = 0.0, τ_l = 10.0, τ_i = 50.0, w = 0.0, dTdt = 0.0, dq_vap_dt = 0.0, Δt = 1.0),
        backend()))
    Γ = Dict(p => [MM2015.coefficients(state(T, p)) for T in Ts] for p in pressures)
    # Γ_l = 1 + L_v ∂x_sl/∂T / c_p, with ∂x_sl/∂T from a central difference of the saturation curve
    problem = state(261.0, 8.0e4)
    k = MM2015.coefficients(problem)
    q_d = MM2015.dry_fraction(problem.basis, problem.state.x_tot)
    x_sl(T) = q_d * MM2015.saturation(backend(), T, 8.0e4, MM2015.Liquid()).r
    c_p = MM2015.cp_m(backend(), problem.state.x_tot, 0.0, 0.0)
    Γ_difference = 1 + MM2015.latent_heat(backend(), 261.0, MM2015.Liquid()) * (x_sl(261.001) - x_sl(260.999)) / 0.002 / c_p
    check(isapprox(k.Γ_l, Γ_difference; rtol = 1e-6), "psychrometric: Γ_l against a difference of the saturation curve")
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (900, 400))
        for (column, (field, title)) in enumerate(((:Γ_l, "liquid Γ_l"), (:Γ_i, "ice Γ_i")))
            axis = CM.Axis(figure[2, column]; title, xlabel = "T [K]", ylabel = column == 1 ? "Γ" : "", yscale = log10, yticks = [1, 1.5, 2, 3, 5, 10, 20])
            for (color, p) in zip(TOKENS.ordinal, pressures)
                values = [getproperty(c, field) for c in Γ[p]]
                CM.lines!(axis, Ts, values; color)
                CM.text!(axis, last(Ts), last(values); text = " $(round(Int, p / 100)) hPa", align = (:left, :center), color = TOKENS.secondary, fontsize = 12)
            end
            CM.xlims!(axis, first(Ts), last(Ts) + 14)
        end
        CM.Legend(figure[1, 1:2], [CM.LineElement(; color, linewidth = 2) for color in TOKENS.ordinal], ["500 hPa", "700 hPa", "900 hPa"]; halign = :left)
        figure
    end
end

# the change over one step against the step, for one regime
function figure_step_change(regime)
    problem = problem_of(regime)
    reference = (; Δts = REFERENCE_STEPS, changes = reference_changes(problem, REFERENCE_STEPS, regime.name))
    return MM2015.plot_step_change(problem, SCHEMES, STEPS; labels = SCHEME_LABELS, reference)
end

# liquid, ice, supersaturations, and temperature over one step, for one regime
function figure_evolution(regime)
    problem = problem_of(regime)
    Δt = regime.Δt
    trajectories = map(scheme -> MM2015.trajectory(scheme, problem, Δt), SCHEMES)
    for (label, trajectory) in zip(SCHEME_LABELS, trajectories)
        final = MM2015.state_at(trajectory, Δt)
        check(isapprox(final.x_l - problem.state.x_liq, trajectory.rates[1] * Δt; rtol = 1e-12, atol = 1e-20) &&
              isapprox(final.x_i - problem.state.x_ice, trajectory.rates[2] * Δt; rtol = 1e-12, atol = 1e-20),
            "evolution of $(regime.name), $label: end state against the rates")
    end
    reference = reference_solution(problem, Δt)
    events = [(s.event, s.t + s.duration) for s in first(trajectories).segments if s.event != MM2015.EndOfStep]
    check(length(events) == length(reference.events), "evolution of $(regime.name): number of events")
    for ((_, t), (_, t_ref)) in zip(events, reference.events)
        check(isapprox(t, t_ref; rtol = 1e-6, atol = 1e-8), "evolution of $(regime.name): event time")
    end
    return MM2015.plot_evolution(trajectories; labels = SCHEME_LABELS)
end

"""`MM2015FixedT` with `Γ_l = Γ_i = α = 1`: the frozen-coefficient step without latent heating."""
struct WithoutLatentHeating end

function MM2015.tendencies(::WithoutLatentHeating, problem::MM2015.MM2015Problem, Δt)
    k = MM2015.coefficients(problem)
    k_1 = MM2015.Coefficients(k.T, k.p, k.x_l, k.x_i, k.x_sl, k.x_si, k.δ, k.δ_i, 1.0, 1.0, 1.0, k.τ_l, k.τ_i, k.A_l, k.below_triple)
    return MM2015.tendencies(MM2015.MM2015FixedT(), k_1, Δt)
end

"""Appendix C (C5–C7) from the start to the end of the step with the phases active at its start: no saturation or exhaustion events."""
struct WithoutEvents end

function MM2015.tendencies(::WithoutEvents, problem::MM2015.MM2015Problem, Δt)
    k = MM2015.coefficients(problem)
    a = MM2015.active_phases(k, k.δ, k.δ_i, k.x_l, k.x_i)
    (; inv_τ, τ, A_δ, A_δi) = MM2015.relaxation(k, a)
    (; I, I_i) = MM2015.frozen_evolution(k.δ, k.δ_i, inv_τ, τ, A_δ, A_δi, Δt)
    (; Δx_l, Δx_i) = MM2015.condensate_increments(k, a, I, I_i)
    return (; liq = Δx_l / Δt, ice = Δx_i / Δt)
end

"""Legend groups of a step-change panel: `phases` as `(label, phase)`, `labels` by line style, and the guides `(element, label)`."""
function step_legend!(position, phases, labels, guides)
    groups = [
        [CM.LineElement(; color = MM2015.phase_color(phase), linewidth = 2.4) for (_, phase) in phases],
        [(s = MM2015.series_style(i); CM.LineElement(; color = TOKENS.ink, linestyle = s[1], linewidth = s[2])) for i in eachindex(labels)],
        first.(guides),
    ]
    names = [first.(phases), collect(labels), last.(guides)]
    return CM.Legend(position, groups, names, [nothing, nothing, nothing]; nbanks = 3, groupgap = 28, rowgap = 1, patchsize = (26, 10),
        halign = :left, gridsvalign = :top)
end

guide_initial() = (CM.LineElement(; color = TOKENS.muted, linestyle = :dashdot, linewidth = 1), "initial rate × Δt")
guide_exhausted() = (CM.LineElement(; color = TOKENS.muted, linewidth = 1), "−q: phase exhausted")
guide_reference() = (CM.MarkerElement(; marker = :circle, color = TOKENS.surface, strokecolor = TOKENS.ink, strokewidth = 0.8, markersize = 7), "reference solution")

# what latent heating and the events change
function figure_limits()
    regimes = (regime_named("between_saturations"), regime_named("warm_supersaturated"))
    variants = (MM2015.MM2015FixedT(), WithoutLatentHeating(), WithoutEvents())
    labels = ("MM2015FixedT", "without latent heating (Γ = 1)", "without events (Appendix C to the end of the step)")
    for regime in regimes
        problem = problem_of(regime)
        k = MM2015.coefficients(problem)
        S_fixed, S_events = MM2015.tendencies(variants[1], problem, 0.1), MM2015.tendencies(WithoutEvents(), problem, 0.1)
        check(all(isapprox(S_fixed[phase], S_events[phase]; rtol = 1e-12) for phase in (:liq, :ice)),
            "limits, $(regime.name): Appendix C without events equals FixedT before the first event")
        S_unit = MM2015.tendencies(WithoutLatentHeating(), problem, 1e-6)
        S_bare = initial_rates(problem) .* (k.Γ_l, k.Γ_i)
        check(all(isapprox.(Tuple(S_unit), S_bare; rtol = 1e-5)), "limits, $(regime.name): initial rates without latent heating are δ/τ")
    end
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (1100, 500), figure_padding = (6, 10, 6, 6))
        for (column, regime) in enumerate(regimes)
            problem = problem_of(regime)
            k = MM2015.coefficients(problem)
            axes = MM2015.plot_step_change!(figure[2, column], problem, variants, STEPS; labels)
            axes.change.title = regime.title
            column > 1 && (axes.change.ylabelvisible = false)
            levels = [(level, phase) for (level, phase, present) in ((k.δ, :liquid, true), (k.δ_i, :ice, regime.q_i > 0)) if present]
            for (level, phase) in levels
                CM.hlines!(axes.change, [level]; color = (MM2015.phase_color(phase), 0.7), linewidth = 1, linestyle = :dashdotdot)
            end
        end
        step_legend!(figure[1, 1:2], [("liquid", :liquid), ("ice", :ice), ("liquid + ice", :total)], labels,
            [guide_initial(), guide_exhausted(), (CM.LineElement(; color = TOKENS.muted, linestyle = :dashdotdot, linewidth = 1), "supersaturation δ, δ_i at the start")])
        CM.rowgap!(figure.layout, 1, 4)
        figure
    end
end

# the same state with and without each forcing
function figure_forcing()
    regime = regime_named("between_saturations")
    forcings = (
        ("no forcing", (;)),
        ("ascent, w = 1 m s⁻¹", (; w = 1.0)),
        ("cooling, dT/dt = −1 mK s⁻¹", (; dTdt = -1e-3)),
        ("moistening, dq_v/dt = 2 mg kg⁻¹ s⁻¹", (; dq_vap_dt = 2e-6)),
    )
    problems = [problem_of(regime; forcing...) for (_, forcing) in forcings]
    drawn = Float64[-1e-4]
    for problem in problems, scheme in SCHEMES, Δt in STEPS
        S = MM2015.tendencies(scheme, problem, Δt)
        append!(drawn, (S[1] * Δt, S[2] * Δt, (S[1] + S[2]) * Δt))
    end
    ylimits = (1.5 * minimum(drawn), 1.5 * maximum(drawn))
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (1000, 780), figure_padding = (6, 10, 6, 6))
        for (n, ((title, _), problem)) in enumerate(zip(forcings, problems))
            row, column = (n - 1) ÷ 2 + 2, (n - 1) % 2 + 1
            reference = (; Δts = REFERENCE_STEPS, changes = reference_changes(problem, REFERENCE_STEPS, "forcing, $title"))
            axes = MM2015.plot_step_change!(figure[row, column], problem, SCHEMES, STEPS; labels = SCHEME_LABELS, reference,
                difference = false, ylimits)
            axes.change.title = title
            column > 1 && (axes.change.ylabelvisible = false; axes.change.yticklabelsvisible = false)
            row == 2 && (axes.change.xlabelvisible = false; axes.change.xticklabelsvisible = false)
        end
        step_legend!(figure[1, 1:2], [("liquid", :liquid), ("ice", :ice), ("liquid + ice", :total)], SCHEME_LABELS,
            [guide_reference(), guide_initial(), guide_exhausted()])
        CM.rowgap!(figure.layout, 1, 4)
        figure
    end
end

"""Largest relative difference of the mean rates `S` from `S_ref`, relative to the larger reference rate; `NaN` when they are equal."""
rate_error(S, S_ref) = (e = maximum(abs.(Tuple(S) .- S_ref)) / maximum(abs, S_ref); iszero(e) ? NaN : e)

# how accurate each scheme is
function figure_accuracy()
    regime = regime_named("between_saturations")
    problem = problem_of(regime)
    model = ParcelReference.Nonlinear(problem)
    Δts = exp10.(range(0.0, log10(3600.0); length = 25))
    references = [ParcelReference.rates(model, reference_solution(problem, Δt), Δt) for Δt in Δts]
    errors = [[rate_error(MM2015.tendencies(scheme, problem, Δt), S_ref) for (Δt, S_ref) in zip(Δts, references)] for scheme in SCHEMES]
    tolerance_regimes = (regime, regime_named("above_triple_point"), regime_named("warm_supersaturated"))
    rtols = exp10.(range(-4.0, -11.0; length = 15))
    tolerance_errors = map(tolerance_regimes) do r
        p = problem_of(r)
        S_ref = ParcelReference.rates(ParcelReference.Nonlinear(p), reference_solution(p, r.Δt), r.Δt)
        [rate_error(MM2015.tendencies(MM2015.MM2015(; rtol), p, r.Δt), S_ref) for rtol in rtols]
    end
    for (r, e) in zip(tolerance_regimes, tolerance_errors)
        e_check = [e[findfirst(==(rtol), rtols)] for rtol in (1e-6, 1e-8, 1e-10)]
        check(e_check[1] > e_check[2] > e_check[3], "accuracy, $(r.name): error falls with the tolerance")
    end
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (1050, 460))
        left = CM.Axis(
            figure[2, 1]; xscale = log10, yscale = log10, xticks = decade_ticks(0:4), yticks = decade_ticks(-12:2:0),
            xlabel = "Δt [s]", ylabel = "relative error of the mean rates", title = "against the reference, $(regime.title)",
        )
        for i in reverse(eachindex(SCHEMES))
            linestyle, linewidth = MM2015.series_style(i)
            CM.lines!(left, Δts, errors[i]; color = TOKENS.ink, linestyle, linewidth)
        end
        right = CM.Axis(
            figure[2, 2]; xscale = log10, yscale = log10, xticks = decade_ticks(-11:-4), yticks = decade_ticks(-12:2:-4), xreversed = true,
            xlabel = "rtol", title = "MM2015 against its tolerance, step of each evolution figure",
        )
        CM.lines!(right, rtols, rtols; color = TOKENS.axis, linewidth = 1)
        CM.text!(right, rtols[3], rtols[3]; text = " error = rtol", align = (:left, :top), color = TOKENS.secondary, fontsize = 12)
        markers = (:circle, :rect, :utriangle)
        for (marker, e) in zip(markers, tolerance_errors)
            CM.scatterlines!(right, rtols, e; color = TOKENS.ink, marker, markersize = 8, strokewidth = 0, linewidth = 1)
        end
        CM.xlims!(right, 2e-4, 2e-12)
        CM.Legend(figure[1, 1], [(s = MM2015.series_style(i); CM.LineElement(; color = TOKENS.ink, linestyle = s[1], linewidth = s[2])) for i in eachindex(SCHEMES)],
            collect(SCHEME_LABELS); halign = :left, valign = :top)
        CM.Legend(figure[1, 2], [CM.MarkerElement(; marker, color = TOKENS.ink, markersize = 8) for marker in markers],
            [r.title for r in tolerance_regimes]; nbanks = 3, orientation = :horizontal, halign = :left, valign = :top)
        figure
    end
end

"""
Smallest positive `t` with `A t + (δ_0 − A)(1 − e^{−t}) = −1`, the exhaustion of a unit condensate at unit relaxation time,
in BigFloat, or `Inf`; and whether the turning point of the left side lies within `1e-12` of the root level, relative to
the size of its terms.
"""
function reference_depletion(δ_0, A)
    setprecision(BigFloat, 256) do
        δ₀, a = big(δ_0), big(A)
        f(t) = a * t + (δ₀ - a) * (1 - exp(-t)) + 1
        turning = (δ₀ < 0) != (a < 0) ? log((a - δ₀) / a) : big(Inf)
        tangent = isfinite(turning) && abs(f(turning)) ≤ big(1e-12) * max(big(1), abs(a * turning), abs(δ₀ - a))
        lo, hi = if δ₀ < 0 && a < 0
            hi = big(1.0)
            while f(hi) > 0
                hi *= 2
            end
            (big(0.0), hi)
        elseif δ₀ < 0
            f(turning) > 0 && return Inf, tangent
            (big(0.0), turning)
        elseif a < 0
            hi = max(turning, big(1.0))
            while f(hi) > 0
                hi *= 2
            end
            (turning, hi)
        else
            return Inf, tangent
        end
        for _ in 1:300
            m = (lo + hi) / 2
            f(m) > 0 ? (lo = m) : (hi = m)
        end
        return Float64(hi), tangent
    end
end

# how accurate the depletion times are
function figure_depletion()
    magnitudes = exp10.(range(-3.0, 3.0; length = 61))
    panels = ((-1, -1, "δ₀ < 0, A < 0"), (-1, 1, "δ₀ < 0, A > 0"), (1, -1, "δ₀ > 0, A < 0"))
    grids = map(panels) do (sδ, sA, _)
        map(Iterators.product(magnitudes, magnitudes)) do (m_δ, m_A)
            δ_0, A = sδ * m_δ, sA * m_A
            t = MM2015.get_t_out_of_q_no_WBF(δ_0, A, 1.0, 1.0, 1.0, 1.0)
            t_ref, tangent = reference_depletion(δ_0, A)
            if isfinite(t) != isfinite(t_ref)
                check(tangent, "depletion: root existence at δ₀ = $δ_0, A = $A")
                return NaN
            end
            isfinite(t) ? abs(t - t_ref) / t_ref / eps() : NaN
        end
    end
    worst = maximum(g -> maximum(filter(!isnan, g); init = 0.0), grids)
    check(worst ≤ 8, "depletion times within 8 eps (worst $worst)")
    ramp = CM.cgrad(["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (1150, 420))
        local plot
        for (column, ((_, _, title), grid)) in enumerate(zip(panels, grids))
            axis = CM.Axis(figure[1, column]; title, xscale = log10, yscale = log10, xlabel = "|δ₀| τ / (τ_c Γ q_c)", ylabel = column == 1 ? "|A_c| τ² / (τ_c Γ q_c)" : "", backgroundcolor = TOKENS.neutral)
            plot = CM.heatmap!(axis, magnitudes, magnitudes, grid; colormap = ramp, colorrange = (0, 8), nan_color = CM.RGBAf(0, 0, 0, 0))
        end
        CM.Colorbar(figure[1, length(panels) + 1], plot; label = "relative error / eps", ticks = 0:2:8)
        CM.Label(figure[2, 1:length(panels)], "gray: the condensate is not exhausted"; color = TOKENS.secondary, fontsize = 12)
        figure
    end
end

SPEED_MARKERS = (:circle, :rect, :utriangle)

# how fast it is
function figure_speed()
    commit = readchomp(`git -C $(repository()) rev-parse --short=7 HEAD`)
    path = joinpath(repository(), "benchmark", "results", "$commit.json")
    isfile(path) || error("the speed figure needs $path: run benchmark/runbenchmarks.jl first")
    report = JSON.parsefile(path)
    results = report["results"]
    names = [string(c.name) for c in ParcelCorpus.CORPUS]
    scheme_keys = ("MM2015", "MM2015FixedT", "MM2015PiecewiseLinear")
    check(all(n -> haskey(results["MM2015_Float64_default_SpecificHumidity"], n), names), "speed: the benchmark covers the corpus")
    allocations = [case["allocations"] for (key, cases) in results if startswith(key, "MM2015") for (_, case) in cases]
    check(!isempty(allocations) && all(iszero, allocations), "speed: zero allocation in every benchmark")
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (1100, 620))
        for (column, FT) in enumerate(("Float64", "Float32"))
            axis = CM.Axis(
                figure[2, column];
                title = "$FT, default backend, specific humidity",
                xscale = log10, xticks = decade_ticks(1:5), xlabel = "time per call [ns]",
                yticks = (eachindex(names), replace.(names, "_" => " ")), yticklabelsvisible = column == 1, yreversed = true,
            )
            for (marker, key) in zip(SPEED_MARKERS, scheme_keys)
                cases = results["$(key)_$(FT)_default_SpecificHumidity"]
                times = [cases[n]["time_ns"] for n in names]
                floors = [cases[n]["floor_ns"] for n in names]
                CM.scatter!(axis, floors, eachindex(names); marker, color = TOKENS.surface, strokecolor = TOKENS.muted, strokewidth = 1.2, markersize = 8)
                CM.scatter!(axis, times, eachindex(names); marker, color = TOKENS.ink, strokewidth = 0, markersize = 8)
            end
        end
        elements = vcat([CM.MarkerElement(; marker, color = TOKENS.ink, markersize = 8) for marker in SPEED_MARKERS],
            CM.MarkerElement(; marker = :circle, color = TOKENS.surface, strokecolor = TOKENS.muted, strokewidth = 1.2, markersize = 8))
        CM.Legend(figure[1, 1:2], elements, [scheme_keys..., "floor of each scheme (hollow)"]; halign = :left)
        figure
    end
end

function main()
    mkpath(joinpath(repository(), "docs", "src", "assets"))
    written = [
        save_asset("regime_map.png", figure_regime_map()),
        save_asset("psychrometric.png", figure_psychrometric()),
        [save_asset("step_change_$(r.name).png", figure_step_change(r)) for r in REGIMES]...,
        [save_asset("evolution_$(r.name).png", figure_evolution(r)) for r in REGIMES]...,
        save_asset("limits.png", figure_limits()),
        save_asset("forcing.png", figure_forcing()),
        save_asset("accuracy.png", figure_accuracy()),
        save_asset("depletion.png", figure_depletion()),
        save_asset("speed.png", figure_speed()),
    ]
    println("wrote ", join(written, ", "))
    return written
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
