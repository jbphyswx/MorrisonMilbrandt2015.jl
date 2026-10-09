module MorrisonMilbrandt2015CairoMakieExt

using CairoMakie: CairoMakie as CM
using MorrisonMilbrandt2015: MorrisonMilbrandt2015 as MM2015

# categorical slots in their fixed order, then text and chrome tokens
const SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
const INK, SECONDARY, MUTED = "#0b0b0b", "#52514e", "#898781"
const SURFACE, GRID, AXIS = "#fcfcfb", "#e1e0d9", "#c3c2b7"

const PHASE_COLORS = (; liquid = SERIES[3], ice = SERIES[1], total = SERIES[2])
const STYLES = ((:solid, 2.4), (:dash, 1.8), (:dot, 1.8), (:dashdot, 1.8))

"""Phase, legend label, and factor from kg kg⁻¹ or K to the panel unit of each quantity of `state_at`."""
const QUANTITIES = (;
    x_l = (:liquid, "liquid", 1e3),
    x_i = (:ice, "ice", 1e3),
    δ = (:liquid, "liquid", 1e3),
    δ_i = (:ice, "ice", 1e3),
    T = (:none, "temperature", 1.0),
)

"""Axis label and quantities of each panel of `plot_evolution`."""
const PANELS = (;
    condensate = ("condensate [g kg⁻¹]", (:x_l, :x_i)),
    supersaturation = ("supersaturation [g kg⁻¹]", (:δ, :δ_i)),
    temperature = ("temperature [K]", (:T,)),
)

function MM2015.figure_theme()
    return CM.Theme(;
        fontsize = 13,
        backgroundcolor = SURFACE,
        textcolor = INK,
        Axis = (;
            backgroundcolor = SURFACE, xgridcolor = GRID, ygridcolor = GRID, xgridwidth = 1, ygridwidth = 1,
            xminorgridvisible = false, yminorgridvisible = false, topspinevisible = false, rightspinevisible = false,
            leftspinecolor = AXIS, bottomspinecolor = AXIS, xtickcolor = AXIS, ytickcolor = AXIS,
            xticklabelcolor = SECONDARY, yticklabelcolor = SECONDARY, xlabelcolor = INK, ylabelcolor = SECONDARY,
            titlefont = :regular, titlealign = :left, titlesize = 13, titlecolor = INK,
        ),
        Legend = (; framevisible = false, labelcolor = INK, orientation = :horizontal, tellheight = true, tellwidth = false),
        Lines = (; linewidth = 2, joinstyle = :round, linecap = :round),
        Scatter = (; markersize = 9, strokecolor = SURFACE, strokewidth = 2),
    )
end

function MM2015.series_color(i::Int)
    i ≤ length(SERIES) || throw(ArgumentError("a figure holds at most $(length(SERIES)) series"))
    return SERIES[i]
end

function MM2015.phase_color(phase::Symbol)
    haskey(PHASE_COLORS, phase) || throw(ArgumentError("phase must be :liquid, :ice, or :total, got :$phase"))
    return PHASE_COLORS[phase]
end

function MM2015.series_style(i::Int)
    i ≤ length(STYLES) || throw(ArgumentError("a figure holds at most $(length(STYLES)) schemes or trajectories"))
    return STYLES[i]
end

scheme_label(scheme) = string(nameof(typeof(scheme)))

"""Ticks at `10^e` for each `e` in `exponents`, labeled as powers of ten."""
decade_ticks(exponents) = (exp10.(exponents), [CM.rich("10", CM.superscript(string(e))) for e in exponents])

"""
Tick magnitudes `(m, e)` for `m 10^e` from `threshold` to `top`, with `m` in `mantissas` and every
`stride`-th power of ten.
"""
function tick_magnitudes(threshold, top, mantissas, stride)
    top < threshold && return Tuple{Int, Int}[]
    e_lo, e_hi = ceil(Int, log10(threshold)), floor(Int, log10(top))
    return [(m, e) for e in (e_lo - 1):stride:e_hi for m in mantissas if threshold ≤ m * exp10(e) ≤ top]
end

tick_label(sign, (m, e)) = CM.rich(sign * (m == 1 ? "" : "$(m)×") * "10", CM.superscript(string(e)))

"""
Ticks of a symmetric logarithmic axis with linear threshold `threshold` from `lo` to `hi`: zero and
signed powers of ten, every second one past eight decades, or 1, 2, and 5 times the powers of ten when
neither side spans two decades.
"""
function symlog_ticks(lo, hi, threshold)
    decades = log10(max(hi, -lo, threshold) / threshold)
    mantissas, stride = decades < 2 ? ((1, 2, 5), 1) : ((1,), decades > 8 ? 2 : 1)
    up, down = tick_magnitudes(threshold, hi, mantissas, stride), reverse(tick_magnitudes(threshold, -lo, mantissas, stride))
    values = vcat([-m * exp10(e) for (m, e) in down], 0.0, [m * exp10(e) for (m, e) in up])
    labels = vcat([tick_label("−", t) for t in down], [CM.rich("0")], [tick_label("", t) for t in up])
    return values, labels
end

"""Times of the events of a trajectory: the ends of its segments that end in an event."""
event_times(trajectory) = [s.t + s.duration for s in trajectory.segments if s.event != MM2015.EndOfStep]

"""Times at which a trajectory is drawn: `samples` even times over the step and the start of each segment."""
sample_times(trajectory, samples) =
    sort!(unique!(vcat(collect(range(zero(trajectory.Δt), trajectory.Δt; length = samples)), [s.t for s in trajectory.segments])))

"""Whether `phase` holds mass or exchanges it in any segment of `trajectories`."""
phase_present(trajectories, phase) = phase === :none || any(trajectories) do trajectory
    any(s -> phase === :liquid ? (s.active.liquid || s.x_l > 0) : (s.active.ice || s.x_i > 0), trajectory.segments)
end

"""Quantities of `panel` drawn for `trajectories`."""
panel_quantities(trajectories, panel) = filter(q -> phase_present(trajectories, first(QUANTITIES[q])), collect(last(PANELS[panel])))

quantity_color(quantity) = (phase = first(QUANTITIES[quantity]); phase === :none ? INK : PHASE_COLORS[phase])

function MM2015.plot_evolution!(
    position, trajectories;
    labels = map(tr -> scheme_label(tr.scheme), trajectories),
    panels = (:condensate, :supersaturation, :temperature),
    samples::Int = 400,
)
    length(labels) == length(trajectories) || throw(ArgumentError("one label per trajectory"))
    layout = CM.GridLayout(position)
    axes = map(enumerate(panels)) do (column, panel)
        CM.Axis(layout[1, column]; ylabel = first(PANELS[panel]), xlabel = "t [s]")
    end
    # the first trajectory is drawn last, on top
    for i in reverse(eachindex(trajectories))
        trajectory = trajectories[i]
        linestyle, linewidth = MM2015.series_style(i)
        times = sample_times(trajectory, samples)
        states = [MM2015.state_at(trajectory, t) for t in times]
        events = event_times(trajectory)
        event_states = [MM2015.state_at(trajectory, t) for t in events]
        for (axis, panel) in zip(axes, panels), quantity in panel_quantities(trajectories, panel)
            factor = last(QUANTITIES[quantity])
            color = quantity_color(quantity)
            CM.lines!(axis, times, [factor * getproperty(s, quantity) for s in states]; color, linestyle, linewidth)
            isempty(events) || CM.scatter!(axis, events, [factor * getproperty(s, quantity) for s in event_states];
                color, markersize = 7, strokecolor = SURFACE, strokewidth = 1.5)
        end
    end
    CM.linkxaxes!(axes...)
    return NamedTuple{Tuple(panels)}(Tuple(axes))
end

function MM2015.plot_evolution(
    trajectories;
    labels = map(tr -> scheme_label(tr.scheme), trajectories),
    panels = (:condensate, :supersaturation, :temperature),
    kwargs...,
)
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (290 * length(panels) + 70, 380), figure_padding = (6, 10, 6, 6))
        MM2015.plot_evolution!(figure[2, 1], trajectories; labels, panels, kwargs...)
        quantities = unique([q for panel in panels for q in panel_quantities(trajectories, panel)])
        phases = unique([(QUANTITIES[q][2], quantity_color(q)) for q in quantities if first(QUANTITIES[q]) !== :none])
        groups = [
            [CM.LineElement(; color, linewidth = 2.4) for (_, color) in phases],
            [(s = MM2015.series_style(i); CM.LineElement(; color = INK, linestyle = s[1], linewidth = s[2])) for i in eachindex(labels)],
            [CM.MarkerElement(; marker = :circle, color = INK, strokecolor = SURFACE, strokewidth = 1.5, markersize = 7)],
        ]
        names = [[label for (label, _) in phases], collect(labels), ["event"]]
        CM.Legend(figure[1, 1], groups, names, [nothing, nothing, nothing]; nbanks = 3, groupgap = 28, rowgap = 1,
            patchsize = (26, 10), halign = :left, gridsvalign = :top)
        CM.rowgap!(figure.layout, 1, 4)
        figure
    end
end

"""Liquid and ice changes over each step in `Δts` from the mean rates of `scheme` for `problem`."""
step_changes(scheme, problem, Δts) = [(S = MM2015.tendencies(scheme, problem, Δt); (S[1] * Δt, S[2] * Δt)) for Δt in Δts]

"""Value of component `c` of a liquid and ice change: 1 liquid, 2 ice, 3 their sum."""
component(change, c) = c == 3 ? change[1] + change[2] : change[c]

"""Everything `plot_step_change` draws and the legend it needs, from the changes of each scheme."""
function step_change_layout(problem, schemes, Δts, labels, reference)
    length(labels) == length(schemes) || throw(ArgumentError("one label per scheme"))
    changes = [step_changes(scheme, problem, Δts) for scheme in schemes]
    moves(c) = any(series -> any(x -> !iszero(x[c]), series), changes)
    liquid, ice = moves(1), moves(2)
    liquid || ice || throw(ArgumentError("neither liquid nor ice changes over the steps"))
    components = [(c, phase, label) for (c, phase, label, shown) in
                  ((1, :liquid, "liquid", liquid), (2, :ice, "ice", ice), (3, :total, "liquid + ice", liquid && ice)) if shown]
    k = MM2015.coefficients(problem)
    a = MM2015.active_phases(k, k.δ, k.δ_i, k.x_l, k.x_i)
    initial = (a.liquid ? k.r_l * k.δ : zero(k.δ), a.ice ? k.r_i * k.δ_i : zero(k.δ))
    losing = [c < 3 && any(series -> any(x -> x[c] < 0, series), changes) for (c, _, _) in components]
    exhausted = [c == 1 ? -k.x_l : -k.x_i for ((c, _, _), lose) in zip(components, losing) if lose]
    drawn = vcat([component(x, c) for series in changes for x in series for (c, _, _) in components], exhausted)
    reference === nothing || append!(drawn, [component(x, c) for x in reference.changes for (c, _, _) in components])
    return (; changes, components, initial, losing, drawn)
end

"""Symmetric logarithmic scale, ticks, and limits of the change axis for the values `drawn`, or for `ylimits`."""
function change_scale(drawn, ylimits)
    lo, hi = ylimits === nothing ? (min(minimum(drawn), 0.0), max(maximum(drawn), 0.0)) : ylimits
    two_signs = lo < 0 < hi
    size = max(-lo, hi)
    threshold = two_signs ? 1e-3 * size : max(1e-3 * size, minimum(abs, filter(!iszero, drawn)))
    scale = CM.Makie.Symlog10(threshold; linscale = two_signs ? 1.0 : 0.25)
    limits = ylimits === nothing ? (lo < 0 ? 1.5 * lo : 0.0, hi > 0 ? 1.5 * hi : 0.0) : ylimits
    return scale, symlog_ticks(1.4 * lo, 1.4 * hi, threshold), limits, two_signs
end

function MM2015.plot_step_change!(
    position, problem::MM2015.MM2015Problem, schemes, Δts;
    labels = map(scheme_label, schemes),
    reference = nothing,
    difference::Bool = true,
    ylimits = nothing,
)
    (; changes, components, initial, losing, drawn) = step_change_layout(problem, schemes, Δts, labels, reference)
    scale, yticks, limits, two_signs = change_scale(drawn, ylimits)
    with_difference = reference !== nothing && difference
    xticks = decade_ticks(floor(Int, log10(minimum(Δts))):ceil(Int, log10(maximum(Δts))))
    layout = CM.GridLayout(position)
    change = CM.Axis(layout[1, 1]; xscale = log10, yscale = scale, xticks, yticks, xlabel = "Δt [s]",
        ylabel = "change over the step [kg kg⁻¹]", xlabelvisible = !with_difference, xticklabelsvisible = !with_difference,
        yminorticksvisible = false)
    CM.ylims!(change, limits...)
    two_signs && CM.hlines!(change, [0.0]; color = AXIS, linewidth = 1)
    k = MM2015.coefficients(problem)
    for ((c, phase, _), lose) in zip(components, losing)
        c < 3 || continue
        color = (PHASE_COLORS[phase], 0.55)
        iszero(initial[c]) || CM.lines!(change, Δts, initial[c] .* Δts; color, linewidth = 1, linestyle = :dashdot)
        lose && CM.hlines!(change, [c == 1 ? -k.x_l : -k.x_i]; color, linewidth = 1)
    end
    # the first scheme is drawn last, on top
    for i in reverse(eachindex(schemes)), (c, phase, _) in components
        linestyle, linewidth = MM2015.series_style(i)
        CM.lines!(change, Δts, [component(x, c) for x in changes[i]]; color = PHASE_COLORS[phase], linestyle, linewidth)
    end
    reference === nothing && return (; change, difference = nothing)
    for (c, phase, _) in components
        CM.scatter!(change, reference.Δts, [component(x, c) for x in reference.changes];
            color = PHASE_COLORS[phase], markersize = 7, strokecolor = INK, strokewidth = 0.8)
    end
    with_difference || return (; change, difference = nothing)
    scale_of(x) = max(abs(x[1]), abs(x[2]))
    differences = map(schemes) do scheme
        map(components) do (c, _, _)
            map(step_changes(scheme, problem, reference.Δts), reference.changes) do x, x_ref
                d = abs(component(x, c) - component(x_ref, c)) / scale_of(x_ref)
                iszero(d) ? NaN : d
            end
        end
    end
    finite = filter(isfinite, reduce(vcat, reduce(vcat, differences)))
    e_lo, e_hi = isempty(finite) ? (-16, 0) : (floor(Int, log10(minimum(finite))), ceil(Int, log10(maximum(finite))))
    e_hi = max(e_hi, e_lo + 1)
    bottom = CM.Axis(layout[2, 1]; xscale = log10, yscale = log10, xticks, xlabel = "Δt [s]",
        ylabel = "relative difference\nfrom the reference", yticks = decade_ticks(e_lo:max(1, cld(e_hi - e_lo, 5)):e_hi),
        yminorticksvisible = false)
    CM.ylims!(bottom, exp10(e_lo), exp10(e_hi))
    for i in reverse(eachindex(schemes)), (j, (_, phase, _)) in enumerate(components)
        linestyle, linewidth = MM2015.series_style(i)
        CM.scatterlines!(bottom, reference.Δts, differences[i][j]; color = PHASE_COLORS[phase], linestyle,
            linewidth = 0.7 * linewidth, markersize = 5, strokewidth = 0)
    end
    CM.linkxaxes!(change, bottom)
    CM.rowsize!(layout, 2, CM.Relative(0.3))
    CM.rowgap!(layout, 1, 6)
    return (; change, difference = bottom)
end

"""Legend groups and labels of `plot_step_change`: phases, schemes, and the guides drawn."""
function step_change_legend(problem, schemes, Δts, labels, reference)
    (; components, initial, losing) = step_change_layout(problem, schemes, Δts, labels, reference)
    guides = [
        (CM.MarkerElement(; marker = :circle, color = SURFACE, strokecolor = INK, strokewidth = 0.8, markersize = 7), "reference solution", reference !== nothing),
        (CM.LineElement(; color = MUTED, linestyle = :dashdot, linewidth = 1), "initial rate × Δt", any(!iszero, initial)),
        (CM.LineElement(; color = MUTED, linewidth = 1), "−q: phase exhausted", any(losing)),
    ]
    groups = [
        [CM.LineElement(; color = PHASE_COLORS[phase], linewidth = 2.4) for (_, phase, _) in components],
        [(s = MM2015.series_style(i); CM.LineElement(; color = INK, linestyle = s[1], linewidth = s[2])) for i in eachindex(schemes)],
        [element for (element, _, shown) in guides if shown],
    ]
    names = [[label for (_, _, label) in components], collect(labels), [label for (_, label, shown) in guides if shown]]
    kept = findall(!isempty, groups)
    return groups[kept], names[kept]
end

function MM2015.plot_step_change(
    problem::MM2015.MM2015Problem, schemes, Δts;
    labels = map(scheme_label, schemes),
    reference = nothing,
)
    return CM.with_theme(MM2015.figure_theme()) do
        figure = CM.Figure(; size = (860, reference === nothing ? 460 : 640), figure_padding = (6, 10, 6, 6))
        MM2015.plot_step_change!(figure[2, 1], problem, schemes, Δts; labels, reference)
        groups, names = step_change_legend(problem, schemes, Δts, labels, reference)
        CM.Legend(figure[1, 1], groups, names, fill(nothing, length(groups)); nbanks = 3, groupgap = 28, rowgap = 1,
            patchsize = (26, 10), halign = :left, gridsvalign = :top)
        CM.rowgap!(figure.layout, 1, 4)
        figure
    end
end

end # module
