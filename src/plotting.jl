"""
    figure_theme() -> Makie.Theme

Theme of the package figures: chart surface, hairline grid and axes, and text in ink. Requires
`using CairoMakie`.
"""
function figure_theme end

"""
    series_color(i) -> color

Color of the `i`-th series of the package figures, from a palette of eight in a fixed order.
Requires `using CairoMakie`.
"""
function series_color end

"""
    phase_color(phase) -> color

Color of liquid (`:liquid`), ice (`:ice`), or their sum (`:total`) in the package figures. Requires
`using CairoMakie`.
"""
function phase_color end

"""
    series_style(i) -> (linestyle, linewidth)

Line style and width of the `i`-th scheme or trajectory of the package figures, for `i` up to 4.
Requires `using CairoMakie`.
"""
function series_style end

"""
    plot_evolution(trajectories; labels, panels = (:condensate, :supersaturation, :temperature)) -> Figure

Liquid and ice (`:condensate`), the supersaturations over liquid and over ice (`:supersaturation`),
and the temperature (`:temperature`) of each [`trajectory`](@ref) of one problem against time over
its step, one panel per entry of `panels`, with a marker at each event. Color marks the phase and
line style the trajectory. A phase that holds no mass and never exchanges any is left out. Requires
`using CairoMakie`.
"""
function plot_evolution end

"""
    plot_evolution!(position, trajectories; labels, panels) -> NamedTuple of axes

[`plot_evolution`](@ref) without its legend into the figure position `position`, with the panels side
by side, returning the axes by panel. Requires `using CairoMakie`.
"""
function plot_evolution! end

"""
    plot_step_change(problem, schemes, Δts; labels, reference = nothing) -> Figure

Change of liquid, ice, and their sum over one step of each of `schemes` for `problem`, against the
step `Δt` for each step in `Δts` [s], on a symmetric logarithmic axis, with two guides: the change
at the rate of the start of the step, and the loss of all of a phase. Color marks the phase and line
style the scheme. `reference = (; Δts, changes)` adds the liquid and ice changes `changes` of a
reference solution at its steps `Δts` as dots, and a panel with the difference of each scheme from
it relative to the larger reference change. Requires `using CairoMakie`.
"""
function plot_step_change end

"""
    plot_step_change!(position, problem, schemes, Δts; labels, reference = nothing, difference = true, ylimits = nothing) -> (; change, difference)

[`plot_step_change`](@ref) without its legend into the figure position `position`, returning the
axes. The difference panel is drawn only when `reference` is given and `difference` holds; without
it, the returned `difference` is `nothing`. `ylimits = (lo, hi)` sets the range of the change axis
[kg kg⁻¹]. Requires `using CairoMakie`.
"""
function plot_step_change! end
