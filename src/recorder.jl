"""Recorder that keeps no segments."""
struct NoRecorder end

@inline record!(::NoRecorder, _, _, _, _, _, _, _, _, _) = nothing

"""
    Segment{FT}

One segment of a step: start time `t` [s], `duration` [s], the `active` phases, the state at its start
(liquid `x_l`, ice `x_i`, temperature `T`, supersaturations `δ` and `δ_i`), and the `event` that ends it.
"""
struct Segment{FT}
    t::FT
    duration::FT
    active::ActivePhases
    x_l::FT
    x_i::FT
    T::FT
    δ::FT
    δ_i::FT
    event::EventKind
end

"""Recorder that keeps every segment."""
struct SegmentRecorder{FT, VS <: AbstractVector{<:Segment{FT}}}
    segments::VS
end

SegmentRecorder{FT}() where {FT} = SegmentRecorder(Segment{FT}[])

record!(r::SegmentRecorder, t, duration, a, x_l, x_i, T, δ, δ_i, kind) =
    push!(r.segments, Segment(t, duration, a, x_l, x_i, T, δ, δ_i, kind))

"""
    Trajectory{FT, S, C}

A step of `scheme` with its segments: the `scheme`, the `context` it evaluates states with (the
[`Coefficients`](@ref) of a frozen-coefficient scheme, the problem of [`MM2015`](@ref)), the step `Δt`
[s], the `segments`, and the mean `rates` of liquid and ice.
"""
struct Trajectory{FT, S, C, VS <: AbstractVector{<:Segment{FT}}}
    scheme::S
    context::C
    Δt::FT
    segments::VS
    rates::NTuple{2, FT}
end

"""
    state_at(trajectory, t) -> (; δ, δ_i, x_l, x_i, T)

Supersaturations over liquid and over ice, liquid, ice, and temperature at time `t` in `[0, Δt]` of a
recorded step.
"""
function state_at(trajectory::Trajectory{FT, <:FrozenCoefficientScheme}, t::Real) where {FT}
    (; scheme, context, segments) = trajectory
    segment = segments[something(findlast(s -> s.t ≤ t, segments), firstindex(segments))]
    δ, δ_i, dx_l, dx_i = advance(scheme, context, segment.active, segment.δ, segment.δ_i, FT(t) - segment.t)
    return (; δ, δ_i, x_l = segment.x_l + dx_l, x_i = segment.x_i + dx_i, T = segment.T)
end
