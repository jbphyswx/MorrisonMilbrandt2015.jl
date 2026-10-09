using Documenter: Documenter
using MorrisonMilbrandt2015: MorrisonMilbrandt2015

Documenter.makedocs(;
    sitename = "MorrisonMilbrandt2015.jl",
    modules = [MorrisonMilbrandt2015],
    authors = "Jordan Benjamin",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = 500_000,
    ),
    pages = [
        "Home" => "index.md",
        "Theory" => [
            "Parcel model" => "theory/parcel_model.md",
            "Appendix C" => "theory/appendix_c.md",
            "Moisture bases" => "theory/moisture_bases.md",
        ],
        "Schemes" => [
            "MM2015FixedT" => "schemes/fixed_T.md",
            "MM2015PiecewiseLinear" => "schemes/piecewise_linear.md",
            "MM2015" => "schemes/t_updating.md",
        ],
        "Numerics" => [
            "Active phases and events" => "numerics/events.md",
            "Integrator" => "numerics/integrator.md",
            "Depletion" => "depletion.md",
            "Performance" => "numerics/performance.md",
        ],
        "Gallery" => "gallery.md",
        "Validation" => "validation.md",
        "API" => "api.md",
        "Internals" => "internals.md",
    ],
    doctest = true,
    checkdocs = :all,
)

Documenter.deploydocs(; repo = "github.com/jbphyswx/MorrisonMilbrandt2015.jl.git", devbranch = "main", push_preview = true)
