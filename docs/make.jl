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
        "Mathematical contract" => "mathematical_contract.md",
        "Equations" => "equations.md",
        "Depletion" => "depletion.md",
        "MM2015FixedT" => "mm2015ep.md",
        "MM2015" => "mm2015.md",
        "API" => "api.md",
    ],
    checkdocs = :all,
    warnonly = [:missing_docs],
)
