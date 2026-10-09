using Aqua: Aqua
using ExplicitImports: ExplicitImports
using MorrisonMilbrandt2015: MorrisonMilbrandt2015
using RootSolvers: RootSolvers
using Test: Test
using Thermodynamics: Thermodynamics

Test.@testset "Aqua.jl" begin
    Aqua.test_all(
        MorrisonMilbrandt2015;
        ambiguities = true,
        persistent_tasks = false,
        stale_deps = true,
    )
end

Test.@testset "Method ambiguities with the extensions loaded" begin
    extensions = map(name -> Base.get_extension(MorrisonMilbrandt2015, name),
        (:MorrisonMilbrandt2015RootSolversExt, :MorrisonMilbrandt2015ThermodynamicsExt))
    Test.@test all(!isnothing, extensions)
    Test.@test isempty(Test.detect_ambiguities(MorrisonMilbrandt2015, extensions...; recursive = true))
end

Test.@testset "ExplicitImports" begin
    ExplicitImports.check_no_implicit_imports(MorrisonMilbrandt2015)
end
