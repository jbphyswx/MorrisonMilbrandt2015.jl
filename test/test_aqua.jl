using Aqua: Aqua
using ExplicitImports: ExplicitImports
using MorrisonMilbrandt2015: MorrisonMilbrandt2015
using Test: Test

Test.@testset "Aqua.jl" begin
    Aqua.test_all(
        MorrisonMilbrandt2015;
        ambiguities = false,
        persistent_tasks = false,
        stale_deps = true,
    )
end

Test.@testset "ExplicitImports" begin
    ExplicitImports.check_no_implicit_imports(MorrisonMilbrandt2015)
end
