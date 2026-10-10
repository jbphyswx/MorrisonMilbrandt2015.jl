using Test: Test
using MorrisonMilbrandt2015: MorrisonMilbrandt2015

include("test_helpers.jl")
include("corpus.jl")
include("reference/parcel_ode.jl")

# Keep each file's helpers and constants in its own namespace.
module AquaTests
include("test_aqua.jl")
end

module ThermodynamicsTests
using ..TestHelpers: TestHelpers
include("test_thermodynamics.jl")
end

module CoefficientsTests
using ..TestHelpers: TestHelpers
using ..ParcelCorpus: ParcelCorpus
using ..ParcelReference: ParcelReference
include("test_coefficients.jl")
end

module ReferenceTests
using ..TestHelpers: TestHelpers
using ..ParcelCorpus: ParcelCorpus
using ..ParcelReference: ParcelReference
include("test_reference.jl")
end

module EquationsTests
using ..ParcelCorpus: ParcelCorpus
include("test_equations.jl")
end

module InterfaceTests
using ..ParcelCorpus: ParcelCorpus
include("test_interface.jl")
end

module RootSolutionsTests
using ..TestHelpers: TestHelpers
include("test_root_solutions.jl")
end

module PiecewiseLinearTests
using ..TestHelpers: TestHelpers
using ..ParcelCorpus: ParcelCorpus
include("test_mm2015piecewiselinear.jl")
end

module FixedTemperatureTests
using ..TestHelpers: TestHelpers
using ..ParcelCorpus: ParcelCorpus
using ..ParcelReference: ParcelReference
include("test_mm2015fixedT.jl")
end

module MM2015Tests
using ..TestHelpers: TestHelpers
using ..ParcelCorpus: ParcelCorpus
using ..ParcelReference: ParcelReference
include("test_mm2015.jl")
end

module InvariantTests
using ..TestHelpers: TestHelpers
using ..ParcelCorpus: ParcelCorpus
include("test_invariants.jl")
end

module InferenceAllocationTests
using ..TestHelpers: TestHelpers
using ..ParcelCorpus: ParcelCorpus
include("test_inference_allocations.jl")
end
