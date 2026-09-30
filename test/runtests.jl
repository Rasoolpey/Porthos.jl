using Porthos
using Test
using SparseArrays
using LinearAlgebra
using ForwardDiff
using Random

const ROOT = normpath(joinpath(@__DIR__, ".."))

include("parity/common.jl")

@testset "Porthos" begin
    @testset "unit" begin
        include("unit/expr.jl")
        include("unit/io.jl")
        include("unit/network.jl")
        include("unit/powerflow.jl")
        include("unit/components.jl")
        include("unit/assembly.jl")
        include("unit/powerfactory.jl")
        include("unit/ph.jl")
        include("unit/roa.jl")
    end
    @testset "parity" begin
        include("parity/p0_pack.jl")
        include("parity/p1_io.jl")
        include("parity/p2_network.jl")
        include("parity/p3_powerflow.jl")
        include("parity/p4_components.jl")
        include("parity/p5_dae.jl")
        include("parity/p6_init.jl")
        include("parity/p7_sim.jl")
        include("parity/p10_ph.jl")
    end
end
