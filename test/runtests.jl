using Porthos
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))

include("parity/common.jl")

@testset "Porthos" begin
    @testset "unit" begin
        include("unit/expr.jl")
        include("unit/io.jl")
        include("unit/network.jl")
        include("unit/powerflow.jl")
    end
    @testset "parity" begin
        include("parity/p0_pack.jl")
        include("parity/p1_io.jl")
        include("parity/p2_network.jl")
        include("parity/p3_powerflow.jl")
    end
end
