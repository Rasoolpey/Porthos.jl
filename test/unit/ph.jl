# Port passivity tools on transfer functions with known answers.

@testset "port passivity tools" begin
    # G(s) = 1/(s + 1): positive real
    lag = PortModel("lag", ["x"], "u", "y", fill(-1.0, 1, 1), [1.0], [1.0], 0.0)
    @test transfer(lag, 0.0) ≈ 1.0
    @test transfer(lag, 1.0im) ≈ 0.5 - 0.5im
    @test isempty(real_part_crossings(lag; sign = 1))
    pc = passivity_certificate(lag; sign = 1)
    @test pc.H0 ≈ 1.0 && !pc.nonpassive && pc.shortage == 0.0
    @test isempty(port_zeros(lag))

    # G(s) = (1 - s)/(1 + s) = 2/(s + 1) - 1: zero at +1, Re G(jw) = (1 - w^2)/(1 + w^2)
    allpass = PortModel("allpass", ["x"], "u", "y", fill(-1.0, 1, 1), [1.0], [2.0], -1.0)
    @test only(port_zeros(allpass)) ≈ 1.0
    @test only(real_part_crossings(allpass; sign = 1)) ≈ 1.0 rtol = 1e-12
    pc = passivity_certificate(allpass; sign = 1)
    @test pc.nonpassive
    @test pc.min_re ≈ -1.0 rtol = 1e-5
    @test pc.shortage ≈ 1.0 rtol = 1e-5
end
