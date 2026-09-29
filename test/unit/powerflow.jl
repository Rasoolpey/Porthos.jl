@testset "power flow" begin
    case = toy_case()
    pf = solve_powerflow(case)
    @test pf.converged
    @test pf.mismatch < 1e-6
    net = pf.spec.net
    @test pf.V[bus_index(net, 1)] == 1.02 && pf.theta[bus_index(net, 1)] == 0.0   # slack
    @test pf.V[bus_index(net, 3)] == 1.01                                          # PV
    S = bus_power(Matrix(ybus_pf(case)), pf.V, pf.theta)
    @test real(S[bus_index(net, 2)]) ≈ -0.5 atol = 1e-6
    @test imag(S[bus_index(net, 2)]) ≈ -0.2 atol = 1e-6
    @test real(S[bus_index(net, 3)]) ≈ 0.3 atol = 1e-6
    # A tighter tolerance takes at least as many iterations and stays on the solution.
    tight = solve_powerflow(case; tol = 1e-12)
    @test tight.converged && tight.iterations >= pf.iterations
    @test maximum(abs, tight.V - pf.V) < 1e-5
    # Not converged within the iteration budget.
    @test !solve_powerflow(case; max_iter = 1).converged

    # skip_pf_solve: voltages come from the Bus table
    skip = toy_case(extra = Dict("config" => Dict("skip_pf_solve" => true)))
    r = solve_powerflow(skip)
    @test r.from_case
    @test r.V == [1.02, 0.98, 1.0]         # buses 1, 2, 3; names default to "<id>"
end
