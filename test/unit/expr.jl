@testset "expression reader" begin
    @test parse_param_expr("2.0 * M_PI * 60.0") === 2.0 * π * 60.0
    @test parse_param_expr("2.0 * M_PI * 60.0") === (2.0 * Float64(π)) * 60.0
    @test parse_param_expr("1.5") === 1.5
    @test parse_param_expr("  -3 ") === -3.0
    @test parse_param_expr("1e-5") === 1e-5
    @test parse_param_expr(".5E+2") === 50.0
    @test parse_param_expr("1 - 2 - 3") === -4.0             # left associative
    @test parse_param_expr("8 / 4 / 2") === 1.0
    @test parse_param_expr("2 + 3 * 4") === 14.0             # precedence
    @test parse_param_expr("(2 + 3) * 4") === 20.0
    @test parse_param_expr("-2 * -3") === 6.0
    @test parse_param_expr("sqrt(4.0) * pi") === 2.0 * π
    @test parse_param_expr("PMIN*R_base", Dict("PMIN" => 0.5, "R_base" => 3.0)) === 1.5
    @test param_value(2) === 2.0
    @test param_value(0.25) === 0.25
    @test param_value("M_PI / 2") === π / 2

    for bad in ("", "2 *", "(1 + 2", "1 + 2)", "x", "exp(1)", "2 ** 3", "1; 2", "`rm`",
                "Base.run(`x`)", "2 3")
        @test_throws Porthos.ParamExprError parse_param_expr(bad)
    end
    @test_throws ArgumentError param_value(true)
    @test_throws ArgumentError param_value(nothing)
end
