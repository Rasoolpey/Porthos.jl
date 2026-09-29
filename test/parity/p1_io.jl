# P1 gate: every JSON under cases/ loads (schema-valid, typed) and round-trips.

const CASES_DIR = joinpath(ROOT, "cases")

function all_case_files(dir)
    out = String[]
    for (root, _, files) in walkdir(dir)
        occursin("schema", relpath(root, dir)) && continue
        for f in files
            endswith(f, ".json") && push!(out, joinpath(root, f))
        end
    end
    return sort!(out)
end

function roundtrip_ok(path)
    x = load_json_input(path)
    mktempdir() do d
        out = joinpath(d, basename(path))
        # Scenario `system` paths are relative: keep the copy next to its system file.
        if x isa Scenario
            cp(x.system_path, joinpath(d, basename(x.system_path)); force = true)
            write_scenario(out, x)
        else
            write_case(out, x)
        end
        y = load_json_input(out)
        same = x isa Scenario ?
            (Porthos.json_identical(x.raw, y.raw) && x.solver == y.solver &&
             x.events == y.events && x.output_dir == y.output_dir) :
            x == y
        return same && Porthos.json_identical(Porthos.read_json(path), Porthos.read_json(out))
    end
end

@testset "P1 inputs load and round-trip" begin
    files = all_case_files(CASES_DIR)
    @test length(files) == 60
    nsys = nscen = 0
    for f in files
        x = load_json_input(f)
        x isa Case ? (nsys += 1) : (nscen += 1)
        @test roundtrip_ok(f)
        if x isa Scenario
            @test isfile(x.system_path)
        end
    end
    @test nsys > 0 && nscen > 0

    # The parameter expressions evaluate to the value C / Python compute.
    for f in files
        x = load_json_input(f)
        x isa Case || continue
        @test x.config.omega_b === 2.0 * π * x.config.fn
    end

    # Contracts: identical to the PHPS file in the pack's commit, all models readable.
    cs = load_contracts()
    @test cs.schema_version == Porthos.CONTRACT_SCHEMA_VERSION
    @test length(cs.entries) == 17
end

# Optional breadth check against every case in a PHPS checkout (read in place).
const PHPS_CASES = joinpath(get(ENV, "PORTHOS_PHPS_DIR", raw"C:\Users\em18736\Documents\PHPS_Opt"),
                            "phps", "cases")
if isdir(PHPS_CASES)
    @testset "P1 all PHPS cases load" begin
        files = all_case_files(PHPS_CASES)
        @test length(files) >= 200
        for f in files
            @test (load_json_input(f); true)
        end
    end
end
