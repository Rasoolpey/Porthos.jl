# P0 gate: the parity-pack loader reads every file and checks every hash.

@testset "P0 parity pack" begin
    m = PACK.manifest
    @test Int(m[:pack_format]) == Porthos.PARITY_PACK_FORMAT
    @test length(string(m[:phps][:commit])) == 40
    # every file listed, readable, and hash-checked (load_parity_pack verified the hashes)
    for (rel, sha) in m[:files]
        @test isfile(joinpath(PACK.dir, string(rel)))
        endswith(string(rel), ".json") && @test Porthos.read_json(joinpath(PACK.dir, string(rel))) !== nothing
    end
    @test [c.name for c in Porthos.pack_cases(PACK)] ==
          ["base", "gfl", "gfl_zif", "vsm", "droop", "voc"]

    # tampering is detected
    mktempdir() do d
        cp(PACK.dir, joinpath(d, "pack"))
        dir = joinpath(d, "pack")
        for (root, _, files) in walkdir(dir), f in files
            chmod(joinpath(root, f), 0o644)
        end
        f = joinpath(dir, "cases.json")
        write(f, read(f, String) * " ")
        @test_throws Porthos.ParityPackError verify_parity_pack(dir)
        cp(PACK.dir, joinpath(d, "pack2"))
        dir2 = joinpath(d, "pack2")
        chmod(dir2, 0o755)
        write(joinpath(dir2, "extra.json"), "{}")
        @test_throws Porthos.ParityPackError verify_parity_pack(dir2)
    end
end
