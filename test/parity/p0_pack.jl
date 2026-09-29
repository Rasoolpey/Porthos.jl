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

    # tampering is detected. The copies are written as new files rather than `cp`'d: `cp`
    # keeps the artifact's read-only permissions, and `chmod` on them failed on Windows
    # (EINVAL, Julia 1.12.7) depending on how the artifact had been extracted.
    function writable_copy(dst)
        for (root, _, files) in walkdir(PACK.dir), f in files
            target = joinpath(dst, relpath(joinpath(root, f), PACK.dir))
            mkpath(dirname(target))
            write(target, read(joinpath(root, f)))
        end
        return dst
    end
    mktempdir() do d
        dir = writable_copy(joinpath(d, "pack"))
        @test verify_parity_pack(dir) !== nothing
        f = joinpath(dir, "cases.json")
        write(f, read(f, String) * " ")
        @test_throws Porthos.ParityPackError verify_parity_pack(dir)
        dir2 = writable_copy(joinpath(d, "pack2"))
        write(joinpath(dir2, "extra.json"), "{}")
        @test_throws Porthos.ParityPackError verify_parity_pack(dir2)
    end
end
