# Bind a generated parity pack as the lazy artifact `parity_pack` in Artifacts.toml.
#
#   julia --project scripts/bind_parity_pack.jl parity/pack v1
#
# Copies the pack into the local artifact store (content-addressed by its git tree hash),
# writes the tarball parity/dist/parity_pack-<version>.tar.gz, and binds it with the
# download URL of the GitHub release `parity-pack-<version>`. Upload the tarball to that
# release (by hand; it is outward-facing) so CI and fresh clones can fetch it.
#
# A new pack is a new baseline: add its reason to parity/README.md.

using Pkg, Pkg.Artifacts, SHA, Tar, p7zip_jll, Porthos

length(ARGS) == 2 || error("usage: bind_parity_pack.jl <pack dir> <version, e.g. v1>")
packdir, version = abspath(ARGS[1]), ARGS[2]

manifest = verify_parity_pack(packdir)          # never bind an inconsistent pack
root = normpath(joinpath(@__DIR__, ".."))
toml = joinpath(root, "Artifacts.toml")
dist = joinpath(root, "parity", "dist")
mkpath(dist)

h = create_artifact() do dir
    for entry in readdir(packdir)
        cp(joinpath(packdir, entry), joinpath(dir, entry))
    end
end
tarball = joinpath(dist, "parity_pack-$version.tar.gz")
tar_sha = archive_artifact(h, tarball)

# On Windows the directory hash from create_artifact can differ from the hash of the tarball
# contents (file modes), and Linux verifies the latter after download. Bind the tarball's
# tree hash, and install the local copy under that hash.
rawtar = tempname() * ".tar"
run(pipeline(`$(p7zip_jll.p7zip()) x -so $tarball`, stdout = rawtar, stderr = devnull))
h_tar = Base.SHA1(Tar.tree_hash(rawtar))
if h_tar != h
    dest = artifact_path(h_tar)
    isdir(dest) || Tar.extract(rawtar, dest)
    remove_artifact(h)
    h = h_tar
end
rm(rawtar)
url = "https://github.com/Rasoolpey/Porthos.jl/releases/download/parity-pack-$version/" *
      basename(tarball)
bind_artifact!(toml, Porthos.PARITY_ARTIFACT, h; download_info = [(url, tar_sha)],
               lazy = true, force = true)

println("bound $(Porthos.PARITY_ARTIFACT) = $(bytes2hex(h.bytes))")
println("PHPS commit  $(manifest[:phps][:commit])")
println("tarball      $tarball")
println("sha256       $tar_sha")
println("upload to    $url")
