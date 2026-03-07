# build/build_sysimage.jl
#
# Creates a custom Julia system image with IgDiscover.jl precompiled.
# This eliminates ~90% of first-run latency.
#
# Usage:
#   julia --project=. build/build_sysimage.jl [--app]
#
# Options:
#   --app    Build a standalone application instead of a sysimage
#
# The sysimage can then be used:
#   julia --sysimage=build/IgDiscover.so --project=. bin/igdiscover.jl run
#
# Prerequisites:
#   julia -e 'using Pkg; Pkg.add("PackageCompiler")'

using PackageCompiler

const PROJECT_DIR = dirname(@__DIR__)
const BUILD_DIR   = joinpath(PROJECT_DIR, "build")
const SNOOP_SCRIPT = joinpath(BUILD_DIR, "snoop_workload.jl")

build_app = "--app" in ARGS

if build_app
    @info "Building standalone application..."
    create_app(
        PROJECT_DIR,
        joinpath(BUILD_DIR, "IgDiscover");
        executables=["igdiscover" => "julia_main"],
        precompile_execution_file=SNOOP_SCRIPT,
        force=true,
        include_lazy_artifacts=true,
    )
    @info "Application built at build/IgDiscover/"
else
    sysimage_path = joinpath(BUILD_DIR, Sys.iswindows() ? "IgDiscover.dll" : "IgDiscover.so")
    @info "Building system image at $sysimage_path..."
    create_sysimage(
        [:IgDiscover];
        sysimage_path=sysimage_path,
        precompile_execution_file=SNOOP_SCRIPT,
        project=PROJECT_DIR,
    )
    @info "System image built. Use with:"
    @info "  julia --sysimage=$sysimage_path --project=. bin/igdiscover.jl run"
end
