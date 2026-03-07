#!/usr/bin/env julia
# IgDiscover.jl command-line interface
#
# Usage:
#   julia --project=. bin/igdiscover.jl init <dir> <database> <reads>
#   julia --project=. bin/igdiscover.jl run [<dir>]
#
# The standalone app (build/IgDiscover/bin/igdiscover) calls IgDiscover.julia_main() directly.

using Pkg
script_dir = @__DIR__
project_dir = dirname(script_dir)
Pkg.activate(project_dir)

using IgDiscover

exit(IgDiscover.julia_main())
