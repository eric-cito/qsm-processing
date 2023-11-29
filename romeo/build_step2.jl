using Pkg
using CUDA

Pkg.add(url="https://github.com/korbinian90/RomeoApp.jl")

Pkg.activate(".")
Pkg.develop(url="https://github.com/korbinian90/CompileMRI.jl")
Pkg.status()


# end
Pkg.build("CompileMRI")

using CompileMRI
compile("./romeo"; apps=["romeo"], force=true)