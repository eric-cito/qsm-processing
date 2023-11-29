using Pkg

Pkg.add("CUDA")
using CUDA

CUDA.set_runtime_version!(v"12.1")
