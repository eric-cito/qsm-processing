set -e

#echo "THIS SCRIPT REQUIRES CUDA 12.1 AND DOES NOT INSTALL IT"

#./julia-1.9.4/bin/julia ] activate .
mkdir romeo -p



pathToJulia=/mnt/nfs/rad/apps/share/versions/julia/julia-1.9.4/bin/julia

export JULIA_DEPOT_PATH=`pwd`/julia-packages

# Uncomment to use on an older node where the above doesn't work
# pathToJulia=/netopt/rhel7/versions/julia/julia-1.9.4/bin/julia
# if ! [ -f $pathToJulia ]; then

#     echo Julia not found at $pathToJulia. Creating a local copy
#     wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.4-linux-x86_64.tar.gz
#     tar zxvf julia-1.9.4-linux-x86_64.tar.gz
#     rm julia-1.9.4-linux-x86_64.tar.gz

#     pathToJulia=./julia-1.9.4/bin/julia

# fi

$pathToJulia ./build_step1.jl
$pathToJulia ./build_step2.jl

# #
# conda activate ../conda-env
# conda clean --all -y
# setenv CONDA_PKGS_DIRS /tmp/reidconda
# conda install -c conda-forge julia -y
# rm -rf /tmp/reidconda

