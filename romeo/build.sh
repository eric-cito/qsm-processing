set -e

#echo "THIS SCRIPT REQUIRES CUDA 11+ AND DOES NOT INSTALL IT"

#./julia-1.9.4/bin/julia ] activate .
mkdir romeo -p

wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.4-linux-x86_64.tar.gz
tar zxvf julia-1.9.4-linux-x86_64.tar.gz
rm julia-1.9.4-linux-x86_64.tar.gz

./julia-1.9.4/bin/julia ./build.jl