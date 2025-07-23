#!/bin/bash
# Note that there is no native apptainer build because it is too problematic.
# Instead, it builds from docker. If the docker container cannot be found then
# You will need to build the docker image yourself (see build-as-docker), altering the
# tag and where it is pushed to to your own docker account, then also alter the
# apptainer definition to match

# set -e

# dir_script="$(dirname "$(readlink -f "$0")")"/
# cd $dir_script

# if [ -f qsm-processing.sif ]; then
#     echo "Already exists at qsm-processing.sif"
#     exit 0
# fi

# # Due to I/O limitations in our own HPC this builds to /dev/shm and is copied back.
# # You can ammend this for your own system as required 
# singularity -v build --ignore-subuid /data/morrison/wip/eric/modular-image-processing-system/qsm-processing-fix.sif apptainer.def
# ##mv /dev/shm/qsm-processing.sif ./qsm-processing.sif

### TESTING THIS CODE
set -e

dir_script="$(dirname "$(readlink -f "$0")")"/
cd "$dir_script"

# Use custom tmp directory to avoid /dev/shm I/O limitations
tmp_build_dir="/data/morrison/wip/eric/tmp-fix"
mkdir -p "$tmp_build_dir"

final_sif_path="/data/morrison/wip/eric/modular-image-processing-system/qsm-processing-fix.sif"

if [ -f "$final_sif_path" ]; then
    echo "Already exists at $final_sif_path"
    exit 0
fi

echo "Building .sif using $tmp_build_dir..."
singularity -v build --ignore-subuid "$tmp_build_dir/qsm-processing-fix.sif" apptainer.def

echo "Copying .sif to final location..."
cp "$tmp_build_dir/qsm-processing-fix.sif" "$final_sif_path"

echo "Build complete: $final_sif_path"