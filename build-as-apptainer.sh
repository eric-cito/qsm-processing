#!/bin/bash
# Note that there is no native apptainer build because it is too problematic.
# Instead, it builds from docker. If the docker container cannot be found then
# You will need to build the docker image yourself (see build-as-docker), altering the
# tag and where it is pushed to to your own docker account, then also alter the
# apptainer definition to match

set -e

dir_script="$(dirname "$(readlink -f "$0")")"/
cd $dir_script

if [ -f qsm-processing-fix.sif ]; then
    echo "Already exists at qsm-processing-fix.sif"
    exit 0
fi

# Due to I/O limitations in our own HPC this builds to /dev/shm and is copied back.
# You can ammend this for your own system as required 
singularity -v build --ignore-subuid /dev/shm/qsm-processing-fix-2.sif apptainer.def
mv /dev/shm/qsm-processing-fix-2.sif ./qsm-processing-fix-2.sif