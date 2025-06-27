#!/bin/bash

# srun --cpus-per-task=4 --time=2:00:00 --mem=64GB --pty /bin/bash

set -e

dir_script="$(dirname "$(readlink -f "$0")")"/
cd $dir_script


# remove __pycache__ etc
rm -rf env .env
find . -type d -name __pycache__ -exec rm -r "{}" \;


if [ -f qsm-processing.sif ]; then
    echo "Already exists at qsm-processing.sif"
    exit 0
fi

singularity -v build --ignore-subuid /dev/shm/qsm-processing.sif apptainer.def
mv /dev/shm/qsm-processing.sif ./qsm-processing.sif