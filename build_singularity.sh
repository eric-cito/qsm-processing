#!/bin/bash
# Run from the top directory of this solution
# Delete any local virtual environments before starting or
# suffer a long build time


# remove __pycache__
find . -type d -name __pycache__ -exec rm -r "{}" \;

rm -rf env oldenv conda-env

set -e
# rm -f ./fix-legui.sif
# singularity build --fakeroot /dev/shm/fix-legui.sif build-as-apptainer/definition.def
# mv /dev/shm/fix-legui.sif ./fix-legui.sif

rm -f ./qsm.sif
# use the /run dir to avoid running out of space on tmp...
export SINGULARITY_TMPDIR=$XDG_RUNTIME_DIR
singularity build --fakeroot /dev/shm/qsm.sif build-singularity.def
mv /dev/shm/qsm.sif ./qsm.sif

