#!/bin/bash
# Run from the top directory of this solution
# Delete any local virtual environments before starting

# remove __pycache__
find . -type d -name __pycache__ -exec rm -r "{}" \;

rm -rf env oldenv conda-env

set -e

FN=qsm-inside-out.sif
rm -f ./$FN
# use the /run dir to avoid running out of space on tmp...
export SINGULARITY_TMPDIR=$XDG_RUNTIME_DIR
singularity build --fakeroot --bind `pwd`:/qsm /dev/shm/$FN build-singularity-inside-out.def
mv /dev/shm/$FN ./$FN

# Now build the env
singularity exec --cleanenv --bind /data $FN ./build_python.sh python3.10 env
chmod -R 777 env/bin/

