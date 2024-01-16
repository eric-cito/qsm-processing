#!/bin/bash

set -e

#singularity run qsm.sif 
singularity exec --cleanenv --bind /data $(dirname $0)/qsm-inside-out.sif $(dirname $0)/main.sh "$@"