#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=gpu
#SBATCH --gpus=2
#SBATCH --output=%x-%j.out

cd /home/jyao3/030_QSM/01_Code/
matlab -nodisplay -r "run_QSM_COSMOSlist; exit"
