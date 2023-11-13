#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --output=%x-%j.out

cd /home/jyao3/030_QSM/01_Code/Manuscript_QSM_comparison/
matlab -nodisplay -r "run_param_optimization('/working/lupolab/jingwen/001_QSM/Test_data/b4489_t12374'); exit"
