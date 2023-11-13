#!/bin/bash

# source to conda
source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh
echo $_CONDA_EXE

# create folders
rm -r HD-bet_subvol
mkdir HD-bet_subvol

# separate subvolumes
# fslsplit $1 Magni -t
# mv Magni0*.nii.gz Magni_subvol

# run HD-BET
conda activate /working/lupolab/jingwen/conda/envs/QSM_DL

if [ "$2" -eq 1 ]
then
    echo "Using GPU"
    hd-bet -i Magni_subvol/ -o HD-bet_subvol/
else
    echo "Using CPU with simplified model"
    hd-bet -i Magni_subvol/ -o HD-bet_subvol/ -device cpu # -mode fast -tta 0
fi

conda deactivate

# combine subvolumes - take the min
fslmerge -t HD-bet_subvol/Magni_mask.nii.gz HD-bet_subvol/Magni*_mask.nii.gz
fslmaths HD-bet_subvol/Magni_mask.nii.gz -Tmin brain_mask_HD.nii.gz

# copy to output folder
# cp brain_mask_HD.nii.gz ../output/brain_mask_HD.nii.gz

