#!/bin/bash

set -e


print_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
    --dicoms DIR        Directory containing DICOM files
    --out DIR           Directory for output files
    --tmp DIR           Temporary directory (default: /tmp/)
    --realimag_ge       Expect real / imaginary dicom pairs from a GE scanner
    --help              Show this help message

Directory Requirements:
    The --dicoms directory must contain two folders:
      - t1
      - qsm
    The qsm directory can contain sub-directories if desired, but this is not necessary.
EOF
}

ParseArgs() {
    realimag_ge=0
    tmp=/tmp/

    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --dicoms) dir_dicoms="$2"; shift ;;
            --out) dir_out="$2"; shift ;;
            --tmp) tmp="$2"; shift ;;
            --realimag_ge) realimag_ge=1 ;;
            --help) print_help; exit 0 ;;
            *) echo "Unknown parameter passed: $1"; print_help; exit 1 ;;
        esac
        shift
    done

    # Validate directories
    [ -d "$dir_dicoms" ] || { echo "Error: --dicoms must be a directory"; print_help; exit 1; }
    [ -e "$dir_out" ] && [ ! -d "$dir_out" ] && { echo "Error: $dir_out exists and is not a directory"; print_help; exit 1; }
    [ -d "$tmp" ] || { echo "Error: --tmp must be a directory. ${tmp} is not a directory"; print_help; exit 1; }

    # Check for required subdirectories in dicoms
    
    [ -d "$dir_dicoms/qsm" ] || { echo "Error: --dicoms must contain a qsm directory"; print_help; exit 1; }
}

ConvertT1Dicoms(){
    if [ -e "$loc_t1" ]; then
        echo "Found $loc_t1 - skipping dicom conversion"
    elif [ -f "$dir_dicoms/t1/t1.nii" ]; then
        echo "Found preconverted T1 at $dir_dicoms/t1/t1.nii — copying"
        cp "$dir_dicoms/t1/t1.nii" "$loc_t1"
    else
        [ -d "$dir_dicoms/t1" ] || { echo "Error: --dicoms must contain a t1 directory"; print_help; exit 1; }
        dcm2niix -o "$dir_anat" -f t1 "$dir_dicoms/t1"
    fi   
}

ConvertQSMDicoms(){
    dcm2niix -o $dir_anat -f qsm_%f_%p_%t_%s $dir_dicoms/qsm

    if [[ "$realimag_ge" -eq 1 ]]; then
        echo "Running ConvertRealAndImaginaryToPhaseAndMag"
        ConvertRealAndImaginaryToPhaseAndMag
    fi
    echo "Running RenamePhaseMagNiftis"
    RenamePhaseMagNiftis  
}

RenamePhaseMagNiftis(){
    local echoCount=$(find "$dir_anat" -maxdepth 1 -type f -name '*_e[0-9].nii' | wc -l)
    for i in $(seq 1 "$echoCount"); do
        mv $dir_anat/qsm_*_e${i}.nii $dir_anat/sub-${subj}_echo-${i}_part-mag_MEGRE.nii
        mv $dir_anat/qsm_*_e${i}.json $dir_anat/sub-${subj}_echo-${i}_part-mag_MEGRE.json

        mv $dir_anat/qsm_*_e${i}_ph.nii $dir_anat/sub-${subj}_echo-${i}_part-phase_MEGRE.nii
        mv $dir_anat/qsm_*_e${i}_ph.json $dir_anat/sub-${subj}_echo-${i}_part-phase_MEGRE.json
    done
}

ConvertRealAndImaginaryToPhaseAndMag(){
    wd=$(pwd)
    cd "$dir_anat"
    local echoCount=$(find "$dir_anat" -maxdepth 1 -type f -name '*_e[0-9].nii' | wc -l)
    for i in $(seq 1 "$echoCount"); do
        # Just point the phase json side car to the mag json sidecar
        ln qsm_*_e${i}.json qsm_*_e${i}_ph.json

        # Convert image/real to phase
        # NB mag will have been created already assuming we are using GE/dcm2niix
        python "$dir_sourceTop"/imaginary-real-to-phase-mag.py qsm_*_e${i}_real.nii \
                                                            qsm_*_e${i}_imaginary.nii \
                                                            qsm_${subj}_e${i}_ph.nii                                           
    done
    cd "$wd"
}


GenerateQSMBrainmask(){

    if [ -e "$loc_qsm_brainmask" ]; then
        echo "QSM brainmask found. Not regenerated"
        return
    fi

    Skullstrip_HDBET_Quick $loc_t1 $loc_t1_brainmask


    if [ -e $loc_t1ToQsmMat ]; then
        echo "Found T1 to QSM transform. Registration skipped"
    else
        cd $tmp
        antsRegistration -d 3 \
            -m MI[ "$loc_qsm_mag_echo1",$loc_t1,1,32,Regular,0.25 ] \
            -c [ 1000x500x250x0,1e-7,5 ] \
            -t Rigid[ 0.1 ] -f 8x4x2x1 -s 4x2x1x0 \
            --use-histogram-matching 1 \
            -z 1 \
            --winsorize-image-intensities [ 0.005, 0.995 ] \
            -x [ , $loc_t1_brainmask ] \
            -o output

        mv output0GenericAffine.mat $loc_t1ToQsmMat
    fi

    antsApplyTransforms --transform $loc_t1ToQsmMat \
                        -i $loc_t1_brainmask \
                        -o $loc_qsm_brainmask \
                        -r "$loc_qsm_mag_echo1" \
                        -u char \
                        -n NearestNeighbor
}

CropQSMToBrainmask(){
    echo "Generating QSM brainmask"
    GenerateQSMBrainmask
    echo "Finished generating QSM brainmask"
    # Crop skulls from QSM
    local echoCount=$(find "$dir_anat" -maxdepth 1 -type f -name '*echo-[0-9]_part-mag_MEGRE.nii' | wc -l)
    for i in $(seq 1 "$echoCount"); do
        python "$dir_sourceTop"/apply-mask.py "$dir_anat/sub-${subj}_echo-${i}_part-mag_MEGRE.nii" $loc_qsm_brainmask "$dir_anat/sub-${subj}_echo-${i}_part-mag_MEGRE.nii" 
        python "$dir_sourceTop"/apply-mask.py "$dir_anat/sub-${subj}_echo-${i}_part-phase_MEGRE.nii" $loc_qsm_brainmask "$dir_anat/sub-${subj}_echo-${i}_part-phase_MEGRE.nii" 
    done
}

current_script="$(realpath "${BASH_SOURCE[0]}")"
dir_sourceTop=$(dirname "$current_script")/
source "$dir_sourceTop"/exe-paths.sh
source "$dir_sourceTop"/hd-bet.sh


ParseArgs "$@"

# Create a temporary directory and add a trap to clean up
tmp=$(mktemp -d "${tmp}/tmpdir-XXXXXX")
trap 'rm -rf "$tmp"' EXIT

dir_bids="${tmp}"/bids
subj=mysubj
dir_anat=$dir_bids/sub-$subj/anat
loc_t1=$dir_anat/t1.nii
loc_qsm_mag_echo1=$dir_anat/sub-${subj}_echo-1_part-mag_MEGRE.nii
loc_t1_brainmask=$dir_out/t1-brainmask.nii.gz
loc_t1ToQsmMat=$dir_out/t1-to-qsm.mat
loc_qsm_brainmask=$dir_out/qsm-brainmask.nii.gz

mkdir -p $dir_anat
mkdir -p $dir_out

echo "Converting T1 dicoms"
ConvertT1Dicoms

echo "Converting QSM dicoms"
ConvertQSMDicoms

echo "Cropping QSM to brainmask"
CropQSMToBrainmask

# Note that this will calculate its own mask. This is important for acquisitions with
# dark slices at the boundaries, like on siemens

qsmxt $dir_bids --premade 'gre' --auto_y

# Copy out results
echo "Copying results to $dir_out"
mv $dir_bids/derivatives/qsmxt-*/sub-*/anat/* "$dir_out"
mv $dir_bids/*/anat/t1.nii "$dir_out"