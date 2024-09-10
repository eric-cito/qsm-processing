source $(dirname "$BASH_SOURCE[0]")/file-or-gz.sh


Skullstrip_HDBET() {
    local loc_in=$(GzFilepathIfOnlyGzFound "$1")
    local loc_mask=$2
    if [ "$#" -ge 3 ]; then
        local loc_brain="$3"
    else
        local loc_brain=""
    fi

    Skullstrip_HDBET_Customisable $loc_in $loc_mask "-device cpu" $loc_brain
}

Skullstrip_HDBET_Quick() {
    local loc_in=$(GzFilepathIfOnlyGzFound "$1")
    local loc_mask=$2
    if [ "$#" -ge 3 ]; then
        local loc_brain="$3"
    else
        local loc_brain=""
    fi

    Skullstrip_HDBET_Customisable $loc_in $loc_mask "-device cpu -mode fast -tta 0" $loc_brain
}


Skullstrip_HDBET_Customisable() {
    local loc_in=$(GzFilepathIfOnlyGzFound "$1")
    local loc_mask=$2
    local hdbetargs=$3
    if [ "$#" -ge 4 ]; then
        local loc_brain="$4"
    fi

    if file_or_gz_exists "$loc_mask" "$loc_brain"; then
        echo $loc_mask and "$loc_brain" already exist. Skullstrip skipped.
        return 0
    fi

    local dir_tmp=`mktemp -d`
    trap "rm -rf $dir_tmp" EXIT

    cd $dir_tmp

    local fn_in=$(basename $loc_in)
    ln -s $loc_in $fn_in


    hd-bet -i $fn_in $hdbetargs

    #local filename_no_suffix="${fn_in%.*}"
    GzSafeMove *_bet_mask.nii.gz $loc_mask
    
    if [ -n "$loc_brain" ]; then
        GzSafeMove *_bet.nii.gz $loc_brain
    fi
}
