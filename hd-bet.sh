source $(dirname "$BASH_SOURCE[0]")/file-or-gz.sh


Skullstrip_HDBET_Quick() {
    local loc_in=$(GzFilepathIfOnlyGzFound "$1")
    local loc_mask=$2
    python qsm/CalcBrainmask.py $loc_in $loc_mask

    if [ "$#" -ge 3 ]; then
        local loc_brain="$3"
        python apply-mask "$loc_in" "$loc_mask" "$loc_brain"
    fi
}
