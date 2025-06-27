#!/bin/bash


file_or_gz_exists() {
    for filename in "$@"; do
        # Check if the file exists
        if [ ! -z "$filename" ] && [ ! -f "$filename" ] && [ ! -f "$filename.gz" ]; then
            return 1  # Failure, at least one file does not exist
        fi
    done
    return 0  # Success, all files exist
}

GzFilepathIfOnlyGzFound() {
    # If the .gz version is found but the original is not, it returns the gz version
    if [ -f "$1" ]; then
        echo "$1"
    elif [ -f "$1.gz" ]; then
        echo "$1.gz"
    else
        echo "$1"
    fi
}

GzSafeMove(){
    # Moves the file if suffixes match
    # If they don't match in terms of .gz then gzip or gunzip is used
    # to write a new file, then the original removed

    loc_from=$1
    loc_to=$2

    if [[ $loc_from == *.gz ]]; then
        if [[ $loc_to == *.gz ]]; then
            mv $loc_from $loc_to
        else
            gunzip -c $loc_from > $loc_to
        fi
    else
        if [[ $loc_to == *.gz ]]; then
            gzip -c $loc_from > $loc_to
        else
            mv $loc_from $loc_to
        fi
    fi
}