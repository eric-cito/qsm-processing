#!/bin/bash

# Read from specified file, or from standard input
infile="${1:-/dev/stdin}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo

done < $infile