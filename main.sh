#!/bin/bash

set -e

current_script="$(realpath "${BASH_SOURCE[0]}")"
dir_sourceTop=$(dirname "$current_script")/

$dir_sourceTop/process-qsm.sh "$@"