#!/bin/bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


cd $SCRIPT_DIR

if [ -z "$XDG_RUNTIME_DIR" ]; then
    export XDG_RUNTIME_DIR="/tmp/"
    echo "XDG_RUNTIME_DIR was not set. Setting it to $XDG_RUNTIME_DIR"
fi

# When julia runs romeo, it can install things if they're not installed yet
# We want them here, not in ~/ because home fills up and may clash with other
# solutions
export JULIA_DEPOT_PATH=$SCRIPT_DIR/romeo/julia-packages

echo `pwd`
source env/bin/activate
python3.10 -m main