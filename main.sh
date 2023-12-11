#!/bin/bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


cd $SCRIPT_DIR

if [ -z "$XDG_RUNTIME_DIR" ]; then
    export XDG_RUNTIME_DIR="/tmp/"
    echo "XDG_RUNTIME_DIR was not set. Setting it to $XDG_RUNTIME_DIR"
fi

echo `pwd`
source env/bin/activate
python3.10 -m main