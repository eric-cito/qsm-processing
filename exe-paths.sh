#!/bin/bash
# This sets PATH so that the correct executables can be called without 
# absolute paths being used. It also activates the local version of python

# You can edit these paths to override use of the installation created
# with install.sh

# Get the absolute path of the script being executed
current_script="$(realpath "${BASH_SOURCE[0]}")"
dir_sourceTop=$(dirname "$current_script")/
dir_hdbet=$dir_sourceTop"HD-BET/HD_BET/"
dir_ants=$dir_sourceTop"ants/bin/"
dir_python_env=$dir_sourceTop"/env/bin"
loc_activate_python=$dir_python_env/activate
loc_python=$dir_python_env/python


# Set path so calls go to internal binaries, not system-installed items
export PATH=$dir_python_env:$dir_hdbet:$dir_ants:$PATH

# Activate the virtual environment    
source $loc_activate_python
