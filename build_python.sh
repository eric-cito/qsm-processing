#!/bin/bash

# Currently builds a python environment that is suitable for this
# Does not build an image, and assumes python is available
# Call from the directory containing this project
# And call the following if needed for python 3.9:
#   conda init
#   conda activate # If needed for python 3.9


set -e

# Uncomment if pip complains it can't find a version
# pip3 install --upgrade pip --user


python3.9 -m venv env
source env/bin/activate

#python3 -m pip install --upgrade pip
python3.9 -m pip install -r requirements.txt

echo 'Installing PythonUtils'
if [ ! -d "PythonUtils" ]; then
    git clone git@git.ucsf.edu:lee-reid/PythonUtils.git
fi

PythonUtils/build.sh

# # HD-BET
# git clone git@git.ucsf.edu:lee-reid/HD-BET-for-python.git HDBET
# bash HDBET/build.sh
# HD-BET
if [ ! -d "HDBET" ]; then
    git clone git@git.ucsf.edu:lee-reid/HD-BET-for-python.git HDBET
fi

HDBET/build.sh