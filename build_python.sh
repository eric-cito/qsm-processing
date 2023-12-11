#!/bin/bash

# Currently builds a python environment that is suitable for this
# Does not build an image, and assumes python is available
# Call from the directory containing this project
# And call the following if needed for python 3.9:
#   conda init
#   conda activate # If needed for python 3.9


# conda create -p `pwd`/conda-env/ -y
# conda activate `pwd`/conda-env
# conda install cuda -c nvidia -y

set -e


if [ -z "$1" ]; then
    python="python3.9"
else
    python="$1"
fi


if [ -z "$2" ]; then
    envName="env"
else
    envName="$2"
fi


# Uncomment if pip complains it can't find a version
# pip3 install --upgrade pip --user


$python -m venv $envName
source env/bin/activate

#python3 -m pip install --upgrade pip
$python -m pip install -r requirements.txt

echo 'Installing PythonUtils'
if [ ! -d "PythonUtils" ]; then
    git clone git@git.ucsf.edu:lee-reid/PythonUtils.git
fi

PythonUtils/build.sh $python
