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

# # HD-BET
# git clone https://github.com/MIC-DKFZ/HD-BET
# mv HD-BET HDBET
# cd HDBET
# touch __init__.py
# echo "folder_with_parameter_files = os.path.join('`pwd`', 'hd-bet_params') # Override by Lee" >> HD_BET/paths.py
# python3.9 -m pip install -e .
# cd ..