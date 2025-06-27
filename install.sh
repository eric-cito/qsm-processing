#!/bin/bash
set -e

install_python() {

    which python
    echo "--------------"
    #apt-get install python3.10 python3.10-venv -y

    # Create a Python virtual environment
    if [ ! -d env ]; then
        # Create a Python virtual environment
        python3 -m venv env
    fi


    # Activate the virtual environment
    chmod -R 700 ./env/bin/
    source $(pwd)/env/bin/activate

    python3 -m pip install -r requirements.txt

}

# Function to check if a command is available
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install unzip using apt-get
install_unzip() {
    if ! command_exists unzip; then
        apt-get install -y unzip
    fi
}

install_hd_bet() {
    # Python must be activated first and requirements installed(!)
    echo "Downloading HD-BET models"

    # Download the model ahead of time
    mkdir models
    wget -O models/0.model https://zenodo.org/record/2540695/files/0.model?download=1 
    #wget -O models/1.model https://zenodo.org/record/2540695/files/1.model?download=1 
    #wget -O models/2.model https://zenodo.org/record/2540695/files/2.model?download=1 
    #wget -O models/3.model https://zenodo.org/record/2540695/files/3.model?download=1 
    #wget -O models/4.model https://zenodo.org/record/2540695/files/4.model?download=1 

    cd ..
    
}

install_ants() {

    if [ -d ants ]; then
        echo "Ants installation found. Delete ants directory and re-run script to reinstall"
        return 0
    fi

    ln -s /opt/ants-2.4.3/ $dir_script/ants
}

# Ensure we are in the dir of this script
dir_script="$(dirname "$(readlink -f "$0")")"/
cd $dir_script

#apt install -y software-properties-common
#add-apt-repository -y 'ppa:deadsnakes/ppa'
#apt-get update -y
#apt-get upgrade -y
#apt-get install dcm2niix wget git -y
install_python
install_hd_bet
install_unzip
install_ants

# Clean some things we will never use
#rm -r $FASTSURFER_HOME/*
#apt-get clean

echo "Install Complete"