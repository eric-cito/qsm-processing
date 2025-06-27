#!/bin/bash

set -e

dir_script="$(dirname "$(readlink -f "$0")")"/
cd $dir_script

docker build . --tag leereid/qsmxt-plus