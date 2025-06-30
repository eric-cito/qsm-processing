#!/bin/bash

set -e

dir_script="$(dirname "$(readlink -f "$0")")"/
cd $dir_script

docker build . --tag leereid/qsmxt-plus

# This will give you access denied unless you are Lee, and logged in
# so it is commented out by default
# docker push leereid/qsmxt-plus:latest
