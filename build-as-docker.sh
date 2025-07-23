#!/bin/bash

set -e

dir_script="$(dirname "$(readlink -f "$0")")"/
cd "$dir_script"

# Build the Docker image using your custom tag
docker build . --tag ericcito/qsmxt-plus

# Uncomment the next line if you later want to push to your DockerHub
# docker push ericcito/qsmxt-plus:latest