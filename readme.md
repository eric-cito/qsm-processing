# QSM_processing 

This is a modified version of QSM-XT which can be containerised for HPC. It includes dicom conversion and extra skull stripping logic.

## Installation

Best done in docker.

### Docker

1. Download from git using `git clone`
1. cd to the directory you have just created
1. `./build-as-docker.sh`

### Apptainer

Alternatively, for HPC, build as apptainer/singularity. This simply pulls from DockerHub. See the build script for apptainer for more details.

1. Download from git using `git clone`
1. cd to the directory you have just created
1. `./build-as-apptainer.sh`


## Running

## Prepare

We work from dicoms. This is because outputs may need to be dicoms.

* Create a dicom directory, with sub-dirs `qsm` and `T1` (case sensitive)
* Put the dicoms for your QSM into `qsm` and dicoms for the t1 in `T1`.
* Create an output directory

These are ideally next to one another. e.g.

```
# Input Dicoms
mkdir -p /my-data/subj1-qsm/dicoms/T1
mkdir /my-data/subj1-qsm/dicoms/qsm

cp -r /some/path/to/dicoms/qsm-wholebrain /my-data/subj1-qsm/dicoms/qsm
cp -r /some/path/to/dicoms/mprage /my-data/subj1-qsm/dicoms/T1

# Output dir
mkdir /my-data/subj1-qsm/results
```

## Start the interactive terminal
You need to run this bound to your data directory.

Start by lanching Docker:
```
docker run --entrypoint /bin/bash -it -v /path/to/your/data:/data/ qsm_qsmxt
```

For example

```
docker run --entrypoint /bin/bash -it -v /my-data/:/data/ qsm_qsmxt
```

## Run it

Run the pipeline specifying the as dicom and output dirs relative to the docker's internal paths (i.e. /data/) 

```
./process-qsm.sh --dicoms /data/rest-of-the-path-to-the-qsm-dicom-folder/ --out /data/your-output-directory
```

e.g.
s
```
./process-qsm.sh --dicoms /data/subj1-qsm/dicoms --out /data/subj1-qsm/results/
```

## Exit and Use Results

Exit the docker session

```
exit
```

Your data are now in the dir you created at the start. Hunt around and try to find the QSM result `*_Chimap.nii`. You can delete everything else.

### Not a registered medical device

This repository and its associated code is not a registered medical device and has not undergone third-party testing or verification of any kind.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Licenses
Copyright 2024 Lee Reid and UCSF
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0