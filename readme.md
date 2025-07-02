# QSM_processing 

This is a modified version of QSM-XT which can be containerised for HPC. It includes dicom conversion and extra skull stripping logic. It is designed to be less fuss - no bids wrangling - and work more reliably for the sequences optimised at UCSF.

We are not the QSM-XT developers. You can find the QSM-XT repo [here](https://github.com/QSMxT/QSMxT).

## Installation

Best done in docker.

### Docker

1. Download from git using `git clone`
1. cd to the directory you have just created
1. `./build-as-docker.sh`

### Apptainer

Alternatively, for HPC, build as apptainer/singularity. 

1. Download from git using `git clone`
1. cd to the directory you have just created
1. `./build-as-apptainer.sh`


Note that tis process simply pulls from DockerHub - changes you make to this codebase will not be reflected in the resulting apptainer container. 
See the build script for apptainer for more details and how to change this process so that you can make changes.


## Preparing your data

We work from dicoms. This is because some parameters are embedded in dicoms and this is supposed to form part of a longer, dicom-oriented, pipeline in our clinic.

* Create a dicom directory, with sub-dirs `qsm` and `t1` (case sensitive)
* Put the dicoms for your QSM into `qsm` and dicoms for the t1 in `t1`. You can use other contrasts - like a T2 or FLAIR image - but the directory should still be named `t1`.
* Create an output directory

These are ideally next to one another. e.g.

```
# Input Dicoms
mkdir -p /my-data/subj1-qsm/dicoms/t1
mkdir /my-data/subj1-qsm/dicoms/qsm

cp -r /some/path/to/dicoms/qsm-wholebrain /my-data/subj1-qsm/dicoms/qsm
cp -r /some/path/to/dicoms/mprage /my-data/subj1-qsm/dicoms/t1

# Output dir
mkdir /my-data/subj1-qsm/results
```

## Run it

### One Step

You need to run this bound to your data directory. Run the pipeline specifying the as dicom and output dirs relative to the docker's internal paths (e.g. /data/) 

```
docker run --entrypoint /bin/bash -it -v /path/to/your/data:/data/ --dicoms /data/rest-of-the-path-to-the-qsm-dicom-folder/ --out /data/your-output-directory
```

This might not display any output as it runs. If this is the case and you are running docker desktop, find it under containers and view the log for a live view. There can be a pause between dicom conversion and QSM-XT steps in this log.

Upon completion, your data appear in the output dir you created at the start. 


### Alternative: Run it interactively

#### Start the interactive terminal
You need to run this bound to your data directory.

Start by launching Docker:
```
docker run --entrypoint /bin/bash -it -v /path/to/your/data:/data/ qsm_qsmxt
```

For example

```
docker run --entrypoint /bin/bash -it -v /my-data/:/data/ qsm_qsmxt
```

#### Run it

Run the pipeline specifying the as dicom and output dirs relative to the docker's internal paths (i.e. /data/) 

```
./process-qsm.sh --dicoms /data/rest-of-the-path-to-the-qsm-dicom-folder/ --out /data/your-output-directory
```

e.g.

```
./process-qsm.sh --dicoms /data/subj1-qsm/dicoms --out /data/subj1-qsm/results/
```

#### Exit and Use Results

Exit the docker session

```
exit
```

Your data are now in the dir you created at the start. 

## Outputs

* `*_Chimap` - Processed QSM images. You typically do not want the singlepass result. See the qsmxt help for more information
* `t1` - Your T1 image converted from dicoms. Note this is not aligned to the QSM image
* `t1-brainmask` - brainmask for the T1 image, created with HD-BET
* `qsm_brainmask` - the brainmask applied to the QSM before it undergoes processing with qsm-xt
* `t1-to-qsm.mat` - an ANTs transform that will move the T1 into QSM space. Apply using AntsApplyTransform as needed

## Not a registered medical device

This repository and its associated code is not a registered medical device and has not undergone third-party testing or verification of any kind.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Licenses
Copyright 2024 Lee Reid and UCSF
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0