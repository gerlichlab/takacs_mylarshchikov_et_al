#!/usr/bin/env bash -l

# Script to run PeakFlow on the CBE cluster using Nextflow
# Load necessary modules and set environment variables
ml build-env/f2022
module --ignore-cache load "nextflow/23.10.1"
export NXF_OPTS='-Xms1g -Xmx8g'
export NXF_TEMP="/scratch/nextflow/"$USER"/nxf_temp"
export NXF_WORK="/scratch/nextflow/"$USER"/nxf_work"
export NXF_ANSI_LOG=false


CURRENT_DIR=$(pwd)

# Run the peakflow pipeline
nohup nextflow -bg run dmitrymyl/peakflow -r main -profile cbe -params-file params.json