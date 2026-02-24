#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --mincpus=8
#SBATCH --time=08:00:00
#SBATCH --mem=10gb
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

export CONTAINER="/groups/gerlich/labinfo/Methods/singularity_images/gerlichlab-replchromconf-jupyter-v0.2_master_test.sif"
export INPUT_PATH=$1
export OUTPUT_PATH=$2
singularity run $CONTAINER /bin/bash calc_scaling.sh $INPUT_PATH $OUTPUT_PATH
